import sys
import tskit, pyslim, msprime
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import false_discovery_control
from scipy.cluster.vq import kmeans
pd.options.mode.copy_on_write = True
# Julia stuff
from juliacall import Main as jl
jl.seval('import Pkg; Pkg.activate(".")')
jl.seval("using GenomicOffsets")
SAMPLE_SIZE = 1000

def centroids(df, nexpected):
    locs = df[df.fraction == 0.0][["x", "y"]]
    K = locs.shape[0] // nexpected
    centroids,_ = kmeans(locs,K)
    return K, centroids

def assign_population(df, centroids):
    points = df[["x", "y"]].to_numpy()
    # Calculate the squared distances from each point to each centroid
    squared_distances = np.sum((points[:, np.newaxis, :] - centroids[np.newaxis, :, :]) ** 2, axis=2)
    # Find the index of the centroid with the minimum distance for each point
    ks = np.argmin(squared_distances, axis=1)
    return ks

def print_output(file, fraction, category, dataset, _type, cor, radj2):
    pd.DataFrame([{
            'file': file,
            'fraction': fraction,
            'category': category,
            'dataset' : dataset,
            'type': _type,
            'statistic': cor.statistic,
            'pvalue': cor.pvalue,
            'radj2': radj2
        }]).to_csv(sys.stdout, index=False, header=False)
    sys.stdout.flush()

def additive_genotype_matrix(smts, individuals, maf):
    # First, we identify sites
    sites = [s for s, locus in enumerate(smts.variants()) if locus.frequencies()[""] < 1-maf]
    Y = np.zeros(shape=(len(individuals), len(sites)))
    # Get list of samples
    samples = list()
    for ind in individuals:
        samples.extend(smts.individual(ind).nodes)
    # Now, create the genotype matrix
    variant = tskit.Variant(smts, samples=samples)
    for i, site_id in enumerate(sites):
        variant.decode(site_id)
        for k in range(Y.shape[0]):
            Y[k, i] = variant.genotypes[2*k]+variant.genotypes[2*k+1]
    return Y

def genomic_offset(model, X, Xpred, candidates=None):
    if len(candidates) == 0:
        return np.zeros(X.shape[0])
    try:
        return np.array(jl.genomic_offset(model, X, Xpred, candidates)).reshape(X.shape[0])
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr, flush=True)
        return np.zeros(X.shape[0])

def adjusted_r2(r2, n, p):
    return 1 - (1 - r2) * (n - 1) / (n - p - 1)

def measure_association(x, y):
    cor = pearsonr(x, y)
    return cor, adjusted_r2(cor.statistic**2, len(x), 1)

infile = sys.argv[1]

ts = tskit.load(infile)
# Sanity check: all individuals coalesce
assert max(t.num_roots for t in ts.trees()) == 1
df = pd.DataFrame.from_dict(ts.metadata["SLiM"]["user_metadata"])
K, centroids = centroids(df, 50)
print(f"Number of clusters: {K}", file=sys.stderr, flush=True)

# Print header
print("file,fraction,category,dataset,statistic,pvalue,r2adj", file=sys.stdout, flush=True)
fractions = np.array([elm for elm in set(df.fraction)])
fractions = fractions[fractions < 1.0]
# Sample 20 fractions with replacement
np.random.seed(123)
fractions = np.random.choice(fractions, 40, replace=True)
for chosen_fraction in fractions:
    try:
        # Sampling genotypes
        sampled_df = df[df.fraction == chosen_fraction]
        print(f"Fraction: {chosen_fraction}", file=sys.stderr, flush=True)
        generations_ago = np.max(df.timePoint) - sampled_df.timePoint.iloc[0]

        ## Simplify the tree sequence
        alive_inds = pyslim.individuals_alive_at(ts, generations_ago)
        keep_nodes = list()
        for i in alive_inds:
            keep_nodes.extend(ts.individual(i).nodes)
        sts = ts.simplify(keep_nodes, keep_input_roots=True)
        # Now, we have to sort alive_inds in the same order as they are in the sampled_df
        alive_inds = pyslim.individuals_alive_at(sts, generations_ago)
        assert len(alive_inds) == sampled_df.shape[0]
        pedigrees = [sts.individual(ind).metadata["pedigree_id"] for ind in alive_inds]
        sampled_df = sampled_df.set_index("ids").loc[pedigrees].reset_index()
        assert np.all(pedigrees == sampled_df.ids)            

        # Sample if the measure is causal or empirical ["causal", "empirical"]
        _type = np.random.choice(["causal", "empirical"], 1)[0]
        # Create genotype matrix
        if _type == "empirical":
            print(f"Empirical genotype matrix", file=sys.stderr, flush=True)
            next_id = pyslim.next_slim_mutation_id(sts)
            smts = msprime.sim_mutations(
                sts,
                rate=1e-9,
                model=msprime.SLiMMutationModel(type=0, next_id=next_id),
                keep=True,
            )
            Y = additive_genotype_matrix(smts, alive_inds, 0.05)
        else:
            Y = additive_genotype_matrix(sts, alive_inds, 0.05)
        ## Take a big sample
        train_index = np.random.choice(range(0, Y.shape[0]), SAMPLE_SIZE, replace=False)
        test_index = np.array([i for i in range(0, Y.shape[0]) if i not in train_index])
        Y = Y[train_index,]
        Y = Y[:, np.std(Y, axis=0) > 0.0]
        print(f"Genotype matrix shape: {Y.shape}", file=sys.stderr, flush=True)
        ## Assign populations to individuals in past, present and future
        sampled_df["population"] = assign_population(sampled_df, centroids)
        df_init = df[df.fraction==0.0]
        df_init["population"] = assign_population(df_init, centroids)
        df_end = df[df.fraction==1.0]
        df_end["population"] = assign_population(df_end, centroids)

        # Prepare environmental data
        X = sampled_df[["pc1", "pc2", "pc3"]].to_numpy()
        # Add two uncorrelated columns
        Xtrain = X[train_index,]
        Xtest = X[test_index,]
        Xpredtest = sampled_df[["pc1Future", "pc2Future", "pc3Future"]].to_numpy()[test_index,]
        Xpredtrain = sampled_df[["pc1Future", "pc2Future", "pc3Future"]].to_numpy()[train_index,]
        if _type == "empirical":
            # Add two uncorrelated columns
            Xtrain = np.concatenate([Xtrain, np.random.normal(0, 1, (Xtrain.shape[0], 2))], axis=1)
            Xtest = np.concatenate([Xtest, np.random.normal(0, 1, (Xtest.shape[0], 2))], axis=1)
            Xpredtrain = np.concatenate([Xpredtrain, np.random.normal(0, 1, (Xpredtrain.shape[0], 2))], axis=1)
            Xpredtest = np.concatenate([Xpredtest, np.random.normal(0, 1, (Xpredtest.shape[0], 2))], axis=1)
        
        neglog = -np.log(sampled_df.futureFitness)
        neglogtrain = neglog.iloc[train_index]
        neglogtest = neglog.iloc[test_index]
        print(f"X shape: {X.shape}", file=sys.stderr, flush=True)
        print(f"Xtrain shape: {Xtrain.shape}", file=sys.stderr, flush=True)
        print(f"Xtest shape: {Xtest.shape}", file=sys.stderr, flush=True)
        print(f"Xpredtest shape: {Xpredtest.shape}", file=sys.stderr, flush=True)
        
        # Start computing genomic offsets
        model = jl.fit(jl.GeometricGO, Y, Xtrain, tw_threshold=1e-5, λ=1e3)
        if _type == "empirical":
            pvalues = np.array( jl.LFMM_Ftest(model, Y, Xtrain))
            qvalues = false_discovery_control(pvalues)
            print(f"Number of significant SNPs with FDR: {np.sum(qvalues < 0.05)}", file=sys.stderr, flush=True)
            candidates = np.array([i+1 for i, val in enumerate(qvalues) if val < 0.05])
        else:
            candidates = np.array([i+1 for i in range(Y.shape[1])])
        trainoffset = genomic_offset(model, Xtrain, Xpredtrain, candidates)
        testoffset = genomic_offset(model, Xtest, Xpredtest, candidates)
        # Train data
        traincor, trainr2 = measure_association(neglogtrain, trainoffset)
        print_output(infile, chosen_fraction, "current", "train", _type, traincor, trainr2)
        # Test data
        testcor, testr2 = measure_association(neglogtest, testoffset)
        print_output(infile, chosen_fraction, "current", "test", _type, testcor, testr2)
        # Only test data & ecotype 1
        isA =  sampled_df.ecotype.iloc[test_index] == "A"
        testcorA, testr2A = measure_association(neglogtest[isA], testoffset[isA])
        print_output(infile, chosen_fraction, "current", "test_ecotypeA", _type, testcorA, testr2A)
        # Only test data & ecotype 2
        isB =  sampled_df.ecotype.iloc[test_index] == "B"
        testcorB, testr2B = measure_association(neglogtest[isB], testoffset[isB])
        print_output(infile, chosen_fraction, "current", "test_ecotypeB", _type, testcorB, testr2B)
        # Now with only env distance
        euclidean = np.sqrt(np.sum((Xpredtest - Xtest)**2, axis=1))
        env_cor, env_r2 = measure_association(euclidean, neglogtest)
        print_output(infile, chosen_fraction, "env_distance", "test", _type, env_cor, env_r2)
        # Compare genomic offset with future env
        sampled_df_test = sampled_df.iloc[test_index,]
        pops_with_inds = [k for k in range(K) if np.sum(sampled_df_test.population==k) > 0]
        mean_offsets = [np.mean(testoffset[sampled_df_test.population==k]) for k in pops_with_inds]
        future_neglog = [
            np.mean(
                -np.log(df_end.currentFitness[df_end.population == k])
            ) if np.sum(df_end.population == k)> 0 else 0.0 for k in pops_with_inds
        ]
        future_cor, future_r2 = measure_association(mean_offsets, future_neglog)
        print_output(infile, chosen_fraction, "future_current_env", "test", _type, future_cor, future_r2)
        # Env distance again
        mean_euclidean = [np.mean(euclidean[sampled_df_test.population==k]) for k in pops_with_inds]
        future_cor_euclidean, future_r2_euclidean = measure_association(mean_euclidean, future_neglog)
        print_output(infile, chosen_fraction, "future_env_distance", "test", _type, future_cor_euclidean, future_r2_euclidean)
        # Compute genomic offset with past env
        Xtrain_past = sampled_df[["population"]].iloc[train_index].merge(
            df_init.groupby('population')[['pc1', 'pc2', 'pc3']].median(),
            left_on = "population", right_index=True
            )[['pc1', 'pc2', 'pc3']].to_numpy()
        Xtest_past = sampled_df[["population"]].iloc[test_index].merge(
            df_init.groupby('population')[['pc1', 'pc2', 'pc3']].median(),
            left_on = "population", right_index=True
            )[['pc1', 'pc2', 'pc3']].to_numpy()
        
        if _type == "empirical":
            # Add two uncorrelated columns
            Xtrain_past = np.concatenate([Xtrain_past, np.random.normal(0, 1, (Xtrain_past.shape[0], 2))], axis=1)
            Xtest_past = np.concatenate([Xtest_past, np.random.normal(0, 1, (Xtest_past.shape[0], 2))], axis=1)
        
        modelpast = jl.fit(jl.GeometricGO, Y, Xtrain_past, tw_threshold=1e-5, λ=1e3)
        if _type == "empirical":
            pvalues = np.array( jl.LFMM_Ftest(modelpast, Y, Xtrain_past))
            qvalues = false_discovery_control(pvalues)
            print(f"Number of significant SNPs with FDR: {np.sum(qvalues < 0.05)}", file=sys.stderr, flush=True)
            candidates = np.array([i+1 for i, val in enumerate(qvalues) if val < 0.05])

        else:
            candidates = np.array([i+1 for i in range(Y.shape[1])])
        trainoffsetpast = genomic_offset(modelpast, Xtrain_past, Xpredtrain, candidates)
        testoffsetpast = genomic_offset(modelpast, Xtest_past, Xpredtest, candidates)

        # Train data
        traincorpast, trainr2past = measure_association(neglogtrain, trainoffsetpast)
        print_output(infile, chosen_fraction, "past_env", "train", _type, traincorpast, trainr2past)
        # Test data
        testcorpast, testr2past = measure_association(neglogtest, testoffsetpast)
        print_output(infile, chosen_fraction, "past_env", "test", _type, testcorpast, testr2past)
        # Only test data & ecotype 1
        testcorpastA, testr2pastA = measure_association(neglogtest[isA], testoffsetpast[isA])
        print_output(infile, chosen_fraction, "past_env", "test_ecotypeA", _type, testcorpastA, testr2pastA)
        # Only test data & ecotype 2
        testcorpastB, testr2pastB = measure_association(neglogtest[isB], testoffsetpast[isB])
        print_output(infile, chosen_fraction, "past_env", "test_ecotypeB", _type, testcorpastB, testr2pastB)
        # Env distance
        euclidean_past = np.sqrt(np.sum((Xtest_past - Xpredtest)**2, axis=1))
        env_cor_past, env_r2_past = measure_association(euclidean_past, neglogtest)
        print_output(infile, chosen_fraction, "past_env_distance", "test", _type, env_cor_past, env_r2_past)
        # Compare genomic offset with future env
        mean_offsets_past = [np.mean(testoffsetpast[sampled_df_test.population==k]) for k in pops_with_inds]
        future_cor_past, future_r2_past = measure_association(mean_offsets_past, future_neglog)
        print_output(infile, chosen_fraction, "past_env_future", "test", _type, future_cor_past, future_r2_past)
        # Env distance again
        mean_euclidean_past = [np.mean(euclidean_past[sampled_df_test.population==k]) for k in pops_with_inds]
        future_cor_euclidean_past, future_r2_euclidean_past = measure_association(mean_euclidean_past, future_neglog)
        print_output(infile, chosen_fraction, "past_env_distance_future", "test", _type, future_cor_euclidean_past, future_r2_euclidean_past)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr, flush=True)       
        continue
    





