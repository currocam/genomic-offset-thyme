import sys
import tskit, pyslim, msprime
import math
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from scipy.cluster.vq import kmeans
pd.options.mode.copy_on_write = True
# Julia stuff
from juliacall import Main as jl
jl.seval('import Pkg; Pkg.activate(".")')
jl.seval("using GenomicOffsets")
SAMPLE_SIZE = 1000

def additive_genotype_matrix(sts):
    Y = sts.genotype_matrix()
    Y2 = np.zeros((Y.shape[0], Y.shape[1] // 2))
    for k in range(Y.shape[1] // 2):
        Y2[:, k] = Y[:, 2 * k] + Y[:, 2 * k + 1]
    return Y2.transpose()

def geometric_high(Y, Xtrain, X, Xpred, candidates=None):
    if candidates is not None:
        candidates = np.array( [i+1 for i, val in enumerate(candidates) if val])
        return jl.genomic_offset(model, X, Xpred, candidates).reshape(X.shape[0])
    return 

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


infile = sys.argv[1]

ts = tskit.load(infile)
# Sanity check: all individuals coalesce
assert max(t.num_roots for t in ts.trees()) == 1
df = pd.DataFrame.from_dict(ts.metadata["SLiM"]["user_metadata"])
K, centroids = centroids(df, 50)
print(f"Number of clusters: {K}", file=sys.stderr, flush=True)

# Print header
print("file,fraction,category,dataset,statatisctic,pvalue", file=sys.stdout, flush=True)
fractions = np.array([elm for elm in set(df.fraction)])
fractions = fractions[fractions < 1.0]
for chosen_fraction in set(fractions):
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
    ## Create causal genotype matrix
    Y = additive_genotype_matrix(sts)
    ## Take a big sample
    train_index = np.random.choice(range(0, Y.shape[0]), SAMPLE_SIZE, replace=False)
    test_index = np.array([i for i in range(0, Y.shape[0]) if i not in train_index])
    Y = Y[train_index,]
    print(f"Causal genotype matrix shape: {Y.shape}", file=sys.stderr, flush=True)
    ## Assign populations to individuals in past, present and future
    sampled_df["population"] = assign_population(sampled_df, centroids)
    df_init = df[df.fraction==0.0]
    df_init["population"] = assign_population(df_init, centroids)
    df_end = df[df.fraction==1.0]
    df_end["population"] = assign_population(df_end, centroids)
    # Compute genomic offset with current env
    X = sampled_df[["pc1", "pc2", "pc3"]].to_numpy()
    Xtrain = X[train_index,]
    Xtest = X[test_index,]
    Xpredtest = sampled_df[["pc1Future", "pc2Future", "pc3Future"]].to_numpy()[test_index,]
    Xpredtrain = sampled_df[["pc1Future", "pc2Future", "pc3Future"]].to_numpy()[train_index,]
    neglog = -np.log(sampled_df.futureFitness)
    neglogtrain = neglog.iloc[train_index]
    neglogtest = neglog.iloc[test_index]
    print(f"X shape: {X.shape}", file=sys.stderr, flush=True)
    print(f"Xtrain shape: {Xtrain.shape}", file=sys.stderr, flush=True)
    print(f"Xtest shape: {Xtest.shape}", file=sys.stderr, flush=True)
    print(f"Xpredtest shape: {Xpredtest.shape}", file=sys.stderr, flush=True)

    model = jl.fit(jl.GeometricGO, Y, Xtrain, tw_threshold=1e-14, λ=1e3)
    trainoffset = np.array(jl.genomic_offset(model, Xtrain, Xpredtrain).reshape(Xtrain.shape[0]))
    testoffset = np.array(jl.genomic_offset(model, Xtest, Xpredtest).reshape(Xtest.shape[0]))
    # Train data
    traincor = pearsonr(neglogtrain, trainoffset)
    pd.DataFrame([{
        'file': infile,
        'fraction': chosen_fraction,
        'category': 'current',
        'dataset' : 'train',
        'statatisctic': traincor.statistic,
        'pvalue': traincor.pvalue
    }]).to_csv(sys.stdout, index=False, header=False)

    # Test data
    testcor = pearsonr(neglogtest, testoffset)
    pd.DataFrame([{
        'file': infile,
        'fraction': chosen_fraction,
        'category': 'current',
        'dataset' : 'test',
        'statatisctic': testcor.statistic,
        'pvalue': testcor.pvalue
    }]).to_csv(sys.stdout, index=False, header=False)

    # Only test data & ecotype 1
    isA =  sampled_df.ecotype.iloc[test_index] == "A"
    testcorA = pearsonr(neglogtest[isA], testoffset[isA])
    pd.DataFrame([{
        'file': infile,
        'fraction': chosen_fraction,
        'category': 'current',
        'dataset' : 'test_ecotypeA',
        'statatisctic': testcorA.statistic,
        'pvalue': testcorA.pvalue
    }]).to_csv(sys.stdout, index=False, header=False)

    # Only test data & ecotype 2
    isB =  sampled_df.ecotype.iloc[test_index] == "B"
    testcorB = pearsonr(neglogtest[isB], testoffset[isB])
    pd.DataFrame([{
        'file': infile,
        'fraction': chosen_fraction,
        'category': 'current',
        'dataset' : 'test_ecotypeB',
        'statatisctic': testcorB.statistic,
        'pvalue': testcorB.pvalue
    }]).to_csv(sys.stdout, index=False, header=False)
    # Compare genomic offset with future env
    sampled_df_test = sampled_df.iloc[test_index,]
    pops_with_inds = [k for k in range(K) if np.sum(sampled_df_test.population==k) > 0]
    mean_offsets = [np.mean(testoffset[sampled_df_test.population==k]) for k in pops_with_inds]
    future_neglog = [
        np.mean(
            -np.log(df_end.currentFitness[df_end.population == k])
        ) if np.sum(df_end.population == k)> 0 else 0.0 for k in pops_with_inds
    ]
    future_cor = pearsonr(mean_offsets, future_neglog)
    pd.DataFrame([{
        'file': infile,
        'fraction': chosen_fraction,
        'category': 'future_current_env',
        'dataset' : 'test',
        'statatisctic': future_cor.statistic,
        'pvalue': future_cor.pvalue
    }]).to_csv(sys.stdout, index=False, header=False)

    # Compute genomic offset with past env
    Xtrain_past = sampled_df[["population"]].iloc[train_index].merge(
        df_init.groupby('population')[['pc1', 'pc2', 'pc3']].median(),
        left_on = "population", right_index=True
        )[['pc1', 'pc2', 'pc3']].to_numpy()
    Xtest_past = sampled_df[["population"]].iloc[test_index].merge(
        df_init.groupby('population')[['pc1', 'pc2', 'pc3']].median(),
        left_on = "population", right_index=True
        )[['pc1', 'pc2', 'pc3']].to_numpy()

    modelpast = jl.fit(jl.GeometricGO, Y, Xtrain_past, tw_threshold=1e-14, λ=1e3)
    trainoffsetpast = np.array(jl.genomic_offset(modelpast, Xtrain_past, Xpredtrain).reshape(Xtrain_past.shape[0]))
    testoffsetpast = np.array(jl.genomic_offset(modelpast, Xtest_past, Xpredtest).reshape(Xtest_past.shape[0]))
    # Train data
    traincorpast = pearsonr(neglogtrain, trainoffsetpast)
    pd.DataFrame([{
        'file': infile,
        'fraction': chosen_fraction,
        'category': 'past_env',
        'dataset' : 'train',
        'statatisctic': traincorpast.statistic,
        'pvalue': traincorpast.pvalue
    }]).to_csv(sys.stdout, index=False, header=False)

    # Test data
    testcorpast = pearsonr(neglogtest, testoffsetpast)
    pd.DataFrame([{
        'file': infile,
        'fraction': chosen_fraction,
        'category': 'past_env',
        'dataset' : 'test',
        'statatisctic': testcorpast.statistic,
        'pvalue': testcorpast.pvalue
    }]).to_csv(sys.stdout, index=False, header=False)

    # Only test data & ecotype 1
    testcorpastA = pearsonr(neglogtest[isA], testoffsetpast[isA])
    pd.DataFrame([{
        'file': infile,
        'fraction': chosen_fraction,
        'category': 'past_env',
        'dataset' : 'test_ecotypeA',
        'statatisctic': testcorpastA.statistic,
        'pvalue': testcorpastA.pvalue
    }]).to_csv(sys.stdout, index=False, header=False)

    # Only test data & ecotype 2
    testcorpastB = pearsonr(neglogtest[isB], testoffsetpast[isB])
    pd.DataFrame([{
        'file': infile,
        'fraction': chosen_fraction,
        'category': 'past_env',
        'dataset' : 'test_ecotypeB',
        'statatisctic': testcorpastB.statistic,
        'pvalue': testcorpastB.pvalue
    }]).to_csv(sys.stdout, index=False, header=False)

    # Compare genomic offset with future env
    mean_offsets_past = [np.mean(testoffsetpast[sampled_df_test.population==k]) for k in pops_with_inds]
    future_cor_past = pearsonr(mean_offsets_past, future_neglog)
    pd.DataFrame([{
        'file': infile,
        'fraction': chosen_fraction,
        'category': 'past_env_future',
        'dataset' : 'test',
        'statatisctic': future_cor_past.statistic,
        'pvalue': future_cor_past.pvalue
    }]).to_csv(sys.stdout, index=False, header=False)
    sys.stdout.flush()

    





