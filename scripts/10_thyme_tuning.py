import sys
import tskit, pyslim, msprime
import pandas as pd
import numpy as np
from scipy.stats import pearsonr, false_discovery_control
# Julia stuff
from juliacall import Main as jl
jl.seval('import Pkg; Pkg.activate(".")')
jl.seval("using GenomicOffsets")

def additive_genotype_matrix(sts):
    Y = sts.genotype_matrix()
    Y2 = np.zeros((Y.shape[0], Y.shape[1] // 2))
    for k in range(Y.shape[1] // 2):
        Y2[:, k] = Y[:, 2 * k] + Y[:, 2 * k + 1]
    return Y2.transpose()

def geometric(Y, Xtrain, X, Xpred, tw, _lambda, empirical, threshold):
    model = jl.fit(jl.GeometricGO, Y, Xtrain, tw_threshold=tw, λ=_lambda)
    if empirical:
        pvalues = false_discovery_control(np.array( jl.LFMM_Ftest(model.model, Y, Xtrain)))
        candidates = np.array( [i+1 for i, pvalue in enumerate(pvalues) if pvalue < threshold])
        if len(candidates) > 0:
            print(f"Found {len(candidates)} significant SNPs", file=sys.stderr, flush=True)
            return jl.genomic_offset(model, X, Xpred, candidates).reshape(X.shape[0]), model.K
        return np.zeros(X.shape[0]), model.K
    return jl.genomic_offset(model, X, Xpred).reshape(X.shape[0]), model.K

N_GENERATIONS = 15
SAMPLE_SIZE = 1000
ENV_NOISY_COLUMNS = 2
KFOLD = 5
# Infile from args
infile = sys.argv[1]
# Outout
res = pd.DataFrame(columns=('file', 'generation_since','causal', 'empirical_mean', 'empirical_std', 'lambda', 'tw', 'causal_K', 'empirical_K', 'significance'))

# Read input
ts = tskit.load(infile)
# Sanity check: all individuals coalesce
assert max(t.num_roots for t in ts.trees()) == 1
df = pd.DataFrame.from_dict(ts.metadata["SLiM"]["user_metadata"])

# Generations
start, end = min(df.timePoint), max(df.timePoint)
elapsed = end - start
print(f"Number of sampled generations: {elapsed}", file=sys.stderr, flush=True)
climate_change_rate = 1 / elapsed
print(f"Climate change rate: {climate_change_rate}", file=sys.stderr, flush=True)
np.random.seed(42)
generations = np.random.choice(range(1, elapsed+1), N_GENERATIONS)
settings =  [(generation, tw, _lambda, significance) for generation in generations for tw in [1e-15, 1e-10, 1e-5, 1e-3] for _lambda in [1e-5, 1e3, 1e10] for significance in [0.1, 0.5, 0.90]]

for generations_ago, tw, _lambda, significance in settings:
    generation_since = end - generations_ago - start
    df_temp = df[df.timePoint == (end - generations_ago)]
    ## Simplify the tree sequence
    alive_inds = pyslim.individuals_alive_at(ts, generations_ago)
    keep_nodes = list()
    for i in alive_inds:
        keep_nodes.extend(ts.individual(i).nodes)
    sts = ts.simplify(keep_nodes, keep_input_roots=True)
    ## Create causal matrices
    Y = additive_genotype_matrix(sts)
    X = df_temp[["pc1", "pc2", "pc3"]].to_numpy()
    Xpred = df_temp[["pc1Future", "pc2Future", "pc3Future"]].to_numpy()
    ## Take a big sample
    indexes = np.random.choice(range(0, X.shape[0]), SAMPLE_SIZE)
    X = X[indexes,]
    Xpred = Xpred[indexes,]
    Y = Y[indexes,]
    df_temp = df_temp.iloc[indexes,]
    neglog = -np.array(np.log(df_temp.futureFitness))
    ## Compute causal
    causal, Kcausal = geometric(
        Y, X, X, Xpred, tw, _lambda,
        False, None
    )
    causal_cor = pearsonr(neglog, causal).statistic

    ## Create empirical matrices
    ## Add uncorrelated environmental variables
    X = np.hstack((X, np.random.normal(0, 1, (X.shape[0], ENV_NOISY_COLUMNS))))
    Xpred = np.hstack((Xpred, np.random.normal(0, 1, (Xpred.shape[0], ENV_NOISY_COLUMNS))))
    ## Assign fold
    df_temp["fold"] = np.random.choice(range(0, KFOLD), X.shape[0])

    next_id = pyslim.next_slim_mutation_id(sts)
    smts = msprime.sim_mutations(
        sts,
        rate=1e-8,
        model=msprime.SLiMMutationModel(type=0, next_id=next_id),
        keep=True,
    )
    Y = additive_genotype_matrix(smts)[indexes,]
    maf = np.array([f if f < 0.5 else 1 - f for f in np.mean(Y, axis=0) / 2 ])
    Y = Y[:, maf > 0.05]
    ## Compute offsets
    empiricals = list()
    for k in range(KFOLD):
        offset, Kempirical = geometric(
            Y[df_temp.fold != k,],
            X[df_temp.fold != k,],
            X[df_temp.fold == k,],
            Xpred[df_temp.fold == k,],
            tw, _lambda,
            True, significance
        )
        stats = pearsonr(neglog[df_temp.fold == k], offset).statistic
        # If stats is NaN, replace with 0
        if np.isnan(stats):
            stats = 0
        empiricals.append(stats)
    entry = pd.DataFrame.from_dict({
        "file": [infile],
        "generation_since": [generation_since],
        "causal": [causal_cor],
        "empirical_mean": [np.nanmean(empiricals)],
        "empirical_std": [np.nanstd(empiricals)],
        "lambda": [_lambda],
        "tw": [tw],
        "causal_K": [Kcausal],
        "empirical_K": [Kempirical],
        "significance": [significance]
        })
    res = pd.concat([res, entry], ignore_index=True)

# Print results to stdout as CSV
res.to_csv(sys.stdout, index=False)