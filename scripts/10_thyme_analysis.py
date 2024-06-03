import sys
import tskit, pyslim, msprime
import math
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
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

def rona(Y, Xtrain, X, Xpred, candidates=None):
    if candidates is not None:
        Y = Y[:, candidates]
    model = jl.fit(jl.RONA, Y, Xtrain)
    return jl.genomic_offset(model, X, Xpred).reshape(X.shape[0])

def rda(Y, Xtrain, X, Xpred, candidates=None):
    if candidates is not None:
        Y = Y[:, candidates]
    model = jl.fit(jl.RDAGO, Y, Xtrain)
    return jl.genomic_offset(model, X, Xpred, weighted=True).reshape(X.shape[0])

def geometric_default(Y, Xtrain, X, Xpred, candidates=None):
    model = jl.fit(jl.GeometricGO, Y, Xtrain, tw_threshold=1e-14)
    if candidates is not None:
        candidates = np.array( [i+1 for i, val in enumerate(candidates) if val])
        return jl.genomic_offset(model, X, Xpred, candidates).reshape(X.shape[0])
    return jl.genomic_offset(model, X, Xpred).reshape(X.shape[0])

def geometric_high(Y, Xtrain, X, Xpred, candidates=None):
    model = jl.fit(jl.GeometricGO, Y, Xtrain, tw_threshold=1e-14, λ=1e3)
    if candidates is not None:
        candidates = np.array( [i+1 for i, val in enumerate(candidates) if val])
        return jl.genomic_offset(model, X, Xpred, candidates).reshape(X.shape[0])
    return jl.genomic_offset(model, X, Xpred).reshape(X.shape[0])

def gradient_forest(Y, Xtrain, X, Xpred, candidates=None):
    if candidates is not None:
        Y = Y[:, candidates]
    model = jl.fit(jl.GradientForestGO, Y, Xtrain)
    return jl.genomic_offset(model, X, Xpred).reshape(X.shape[0])

if __name__ == "__main__":
    N_GENERATIONS = 10
    SAMPLE_SIZE = 700
    KFOLD = 5
    SETTINGS = [
        (rda, "rda"),
        (rona, "rona"),
        (geometric_default, "geometric_default"),
        (geometric_high, "geometric_high"),
        (gradient_forest, "gradient_forest"),
    ]

    # Infile from args
    infile = sys.argv[1]
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
    for generations_ago in generations:
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
        ## Assign fold
        df_temp["fold"] = np.random.choice(range(0, KFOLD), X.shape[0])

        ## Create empirical matrices
        next_id = pyslim.next_slim_mutation_id(sts)
        smts = msprime.sim_mutations(
            sts,
            rate=1e-8,
            model=msprime.SLiMMutationModel(type=0, next_id=next_id),
            keep=True,
        )
        Y2 = additive_genotype_matrix(smts)[indexes,]
        Y2 = Y2[:, np.std(Y2, axis=0) > 0.0]
        # Identify candidate SNPs
        for fn, fn_name in SETTINGS:
            print(f"Method: {fn_name}", file=sys.stderr)
            ## Compute
            causal2empirical, causal2neglog, empirical2neglog = list(), list(), list()
            for k in range(KFOLD):
                causal = fn(
                        Y[df_temp.fold != k,],
                        X[df_temp.fold != k,],
                        X[df_temp.fold == k,],
                        Xpred[df_temp.fold == k,],
                    )
                # GEA candidates
                lfmm = jl.fit(jl.GeometricGO, Y2[df_temp.fold != k,], X[df_temp.fold != k,], 1)
                pvalues = np.array( jl.LFMM_Ftest(lfmm, Y2[df_temp.fold != k,], X[df_temp.fold != k,]))
                candidates = pvalues < 0.01
                if sum(candidates) > 0:
                    empirical = fn(
                        Y2[df_temp.fold != k,],
                        X[df_temp.fold != k,],
                        X[df_temp.fold == k,],
                        Xpred[df_temp.fold == k,],
                        candidates=candidates,
                    )
                else:
                    empirical = np.zeros(Xpred[df_temp.fold == k,].shape[0])
                # Cors
                causal2empirical.append(pearsonr(causal, empirical).statistic)
                causal2neglog.append(pearsonr(causal, neglog[df_temp.fold == k]).statistic)
                empirical2neglog.append(pearsonr(empirical, neglog[df_temp.fold == k]).statistic)
            # Print
            causal2empirical_mean = np.mean(causal2empirical)
            causal2empirical_std = np.std(causal2empirical)
            causal2neglog_mean = np.mean(causal2neglog)
            causal2neglog_std = np.std(causal2neglog)
            empirical2neglog_mean = np.mean(empirical2neglog)
            empirical2neglog_std = np.std(empirical2neglog)        
            print(
                f"{infile},{fn_name},{climate_change_rate},{generation_since},{causal2empirical_mean},{causal2empirical_std},{causal2neglog_mean},{causal2neglog_std},{empirical2neglog_mean},{empirical2neglog_std}", flush=True
            )