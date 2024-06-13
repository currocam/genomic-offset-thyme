using RCall, GenomicOffsets, Statistics, DataFrames, CSV, Dates, Random, MultipleTesting, LinearAlgebra
import StatsBase: corspearman
R"library(tidyverse)"

function logging(x...)
    println(stderr, Dates.now(), " ", join(x, " ")...)
    flush(stderr)
end

function fix_loci(loci, removed)
    new = Int[]
    for locus in loci
        if locus in removed
            continue
        end
        push!(new, locus - count(x -> x < locus, removed))
    end
    new
end

function measure_time(func, args...)
    result = func(args...)
    elapsed_time = @elapsed func(args...)
    return result, elapsed_time
end

function find_candidates(Y, X, significance)
    mx = mean(X, dims=1)
    sdx = std(X, dims=1)
    X = (X .- mx) ./ sdx
    Y = Y .- mean(Y, dims=1)
    eigenvalues = eigvals(Y * Y' / (size(Y, 1) - 1))
    _, pvalues = TracyWidom(eigenvalues)
    K = max(findfirst(pvalues .> 1e-5) - 1, 1)
    model = RidgeLFMM(Y, X, K; center=false)
    pvalues = LFMM_Ftest(model, Y, X; genomic_control=true, center=false)
    qvalues = MultipleTesting.adjust(pvalues, BenjaminiHochberg())
    findall(qvalues .< significance)
end

function empirical_offset(::Type{GO}, Y, X, Xpred, significance) where {GO<:GenomicOffsets.AbstractGO}
    candidates = find_candidates(Y, X, significance)
    genomic_offset(GenomicOffsets.fit(GO, Y[:, candidates], X), X, Xpred)
end

function empirical_offset(::Type{GeometricGO}, Y, X, Xpred, significance)
    candidates = find_candidates(Y, X, significance)
    genomic_offset(GenomicOffsets.fit(GeometricGO, Y, X), X, Xpred, candidates)
end

function causal_offset(::Type{GO}, Y, X, Xpred, causal_loci) where {GO<:GenomicOffsets.AbstractGO}
    genomic_offset(GenomicOffsets.fit(GO, Y[:, causal_loci], X), X, Xpred)
end

function causal_offset(::Type{GeometricGO}, Y, X, Xpred, causal_loci)
    genomic_offset(GenomicOffsets.fit(GeometricGO, Y, X), X, Xpred, causal_loci)
end

function handle_file(infile)
    logging("Processing $infile")
    # Use R code from other scripts for simplicity
    R"""
    simulation <- read_rds($infile)
    shifted_fitness <- simulation[["Future fitness"]]
    causal_loci <- c(simulation[["Index QTLs 1"]], simulation[["Index QTLs 2"]])
    Y <- simulation$Genotype
    X <- matrix(c(
      simulation[["Current env 1"]], simulation[["Current env 2"]],
      rnorm(100), rnorm(100)
    ),ncol=4)
    Xstar <- matrix(c(
      simulation[["Future env 1"]], simulation[["Future env 2"]],
      rnorm(100), rnorm(100)
    ),ncol=4)
    minuslogfitness <- -log(shifted_fitness)
    """
    @rget Y
    @rget X
    @rget Xstar
    @rget minuslogfitness
    @rget causal_loci
    # Find columns with zero variance
    # Fix the causal_loci to be the correct indices after removing zero variance columns

    sds = vec(std(Y, dims=1))
    zero_var = findall(sds .== 0)
    causal_loci = fix_loci(causal_loci, zero_var)
    Y = Y[:, vec(sds .> 0)]
    logging("Loading data with $(size(Y, 2)) loci")

    df = DataFrame(
        file=String[], method=String[], FDR=Float64[],
        causal=Float64[], causalruntime=Float64[],
        empirical=Float64[], empiricalruntime=Float64[],
        minconf95=Float64[], maxconf95=Float64[], bootruntime=Float64[]
    ) 
    # Add if needed
    significances = [0.05]
    methods = [RONA, RDAGO, GeometricGO, GradientForestGO]
    names = ["RONA", "RDAGO", "GeometricGO", "GradientForestGO"]
    for (method, name) in zip(methods, names)
        logging("Running $name...")
        causal, causalruntime = measure_time(() -> causal_offset(method, Y, X, Xstar, causal_loci) |> vec)
        for significance in significances
            empiricaloffset, empiricalruntime = measure_time(() -> empirical_offset(method, Y, X, Xstar, significance) |> vec)
            boots, bootruntime = measure_time(() -> bootstrap_with_candidates(
                method, Y, X, Xstar, 100, candidates_threshold=significance, tw_threshold=1e-5)
            )
            bootscorrelation = [corspearman(boot, minuslogfitness) for boot in eachcol(boots)]
            push!(df, Dict(
                :file => infile, :method => name, :FDR => significance,
                :causal => corspearman(causal, minuslogfitness),
                :causalruntime => causalruntime,
                :empirical => corspearman(empiricaloffset, minuslogfitness),
                :empiricalruntime => empiricalruntime,
                :minconf95 => quantile(bootscorrelation, 0.025),
                :maxconf95 => quantile(bootscorrelation, 0.975),
                :bootruntime => bootruntime
            )
            )
        end
    end
    # # Gradient Forest GO
    # logging("Running GF...")
    # # We are going to discard extrapolated individuals
    # # Find individuals in Xstar that are not in range of X
    # out_of_range_individuals = [i for i in 1:size(Xstar, 1) if any(Xstar[i, :] .< vec(minimum(X, dims=1)) .|| Xstar[i, :] .> vec(maximum(X, dims=1)))]
    # inrange = setdiff(1:size(X, 1), out_of_range_individuals)
    # Y = Y[inrange, :]
    # X = X[inrange, :]
    # Xstar = Xstar[inrange, :]
    # minuslogfitness = minuslogfitness[inrange]
    # causal, causalruntime = measure_time(() -> genomic_offset(GenomicOffsets.fit(GradientForestGO, Y[:, causal_loci], X), X, Xstar) |> vec)
    # for significance in significances
    #     empiricaloffset, empiricalruntime = measure_time(() -> empirical_offset(GradientForestGO, Y, X, Xstar, significance) |> vec)
    #     boots, bootruntime = measure_time(() -> bootstrap_with_candidates(
    #         GradientForestGO, Xoshiro(1000), Y, X, Xstar, 100, candidates_threshold=significance)
    #     )
    #     bootscorrelation = [cor(boot, minuslogfitness) for boot in eachcol(boots)]
    #     push!(df, Dict(
    #         :file => infile, :method => "GradientForestGO", :FDR => significance,
    #         :causal => cor(causal, minuslogfitness),
    #         :causalruntime => causalruntime,
    #         :empirical => cor(empiricaloffset, minuslogfitness),
    #         :empiricalruntime => empiricalruntime,
    #         :minconf95 => quantile(bootscorrelation, 0.025),
    #         :maxconf95 => quantile(bootscorrelation, 0.975),
    #         :bootruntime => bootruntime
    #     )
    #     )
    # end
    df
end

function run(infiles, outfile)
    dfs = map(handle_file, infiles)
    df = vcat(dfs...)
    CSV.write(outfile, df; compress=true)
end

infiles = ARGS[1:end-1]
outfile = ARGS[end]

run(infiles, outfile)