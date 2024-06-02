using RCall, GenomicOffsets, Statistics, DataFrames, CSV, Dates, Random, MultipleTesting, LinearAlgebra
R"library(tidyverse)"

function logging(x...)
    println(Dates.now(), " ", join(x, " ")...)
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
    K = max(findfirst(pvalues .> 0.01) - 1, 1)
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

function handle_file(infile, outfile)
    FDR = 0.1
    logging("Processing $infile")
    # Use R code from other scripts for simplicity
    R"""
    # For later
    add_entries <- function(outfile, infile, causal, empirical, boots, minuslogfitness, method) {
        boots <- data.frame(boots)
        colnames(boots) <- paste0("boot", 1:ncol(boots))
        df <- boots |>
            mutate(
                causal = causal,
                empirical = empirical,
                negative_log_fitness = minuslogfitness,
                method = method
                )
        if (file.exists(outfile)) {
            df <-  bind_rows(df, read_csv(outfile))
        }
        write_csv(df, outfile)
    }

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

    methods = [RONA, GeometricGO, RDAGO, GradientForestGO]
    names = ["RONA", "GeometricGO", "RDAGO", "GradientForestGO"]
    for (method, name) in zip(methods, names)
        logging("Running $name...")
        causal = genomic_offset(GenomicOffsets.fit(method, Y[:, causal_loci], X), X, Xstar) |> vec
        empirical = empirical_offset(method, Y, X, Xstar, FDR) |> vec
        boots = bootstrap_with_candidates(
            method, Xoshiro(1000), Y, X, Xstar, 500,
            candidates_threshold=FDR, tw_threshold=1e-5
            )
        @rput causal
        @rput empirical
        @rput boots

        R"add_entries($outfile, $infile, causal, empirical, boots, minuslogfitness, $name)"
    end
end
function run(infiles, outfile)
    for infile in infiles
        handle_file(infile, outfile)
    end
end

infiles = ARGS[1:end-1]
outfile = ARGS[end]

run(infiles, outfile)