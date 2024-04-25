using RCall, GenomicOffsets, Statistics, DataFrames, MultipleTesting, CSV, Dates, Distributed
R"library(tidyverse)"

function log(x...)
    println(Dates.now(), " ", join(x, " ")...)
    flush(stderr)
end

function empirical_genomic_offset(Y::Matrix{T1}, X::Matrix{T2}, Xstar::Matrix{T2}; scale=true, λ=0.0001, FTest_FDR=0.1, TW_threshold=0.01) where {T1<:Real,T2<:Real}
    if scale
        mx = mean(X, dims=1)
        sdx = std(X, dims=1)
        X = (X .- mx) ./ sdx
        Xstar = (Xstar .- mx) ./ sdx
    end
    Y = Y .- mean(Y, dims=1)
    # Undefined offset matrix
    K = latent_factors_tracy_widom(Y, TW_threshold, false)
    model = RidgeLFMM(Y, X, K, λ; center=false)
    pvalues = LFMM_Ftest(model, Y, X) |> vec
    qvalues = MultipleTesting.adjust(pvalues, BenjaminiHochberg())
    putative_loci = findall(qvalues .< FTest_FDR)
    Bt = model.Bt[:, putative_loci]
    offsets = geometric_genomic_offset(Bt, X, Xstar)
    return (offsets=offsets, K=K, putative_loci=putative_loci)
end

function causal_genomic_offset(Y::Matrix{T1}, X::Matrix{T2}, Xstar::Matrix{T2}, loci::Array{Int}; scale=true, λ=0.0001, TW_threshold=0.01) where {T1<:Real,T2<:Real}
    if scale
        mx = mean(X, dims=1)
        sdx = std(X, dims=1)
        X = (X .- mx) ./ sdx
        Xstar = (Xstar .- mx) ./ sdx
    end
    Y = Y .- mean(Y, dims=1)
    Y = Y[:, loci]
    K = latent_factors_tracy_widom(Y, TW_threshold, false)
    model = RidgeLFMM(Y, X, K, λ; center=false)
    Bt = model.Bt
    offsets = geometric_genomic_offset(Bt, X, Xstar)
    return (offsets=offsets, K=K, loci=loci)
end

function bootstrap(Y::Matrix{T1}, X::Matrix{T2}, Xstar::Matrix{T2}, nboot::Int; scale=true, λ=0.0001, FTest_FDR=0.1, TW_threshold=0.01) where {T1<:Real,T2<:Real}
    if scale
        mx = mean(X, dims=1)
        sdx = std(X, dims=1)
        X = (X .- mx) ./ sdx
        Xstar = (Xstar .- mx) ./ sdx
    end
    Y = Y .- mean(Y, dims=1)
    # Undefined offset matrix
    offsets = Array{Float64}(undef, size(X, 1), nboot)
    latent_factors = Array{Int}(undef, nboot)
    nputativeloci = Array{Int}(undef, nboot) 

    for index in 1:nboot
        idx = rand(1:size(Y, 2), size(Y, 2))
        Ytemp = Y[:, idx]
        K = latent_factors_tracy_widom(Ytemp, TW_threshold, false)
        latent_factors[index] = K
        model = RidgeLFMM(Ytemp, X, K, λ; center=false)
        pvalues = LFMM_Ftest(model, Ytemp, X) |> vec
        qvalues =  MultipleTesting.adjust(pvalues, BenjaminiHochberg())
        putative_loci = findall(qvalues .< FTest_FDR)
        nputativeloci[index] = length(putative_loci)
        Bt = model.Bt[:, putative_loci]
        offsets[:, index] = geometric_genomic_offset(Bt, X, Xstar)
    end
    return (offsets=offsets, latent_factors=latent_factors, nputativeloci=nputativeloci)

end


function grid_conditions(infile)
    R"""
    #thresholds <- c(0.01, 0.05, 0.1)
    thresholds <- c(0.01)
    lambdas <- c(0.00001)
    conditions <- expand_grid(FTest_FDR = thresholds, TW_threshold = thresholds, lambda = lambdas) |>
        mutate(file = $infile) |>
        separate_wider_delim(file, "_", names = c("slim", "seed", NA, NA, "QTLs")) |>
        mutate(
            slim = basename(slim),
            QTLs = str_remove_all(QTLs, ".Rds") |> as.integer(),
            #m3.1 non local, m3.2 local, other NA
            scenario = case_when(
                grepl("m3.1", slim) ~ "Non local adaptation",
                grepl("m3.2", slim) ~ "Local adaptation",
                TRUE ~ NA_character_
            )
        )
    """
    @rget conditions
    conditions[!, :empirical_K] .= NaN
    conditions[!, :causal_K] .= NaN
    conditions[!, :minconf_K] .= NaN
    conditions[!, :maxconf_K] .= NaN
    conditions[!, :empirical_nloci] .= NaN
    conditions[!, :causal_nloci] .= NaN
    conditions[!, :minconf_nloci] .= NaN
    conditions[!, :maxconf_nloci] .= NaN
    conditions[!, :empirical_cor_fitness] .= NaN
    conditions[!, :causal_cor_fitness] .= NaN
    conditions[!, :minconf_cor_fitness] .= NaN
    conditions[!, :maxconf_cor_fitness] .= NaN
    conditions
end

function handle_row!(row, Y, X, Xstar, minuslogfitness, causal_loci)
    # Compute offsets
    causal_offset = causal_genomic_offset(Y, X, Xstar, causal_loci; scale=true, λ=row.lambda, TW_threshold=row.TW_threshold)
    empirical_offset = empirical_genomic_offset(Y, X, Xstar; scale=true, λ=row.lambda, FTest_FDR=row.FTest_FDR, TW_threshold=row.TW_threshold)
    boots = bootstrap(Y, X, Xstar, 100; λ=row.lambda, scale=true, FTest_FDR=row.FTest_FDR, TW_threshold=row.TW_threshold)
    # Latent factors
    row.empirical_K = empirical_offset.K
    row.causal_K = causal_offset.K
    row.minconf_K = quantile(boots.latent_factors, 0.025)
    row.maxconf_K = quantile(boots.latent_factors, 0.975)
    # Loci
    row.empirical_nloci = length(empirical_offset.putative_loci)
    row.causal_nloci = length(causal_offset.loci)
    row.minconf_nloci = quantile(boots.nputativeloci, 0.025)
    row.maxconf_nloci = quantile(boots.nputativeloci, 0.975)
    # Correlations
    row.empirical_cor_fitness = cor(empirical_offset.offsets, minuslogfitness)
    row.causal_cor_fitness = cor(causal_offset.offsets, minuslogfitness)
    boot_cors = [cor(col, minuslogfitness) for col in eachcol(boots.offsets)]
    row.minconf_cor_fitness = quantile(boot_cors, 0.025)
    row.maxconf_cor_fitness = quantile(boot_cors, 0.975)
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

function handle_file(infile)
    log("Processing $infile")
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
    log("Loading data with $(size(Y, 2)) loci")
    df = grid_conditions(infile)
    log("Testing $(size(df, 1)) conditions")

    for row in eachrow(df)
        log("Processing condition $(row.FTest_FDR) $(row.TW_threshold) $(row.lambda)")
        handle_row!(row, Y, X, Xstar, minuslogfitness, causal_loci)
    end
    df
end

function run(infiles, outfile)
    dfs = map(handle_file, infiles);
    df = vcat(dfs...)
    CSV.write(outfile, df; compress=true)
end

infiles = ARGS[1:end-1]
outfile = ARGS[end]

run(infiles, outfile)