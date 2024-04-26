using RCall, GenomicOffsets, Statistics, DataFrames, MultipleTesting, CSV, Dates

R"library(tidyverse)"

function logging(x...)
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

function handle_file(infile, outfile)
    @rput infile
    @rput outfile
    logging("Processing $infile")
    # Use R code from other scripts for simplicity
    R"""
    simulation <- read_rds(infile)
    shifted_fitness <- simulation[["Future fitness"]]
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
    # Find columns with zero variance
    # Fix the causal_loci to be the correct indices after removing zero variance columns

    sds = vec(std(Y, dims=1))
    Y = Y[:, vec(sds .> 0)]
    logging("Loading data with $(size(Y, 2)) loci")
    logging("Running empirical_genomic_offset")
    empirical = empirical_genomic_offset(Y, X, Xstar; scale=true, λ=0.0001, FTest_FDR=0.1, TW_threshold=0.01)
    logging("Running bootstrap")    
    boots = bootstrap(Y, X, Xstar, 200; scale=true, λ=0.0001, FTest_FDR=0.1, TW_threshold=0.01).offsets
    @rput empirical
    @rput boots

    R"""
    boots <- data.frame(boots)
    colnames(boots) <- paste0("boot", 1:ncol(boots))
    df <- boots |>
        mutate(
            empirical = empirical$offsets,
            negative_log_fitness = minuslogfitness,
            ) |>
        pivot_longer(-negative_log_fitness, names_to = "boot_number", values_to = "offset") |>
        mutate(file = infile) |>
        separate_wider_delim(file, "_", names = c("slim", "seed", NA, NA, "QTLs")) |>
        mutate(
          slim = basename(slim),
          QTLs = str_remove_all(QTLs, ".Rds") |> as.integer(),
        )
    if (file.exists(outfile)) {
        df <-  bind_rows(df, read_csv(outfile))
    }
    write_csv(df, outfile)
    """
    
end
function run(infiles, outfile)
    for infile in infiles
        handle_file(infile, outfile)
    end
end

infiles = ARGS[1:end-1]
outfile = ARGS[end]

run(infiles, outfile)