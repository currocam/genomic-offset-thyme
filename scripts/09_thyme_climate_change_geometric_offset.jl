using GenomicOffsets, RCall, LinearAlgebra, Statistics, MultipleTesting
using Clustering, DataFrames, CSV, Random

R"library(tidyverse)"

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

function empirical_genomic_offset(Y::Matrix{T1}, X::Matrix{T2}, Xstar::Matrix{T2}; scale=true, λ=0.0001) where {T1<:Real,T2<:Real}
    if scale
        mx = mean(X, dims=1)
        sdx = std(X, dims=1)
        X = (X .- mx) ./ sdx
        Xstar = (Xstar .- mx) ./ sdx
    end
    K = latent_factors_tracy_widom(Y, 0.05, true)
    Y = Y .- mean(Y, dims=1)
    model = RidgeLFMM(Y, X, K, λ; center=false)
    pvalues = LFMM_Ftest(model, Y, X) |> vec
    qvalues = adjust(pvalues, Bonferroni())
    putative_loci = findall(qvalues .< 0.05)
    Bt = model.Bt[:, putative_loci]
    offsets = geometric_genomic_offset(Bt, X, Xstar)
    return (offsets=offsets, K=K, putative_loci=putative_loci)
end

function causal_genomic_offset(Y::Matrix{T1}, X::Matrix{T2}, Xstar::Matrix{T2}, loci; scale=true, λ=0.0001) where {T1<:Real,T2<:Real}
    if scale
        mx = mean(X, dims=1)
        sdx = std(X, dims=1)
        X = (X .- mx) ./ sdx
        Xstar = (Xstar .- mx) ./ sdx
    end
    Y = Y .- mean(Y, dims=1)
    Y = Y[:, loci]
    K = latent_factors_tracy_widom(Y, 0.05, true)
    model = RidgeLFMM(Y, X, K, λ; center=false)
    Bt = model.Bt
    offsets = geometric_genomic_offset(Bt, X, Xstar)
    return (offsets=offsets, K=K, loci=loci)
end

function find_centers_geo(data, K)
    before_geo = [data[:LocationX_T0] data[:LocationY_T0]]
    clusters = kmeans(before_geo', K)
    pop_index = findall(clusters.counts .> 4)
    centers = clusters.centers[:,pop_index]
    centers
end

function assign_centers(geo, centers)
    pops = Array{Int}(undef, size(geo, 1))
    for (i, point) in enumerate(eachrow(geo))
        dist = [norm(point - center) for center in eachcol(centers)]
        closest_center = argmin(dist)
        pops[i] = closest_center
    end
    pops

end

function handle_row(infile, sample_size, seed)
    R"data <- readr::read_rds($infile)"
    @rget data
    
    stp = Int(data[:SamplingTimepoint])    
    ftp = Int(data[:LastTimepoint])
    @rput stp
    @rput ftp
    Y = data[:Genotype]
    
    X = [data[Symbol("Env1_T$stp")] data[Symbol("Env2_T$stp")]]
    Xstar = [data[Symbol("PredEnv1_T$stp")] data[Symbol("PredEnv2_T$stp")]]    

    R"""
    set.seed($seed)
    sample_individuals <- function(data, n){
    if (n > 0.5 * nrow(data$Genotype)) {
        sample(1:nrow(data$Genotype), n)
    }
    with(data, {
        breaks <- c(0, 4, 8, 12)
        tibble(x = data[[paste0("LocationX_T", stp)]], y = data[[paste0("LocationY_T", stp)]], row = 1:nrow(data$Genotype)) |>
        # Bin the locations
        mutate(x = cut(x, breaks = breaks), y = cut(y, breaks = breaks)) |>
        group_by(x, y) |>
        sample_n(ceiling(n/9), replace = TRUE) |>
        # Shuffle the rows
        ungroup() |>
        sample_n(n) |>
        pull("row")
        })
    }
    """
    sample_index = Int.(R"sample_index <- sample_individuals(data, $sample_size)")
    Y = Y[sample_index, :]
    X = X[sample_index, :]
    Xstar = Xstar[sample_index, :]
    pred_neglogTs = -log.(data[Symbol("PredFitness_T$stp")])
    causal_loci = Int.(data[:QTLs])

    # Remove loci with zero variance
    sds = vec(std(Y, dims=1))
    zero_var = findall(sds .== 0)
    causal_loci = fix_loci(causal_loci, zero_var)
    Y = Y[:, vec(sds .> 0)]

    # Compute metrics
    empirical_go = empirical_genomic_offset(Y, X, Xstar).offsets
    empirical_cor_pred = cor(
        pred_neglogTs[sample_index],
        empirical_go
        )
    causal_go = causal_genomic_offset(Y, X, Xstar, causal_loci).offsets
    causal_cor_pred = cor(
        pred_neglogTs[sample_index],
        causal_go
        )
    
    ## Compute populations
    ## First, we have to identify populations (so we can compare the mean go of the population with their mean fitness)
    centers = find_centers_geo(data, 6*6)
    # Assign each sampled individual to a population
    now_geo = [data[Symbol("LocationX_T$stp")] data[Symbol("LocationY_T$stp")]]
    current_pops = assign_centers(now_geo, centers)
    
    # Now, we can compute the mean fitness and mean go of each population
    K = size(centers, 2)
    mean_empirical_go = zeros(K); mean_causal_go = zeros(K); mean_predfit = zeros(K)

    for i in 1:K
        mean_empirical_go[i] = mean(empirical_go[current_pops[sample_index] .== i])
        mean_causal_go[i] = mean(causal_go[current_pops[sample_index] .== i])
    end

    for i in 1:K
        if sum(current_pops .== i) == 0
            continue
        end
        mean_predfit[i] = mean(pred_neglogTs[current_pops .== i])
    end
    
    # We can know see if averaring by population gives us a better correlation
    indexes = findall(.!isnan.(mean_empirical_go))
    mean_empirical_cor_pred = cor(mean_empirical_go[indexes],mean_predfit[indexes])
    indexes = findall(.!isnan.(mean_causal_go))
    mean_causal_cor_pred = cor(mean_causal_go[indexes],mean_predfit[indexes])

    # Now, we can compare with the actual fitness after climate change
    fut_geo = [data[Symbol("LocationX_T$ftp")] data[Symbol("LocationY_T$ftp")]]
    future_pops = assign_centers(fut_geo, centers)

    mean_fut_predfit = zeros(K)
    futneglog = -log.(data[Symbol("Fitness_T$ftp")])
    for i in 1:K
        if sum(current_pops .== i) == 0
            continue
        end
        mean_fut_predfit[i] = mean(futneglog[future_pops .== i])
    end

    indexes = findall(.!isnan.(mean_empirical_go))
    mean_empirical_cor_actual = cor(mean_empirical_go[indexes],mean_fut_predfit[indexes])
    indexes = findall(.!isnan.(mean_causal_go))
    mean_causal_cor_actual = cor(mean_causal_go[indexes],mean_fut_predfit[indexes])

    # Create a dataframe with the results
    DataFrame(
        empirical_cor_predTs=empirical_cor_pred,
        causal_cor_predTs=causal_cor_pred,
        mean_empirical_cor_pred=mean_empirical_cor_pred,
        mean_causal_cor_pred=mean_causal_cor_pred,
        mean_empirical_cor_actual=mean_empirical_cor_actual,
        mean_causal_cor_actual=mean_causal_cor_actual,
        climate_change_rate = data[:ClimateChangeRate],
        sample_size = sample_size,
        sample_seed = seed,
        file = infile
    )
end

function run(infiles, outfile)
    df = DataFrame()
    sample_sizes = [100, 200, 500, 1000]
    seeds = [1000, 2000, 3000]

    combs = Base.product(infiles, sample_sizes, seeds)

    for (infile, n, seed) in combs
        println("Processing $infile with $n samples and seed $seed")
        try
            row = handle_row(infile, n, seed)
            df = vcat(df, row)
        catch e
            println("Error processing $infile with $n samples and seed $seed: $e")
            continue
        end
    end
    CSV.write(outfile, df; compress=true)
end

infiles = ARGS[1:end-1]
outfile = ARGS[end]

run(infiles, outfile)