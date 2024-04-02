### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 48675696-8dd2-40c6-a345-f5ca8b69e916
begin
  import Pkg
  # activate the shared project environment
  Pkg.activate(Base.current_project())
  Pkg.instantiate()
  using CSV, DataFrames, VCFTools, Distributions, StatsPlots, RCall, Random
  Random.seed!(555)
end

# ╔═╡ 522e50d4-cdc4-11ee-1932-df7aa47c76a3
md"""
Author: @currocam

This notebook explores the effects of misidentified SNPs when applying genomic offsets. 

## Input data
"""

# ╔═╡ 28a6f773-fb47-4844-a5b6-ceb9e977b6c7
project_path = splitdir(@__DIR__)[1]

# ╔═╡ 55103a97-805f-4fab-b73b-9159b7d3d27a
outdir = "$project_path/steps/misidentified_snps"

# ╔═╡ e90bc337-dedb-49ac-b86b-d1cc9bce79aa
md"""
First, we define a few functions to (1) parse data, (2) compute actual shift, (3) correlation between that -log shift and the geometric genomic gap. 
"""

# ╔═╡ d7f0ccad-b9e3-49d3-a29d-424a0912bfa6
function parse_slim(basename::String)
  file = open("$basename.txt", "r")
  @assert readline(file) == "#QTL"
  QTLPositions = parse.(Int, split(readline(file)[2:end]))
  @assert readline(file) == "#Populations"
  populations = parse.(Int, split(readline(file)[2:end]))
  @assert readline(file) == "#Trait"
  trait = parse.(Float64, split(readline(file)[2:end]))
  close(file)
  individuals = CSV.read("$basename.txt", DataFrame; comment="#")
  individuals = individuals[populations, 1:end]
  individuals.trait = trait
  locations_raw = read(`bash -c "grep -v '#' < $basename.vcf | cut -f 2"`, String)
  locations = parse.(Int, split(locations_raw))
  QTLPositions_idx = [findfirst(x -> qtl == x, locations) for qtl in QTLPositions]
  A = convert_gt(
    Int8, "$basename.vcf"; model=:additive,
    impute=false, center=false, scale=false
  )
  return (genotypes=A, inds=individuals, qtls=[i for i in QTLPositions_idx if !isnothing(i)])
end

# ╔═╡ 5ff3956c-8575-4959-ada0-f009abfc3269
R"library(LEA)"

# ╔═╡ 5ad983f0-8890-479d-ba1e-4b34332ed511
function genetic_gap(Y, X, Xstar, K)
  @rput Y
  @rput X
  @rput Xstar
  @rput K
  rcopy(R"gap <- genetic.gap(Y, X, pred.env = Xstar,  K = $K)$offset")
end

# ╔═╡ c26544be-3606-4643-85bc-23ca06aaf6e4
function fitness_shift(individuals)
  fitness = pdf.(Normal(0, 1), individuals.trait - individuals.optimum)
  new_fitness = pdf.(Normal(0, 1), individuals.trait .- (-10))
  -log.(fitness - new_fitness)
end

# ╔═╡ 216781db-3da9-4cd7-85fa-69fe06e4c399
function pearson_shift_gap(individuals, Y, X, Xstar, K)
  shifts = fitness_shift(individuals)
  gap = genetic_gap(Y, X, Xstar, K)
  cor(shifts, gap)
end

# ╔═╡ 6978f549-52df-42e3-be51-0e85cf7ca17b
md"""
Now, we can compute the GO under different conditions. First, let us consider only using causal ones. 

## Only causal SNPs, one causal predictor
"""

# ╔═╡ 86258821-ba90-4abd-a47a-4c529577e928
basenames = ["$project_path/steps/slim/m1.1_s$seed" for seed in 500:549]

# ╔═╡ b1536a7f-cb90-46bc-9b34-2549175bc464
function causal_snps_one_causal_predictor(basenames)
  nQTLs, pearsons = Int[], Float64[]
  sizehint!(nQTLs, length(basenames))
  sizehint!(pearsons, length(basenames))
  for basename in basenames
    Y, inds, qtls = parse_slim(basename)
    Y = Y[:, qtls] # Remove non causal SNPs from the genotypes
    X = Matrix(inds[:, 1:1]) .* 2 + randn(size(inds, 1), 1)
    Xstar = -10 .* ones(size(X)) + randn(size(X))
    push!(pearsons, pearson_shift_gap(inds, Y, X, Xstar, 1))
    push!(nQTLs, length(qtls))
  end
  DataFrame("nQTLs" => nQTLs, "PearsonCorrelation" => pearsons)
end

# ╔═╡ 75db95b1-2026-4749-824a-cc53607f45e8
let
  println("Running...")
  data = causal_snps_one_causal_predictor(basenames)
  CSV.write("$outdir/only_causal_snps_one_causal_pred.csv", data)
end

# ╔═╡ 77253601-82c8-4fd3-a3a7-e516c3906f99
md"""
## Only causal SNPs, one causal predictor, one confounded
Now, let's do the same again, but including an confounded covariate:
"""

# ╔═╡ 9169d228-4959-4aa5-816f-c9efcdc39bec
function causal_snps_two_predictors_confounded(basenames)
  nQTLs, pearsons = Int[], Float64[]
  sizehint!(nQTLs, length(basenames))
  sizehint!(pearsons, length(basenames))
  for basename in basenames
    Y, inds, qtls = parse_slim(basename)
    Y = Y[:, qtls] # Remove non causal SNPs from the genotypes
    X = Matrix(inds[:, 1:2]) .* 2 + randn(size(inds, 1), 2)
    Xstar = -10 .* ones(size(X)) + randn(size(X))
    push!(pearsons, pearson_shift_gap(inds, Y, X, Xstar, 1))
    push!(nQTLs, length(qtls))
  end
  DataFrame("nQTLs" => nQTLs, "PearsonCorrelation" => pearsons)
end

# ╔═╡ 010e8fca-9832-4a61-b0b7-0a91de1c099a
let
  println("Running...")
  data = causal_snps_two_predictors_confounded(basenames)
  CSV.write("$outdir/only_causal_snps_causal_confounded_pred.csv", data)
end

# ╔═╡ e44c8e18-fc3b-451a-9a20-b4d24e12cae9
md"""
## All SNPs, one causal predictor
"""

# ╔═╡ a4034363-2b87-444a-b699-e1d32b16dd32
function all_snps_one_causal_predictor(basenames)
  nQTLs, nloci, pearsons = Int[], Int[], Float64[]
  sizehint!(nQTLs, length(basenames))
  sizehint!(nloci, length(basenames))
  sizehint!(pearsons, length(basenames))
  for basename in basenames
    Y, inds, qtls = parse_slim(basename)
    X = Matrix(inds[:, 1:1]) .* 2 + randn(size(inds, 1), 1)
    Xstar = -10 .* ones(size(X)) + randn(size(X))
    push!(pearsons, pearson_shift_gap(inds, Y, X, Xstar, 25))
    push!(nQTLs, length(qtls))
    push!(nloci, size(Y, 2))
  end
  DataFrame("nQTLs" => nQTLs, "nloci" => nloci, "PearsonCorrelation" => pearsons)
end

# ╔═╡ 5089c54e-bf31-4aae-ae1a-b83b85217d12
let
  println("Running...")
  data = all_snps_one_causal_predictor(basenames)
  CSV.write("$outdir/all_snps_one_causal_pred.csv", data)
end

# ╔═╡ 779e5475-0869-40ff-b6b5-dfe84a1ef0e6
md"""
## All SNPs, one causal predictor, one confounded
"""

# ╔═╡ 0a780eeb-fbfc-451d-babc-a3a3ef33b847
function all_snps_two_predictors_confounded(basenames)
  nQTLs, nloci, pearsons = Int[], Int[], Float64[]
  sizehint!(nQTLs, length(basenames))
  sizehint!(nloci, length(basenames))
  sizehint!(pearsons, length(basenames))
  for basename in basenames
    Y, inds, qtls = parse_slim(basename)
    X = Matrix(inds[:, 1:2]) .* 2 + randn(size(inds, 1), 2)
    Xstar = -10 .* ones(size(X)) + randn(size(X))
    push!(pearsons, pearson_shift_gap(inds, Y, X, Xstar, 25))
    push!(nQTLs, length(qtls))
    push!(nloci, size(Y, 2))
  end
  DataFrame("nQTLs" => nQTLs, "nloci" => nloci, "PearsonCorrelation" => pearsons)
end

# ╔═╡ edf17bbf-3d31-429d-ab74-889f59030c64
let
  println("Running...")
  data = all_snps_two_predictors_confounded(basenames)
  CSV.write("$outdir/all_snps_two_predictors_confounded.csv", data)
end

# ╔═╡ 2c1c8a8c-cb52-4e42-b3a4-bfb39b594620
md"""
## All SNPs, one causal predictor, one uncorrelated
"""

# ╔═╡ c98caa23-83ee-4a35-9ceb-0d94302b77f2
function all_snps_two_predictors_uncorrelated(basenames)
  nQTLs, nloci, pearsons = Int[], Int[], Float64[]
  sizehint!(nQTLs, length(basenames))
  sizehint!(nloci, length(basenames))
  sizehint!(pearsons, length(basenames))
  for basename in basenames
    Y, inds, qtls = parse_slim(basename)
    X = [inds[:, 1] randn(size(inds[:, 1]))] .* 2 + randn(size(inds, 1), 2)
    Xstar = -10 .* ones(size(X)) + randn(size(X))
    push!(pearsons, pearson_shift_gap(inds, Y, X, Xstar, 25))
    push!(nQTLs, length(qtls))
    push!(nloci, size(Y, 2))
  end
  DataFrame("nQTLs" => nQTLs, "nloci" => nloci, "PearsonCorrelation" => pearsons)
end

# ╔═╡ bd9ab283-7148-4c51-8d92-07e6f28217fe
let
  println("Running...")
  data = all_snps_two_predictors_uncorrelated(basenames)
  CSV.write("$outdir/all_snps_two_predictors_uncorrelated.csv", data)
end

# ╔═╡ 1dd09742-5318-444a-a5ab-077104269d20
md"""
## Misidentified
Finally, we are going to assess the impact of having fewer loci than expected. 
"""

# ╔═╡ 0d0ca699-3dbd-4c29-afb1-aecab596f2ae
function sample_two_predictors_uncorrelated(basenames, k)
  nQTLs, nloci, allQTLs, pearsons = Int[], Int[], Int[], Float64[]
  sizehint!(nQTLs, length(basenames) * k)
  sizehint!(nloci, length(basenames) * k)
  sizehint!(nloci, length(basenames) * k)
  sizehint!(pearsons, length(basenames) * k)
  for basename in basenames
    Y, inds, qtls = parse_slim(basename)
    for _ in 1:k
      X = [inds[:, 1] randn(size(inds[:, 1]))] .* 2 + randn(size(inds, 1), 2)
      Xstar = -10 .* ones(size(X)) + randn(size(X))
      loci = sample(1:size(Y, 2), sample(50:300))
      push!(pearsons, pearson_shift_gap(inds, Y[:, loci], X, Xstar, 25))
      push!(nQTLs, sum([i in loci for i in qtls]))
      push!(nloci, length(loci))
      push!(allQTLs, length(qtls))
    end
  end
  DataFrame(
    "nQTLs" => nQTLs, "nloci" => nloci,
    "allQTLs" => allQTLs,
    "PearsonCorrelation" => pearsons
  )
end

# ╔═╡ 6aec87ce-e4d9-426c-88ee-970a44026317
let
  println("Running...")
  data = sample_two_predictors_uncorrelated(basenames, 100)
  CSV.write("$outdir/sample_snps_two_predictors_uncorrelated.csv", data)
end

# ╔═╡ d9eecb4d-3d5d-4dc3-ae77-79e978732686
md"""
## Sample of QTLs with unlinked neutral variation
"""

# ╔═╡ 58e8d646-5a82-4a29-b2d1-63e350acff20
function sample_two_predictors_uncorrelated_unlinked(basenames, k)
  nQTLs, nloci, allQTLs, pearsons = Int[], Int[], Int[], Float64[]
  sizehint!(nQTLs, length(basenames) * k)
  sizehint!(nloci, length(basenames) * k)
  sizehint!(nloci, length(basenames) * k)
  sizehint!(pearsons, length(basenames) * k)
  for basename in basenames
    genotypes, inds, qtls = parse_slim(basename)
    for _ in 1:k
      X = [inds[:, 1] randn(size(inds[:, 1]))] .* 2 + randn(size(inds, 1), 2)
      Xstar = -10 .* ones(size(X)) + randn(size(X))
      nUnlinked = sample(25:300)
      frecuencies = rand(Uniform(0.2, 0.8), nUnlinked)
      Y = mapreduce(p -> rand(Binomial(2, p), size(X, 1)), hcat, frecuencies)
      push!(allQTLs, length(qtls))
      selected = max(1, Int(ceil(length(qtls) / 2)))
      current_qtls = ifelse(
        rand(Bernoulli(0.85)), sample(qtls, sample(1:selected)), Int[]
      )
      if !isempty(current_qtls)
        Y = hcat(Y, genotypes[:, current_qtls])
      end
      push!(pearsons, pearson_shift_gap(inds, Y, X, Xstar, 1))
      push!(nQTLs, length(current_qtls))
      push!(nloci, size(Y, 2))
    end
  end
  DataFrame(
    "nQTLs" => nQTLs, "nloci" => nloci,
    "allQTLs" => allQTLs,
    "PearsonCorrelation" => pearsons
  )
end

# ╔═╡ 4fa85157-88da-4b97-936c-e2aa519f3b81
let
  println("Running...")
  data = sample_two_predictors_uncorrelated_unlinked(basenames, 500)
  CSV.write("$outdir/sample_snps_two_predictors_uncorrelated_unlinked.csv", data)
end

# ╔═╡ Cell order:
# ╟─522e50d4-cdc4-11ee-1932-df7aa47c76a3
# ╠═48675696-8dd2-40c6-a345-f5ca8b69e916
# ╠═28a6f773-fb47-4844-a5b6-ceb9e977b6c7
# ╠═55103a97-805f-4fab-b73b-9159b7d3d27a
# ╠═e90bc337-dedb-49ac-b86b-d1cc9bce79aa
# ╠═d7f0ccad-b9e3-49d3-a29d-424a0912bfa6
# ╠═5ff3956c-8575-4959-ada0-f009abfc3269
# ╠═5ad983f0-8890-479d-ba1e-4b34332ed511
# ╠═c26544be-3606-4643-85bc-23ca06aaf6e4
# ╠═216781db-3da9-4cd7-85fa-69fe06e4c399
# ╠═6978f549-52df-42e3-be51-0e85cf7ca17b
# ╠═86258821-ba90-4abd-a47a-4c529577e928
# ╠═b1536a7f-cb90-46bc-9b34-2549175bc464
# ╠═75db95b1-2026-4749-824a-cc53607f45e8
# ╠═77253601-82c8-4fd3-a3a7-e516c3906f99
# ╠═9169d228-4959-4aa5-816f-c9efcdc39bec
# ╠═010e8fca-9832-4a61-b0b7-0a91de1c099a
# ╠═e44c8e18-fc3b-451a-9a20-b4d24e12cae9
# ╠═a4034363-2b87-444a-b699-e1d32b16dd32
# ╠═5089c54e-bf31-4aae-ae1a-b83b85217d12
# ╠═779e5475-0869-40ff-b6b5-dfe84a1ef0e6
# ╠═0a780eeb-fbfc-451d-babc-a3a3ef33b847
# ╠═edf17bbf-3d31-429d-ab74-889f59030c64
# ╠═2c1c8a8c-cb52-4e42-b3a4-bfb39b594620
# ╠═c98caa23-83ee-4a35-9ceb-0d94302b77f2
# ╠═bd9ab283-7148-4c51-8d92-07e6f28217fe
# ╠═1dd09742-5318-444a-a5ab-077104269d20
# ╠═0d0ca699-3dbd-4c29-afb1-aecab596f2ae
# ╠═6aec87ce-e4d9-426c-88ee-970a44026317
# ╠═d9eecb4d-3d5d-4dc3-ae77-79e978732686
# ╠═58e8d646-5a82-4a29-b2d1-63e350acff20
# ╠═4fa85157-88da-4b97-936c-e2aa519f3b81
