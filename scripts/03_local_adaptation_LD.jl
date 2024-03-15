### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 76c2c384-e215-11ee-1723-231ec963f60e
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	using DelimitedFiles, RCall, Plots, StatsPlots, Statistics
end

# ╔═╡ 31562294-e823-4b5c-80f2-ad247f54fdda
md"""
When computing the impact of adding noise to our GO statistics, we got very good result even when the % of missing QTLs was very high. This notebook checks whether this was because QTLs were in high LD with other QTLs. 
"""

# ╔═╡ 25f20e2b-48cf-4e4b-92f9-a8631637fb1d
function r2_between_qtls(model)
	mat = readdlm("../steps/pairwise_ld/$model.ld")
	rds = "../steps/slim/$model.Rds"
	@rput rds
	QTLs = rcopy(
	R"""
	meta <- readr::read_rds(rds)
	list(meta[["Index QTLs 1"]], meta[["Index QTLs 2"]])
	"""
	)
	mat[QTLs[1], QTLs[2]]
end

# ╔═╡ 064dd7fc-774d-4162-b10c-91fd385c6c1a
highly = r2_between_qtls("m3.2_s100_nQTL1s100_nQTL2s_100")

# ╔═╡ d63597ab-1cc6-47ae-b46b-94fb12d0c8e1
medium = r2_between_qtls("m3.2_s100_nQTL1s50_nQTL2s_50")

# ╔═╡ d897ed36-f6ec-4924-b02a-66eb0dc80559
low = r2_between_qtls("m3.2_s100_nQTL1s20_nQTL2s_20")

# ╔═╡ e58fa5cc-d5d1-4795-b24a-e415a2f6df64
let
	histogram(vec(highly), bins = 0:0.01:1, labels = "Highly")
	histogram!(vec(medium), bins = 0:0.01:1, labels = "Medium")
	histogram!(vec(low), bins = 0:0.01:1, labels = "Low")
	xlabel!("Linkage disequilibrium between QTL1 and QTL2 (R2)")
	ylabel!("Frecuency")
end

# ╔═╡ dc42c409-f25b-4a40-8641-fc516a632571
function metric_across_simulations(f)
	metric = x -> f(filter(!isnan, vec(x)))
	seeds = 100:109
	res = zeros((length(seeds), 3))
	for (i, seed) in enumerate(seeds)
		for (j, n) in enumerate([20, 50, 100])
			res[i, j] = metric(
				r2_between_qtls("m3.2_s$(seed)_nQTL1s$(n)_nQTL2s_$(n)")
			)
		end
	end
	res
end

# ╔═╡ 46174662-1e03-4438-8d58-075a9f402996
mean_ld = metric_across_simulations(mean)

# ╔═╡ 6aa13449-7edf-4106-88e7-115f51269a4e
median_ld = metric_across_simulations(median)

# ╔═╡ 1dc6405d-167a-4906-bb40-0aba64bf463f
let
	scatter(median_ld[:, 1], median_ld[:, 2])
end

# ╔═╡ 27ea6cc1-ba97-4f72-abcb-21822004084a
[mean(filter(!isnan, vec(x))) for x in [highly, medium, low]]

# ╔═╡ e2d7a89c-68eb-4179-91ec-70d01a883aff


# ╔═╡ Cell order:
# ╠═76c2c384-e215-11ee-1723-231ec963f60e
# ╠═31562294-e823-4b5c-80f2-ad247f54fdda
# ╠═25f20e2b-48cf-4e4b-92f9-a8631637fb1d
# ╠═064dd7fc-774d-4162-b10c-91fd385c6c1a
# ╠═d63597ab-1cc6-47ae-b46b-94fb12d0c8e1
# ╠═d897ed36-f6ec-4924-b02a-66eb0dc80559
# ╠═e58fa5cc-d5d1-4795-b24a-e415a2f6df64
# ╠═dc42c409-f25b-4a40-8641-fc516a632571
# ╠═46174662-1e03-4438-8d58-075a9f402996
# ╠═6aa13449-7edf-4106-88e7-115f51269a4e
# ╠═1dc6405d-167a-4906-bb40-0aba64bf463f
# ╠═27ea6cc1-ba97-4f72-abcb-21822004084a
# ╠═e2d7a89c-68eb-4179-91ec-70d01a883aff
