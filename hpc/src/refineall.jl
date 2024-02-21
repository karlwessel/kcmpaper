### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 7739ece6-0fc7-11ee-29a2-f71a8ac1f35f
using DrWatson

# ╔═╡ 9ee9d0b9-1197-4dd7-ac67-98c5e2d595a8
@quickactivate

# ╔═╡ fe23e367-90d5-4414-9ef1-c45fe838efc0
using DataFrames

# ╔═╡ e8fa90a6-56dd-4ffd-9b25-f000f8eb4f36
using JDF

# ╔═╡ 32bf2eb9-2f8f-4861-a92b-da87bbe9d12e
using Statistics

# ╔═╡ e3e629f1-ab5e-4758-9728-0b7ffc902ecb
using PlutoLinks

# ╔═╡ 4a25bb5a-c549-4505-b9d4-e99706228d8e
using Memoize

# ╔═╡ 34e90af5-dff3-467d-801b-289a7d70a6d5
tools = @ingredients(srcdir("refinetools.jl"))

# ╔═╡ daf3a9d6-cd73-4175-a582-dbae1b191e2d
model = @ingredients(srcdir("model.jl"))

# ╔═╡ 45b6d56f-0caa-421d-ba71-a68088f05288
datapath = if isdir("/scratch/users")
	joinpath("/scratch/users", ENV["USER"], "raw", "newdisorder2")
else
	datadir("raw", "newdisorder2")
end

# ╔═╡ ae615072-48d6-4701-985d-8e34d3f721b6
println("""
# One particle measures
""")


# ╔═╡ 7831803e-dde8-4bb9-8314-7d3d7c497459
println("Scaling...")

# ╔═╡ fea1ec56-a5e7-4db2-b21b-e69c0c49cf1c
occdist(x) = min.(x, 1.0 .- x)

# ╔═╡ 7a2b67c8-40a4-4c61-9c93-02107d425e9f
function gapratio(x)
	diffs = diff(sort(x))
	ratios = diffs[begin+1:end] ./ diffs[begin:end-1]
	rev = ratios .> 1
	ratios[rev] .= 1 ./ ratios[rev]
	ratios[isfinite.(ratios)]
end

# ╔═╡ d0919127-ad60-4995-8200-e7a980e1c115
avrggapratio(x) = mean(gapratio(x))

# ╔═╡ 623ad489-0736-4f2e-91ab-09eb07c13ee1
avrgoccdist(a) = mean(mean.(occdist.(a)))

# ╔═╡ d6a975b1-2b25-46f0-8edf-a4365fec776c
function refineoneparticlemeasures(path)
	dffiles = tools.collectresultfiles(path, [])
	df = tools.loadresults(dffiles.path)

	df.avrgdistance = avrgoccdist.(df.densities)
	df.avrggapratio = avrggapratio.(df.energy)
	combine(groupby(df, [:L, :N, :J, :V, :W, :unconstrained]), 
		:avrgdistance => mean => :avrgdistance, 
		:avrgdistance => tools.uncertainty => :stddistance,
		:avrggapratio => mean => :avrggapratio, 
		:avrggapratio => tools.uncertainty => :stdgapratio,
		:seed => (x -> length(x)) => :numsamples)
end

# ╔═╡ 0fa91f7e-ed9f-4db0-8aeb-a8b7d5eae7a9
dfscaling = refineoneparticlemeasures(joinpath(datapath, "scaling"))

# ╔═╡ 1792af0c-0708-45a3-a4fb-4883136f0be0
JDF.save(datadir("newdisorder", "scaling", "avrgoccdistcombined.jdf"), dfscaling)

# ╔═╡ 0bb40560-03dc-4476-9849-09dfcb5a2d11
println("WV map...")

# ╔═╡ 7b0a92cb-f4ad-476b-b477-1f3021ba3f56
dfwvmap = refineoneparticlemeasures(joinpath(datapath, "wvmap"))

# ╔═╡ f5df635d-cb1b-447b-92e2-66e5e761578c
JDF.save(datadir("newdisorder", "wvmap", "avrgoccdistcombined.jdf"), dfwvmap)

# ╔═╡ 85228f12-2d65-47f1-aef9-0ab2b06b0a98
println("""
# Time evolutions
""")

# ╔═╡ cb7eba67-4f00-4ffe-9ae0-003e257d5c43
println("Paper figures...")

# ╔═╡ fe89ab10-ba0d-4071-a643-1adcbe7eb413
dfrepsrep = let
	dffiles = tools.collectresultfiles(joinpath(datapath, "paperfigures"), ["L"])
	df = tools.loadresults(dffiles.path)

	#timeaverage
	df.timeavrgcorr = tools.blocktimeavrg.(df.ts, df.autocorr)
	
	df = tools.unpack(df, [:timeavrgcorr, :autocorr], :state)
	df = tools.classifystates(df)

	# disorder average
	df = combine(groupby(df, [:L, :N, :J, :V, :W, :ts, :holes, :vdimers, :hdimers, :representator, :instance, :unconstrained, :energy]), 
		:timeavrgcorr => (x -> [mean(x)]) => :avrgtimeavrgcorr,
		:timeavrgcorr => (x -> [tools.uncertainty(x)]) => :stdtimeavrgcorr, 
		:autocorr => (x -> [mean(x)]) => :avrgcorr,
		:autocorr => (x -> [tools.uncertainty(x)]) => :stdcorr,
		:seed => (x -> length(unique(x))) => :numsamples,
		:seed => (x -> length(x) / length(unique(x))) => :periodicity)

	df = flatten(df, [:ts, :avrgtimeavrgcorr, :avrgcorr, 
		:stdtimeavrgcorr, :stdcorr])
	df.ts = identity.(df.ts) # convert Any column to float	
	df
end

# ╔═╡ b52b83e4-6102-491d-9622-5ec9898d5a19
JDF.save(datadir("newdisorder", "paperfigures", "autocorrelation.jdf"), dfrepsrep)

# ╔═╡ b9a51118-c7ed-449d-a58c-d48f5c601828
md"""
# Inf time vals
"""

# ╔═╡ 71bcf354-ab2c-4420-84d5-f87e32d9ddc1
dfinfvals = let
	df = collect_results!(joinpath(datapath, "inftimevals"); subfolders=true)
	df.W = round.(df.W, sigdigits=3)
	df
end

# ╔═╡ 468cc48d-d1a9-4f3a-be45-699285e8e8bd
dfinftimevals12 = let	
	df = dfinfvals
	df = df[df.L .== 12, :]
	df = tools.unpack(df, [:autocorr], :state)
	tools.classifystates(df)
end

# ╔═╡ e2f5a6c8-7386-4045-8755-64de54354eaa
dfinftimesavrg = let
	df = dfinfvals
	combine(groupby(df, [:L, :unconstrained, :W, :J, :V, :N]), :autocorr => (x -> mean(vcat(x...))) => :avrgcorr, :autocorr => (x -> tools.uncertainty(vcat(x...))) => :stdcorr)
end

# ╔═╡ 3d18ae11-6c7f-4497-9af4-4e8db73e967a
dfinftiemavrg12 = combine(groupby(dfinftimevals12, [:L, :N, :J, :V, :W, :holes, :vdimers, :hdimers, :representator, :instance, :unconstrained, :energy]),
		:autocorr => (x -> [mean(x)]) => :avrgcorr,
		:autocorr => (x -> [tools.uncertainty(x)]) => :stdcorr,
		:seed => (x -> length(unique(x))) => :numsamples,
		:seed => (x -> length(x) / length(unique(x))) => :periodicity)

# ╔═╡ 3117d495-39d4-4b35-b0d3-4dc8cacff1a4
JDF.save(datadir("newdisorder", "paperfigures", "inftimevals12.jdf"), dfinftiemavrg12)

# ╔═╡ 27857f11-accb-4171-8851-06967ccc32c9
JDF.save(datadir("newdisorder", "paperfigures", "inftimevals.jdf"), dfinftimesavrg)

# ╔═╡ f72a6b72-ded6-40af-ba91-df5df5d1389f
println("Meta stable...")

# ╔═╡ fafc8ca1-30cb-4e34-8889-d092bb1ca5dd
df1000 = let
	dffiles = tools.collectresultfiles(joinpath(datapath, "metastable"), [])
	df = tools.loadresults(dffiles.path)

	# Time average correlation
	df = tools.unpack(df, [:autocorr], :state)
	df.ts = mean.(df.ts)
	df.timeavrgcorr = mean.(df.autocorr)

	dfbystate = tools.classifystates(df)

	# disorder average correlation
	combine(groupby(dfbystate, [:representator, :L, :N, :J, :V, :unconstrained, :W, :ts, :energy]),
	:timeavrgcorr => mean => :avrgcorr,
	:timeavrgcorr => tools.uncertainty => :uncertcorr)
end

# ╔═╡ c1a325f1-66df-4e7d-ae7e-4772122b1dbe
JDF.save(datadir("newdisorder", "paperfigures", "metastable.jdf"), df1000)

# ╔═╡ 5e1b61e1-110d-4c8c-a385-03f2b0aaea37
println("""
# Distribution of eigenentropies
""")

# ╔═╡ 05cb1ee8-2d8d-48fc-8982-dcd694831f29
dfeigen = let 
	df = collect_results(joinpath(datapath, "eigenentropydistribution");
		subfolders = true)
	df = select!(df, Not(:path))
	
	df = unique(flatten(df, [:energy, :entropy]))
	df.energy = identity.(df.energy) # convert Any column to float	
	df.entropy = identity.(df.entropy) # convert Any column to float	
	df
end

# ╔═╡ de805c90-95ea-49f1-b7d2-4f26c0c32236
JDF.save(datadir("newdisorder", "entropydistribution", "entropyeigen.jdf"), dfeigen)

# ╔═╡ 4efd6290-5906-4bae-a297-222a17635b5b
println("""
# Evolution of entropy for homogeneous state 1110 1110...
""")

# ╔═╡ ab7284cd-6e46-45da-a75b-d76ad3c5bde4
dfdisorderavrg = let
	dffiles = tools.collectresultfiles(joinpath(datapath, "logincrease"), ["ts"])
	dffiles = dffiles[dffiles.ts .== "rough", :]
	
	df = tools.loadresults(dffiles.path)
	
	df = combine(groupby(df, [:L, :N, :V, :J, :W, :ts, :unconstrained, :state]), 
		:entropy => (x -> [mean(x)]) => :avrgentropy,
		:entropy => (x -> [tools.uncertainty(x)]) => :uncertentropy,
		:entropy => length => :numsamples)
	flatten(df, [:ts, :avrgentropy, :uncertentropy])
end

# ╔═╡ 2ef703a7-f59e-423a-aaf1-7d99eaf3aa42
JDF.save(datadir("newdisorder", "logincrease", "11101110etc.jdf"), dfdisorderavrg)

# ╔═╡ Cell order:
# ╠═7739ece6-0fc7-11ee-29a2-f71a8ac1f35f
# ╠═9ee9d0b9-1197-4dd7-ac67-98c5e2d595a8
# ╠═fe23e367-90d5-4414-9ef1-c45fe838efc0
# ╠═e8fa90a6-56dd-4ffd-9b25-f000f8eb4f36
# ╠═32bf2eb9-2f8f-4861-a92b-da87bbe9d12e
# ╠═e3e629f1-ab5e-4758-9728-0b7ffc902ecb
# ╠═4a25bb5a-c549-4505-b9d4-e99706228d8e
# ╠═34e90af5-dff3-467d-801b-289a7d70a6d5
# ╠═daf3a9d6-cd73-4175-a582-dbae1b191e2d
# ╠═45b6d56f-0caa-421d-ba71-a68088f05288
# ╠═ae615072-48d6-4701-985d-8e34d3f721b6
# ╠═d6a975b1-2b25-46f0-8edf-a4365fec776c
# ╠═7831803e-dde8-4bb9-8314-7d3d7c497459
# ╠═fea1ec56-a5e7-4db2-b21b-e69c0c49cf1c
# ╠═7a2b67c8-40a4-4c61-9c93-02107d425e9f
# ╠═d0919127-ad60-4995-8200-e7a980e1c115
# ╠═623ad489-0736-4f2e-91ab-09eb07c13ee1
# ╠═0fa91f7e-ed9f-4db0-8aeb-a8b7d5eae7a9
# ╠═1792af0c-0708-45a3-a4fb-4883136f0be0
# ╠═0bb40560-03dc-4476-9849-09dfcb5a2d11
# ╠═7b0a92cb-f4ad-476b-b477-1f3021ba3f56
# ╠═f5df635d-cb1b-447b-92e2-66e5e761578c
# ╠═85228f12-2d65-47f1-aef9-0ab2b06b0a98
# ╠═cb7eba67-4f00-4ffe-9ae0-003e257d5c43
# ╠═fe89ab10-ba0d-4071-a643-1adcbe7eb413
# ╠═b52b83e4-6102-491d-9622-5ec9898d5a19
# ╠═b9a51118-c7ed-449d-a58c-d48f5c601828
# ╠═71bcf354-ab2c-4420-84d5-f87e32d9ddc1
# ╠═468cc48d-d1a9-4f3a-be45-699285e8e8bd
# ╠═e2f5a6c8-7386-4045-8755-64de54354eaa
# ╠═3d18ae11-6c7f-4497-9af4-4e8db73e967a
# ╠═3117d495-39d4-4b35-b0d3-4dc8cacff1a4
# ╠═27857f11-accb-4171-8851-06967ccc32c9
# ╠═f72a6b72-ded6-40af-ba91-df5df5d1389f
# ╠═fafc8ca1-30cb-4e34-8889-d092bb1ca5dd
# ╠═c1a325f1-66df-4e7d-ae7e-4772122b1dbe
# ╠═5e1b61e1-110d-4c8c-a385-03f2b0aaea37
# ╠═05cb1ee8-2d8d-48fc-8982-dcd694831f29
# ╠═de805c90-95ea-49f1-b7d2-4f26c0c32236
# ╠═4efd6290-5906-4bae-a297-222a17635b5b
# ╠═ab7284cd-6e46-45da-a75b-d76ad3c5bde4
# ╠═2ef703a7-f59e-423a-aaf1-7d99eaf3aa42
