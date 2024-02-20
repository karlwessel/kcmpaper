### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ f88bfa04-5912-40e3-84e4-aa2a419284b0
using DrWatson

# ╔═╡ f17cfd29-e7ec-4a26-bb48-bed8a18e0e3e
@quickactivate

# ╔═╡ a3e411ad-f67d-4fd2-9125-c5946de05861
using AlgebraOfGraphics, CairoMakie

# ╔═╡ a27db53a-e8e5-46b8-8722-fb69dbc19760
using LaTeXStrings

# ╔═╡ 341fbb46-3cdb-4ec3-a596-90f430a546b3
using JDF

# ╔═╡ 48e1341d-de53-4fa7-975b-37949d56f91d
using DataFrames

# ╔═╡ 1ad18b54-4844-4135-b03a-4b76d3c138e8
using PlutoLinks

# ╔═╡ e41902d0-14ab-479f-8ad1-e65300352f0b
using Statistics

# ╔═╡ c2431bca-79da-4d82-b607-f64df97a4701
using PrettyTables

# ╔═╡ c9dfa4e1-7c2b-48c2-995e-aada20b804bf
pt = @ingredients("../src/plottingtools.jl")

# ╔═╡ 1a9c9d45-9f40-4b6c-b1c9-e1fe11275065
datapath = datadir("newdisorder2")

# ╔═╡ d8a7a0a2-79be-4ba8-bb9b-b7808ad55626
md"""
# Model
"""

# ╔═╡ 2b000d71-438a-4065-819e-7f54f13d8ca5
md"""
# $C(t)$, $W=0$, vary $V$, (a) constraint, (b) no constraint (2x6 panels)
"""

# ╔═╡ 76aedbe9-861e-45e5-a4cd-fdd39b15bd21
dfpaperfigures = let
	df = DataFrame(JDF.load(joinpath(datapath, "paperfigures",  "autocorrelation.jdf")))
	df.filling = df.N ./ df.L
	df
end

# ╔═╡ 20eb927a-37ee-4029-8cd7-cdcfb5963dff
let
	df = dfpaperfigures
	
	df = df[(df.filling .== 0.75) .& (df.L .== 12) .& (df.instance .< 5) .& (df.W .== 0), :]
	df = sort(df, [:representator, :ts])
	intsort = renamer([x => "V/J = $x" for x in sort(unique(df.V))])
	select!(df, Not(:numsamples))
	xlabel = pt.timelabel
	ylabel = pt.autocorrlabel
	dissort = sorter(["W = $x" for x in sort(unique(df.W))])

	datapoints = (data(df)
		* pt.statebydimersformat 
		* visual(Lines, markersize=3) 
		* mapping(:ts => xlabel, :avrgtimeavrgcorr => ylabel, group = :representator))
	label = (data(DataFrame(
			V = [1,4,16,1,4,16],
			x = ones(6)*2e7,
			y = ones(6)*0.65,
			unconstrained=[false, false, false, true, true, true], 
			label=["($x)" for x in 'a':'f'])) *
		mapping(:x => xlabel, :y => ylabel, text = :label => verbatim) *
		visual(Makie.Text))

	pt.figuredraw("paperfigurescorrect", "evolutionsmatrixpaper", (datapoints + label) 
		* mapping(; col = pt.constraintfacet, row = pt.interactionfacet);
		ignore = [:representator, :instance, :avrgcorr, :stdcorr, :timestdcorr, :vdimers, :hdimers, :holes, :periodicity, :stdtimeavrgcorr, :N],
		addtitle = false,
		axis=(;xscale=log10, limits=((1e-1, 1e9), nothing),
			width=150, height=70),
		legend=(;position=:right,))
end

# ╔═╡ 21d90b72-27fb-4287-a25f-cc96f9ad2db9
let
	df = dfpaperfigures
	combine(groupby(df, [:L, :N, :W, :J, :V, :unconstrained, :representator]), :ts => (x -> mean(diff(log10.(x)))) => :avrgstep, :ts => (x -> std(diff(log10.(x)))) => :stdstep, :ts => minimum => :t0)
end

# ╔═╡ cd5b2765-9485-480a-8d76-59a1c03b90e5
md"""
# Sketch of initial states (group)
"""

# ╔═╡ a02e7455-6807-40e2-a21b-fd4f5cfa0dfd
function plotstate(f, ψ; ms=5, nodesize=12, edgewidth=2, edgealpha=1, width=250, height=40, arrowsize=4, kwargs...)
	ψ = pt.strtoarray(ψ)
	L = length(ψ)
	ax = Axis(f, width=width, height=height)
	hidexdecorations!(ax)
	hideydecorations!(ax)
	hidespines!(ax)
	pt.plotlattice!(ax, L; inverse=false, alpha=edgealpha, linewidth=edgewidth, nodesize=nodesize)
	ilimsx, ilimsy = xlims!, ylims!
	ilimsy(ax, 0.2*1.66, 2.8*1.66)
	ilimsx(ax, 0, L+1)
	pt.plotstate!(ax, ψ; color=:tomato2, inverse=false, ms=ms, markershape=:hexagon, kwargs...)

	nc, nnc = pt.contraintedges(ψ)
	ni, nni = pt.interactionedges(ψ)
	pt.plotinteractions!(ax, nc .& ni, nnc .& nni; inverse=false, lw=arrowsize)
	ax
end

# ╔═╡ 314e697a-707d-47ef-959b-612476f0e718
function plotreps(df)
	cols = [:hdimers]
	rows = [:energy]
	f = Figure()
	layout = f[1,1] = GridLayout()
	for (i, (k, dfi)) in enumerate(pairs(groupby(df, rows)))

		#layoutcolumn = layout[2, i] = GridLayout()		
		for (j, (k, subdf)) in enumerate(pairs(groupby(dfi, cols)))
			state = first(subdf.representator)
			state = state[3:end]*state[1:2]
			plotstate(layout[i+1, k.hdimers+2], state, tiny=true,ms=17, width=130*12/20, height=25,
				color=:tomato2, nodesize=13*12/20, edgewidth=1, arrowsize=4*12/20)
		end
		#rowgap!(layoutcolumn, 5)
	end
	
	Label(layout[0, 2:end], "horizontal dimers", halign=:center)
	Label(layout[1:end, 0], "energy", rotation=0.5π)
	styles = [:dash, :solid, :dot]
	colors = Makie.wong_colors()
	for (i, (k, dfi)) in enumerate(pairs(groupby(df, rows)))
		layoutheader = layout[i+1, 1] = GridLayout()
		Label(layoutheader[1, 1], LaTeXString(k.energy == 0 ? "\$E = $(Int(k.energy))\$" : "\$E = $(Int(k.energy)) V\$"))

		# add legend line
		ax = Axis(layoutheader[1, 2], width=k.energy == 0 ? 24 : 25*12/20, height=20*12/20)
		hidedecorations!(ax)
		hidespines!(ax)
		lines!(ax, [0, 1], [1,1], color=colors[i],linewidth=5*12/20)
		colgap!(layoutheader, 5*12/20)
	end

	dimmap = Dict(0 => "none", 1 => "one", 2 => "two")

	for (i, (k, dfi)) in enumerate(pairs(groupby(df, cols)))
		layoutheader = layout[1, i+1] = GridLayout()
		Label(layoutheader[1, 1], dimmap[k.hdimers])

		#= add legend line
		ax = Axis(layoutheader[1, 2], width=110*12/20, height=20*12/20)
		hidedecorations!(ax)
		hidespines!(ax)
		lines!(ax, [0, 1], [1,1], linestyle = styles[i], color=:black,  linewidth=5*12/20)
		colgap!(layoutheader, 5*12/20)=#
	end
	
	colgap!(layout, 10*12/20)
	rowgap!(layout, 10*12/20)
	resize_to_layout!(f)
	f
end

# ╔═╡ 359ddbf4-2851-45f8-85ca-062bebcd36d4
md"""
# Autocorrelations $V$, $W$ (9 panel figure) - with constraint
"""

# ╔═╡ 21bee613-0ea5-41aa-ac43-0974507b741e
function weightedmean(x, w)
	sum(x .* w) / sum(w)
end

# ╔═╡ 5add8b59-62f4-4f2a-84df-f24970ccd0ac
statebydimersformattight = let
	s = renamer([x => "$x" for x in [1, 0, 2]])
	mapping(color = :energy => (x -> LaTeXString(x == 0 ? "\$0\$" : "\$$(Int(x)) V\$")) => L"energy $E$", linestyle = :hdimers => s => "dimers")
end

# ╔═╡ f46661eb-506f-499f-aaae-acef7fb78903
md"""
# Autocorrelations $V$, $W$ (9 panel figure) - no constraint
"""

# ╔═╡ 9a436493-3ea6-4421-bcc3-3878e7ea4072
ggplot_theme = Theme(
	colgap = 0,
	fonts = (; regular = "Computer Modern", bold = "Computer Modern"),
	rowgap = 0,
	fontsize = 12,
	figure_padding = 0,
	Legend = (labelsize=12, titlesize=9, rowgap=-9,
		padding=0, framevisible=false, orientation = :horizontal, nbanks=4, titlegap=-6, groupgap=5, patchsize=(15.0f0, 20.0f0)),
	Axis = (
		xgridvisible=false,
		ygridvisible=false,
			width=80, height=70,
	)
)

# ╔═╡ 64658026-f0f1-408d-821c-5bd13981c9c8
md"""
# Phase diagrams (a) constraint (b) no constraint - occupation distance
"""

# ╔═╡ 68b21fa5-38b8-4f10-84cd-2654e4c8b835
paper_theme = Theme(
	fonts = (; regular = "Computer Modern", bold = "Computer Modern"), fontsize = 12,
	markersize = 8,
	colgap = 5,
	rowgap = 5,
	figure_padding = 0,
)

# ╔═╡ ba537bff-640f-4a4a-92d1-7416959d1c53
with_theme(paper_theme) do
	f = pt.plotstatemakie("1101101011", tiny=true,
				color=:tomato2, dgewidth=1, markersize=55, width=250)
	save(projectdir("figures", "paperfigurescorrect", "model.pdf"), f)
	save(projectdir("figures", "paperfigurescorrect", "model.png"), f)
	f
end

# ╔═╡ e05d6bbd-c96f-486c-bb77-1369a1edb526
with_theme(paper_theme) do	
	df = dfpaperfigures
	
	df = df[(df.L .== 12), :]
	f = plotreps(sort(df, :energy))
	save(projectdir("figures", "paperfigurescorrect", "legendbyenergy.png"), f)
	save(projectdir("figures", "paperfigurescorrect", "legendbyenergy.pdf"), f)
	f
end

# ╔═╡ a331cca6-ca56-4e63-8d2a-9b7b4c1c33c4
dfwv = let
	df = DataFrame(JDF.load(joinpath(datapath, "wvmap", "avrgoccdistcombined.jdf")))	
	df
end

# ╔═╡ a73b3081-4900-499a-a5fd-f7df5ce931a2
with_theme(paper_theme) do
	df = copy(dfwv)
	df.filling = df.N ./ df.L
	
	
	df = select!(df, Not(:stdgapratio))
	df = select!(df, Not(:avrggapratio))
	df.avrgdistance = df.avrgdistance ./ (1 .- df.filling)

	ylabel = L"interaction strength $V/J$"

	heatmap = (data(df)
		* mapping(pt.normdisorderaxis, :V => ylabel, :avrgdistance => L"occup. dist. $\overline{\delta n_i} / (1-\nu)$")
		* (visual(Heatmap, colormap=cgrad(:balance; rev=false), colorrange=(0, 1))+visual(Contour)))

	label = (data(DataFrame(
			unconstrained = [false, true],
			x = ones(2)*0.2,
			y = ones(2)*0.1, 
			label=["($x)" for x in 'a':'b'])) *
		mapping(:x => pt.normdisorderlabel, :y => ylabel, text = :label => verbatim) *
		visual(Makie.Text, color = :white))
	
	pt.figuredraw("paperfigurescorrect", "occdistbyinteractiondisorder", 
		(heatmap + label) * mapping(col = pt.constraintfacet),
		ignore = [:stddistance, :N],
		addtitle = false,
		colorbar=(ticks=[0, 1], width=8),
		axis=(width=100, height=120, xscale=log10, yscale=log10),
		facet=(linkxaxes=:minimal,))
end

# ╔═╡ 067a8e00-b357-4d90-8ac0-edb22f71e118
extrema(dfwv.stddistance)

# ╔═╡ 3e24bc9a-7c22-4d23-95c9-1abf18fb7f6b
extrema(dfwv.stddistance) ./ 0.25

# ╔═╡ 79cd99f9-735b-4251-97c9-38a81c00f88d
md"""
# L-dependence (fixed $V$)

For interaction $V=4$ and filling $\sigma = 0.75$:
"""

# ╔═╡ dc629992-e7b3-434e-baf8-622de582b35e
md"""
## Occupation distance vs $W$
"""

# ╔═╡ 8a2a09d5-d537-4317-b867-01754cc26ff8
dfscaling = let
	df = DataFrame(JDF.load(joinpath(datapath, "scaling", "avrgoccdistcombined.jdf")))
	df
end

# ╔═╡ e1d2bf2c-46c8-4e4a-a739-6bc9ee97d08e
md"""
## nr of metastables vs $W$
"""

# ╔═╡ 17b7af59-2b9c-41d6-bc14-9e7a2f44476a
df1000 = DataFrame(JDF.load(joinpath(datapath, "paperfigures", "metastable.jdf")))

# ╔═╡ 74dbd88f-ae90-40cf-81f7-cda683426e42
function countrelaxed(df, ϵ)
	df = combine(groupby(df, [:L, :N, :J, :V, :W, :unconstrained]), :avrgcorr => (x -> count(x .> ϵ)) => :nummeta,
	:avrgcorr => length => :numstates)
	df.ϵ .= ϵ
	df
end

# ╔═╡ 00b6b3b4-184f-40b6-be4a-ca13701b7876
dfcounted = let
	df = df1000
	dfs = vcat([countrelaxed(df,  e) for e in [0.15]]...)
end

# ╔═╡ 415278fe-7b23-4a1b-82f7-73cbd48b8c70
nogrid_theme = Theme(
	fonts = (; regular = "Computer Modern", bold = "Computer Modern"), fontsize = 12,
	markersize = 8, 
	rowgap = 5,
	colgap = 5,
	figure_padding = 0,
	Legend = (framevisible=false, rowgap = -5),
	Axis = (
		xgridvisible=false,
		ygridvisible=false,
	)
)

# ╔═╡ 8b029e63-3905-4067-b487-446530df79ca
with_theme(nogrid_theme) do
	df2 = sort(dfscaling, :W)
	df2.normalV = df2.V ./ df2.J
	df2.filling = df2.N ./ df2.L
	df2.λ = df2.J
	df = df2
	df = select!(df, Not(:stdgapratio))
	df = select!(df, Not(:avrggapratio))
	xlabel = pt.normdisorderlabel
	ylabel = pt.occdistlabel

	dataplot = (mapping(color=pt.siteslegend, marker=pt.siteslegend) 
		* (data(df) 
			* ((visual(Lines) + visual(Scatter)) 
				* mapping(pt.normdisorderaxis, pt.occdistaxis))))

	label = (data(DataFrame(
			unconstrained = [false, true],
			x = ones(2)*0.02,
			y = ones(2)*0.0, 
			label=["(a) constrained", "(b) unconstrained"])) *
		mapping(:x => xlabel, :y => ylabel, text = :label => verbatim) *
		visual(Makie.Text))

	hline = (data((d=[0.25],)) * mapping(:d => xlabel) * visual(HLines, linestyle=:dash))
	layers = (dataplot + label) * mapping(row = pt.constraintfacethidden())
	
	f = Figure()
	p = pt.draw!(f[1, 1], hline + layers, 
		axis=(;width=350*0.67, height=155*0.67, xscale=log10,))
	legend!(f[1,1], p; halign=:right, valign=:top, margin=(10, 4, 10, 2))
	resize_to_layout!(f)
	pt.savefig("paperfigurescorrect", "occdistbydisorderpaper", f,
		layers, ignore = [:numsamples, :stddistance, :N])
	
	f
end

# ╔═╡ 32cd43d3-1202-4044-8066-cc37a69f8634
with_theme(nogrid_theme) do
	df2 = sort(dfcounted, :W)
	df2.filling = df2.N ./ df2.L
	df = df2[(df2.filling .== 0.75) .& (df2.V .== 4), :]	
	df.ratio = df.nummeta ./ df.numstates

	ylabel = LaTeXString("fraction of metastable states")
	xlabel = pt.normdisorderlabel

	dataplot = (mapping(color=pt.siteslegend, marker=pt.siteslegend) 
		* (data(df) 
			* ((visual(Lines) + visual(Scatter)) 
				* mapping(pt.normdisorderaxis, :ratio => ylabel))))

	label = (data(DataFrame(
			unconstrained = [false, true],
			x = ones(2)*0.02,
			y = ones(2)*1, 
			label=["(a) constrained", "(b) unconstrained"])) *
		mapping(:x => xlabel, :y => ylabel, text = :label => verbatim) *
		visual(Makie.Text, align=(:left, :top)))
	layers = (dataplot + label) * mapping(row = pt.constraintfacethidden())

	f = Figure()
	p = draw!(f[1,1], layers,
		axis=(;width=350*0.67, height=155*0.67, xscale=log10,))
	legend!(f[1, 1], p; halign=:left, valign=:bottom, margin=(5, 30, 20, 10))
	resize_to_layout!(f)
	pt.savefig("paperfigurescorrect", "metastablebydisorderpaper", f,
		layers, 
		ignore = [:N, :nummeta, :numstates])
	f
end

# ╔═╡ 15f7f6b9-5307-477a-97a1-1a12bd7ccdab
md"""
# Inftimevals
"""

# ╔═╡ 29f9380c-3478-4dc8-89fa-01e0326948b6
dfinf = DataFrame(JDF.load(joinpath(datapath, "paperfigures", "inftimevals12.jdf")))

# ╔═╡ 12d01a3c-b378-4374-af4d-a191d58cf5e4
dfinfavrg = let
	df = DataFrame(JDF.load(joinpath(datapath, "paperfigures", "inftimevals.jdf")))
	df.N = Int.(0.75*df.L)
	df.J .= 1
	df.V .= 4
	df
end

# ╔═╡ e53058bc-61a1-4dc4-a56a-e32280a19f84
combine(groupby(dfinf[dfinf.avrgcorr .< 0.01, :], [:L, :N, :unconstrained, :J, :V]), :W => maximum)

# ╔═╡ 4abc9099-56ec-4d16-98d1-ee79dc3e151c
combine(groupby(dfinfavrg[dfinfavrg.avrgcorr .< 0.5, :], [:L, :N, :unconstrained, :J, :V]), :W => maximum)

# ╔═╡ 88a9b85a-8e1b-4082-84d1-d7e3b47f58ba
dfthresh = let
	df = dfinfavrg
	df = df[df.L .== 20, :]
	t=0.2
	vcat([combine(groupby(df[df.avrgcorr .< t, :], [:L, :N, :unconstrained, :J, :V]), :W => maximum => :Wₜ, :W => (x -> t) => :threshold) for t in [0.05, 0.1, 0.15, 0.2, 0.5]]...)
end

# ╔═╡ 8202fdab-d1f8-4117-9506-21a1a17ed527
pretty_table(dfthresh; tf=tf_markdown)

# ╔═╡ 30b9ea90-5dad-4885-8c7c-f344ee870678
let
	df = sort(dfinf, :W)
	df = df[df.unconstrained .== false, :]
	df = df[df.L .== 12, :]

	df.totalsamples = df.numsamples .* df.periodicity
	df = combine(groupby(df, [:L, :N, :J, :V, :W, :hdimers, :unconstrained, :energy]), [:avrgcorr, :totalsamples] => weightedmean => :avrgcorr, [:hdimers, :energy] => ((x, y) -> "$(only(unique(x))), $(only(unique(y)))") => :id)
end

# ╔═╡ b727eeb7-a2b9-45c5-948b-a206242b8504
byenergylegend = :energy => (x -> x == 0 ? LaTeXString("\$E = $(Int(x))\$") : LaTeXString("\$E = $(Int(x))V\$")) => ""

# ╔═╡ b910b791-e1a8-4e7b-a52a-aac43ec8a74c
with_theme(ggplot_theme) do
	df = dfpaperfigures
	
	df = df[(df.filling .== 0.75) .& (df.L .== 12) .& in.(df.V, Ref([1, 4, 16])) .& (df.unconstrained .== false), :]

	df.totalsamples = df.numsamples .* df.periodicity
	df = combine(groupby(df, [:numsamples, :filling, :L, :N, :J, :V, :W, :ts, :unconstrained, :energy]), [:avrgtimeavrgcorr, :totalsamples] => weightedmean => :avrgtimeavrgcorr)

	df = sort(df, [:ts])
	disorders = sort(unique(df.W))
	
	xlabel = pt.timelabel
	ylabel = pt.autocorrlabel

	datapoints = (data(df)
		* mapping(color=byenergylegend)
		* (visual(Lines) 
		* mapping(:ts => xlabel, :avrgtimeavrgcorr => ylabel)
		#+ visual(Errorbars)*mapping(:ts=>xlabel, :timeavrgcorr => ylabel, :timestdcorr => ylabel)
		))
	label = (data(DataFrame(
			V = repeat([1, 4, 16], 3),
			x = ones(9)*0.2,
			y = ones(9)*-0.2,
			W=repeat(disorders', 3)[1:end], 
			label=["($x)" for x in 'a':'i'])) *
		mapping(:x => xlabel, :y => ylabel, text = :label => verbatim) *
		visual(Makie.Text))

	f = Figure()
	p = pt.draw!(f[1, 1], (datapoints + label) 
		* mapping(; col = pt.disorderfacet, row = pt.interactionfacet), 
		axis=(;xscale=log10, limits=((1e-1, 1e9), nothing), xticks=([1, 1e4, 1e8], [L"10^0", L"10^4", L"10^8"])))
	legend!(f[1,1], p; halign=:right, valign=:bottom, margin=(5, 7, 5, 10))
	resize_to_layout!(f)
	pt.savefig("paperfigurescorrect", "evolutionsmatrix3by3", f,
		datapoints* mapping(; col = pt.disorderfacet, row = pt.interactionfacet), ignore = [:instance, :representator, :avrgcorr, :stdcorr, :timestdcorr, :vdimers, :hdimers, :holes, :periodicity, :stdtimeavrgcorr, :N],)
	f
end

# ╔═╡ 09fdf42b-3129-45bb-9074-6ff969765774
with_theme(ggplot_theme) do
	df = dfpaperfigures
	
	df = df[(df.filling .== 0.75) .& (df.L .== 12) .& (df.instance .< 5) .& in.(df.V, Ref([1, 4, 16])) .& (df.unconstrained .== true), :]

	df.totalsamples = df.numsamples .* df.periodicity
	df = combine(groupby(df, [:numsamples, :filling, :L, :N, :J, :V, :W, :ts, :unconstrained, :energy]), [:avrgtimeavrgcorr, :totalsamples] => weightedmean => :avrgtimeavrgcorr)
	
	df = sort(df, [:ts])
	disorders = sort(unique(df.W))
	
	xlabel = pt.timelabel
	ylabel = pt.autocorrlabel

	datapoints = (data(df)
		* mapping(color=byenergylegend)
		* visual(Lines) 
		* mapping(:ts => xlabel, :avrgtimeavrgcorr => ylabel))
	
	label = (data(DataFrame(
			V = repeat([1, 4, 16], 3),
			x = ones(9)*0.15,
			y = ones(9)*-0.2,
			W=repeat(disorders', 3)[1:end], 
			label=["($x)" for x in 'a':'i'])) *
		mapping(:x => xlabel, :y => ylabel, text = :label => verbatim) *
		visual(Makie.Text))

	

	f = Figure()
	p = pt.draw!(f[1, 1], (datapoints + label) 
		* mapping(; col = pt.disorderfacet, row = pt.interactionfacet), 
		axis=(;xscale=log10, limits=((1e-1, 1e9), nothing), xticks=([1, 1e4, 1e8], [L"10^0", L"10^4", L"10^8"])))
	legend!(f[1,1], p; halign=:left, valign=:top, margin=(17, 5, 5, 5),)
	resize_to_layout!(f)
	pt.savefig("paperfigurescorrect", "evolutionsmatrix3by3", f,
		datapoints* mapping(; col = pt.disorderfacet, row = pt.interactionfacet), ignore = [:instance, :representator, :avrgcorr, :stdcorr, :timestdcorr, :vdimers, :hdimers, :holes, :periodicity, :stdtimeavrgcorr, :N],)
	f
end

# ╔═╡ 878ab1a7-389d-4c48-9171-3254c2bd77bc
with_theme(nogrid_theme) do
	df = sort(dfinf, :W)
	df = df[df.unconstrained .== false, :]
	df = df[df.L .== 12, :]

	df.totalsamples = df.numsamples .* df.periodicity
	df = combine(groupby(df, [:L, :N, :J, :V, :W, :unconstrained, :energy]), [:avrgcorr, :totalsamples] => weightedmean => :avrgcorr)
	
	dfavrg = sort(dfinfavrg, :W)	
	df = sort(df, :W)	
	labels = [L"(a) constrained, $L=12$", "(b) constrained", "(c) unconstrained"]
	hidefacet = renamer([a => "" for a in labels])
	df.label .= L"(a) constrained, $L=12$"
	dfavrg.label .= "(b) constrained"
	dfavrg[dfavrg.unconstrained .== true, :label] .= "(c) unconstrained"
	ylabel = L"diagonal ensemble expectation value $c_\text{diag}$"
	perstate = (data(df) *
		(mapping(pt.normdisorderaxis, :avrgcorr => ylabel) * visual(Lines))*
		mapping(color=byenergylegend) *
		mapping(group=:energy))
	avrgs = (data(dfavrg) * 
		(mapping(pt.normdisorderaxis, :avrgcorr => ylabel) *
		visual(Lines, color=:black)) *
		mapping(linestyle=:L => (x -> LaTeXString("\$L = $x\$")) => "avrg."))
	
	datapoints = (perstate + avrgs)

	label = (data(DataFrame(
			x = ones(3)*1e-8,
			y = ones(3)*-0.20,
			label=labels)) *
		mapping(:x => pt.normdisorderlabel, :y => ylabel, text = :label => verbatim) *
		visual(Makie.Text))

	f = Figure()
	p = pt.draw!(f[1, 1], (datapoints + label) 
		* mapping(; row = :label => hidefacet), 
		axis=(;xscale=log10,
			xticks=(10. .^ (-7:2:3), [LaTeXString("\$10^{$x}\$") for x in -7:2:3]),
			width=350*0.7, height=155*0.7))
	legend!(f[1,1], p; halign=:left, valign=:top, margin=(20, 5, 20, 5), rowgap=-4,
		padding=2, titlevisible=false, framevisible=false, orientation = :vertical, nbanks=1, titlegap=0, groupgap=160)
	resize_to_layout!(f)
	pt.savefig("paperfigurescorrect", "inftimevals", f,
		datapoints* mapping(; row = pt.constraintfacet), ignore = [:instance, :representator, :avrgcorr, :stdcorr, :timestdcorr, :vdimers, :hdimers, :holes, :periodicity, :stdtimeavrgcorr, :numsamples],)
	f
end

# ╔═╡ c63ed763-3088-45ff-9f1b-0ebd2d56dad2
let
	df = sort(dfinf, :W)
	df.W = round.(df.W, sigdigits=2)
	df = df[df.L .== 12, :]
	unique(df.numsamples)
end

# ╔═╡ 496a1500-2a3e-4443-9c16-b3c23e823aa4
let
	df = sort(dfinf, :W)
	df = df[df.L .== 12, :]
	df[df.numsamples .== 19, :]
end

# ╔═╡ 9c8557a8-7c0a-4140-93bf-9f29257fafe8
md"""
# Eigenentropies
"""

# ╔═╡ 15e3afdf-cd8d-4a57-827e-288c2310fd5f
dfeigen = DataFrame(JDF.load(projectdir(datapath, "entropydistribution", "entropyeigen.jdf")))

# ╔═╡ d6c1bb84-ecc9-4cdd-b9f2-d58cb8d9b58e
histogram

# ╔═╡ 5f891301-0c80-4297-9934-9827a5950058
let
	df = dfeigen
	df = df[(df.V .== 4) .& in.(df.W, Ref([2, 10, 20])) .& (df.L .== 16), :]
	combine(groupby(df, [:unconstrained, :W]), :seed => (x -> length(unique(x))) => :numsamples)
end

# ╔═╡ 87e44290-5f85-4bb6-9446-ef35645506ee
function disorderfacet(df) 
	r = renamer([w => "W = $(Int(w)) J" for w in sort(unique(df.W))])
	:W => r
end

# ╔═╡ 1a805a5d-eae0-4f2c-9767-3ffb0f4c76bc
with_theme(nogrid_theme) do
	df = dfeigen
	df = df[(df.V .== 4) .& in.(df.W, Ref([2, 10, 20])) .& (df.L .== 16), :]
	df.normentropies = round.(df.entropy ./ df.L, digits = 4)
	xlabel = L"normalized halfchain entanglement entropy $S_\text{vN}/L$"
	ylabel = L"P(S_\text{vN})"
	label = (data(DataFrame(
			unconstrained = [false, false, false, true, true, true],
			W = [2, 10, 20, 2, 10, 20],
			x = ones(6)*0.015,
			y = ones(6)*15, 
			label=["($x)" for x in 'a':'f'])) *
		mapping(:x => xlabel, :y => ylabel, text = :label => verbatim) *
		visual(Makie.Text, color = :black))

	datam = (data(df) *
		mapping(:normentropies => xlabel) *
		histogram(bins=64, normalization=:pdf))
	pt.figuredraw("paperfigurescorrect", "entropiesd", (datam + label) *
		mapping(row = disorderfacet(df), col = pt.constraintfacet),
		ignore = [:energy, :data, :entropy, :seed, :entropies],
		addtitle = false,
		axis=(width=120, height=70, limits=(nothing, (-1, 21)),
			ylabel=ylabel),
		facet=(linkyaxes = :minimal,))
end

# ╔═╡ 20d09713-4851-4346-af95-1454db946dc2
let
	df = dfeigen
	df = df[(df.V .== 4) .& in.(df.W, Ref([2, 10, 20])) .& (df.L .== 16) , :]
	df.normentropies = round.(df.entropy, digits = 4)
	xlabel = L"halfchain entanglement entropy $S_\text{vN}$"
	ylabel = L"P(S_\text{vN})"
	label = (data(DataFrame(
			unconstrained = [false, false, false, true, true, true],
			W = [2, 10, 20, 2, 10, 20],
			x = ones(6)*0.015,
			y = ones(6)*15, 
			label=["($x)" for x in 'a':'f'])) *
		mapping(:x => xlabel, :y => ylabel, text = :label => verbatim) *
		visual(Makie.Text, color = :black))

	datap = (data(df) *
		mapping(:normentropies => xlabel) *
		histogram(bins=64, normalization=:pdf))
	draw((datap) *
		mapping(row = disorderfacet(df), col = pt.constraintfacet) + data((x=[0.7],))*mapping(:x => xlabel)*visual(VLines, linestyle=:dot, alpha=0.5),
		axis=(width=170, height=100, limits=(nothing, (0, 20/16)),
			ylabel=ylabel),
		facet=(linkyaxes = :minimal,))
end

# ╔═╡ Cell order:
# ╠═f88bfa04-5912-40e3-84e4-aa2a419284b0
# ╠═f17cfd29-e7ec-4a26-bb48-bed8a18e0e3e
# ╠═a3e411ad-f67d-4fd2-9125-c5946de05861
# ╠═a27db53a-e8e5-46b8-8722-fb69dbc19760
# ╠═341fbb46-3cdb-4ec3-a596-90f430a546b3
# ╠═48e1341d-de53-4fa7-975b-37949d56f91d
# ╠═1ad18b54-4844-4135-b03a-4b76d3c138e8
# ╠═e41902d0-14ab-479f-8ad1-e65300352f0b
# ╠═c9dfa4e1-7c2b-48c2-995e-aada20b804bf
# ╠═1a9c9d45-9f40-4b6c-b1c9-e1fe11275065
# ╠═d8a7a0a2-79be-4ba8-bb9b-b7808ad55626
# ╠═ba537bff-640f-4a4a-92d1-7416959d1c53
# ╠═2b000d71-438a-4065-819e-7f54f13d8ca5
# ╠═20eb927a-37ee-4029-8cd7-cdcfb5963dff
# ╠═76aedbe9-861e-45e5-a4cd-fdd39b15bd21
# ╠═21d90b72-27fb-4287-a25f-cc96f9ad2db9
# ╠═cd5b2765-9485-480a-8d76-59a1c03b90e5
# ╠═a02e7455-6807-40e2-a21b-fd4f5cfa0dfd
# ╠═314e697a-707d-47ef-959b-612476f0e718
# ╠═e05d6bbd-c96f-486c-bb77-1369a1edb526
# ╠═359ddbf4-2851-45f8-85ca-062bebcd36d4
# ╠═b910b791-e1a8-4e7b-a52a-aac43ec8a74c
# ╠═21bee613-0ea5-41aa-ac43-0974507b741e
# ╠═5add8b59-62f4-4f2a-84df-f24970ccd0ac
# ╠═f46661eb-506f-499f-aaae-acef7fb78903
# ╠═9a436493-3ea6-4421-bcc3-3878e7ea4072
# ╠═09fdf42b-3129-45bb-9074-6ff969765774
# ╠═64658026-f0f1-408d-821c-5bd13981c9c8
# ╠═68b21fa5-38b8-4f10-84cd-2654e4c8b835
# ╠═a73b3081-4900-499a-a5fd-f7df5ce931a2
# ╠═067a8e00-b357-4d90-8ac0-edb22f71e118
# ╠═3e24bc9a-7c22-4d23-95c9-1abf18fb7f6b
# ╠═a331cca6-ca56-4e63-8d2a-9b7b4c1c33c4
# ╠═79cd99f9-735b-4251-97c9-38a81c00f88d
# ╠═dc629992-e7b3-434e-baf8-622de582b35e
# ╠═8b029e63-3905-4067-b487-446530df79ca
# ╠═8a2a09d5-d537-4317-b867-01754cc26ff8
# ╠═e1d2bf2c-46c8-4e4a-a739-6bc9ee97d08e
# ╠═17b7af59-2b9c-41d6-bc14-9e7a2f44476a
# ╠═00b6b3b4-184f-40b6-be4a-ca13701b7876
# ╠═74dbd88f-ae90-40cf-81f7-cda683426e42
# ╠═415278fe-7b23-4a1b-82f7-73cbd48b8c70
# ╠═32cd43d3-1202-4044-8066-cc37a69f8634
# ╠═15f7f6b9-5307-477a-97a1-1a12bd7ccdab
# ╠═29f9380c-3478-4dc8-89fa-01e0326948b6
# ╠═12d01a3c-b378-4374-af4d-a191d58cf5e4
# ╠═e53058bc-61a1-4dc4-a56a-e32280a19f84
# ╠═4abc9099-56ec-4d16-98d1-ee79dc3e151c
# ╠═88a9b85a-8e1b-4082-84d1-d7e3b47f58ba
# ╠═c2431bca-79da-4d82-b607-f64df97a4701
# ╠═8202fdab-d1f8-4117-9506-21a1a17ed527
# ╠═30b9ea90-5dad-4885-8c7c-f344ee870678
# ╠═b727eeb7-a2b9-45c5-948b-a206242b8504
# ╠═878ab1a7-389d-4c48-9171-3254c2bd77bc
# ╠═c63ed763-3088-45ff-9f1b-0ebd2d56dad2
# ╠═496a1500-2a3e-4443-9c16-b3c23e823aa4
# ╠═9c8557a8-7c0a-4140-93bf-9f29257fafe8
# ╠═15e3afdf-cd8d-4a57-827e-288c2310fd5f
# ╠═d6c1bb84-ecc9-4cdd-b9f2-d58cb8d9b58e
# ╠═1a805a5d-eae0-4f2c-9767-3ffb0f4c76bc
# ╠═5f891301-0c80-4297-9934-9827a5950058
# ╠═20d09713-4851-4346-af95-1454db946dc2
# ╠═87e44290-5f85-4bb6-9446-ef35645506ee
