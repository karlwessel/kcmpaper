### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 0473271b-b8ac-4260-ba8b-dcc044a4584d
using DrWatson

# ╔═╡ 47c268dc-b433-407a-b22b-89f1e32e9321
@quickactivate

# ╔═╡ 735764ef-9b8e-4c4e-96f9-1ceb63d15a0a
using AlgebraOfGraphics, CairoMakie

# ╔═╡ a98ced4c-1d7d-4909-bd30-5d143c12abdb
using CSV

# ╔═╡ 6d464ced-bf05-4057-8a8c-8f9720768ca5
using LaTeXStrings

# ╔═╡ ca87392d-9633-4828-9a91-f904da80c4f8
using DataFrames

# ╔═╡ 2a5f3965-ac1f-4e41-9236-7a1bacce2944
set_theme!(fonts = (; regular = "Computer Modern", bold = "Computer Modern"))

# ╔═╡ b48c9874-f8a9-41f5-9eca-565ffc555398
md"""
# Legend entries and facets
"""

# ╔═╡ 9abb9cc8-015b-4036-8678-c15293754da7
interactionfacet = let
	intsort = renamer([x => x == 1 ? "V = J" : "V = $x J" for x in [1, 4, 16]])
	:V => intsort
end

# ╔═╡ 0a7f393a-7739-4b55-b7d6-3e2ccb49b8ddlet
constraintfacet = :unconstrained => (x -> x ? LaTeXString("unconstrained") : LaTeXString("constrained"))

function constraintfacethidden(labels = [false, true])
        hide = renamer([a => "" for a in labels])
        :unconstrained =>hide
end

# ╔═╡ 7adc5c8d-d723-4932-8b97-cc2a86b49bef
autocorrlabel = L"time averaged autocorrelation $\bar{c}$"

# ╔═╡ 9263b1bc-cadb-4d43-9ccc-47e426145c01
timelabel = L"time $Jt$"

# ╔═╡ 0b6bd6bd-eebe-4b02-b291-97f31f5b5d8b
statebydimersformat = let	
	mapping(color = :energy => (x -> LaTeXString("\$e/V = $x\$")) => "energy")
end

# ╔═╡ ce3e4422-19ba-4c52-b14f-b9fe9d307f93
disorderfacet = let
        wsort = renamer([x => x == 0 ? "W = 0" : x == 1 ? "W = J" : "W = $(Int(x)) J" for x in [0, 1, 2, 10, 20]])
	:W => wsort
end

# ╔═╡ f131f150-ef5d-4bd9-9511-01391648b982
siteslegend = :L=> (x -> LaTeXString("\$L = $x\$")) => ""

# ╔═╡ e62552c3-7ec6-43be-b77d-f595e33b6fdb
occdistlabel = L"average occupation distance $\overline{\delta n_i}$"

# ╔═╡ c32b1f1c-3de7-4cb1-ba7b-92062f51e536
occdistaxis = :avrgdistance => occdistlabel

# ╔═╡ 843efd9a-7724-4830-95a7-d888812bcdf2
normdisorderlabel = L"logarithmic disorder strength $W/J$"

# ╔═╡ a2e60818-4f15-4841-8d17-8f4d3f447a3d
normdisorderaxis = :W => normdisorderlabel

# ╔═╡ c14fc29f-e686-4a6c-8454-990fcacd0cac
md"""
# Plot and save figures togehter with their dataset
"""

# ╔═╡ 0670c838-4b59-4473-8d0d-00d5bd00de0c
getdata(layer) = [layer.data]

# ╔═╡ fd67b127-8d6b-4868-a42f-eac31611c203
getdata(layers::Layers) = vcat(getdata.(layers)...)

# ╔═╡ e20d3971-3c6a-4765-b894-1441df13331b
usedsymbol(x) = x isa Pair ? first(x) : x

# ╔═╡ 62e53ecf-744a-4664-bd1e-8f680d3915d0
function unused(layer)
	cols = Symbol.(names(layer.data))
	namedcols = [usedsymbol(x) for x in values(layer.named)]
	setdiff(cols, usedsymbol.(layer.positional), namedcols)
end

# ╔═╡ a918682d-ab61-4ccc-bdde-3cfa56f8ebd7
function unused(layers::Layers)
	union([unused(l) for l in layers]...)
end

# ╔═╡ 66ed823b-9d5e-47c9-94f1-769dc5211e63
function unusedvals(layer; ignore = [])
	undrawn = setdiff(unused(layer), ignore)
	data = first(unique(getdata(layer)))
	undrawnvals = [x => unique(getproperty(data, Symbol(x))) for x in undrawn]
	undrawnvals = Dict([p => length(v) == 1 ? first(v) : "many" for (p, v) in undrawnvals])
	undrawnvals, data
end

# ╔═╡ 8fd77707-01aa-4e32-8d81-e94cbc18c52b
function savefig(path, name, p, layer; ignore = [])
	undrawnvals, data = unusedvals(layer; ignore = ignore)
	mkpath(projectdir("figures", path))
	save(projectdir("figures", path, savename(name, undrawnvals, "pdf")), p)
	save(projectdir("figures", path, savename(name, undrawnvals, "png")), p)
	mkpath(datadir("figures", path))
	CSV.write(datadir("figures", path, savename(name, undrawnvals, "csv")), data)
end

# ╔═╡ fe86b1e8-f5eb-4081-97b7-f5a61b9664a4
function figuredraw(path, name, layer; ignore = [], addtitle=true, kwargs...)
	undrawnvals, _ = unusedvals(layer; ignore = ignore)
	p = draw(layer; kwargs...)
	if addtitle
		Label(p.figure[0, :], savename(undrawnvals; connector = ", "), fontsize=20)
		resize_to_layout!(p.figure)
	end
	savefig(path, name, p, layer; ignore = ignore)
	p
end

# ╔═╡ 7585ce3f-53fa-431f-be83-180bf713ac9a
# ╠═╡ disabled = true
#=╠═╡
figuredraw("test", "test", 
	data(DataFrame(x = rand(4), y=ones(4)))*mapping(:x)*histogram())
  ╠═╡ =#

# ╔═╡ ff79e886-5192-4c73-bb47-c883276b923b
md"""
## Plot states
"""

# ╔═╡ df80a991-9bac-4e82-b538-e42219cb0093
function plotstatemakie(ψ; tiny=false, inverse=false, width=250, resolution=nothing, kwargs...)
	xdecoheight = 45
	L = length(ψ)

	aspect = 1.66
	height = width*aspect*3/(L+1) + xdecoheight
	axisheight = height - xdecoheight
	axiswidth = width
	ms = axisheight / 2 - 6
	rwidth, rheight = width, height
	if inverse
		aspect = 1/aspect
		rwidth, rheight = rheight, rwidth
	end
	if isnothing(resolution)
		resolution = (rwidth, rheight)
	end
	
	f = Figure(;resolution=resolution)
	ax = plotstatemakie(f[1,1], ψ; tiny=tiny, inverse=inverse, width=width, kwargs...)
	if tiny
		ax.scene
	else
		f
	end
end

# ╔═╡ 46dbdeca-877b-4ad2-b648-219926b45361
shorten(x, f) = [x[2]*(1-f)+f*x[1], x[1]*(1-f)+f*x[2]]

# ╔═╡ fbeee4e5-d4ae-473c-b0da-db01db7b3c0a
applyinteraction(i, j, b) = b[i] != b[mod1(j, length(b))]

# ╔═╡ ace3f01f-a663-4a0e-8b50-4f5cee2c560b
function interactionedges(ψ)
	st = ψ .> 0
	nearest = [applyinteraction(i, i+1, st) for i in 1:(length(st)-1)]
	nnearest = [applyinteraction(i, i+2, st) for i in 1:(length(st)-2)]
	return nearest, nnearest
end

# ╔═╡ b1908130-1f3e-43a1-8252-fa13d36bb7da
function applyconstraint(i, j, b)
	N = length(b)
	if j == i + 1
		return !b[mod1(i-1, N)] || !b[mod1(j+1, N)]
	elseif j == i + 2
		return !b[mod1(i+1, N)]
	else
		return 0
	end
end

# ╔═╡ d6768f2e-b6d8-42a2-bc1c-38c5c6a30130
function contraintedges(ψ)
	st = ψ .> 0
	nearest = [false;
		[applyconstraint(i, i+1, st) for i in 2:(length(st) - 2)];
		false]
	nnearest = [applyconstraint(i, i+2, st) for i in 1:(length(st) - 2)]
	return nearest, nnearest
end

# ╔═╡ 498b8e17-1f21-4830-a497-65e07f21a649
function iscatter!(p, x, y; inverse=false, rotations=0, kwargs...)
	if inverse
		x, y = y, x
		rotations += π/2
	end
	scatter!(p, x, y; rotations=rotations, kwargs...)
end

# ╔═╡ 26b1ea31-29a2-4369-956c-53966900b629
function iarrows!(p, x, y, u, v; inverse=false, kwargs...)
	if inverse
		x, y = y, x
		u, v = v, u
	end
	arrows!(p, x, y, u, v; kwargs...)
end

# ╔═╡ 07d29080-cb78-45ba-9ed3-84b425f4178a
function arrow!(ax, x, y; lw=5, kwargs...)
	u = diff(x)
	v = diff(y)
	x = x[1:end-1]
	y = y[1:end-1]
	iarrows!(ax, x, y, u, v; arrowsize=lw*3, linewidth=lw, kwargs...)
end

# ╔═╡ 54602680-9bdf-4cf1-8777-c7a07d716fe6
function plotinteractions!(ax, nearest, nextnearest; kwargs...)
	for i in findall(isone.(nearest))
		x = shorten([i, (i+1)], 0.8)
		y = shorten(mod1.([i, (i+1)], 2), 0.8) * 1.66
		arrow!(ax, x, y; color=:steelblue, kwargs...)
		arrow!(ax, reverse(x), reverse(y); color=:steelblue, kwargs...)
	end
	for i in findall(isone.(nextnearest))
		x = shorten([i, i+2], 0.8)
		y = mod1.([i, i+2], 2) * 1.66
		arrow!(ax, x, y; color=:steelblue, kwargs...)
		arrow!(ax, reverse(x), reverse(y); color=:steelblue, kwargs...)
	end
end

# ╔═╡ 3eecf733-d4c9-43a4-8850-8088c1e3b817
function ilines!(p, x, y; inverse=false, kwargs...)
	if inverse
		x, y = y, x
	end
	lines!(p, x, y; kwargs...)
end

# ╔═╡ 3b922a5c-a896-4198-aa3e-921fcbd57e86
function plotnodes!(ax, nodes; nodesize=8, kwargs...)
	x = findall(isone.(nodes))
	y = mod1.(x, 2)
	iscatter!(ax, x, y*1.66; markersize=nodesize-4, rotations=π/2, strokewidth=1, kwargs...)
end

# ╔═╡ a03bf078-64c5-425e-9efa-36a3021b7d7f
function plotstate!(ax, ψ; inverse=false, ms=30, xlims=(0, length(ψ)+1), color=:tomato2, kwargs...) 
	L = length(ψ)
	plotnodes!(ax, ψ .> 0; inverse=inverse, marker='⬣', color=collect(zip(fill(color, L), ψ[ψ .> 0])), nodesize=ms, xlims=xlims, label="", strokecolor=collect(zip(fill(:black, L), ψ[ψ .> 0])), kwargs...)
	iscatter!(ax, 1:2:L+1, ones(Int(L/2)+1)*3*1.66;inverse=inverse, marker='⬣', rotations=π/2, markersize=ms-4, color=:black)
	iscatter!(ax, 0:2:L, zeros(Int(L/2)+1);inverse=inverse, marker='⬣', rotations=π/2, markersize=ms-4, color=:black)
	iscatter!(ax, [0, L+1], [2,1]*1.66;inverse=inverse, marker='⬣', rotations=π/2, color=(color, 0.5), markersize=ms-4, strokewidth=1, strokecolor=(:black, 0.5), kwargs...)
end

# ╔═╡ a5d1d184-1140-4ff8-b284-f5f2a5b57743
function plotedges!(ax, nearest, nextnearest; color=:black, kwargs...)
	for i in findall(isone.(nearest))
		x = [i, i+1]
		y = mod1.(x, 2)
		ilines!(ax, x, y*1.66; color=color, label="", kwargs...)
	end
	for i in findall(isone.(nextnearest))
		x = [i, i+2]
		y = mod1.(x, 2)
		ilines!(ax, x, y*1.66; color=color, label="", kwargs...)
	end
end

# ╔═╡ 82872e47-5e85-4b3a-ba65-f055af18e398
plotedges!(ax, L::Int; kwargs...) = plotedges!(ax, ones(L-1), ones(L-2); kwargs...)

# ╔═╡ 0ab0d142-6460-4fbb-9ddf-3198e9d6fbee
function plotlattice!(ax, L; inverse=false, alpha=1, nodesize=8, kwargs...)
	plotedges!(ax, L; alpha=alpha, inverse=inverse, kwargs...)
	plotnodes!(ax, ones(L); inverse=inverse, color=:black, label=:none, nodesize=nodesize)
end

# ╔═╡ e56fe691-3fa5-4835-8136-6eeea1ca3e25
function plotstatemakie(f, ψ; ms=nothing, inverse=false, tiny=false, nodesize=12, edgewidth=2, edgealpha=1, width=250, xticks=nothing, xlabel="", arrowsize=4, kwargs...)
	xdecoheight = 45
	L = length(ψ)

	aspect = 1.66
	height = width*aspect*3/(L+1) + xdecoheight
	axisheight = height - xdecoheight
	axiswidth = width
	if isnothing(ms)
		ms = axisheight / 2 - 6
	end
	if isnothing(xticks)
		xticks=1:1:L
	end
	if inverse
		aspect = 1/aspect
		width, height = height, width
		axiswidth, axisheight = axisheight, axiswidth
		ticks = (yticks=xticks, ylabel=xlabel)
	else
		ticks = (xticks=xticks, xlabel=xlabel)
	end
	ax = Axis(f, aspect=DataAspect(), height=axisheight; ticks...)
	if tiny || !inverse
		hideydecorations!(ax)
	end
	if tiny || inverse
		hidexdecorations!(ax)
	end
	hidespines!(ax)
	plotlattice!(ax, L; inverse=inverse, alpha=edgealpha, linewidth=edgewidth, nodesize=nodesize)
	if inverse
		ilimsx, ilimsy = ylims!, xlims!
	else
		ilimsx, ilimsy = xlims!, ylims!
	end
	if tiny
		ilimsy(ax, 0.2*1.66, 2.8*1.66)
	else
		ilimsy(ax, 0.0, 3*1.66)
	end
	ilimsx(ax, 0, L+1)
	plotstate!(ax, ψ; color=:tomato2, inverse=inverse, ms=ms, markershape=:hexagon, kwargs...)

	nc, nnc = contraintedges(ψ)
	ni, nni = interactionedges(ψ)
	plotinteractions!(ax, nc .& ni, nnc .& nni; inverse=inverse, lw=arrowsize)
	ax
end

# ╔═╡ 962cc8a8-6ab3-4faf-b6c0-4bab765f4aaa
strtoarray(s) = BitArray([c == '1' ? 1 : 0 for c in filter(x -> ~isspace(x), s)])

# ╔═╡ d9480c26-43fd-4af4-a881-a2365c84ae11
plotstatemakie(s::String; kwargs...) = plotstatemakie(strtoarray(s); kwargs...)

# ╔═╡ 1ebf68a1-9cae-44e5-9b35-b88ad11c10c2
plotstatemakie(f, s::String; kwargs...) = plotstatemakie(f, strtoarray(s); kwargs...)

# ╔═╡ Cell order:
# ╠═0473271b-b8ac-4260-ba8b-dcc044a4584d
# ╠═47c268dc-b433-407a-b22b-89f1e32e9321
# ╠═735764ef-9b8e-4c4e-96f9-1ceb63d15a0a
# ╠═a98ced4c-1d7d-4909-bd30-5d143c12abdb
# ╠═6d464ced-bf05-4057-8a8c-8f9720768ca5
# ╠═2a5f3965-ac1f-4e41-9236-7a1bacce2944
# ╠═b48c9874-f8a9-41f5-9eca-565ffc555398
# ╠═9abb9cc8-015b-4036-8678-c15293754da7
# ╠═0a7f393a-7739-4b55-b7d6-3e2ccb49b8dd
# ╠═7adc5c8d-d723-4932-8b97-cc2a86b49bef
# ╠═9263b1bc-cadb-4d43-9ccc-47e426145c01
# ╠═0b6bd6bd-eebe-4b02-b291-97f31f5b5d8b
# ╠═ce3e4422-19ba-4c52-b14f-b9fe9d307f93
# ╠═c32b1f1c-3de7-4cb1-ba7b-92062f51e536
# ╠═a2e60818-4f15-4841-8d17-8f4d3f447a3d
# ╠═f131f150-ef5d-4bd9-9511-01391648b982
# ╠═e62552c3-7ec6-43be-b77d-f595e33b6fdb
# ╠═843efd9a-7724-4830-95a7-d888812bcdf2
# ╠═c14fc29f-e686-4a6c-8454-990fcacd0cac
# ╠═0670c838-4b59-4473-8d0d-00d5bd00de0c
# ╠═fd67b127-8d6b-4868-a42f-eac31611c203
# ╠═e20d3971-3c6a-4765-b894-1441df13331b
# ╠═62e53ecf-744a-4664-bd1e-8f680d3915d0
# ╠═a918682d-ab61-4ccc-bdde-3cfa56f8ebd7
# ╠═66ed823b-9d5e-47c9-94f1-769dc5211e63
# ╠═8fd77707-01aa-4e32-8d81-e94cbc18c52b
# ╠═fe86b1e8-f5eb-4081-97b7-f5a61b9664a4
# ╠═ca87392d-9633-4828-9a91-f904da80c4f8
# ╠═7585ce3f-53fa-431f-be83-180bf713ac9a
# ╠═ff79e886-5192-4c73-bb47-c883276b923b
# ╠═df80a991-9bac-4e82-b538-e42219cb0093
# ╠═d9480c26-43fd-4af4-a881-a2365c84ae11
# ╠═07d29080-cb78-45ba-9ed3-84b425f4178a
# ╠═46dbdeca-877b-4ad2-b648-219926b45361
# ╠═54602680-9bdf-4cf1-8777-c7a07d716fe6
# ╠═fbeee4e5-d4ae-473c-b0da-db01db7b3c0a
# ╠═ace3f01f-a663-4a0e-8b50-4f5cee2c560b
# ╠═b1908130-1f3e-43a1-8252-fa13d36bb7da
# ╠═d6768f2e-b6d8-42a2-bc1c-38c5c6a30130
# ╠═a03bf078-64c5-425e-9efa-36a3021b7d7f
# ╠═498b8e17-1f21-4830-a497-65e07f21a649
# ╠═26b1ea31-29a2-4369-956c-53966900b629
# ╠═3eecf733-d4c9-43a4-8850-8088c1e3b817
# ╠═3b922a5c-a896-4198-aa3e-921fcbd57e86
# ╠═a5d1d184-1140-4ff8-b284-f5f2a5b57743
# ╠═82872e47-5e85-4b3a-ba65-f055af18e398
# ╠═0ab0d142-6460-4fbb-9ddf-3198e9d6fbee
# ╠═e56fe691-3fa5-4835-8136-6eeea1ca3e25
# ╠═962cc8a8-6ab3-4faf-b6c0-4bab765f4aaa
# ╠═1ebf68a1-9cae-44e5-9b35-b88ad11c10c2
