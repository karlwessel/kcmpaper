### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ cb14c839-4e5b-4f87-a97b-d8dd157cf443
using DrWatson

# ╔═╡ 079b8f84-173e-400e-bff7-846c15ed2828
@quickactivate

# ╔═╡ 3d5efb00-5776-4831-9e36-affb8d62c542
using DataFrames

# ╔═╡ 08981398-4465-4691-8434-e52a239429e5
using Statistics

# ╔═╡ e4e65bde-9ce7-49df-a701-4e4054fe8df8
using Memoize

# ╔═╡ a8af5681-cb31-4abd-980e-57f7fec60559
using PlutoLinks

# ╔═╡ e31b01b6-4f79-4442-a5dd-551e23faff8c
model = @ingredients(srcdir("model.jl"))

# ╔═╡ d8fdc06a-9290-4d6b-b20a-d002ec2ee069
md"""
# Filter result files before loading
"""

# ╔═╡ 780adc73-aecd-4389-88cc-6785aa5c95ef
parseboolorfloat(x) = 
	if '.' in x
		try
			parse(Float64, x)
		catch
			x
		end
	elseif x == "true"
		true
	elseif x == "false"
		false
	else
		try
			parse(Int64, x)
		catch
			x
		end
	end

# ╔═╡ 1425755b-da8a-4863-967d-1f97f69eec9a
parseboolorfloat("true")

# ╔═╡ 33695655-4340-4e2f-be2e-d3db11502e5a
"""
parsepath(path, keys)

Parse specific configuration keys from file names of result files.
"""
function parsepath(p, keys)
	matches = [k => match(Regex(".*$k=([^\\.^_]+)"), p) for k in keys]
	[k => (isnothing(m) ? missing : parseboolorfloat(m[1])) for (k, m) in matches]
end

# ╔═╡ 111d2e32-5b79-4510-9f05-75bb492f7921
"""
collectresultfiles(path, keys)

Collect all result files in the passed directory into a dataframe that also contains the configuration of the results parsed from the filenames.

This does not load the data in the result files yet and allows to prefilter results based on their configuration.
"""
function collectresultfiles(path, keys)
   files = []
   for (dir, _, f) in walkdir(path)
       for x in f
           push!(files, (dir, x))
       end
   end
	
	df = DataFrame([Dict(parsepath(f, keys)) for (dir, f) in files])
	df.path = [joinpath(x...) for x in files]
	df
end

# ╔═╡ bf83b424-3223-4848-a124-adf6ba569887
"""
loadresults(files)

Load the data in the passed list of files into one big DataFrame, similar to DrWatsons `collectresults` method.
"""
function loadresults(files)
	res = []
	for f in files
		try
			push!(res, DrWatson.to_data_row(wload(f), f))
		catch e
			@warn("Could not load file $f: $e")
		end
	end
	vcat(res...)
	# vcat([DrWatson.to_data_row(wload(f), f) for f in files]...)
end

# ╔═╡ d225a599-31ca-4902-8784-907124dc56a1
md"""
# Flatten rows and add index as new column
"""

# ╔═╡ 02a6f357-0b56-4654-be88-45643013968b
groupon(df, columns) = groupby(df, setdiff(Symbol.(names(df)), columns))

# ╔═╡ 56fa1c70-9475-461b-8c23-b77a16071cf0
"""
unpack(df, columns, indexname)

Same as DataFrames `flatten` but adds a new column containing the indices of array that was flattened.
"""
function unpack(df, columns, indexname)
	grouped = groupon(df, columns)
	maps = [[k => first => k for k in columns];
		first(columns) => (x -> eachindex(first(x))) => indexname]
	combine(grouped, maps)
end

# ╔═╡ bf3bfedd-8560-4639-adf9-689449fec3d2
"""
unpack(df, columns)

Same as `flatten` but maybe faster or less memory consuming?
"""
function unpack(df, columns)
	grouped = groupon(df, columns)
	maps = [k => first => k for k in columns]
	combine(grouped, maps)
end

# ╔═╡ 3015f3dd-921d-4dfb-bd74-38f472ea8d3c
md"""
# Disorder uncertainty
"""

# ╔═╡ f37b0cd2-31d3-4d81-953a-1981babeccb3
"""
uncertainty(\$X\$)

Returns the uncertainty in the mean value of \$X\$ as 

\$\$\\sqrt{\\frac{VAR(X)}{N-1}}\$\$

where \$VAR(X)\$ is the variance of the data in \$X\$ and \$N\$ is the number of samples in \$X\$.
"""
uncertainty(x) = std(x) / sqrt(length(x) - 1)

# ╔═╡ cc29033a-66e9-4130-8dcb-bef3c8d05987
md"""
# Time average
"""

# ╔═╡ 6d7b3dd9-cb45-4a20-b3d0-06e3c4e4c5af
norm(x, dts, norms) = cumsum(x .* dts) ./ norms

# ╔═╡ cefea845-b441-407a-9e0f-8aa22815c72e
"""
Time average 

\$\\bar y_i = \\frac{1}{t_i - t_0} \\sum_{j = 1}^i y_i(t_i - t_{i-1})\$ 

of samples \$y_i\$ at times \$t_i\$ and \$t_0 = 0\$.

"""
function blocktimeavrg(ts, blocks)
	dts = diff([0; ts])
	norms = ts
	norm.(blocks, Ref(dts), Ref(norms))
end

# ╔═╡ 41a66486-b4f2-4622-b363-95b95e8006b7
md"""
# Classify states
"""

# ╔═╡ 8a9e836a-23cf-483e-bf7a-27789c442593
function convolve(x, k)
	K = length(k)
	r = 0
	for i in eachindex(x)
		v = circshift(x, i)
		r += all(k .== v[begin:K])
	end
	r
end

# ╔═╡ 8943cb4d-cc4e-4e50-88c5-1508e1a9f7ed
@memoize arrtostr(a) = String([x == 0 ? '0' : '1' for x in a])

# ╔═╡ 21646d2b-9705-41eb-a987-8f1e3c89e932
@memoize holecount(x) = convolve(x, [1,1,0,1,1])

# ╔═╡ 573ac072-cec4-4738-94f8-b3198d5aa94c
@memoize hdimercount(x) = convolve(x, [0,1,0]) + convolve(x, [0,0,0])

# ╔═╡ 4a0b5417-5e16-4f56-a909-a3a96d79b50d
@memoize vdimercount(x) = convolve(x, [0,0])

# ╔═╡ b839fcbf-4fcf-4d5b-86cb-07f7905becfd
function energies(L, N)
	basis = model.basis(L, N)
	H = model.hamiltoniandisordernew(basis, 1, 0, 1, J=1, unconstrained=false)
	reps = filter(model.isrepresentator, basis)
	vs = model.basisvector.(reps)
	energy = [v'*(H*v) for v in vs]
	DataFrame(representator = reps, energy = energy, L=L, N=N)
end

# ╔═╡ 1ac41f54-8ad5-4d45-bb8a-0c68740bce86
function classifystates(df)
	dfs = []
	for (k, dfsub) in pairs(groupby(df, [:L, :N]))
		reps = model.representator.(model.basis(k.L, k.N)[dfsub.state])
		dfsub.representator = reps
		dfenerg = energies(k.L, k.N)
		dfsub = leftjoin(dfsub, dfenerg, on=[:L, :N, :representator])
		push!(dfs, dfsub)
	end
	df = vcat(dfs...)
	df.holes = holecount.(df.representator)
	df.vdimers = vdimercount.(df.representator)
	df.hdimers = hdimercount.(df.representator)

	df.representator = arrtostr.(df.representator)
	df = transform(groupby(df, [:L, :N, :holes, :hdimers]), :representator => (x -> indexin(x, unique(x))) => :instance, :state => (x -> indexin(x, unique(x))) => :stateinstance)
	df
end

# ╔═╡ 8699013b-041d-4bea-ad66-fbf41938a24e


# ╔═╡ Cell order:
# ╠═cb14c839-4e5b-4f87-a97b-d8dd157cf443
# ╠═079b8f84-173e-400e-bff7-846c15ed2828
# ╠═3d5efb00-5776-4831-9e36-affb8d62c542
# ╠═08981398-4465-4691-8434-e52a239429e5
# ╠═e4e65bde-9ce7-49df-a701-4e4054fe8df8
# ╠═a8af5681-cb31-4abd-980e-57f7fec60559
# ╠═e31b01b6-4f79-4442-a5dd-551e23faff8c
# ╟─d8fdc06a-9290-4d6b-b20a-d002ec2ee069
# ╠═1425755b-da8a-4863-967d-1f97f69eec9a
# ╠═780adc73-aecd-4389-88cc-6785aa5c95ef
# ╠═33695655-4340-4e2f-be2e-d3db11502e5a
# ╠═111d2e32-5b79-4510-9f05-75bb492f7921
# ╟─bf83b424-3223-4848-a124-adf6ba569887
# ╟─d225a599-31ca-4902-8784-907124dc56a1
# ╠═02a6f357-0b56-4654-be88-45643013968b
# ╟─56fa1c70-9475-461b-8c23-b77a16071cf0
# ╟─bf3bfedd-8560-4639-adf9-689449fec3d2
# ╟─3015f3dd-921d-4dfb-bd74-38f472ea8d3c
# ╟─f37b0cd2-31d3-4d81-953a-1981babeccb3
# ╠═cc29033a-66e9-4130-8dcb-bef3c8d05987
# ╠═6d7b3dd9-cb45-4a20-b3d0-06e3c4e4c5af
# ╠═cefea845-b441-407a-9e0f-8aa22815c72e
# ╠═41a66486-b4f2-4622-b363-95b95e8006b7
# ╠═8a9e836a-23cf-483e-bf7a-27789c442593
# ╠═8943cb4d-cc4e-4e50-88c5-1508e1a9f7ed
# ╠═21646d2b-9705-41eb-a987-8f1e3c89e932
# ╠═573ac072-cec4-4738-94f8-b3198d5aa94c
# ╠═4a0b5417-5e16-4f56-a909-a3a96d79b50d
# ╠═b839fcbf-4fcf-4d5b-86cb-07f7905becfd
# ╠═1ac41f54-8ad5-4d45-bb8a-0c68740bce86
# ╠═8699013b-041d-4bea-ad66-fbf41938a24e
