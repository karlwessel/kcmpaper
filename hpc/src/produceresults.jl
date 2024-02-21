using DrWatson
@quickactivate

using PlutoLinks
using LinearAlgebra

Model = @ingredients("model.jl")

occupation(L, N, i) = diag(Model.tobasis(Model.occupation(L, i), N))

function evolvestate(state, T, et)
	v = Model.basisvector(state)
	vev = T*(et .* (T'*v))
end

expvaldiag(d, ψs) = sum(d .* abs2.(ψs), dims=1)

function autocorr(b0, ψs, n)
	L = length(b0)
	occupiedsites = (1:L)[b0 .> 0]
	N = length(occupiedsites)
	@assert N > 0
	
	C = sum(n[b0])	
	
	ϕ = N/L
	c = expvaldiag(C, ψs)
	c = (c .- L * ϕ^2) / (L*ϕ*(1-ϕ))
	return c
end

function loadorproduce(args...; kwargs...)
    try
        produce_or_load(args...; kwargs...)
    catch e
        @warn("Error while producing $(args...). Recreating result.")
        if :force in keys(kwargs)
            delete!(:force, kwargs)
        end
        produce_or_load(args...; force=true, kwargs...)
    end
end

function producedecomposition(c)
    @unpack L, N, J, V, W, seed, unconstrained = c
    base = Model.basis(L, N)
    H = Model.hamiltoniandisordernew(base, V, W, seed; J=J, unconstrained=unconstrained)
    eigdecomp = eigen(Symmetric(Matrix(H)))
end

function producedecompositionold(c)
    @unpack L, N, J, V, W, seed, unconstrained = c
    base = Model.basis(L, N)
    H = Model.hamiltoniandisorder(base, V, W, seed; J=J, unconstrained=unconstrained)
    eigdecomp = eigen(Symmetric(Matrix(H)))
end

function producedecompositionfrompot(c)
    @unpack L, N, J, V, unconstrained, disorder = c
    base = Model.basis(L, N)
    H = Model.hamiltoniandisorderfrompot(base, V, disorder; J=J, unconstrained=unconstrained)
    eigdecomp = eigen(Symmetric(Matrix(H)))
end

function produceautocorrs(c; decomp = producedecomposition)
    @unpack L, N, ts = c
	
	if ts == "detail"
		ts = 10. .^ range(-1, 9, length=250) / 2
	elseif ts == "rough"
		ts = 10. .^ range(0, 8, length=50) / 2
	elseif ts == "lifetime"
		ts = 10. .^ range(-1, 15, length=250) / 2
	elseif ts == "meta"
		ts = collect(range(900, 1100, length=20)) / 2
	elseif ts == "inftime"
		ts = collect(range(0.9e15, 1.1e15, length=100)) 
	end

	basis = Model.basis(L, N)
	D, T = decomp(c)

	et = exp.(-im*D*ts')
	ns = [occupation(L, N, k) for k in 1:L]
	
	corrs = [autocorr(b0, evolvestate(b0, T, et), ns) for b0 in basis]
	[reshape(x, length(x)) for x in corrs], ts
end

function autocorrop(b0, n)
	L = length(b0)
	occupiedsites = (1:L)[b0 .> 0]
	N = length(occupiedsites)
	@assert N > 0
	ϕ = N / L
	C = Diagonal(sum(n[b0]))
	C / (L*ϕ*(1-ϕ)) - (I*(ϕ / (1-ϕ)))
end

diagensop(O, T) = abs2.(T) * diag(O)

inftimeval(O, v, Td) = only(expvaldiag(diagensop(O, Td), Td * v))

inftimeval_(b0, Td, ns) = inftimeval(autocorrop(b0, ns), Model.basisvector(b0), Td)


function inftimevals(c; decomp = producedecomposition)
        @unpack L, N = c
	D, T = decomp(c)
	Td = T'
	ns = [occupation(L, N, k) for k in 1:L]
	[inftimeval_(b0, Td, ns) for b0 in Model.basis(L, N)]
end

function producedensityevolutions(c; decomp = producedecomposition)
    @unpack L, N, ts = c
	
	if ts == "detail"
		ts = 10. .^ range(-1, 9, length=250) / 2
	elseif ts == "rough"
		ts = 10. .^ range(0, 8, length=50) / 2
	elseif ts == "lifetime"
		ts = 10. .^ range(-1, 15, length=250) / 2
	elseif ts == "meta"
		ts = collect(range(900, 1100, length=20)) / 2
	end

	basis = Model.basis(L, N)
	D, T = decomp(c)

	et = exp.(-im*D*ts')
	ns = [occupation(L, N, k) for k in 1:L]
	evols = [evolvestate(b0, T, et) for b0 in basis]
	dens = [expvaldiag(O, evol) for O in ns, evol in evols]
	dens, ts
end

# --- densities
function produceeigendensities(c; decomp = producedecomposition)
    @unpack L, N = c
    
    energies, Te = decomp(c)

	res = []
	
	for site in 1:L
		O =  occupation(L, N, site)
		dimersfast = expvaldiag(O, Te)
		
		push!(res, dimersfast[:])
	end
	
	return res, energies
end

# --- Entropies
function entropy(C)
	eval = svd!(C, full=false).S .^ 2
	D = Diagonal(eval[eval .!= 0])
	-tr(D*log(D))
end

function tomat(v, basis, idx, idx2, Mh)
	C = zeros(ComplexF64, Mh, Mh)
	for (i, (b, n, m)) in enumerate(zip(basis, idx, idx2))
		C[n,m] = v[i]
	end
	C
end

function entropy(v, basis)
	L = length(first(basis))
	hL = Int(L/2)

	halfbase = unique([b.data[1:hL] for b in basis])
	idx = indexin([b.data[1:hL] for b in basis], halfbase)
	idx2 = indexin([b.data[hL+1:end] for b in basis], halfbase)
	Mh = length(halfbase)

	[entropy(tomat(s, basis, idx, idx2, Mh)) for s in eachcol(v)]
end


function produceeigenentropies(c; decomp = producedecomposition)
    @unpack L, N, V, W, J, unconstrained = c
    
	basis = Model.basis(L, N)
	energies, T = decomp(c)
	
	entropies = entropy(T, basis)
	reshape(entropies, length(entropies)), energies
end

strtoarr(s) = BitVector([x == '1' ? 1 : 0 for x in s if x != '_'])

function produceentropyevolutions(c; decomp = producedecomposition)
    @unpack L, N, state, V, ts, W, J, unconstrained = c
    state = strtoarr(state)
	@assert L == length(state)
	@assert N == sum(state)
	
	if ts == "detail"
		all = range(-1, 9, length=250) / 2
		ts = 10. .^ all
	elseif ts == "rough"
		ts = unique(10. .^ [range(-1, 2, length=14); 3:16]) / 2
	elseif ts == "veryrough"
		ts = 10. .^ (-1:1:16) / 2
	end
	
	basis = Model.basis(L, N)	
	D, T = decomp(c)
	
	et = exp.(-im*D*ts')
	entropies = entropy(evolvestate(state, T, et), basis)
	reshape(entropies, length(entropies)), ts
end

function produceresult(c)
    @unpack data = c
    
    res = Dict{String, Any}(c)	
    
    if data == "eigendensities"
        densities, energies = produceeigendensities(c)
        res["energy"] = energies
	    res["densities"] = densities
	elseif data == "autocorrelationevolutions"
        autocorr, ts = produceautocorrs(c)
	    res["autocorr"] = autocorr
	    res["ts"] = ts
	elseif data == "eigenentropies"
        entropies, energies = produceeigenentropies(c)
        res["energy"] = energies
	    res["entropy"] = entropies
	elseif data == "entropyevolutions"
        entropies, ts = produceentropyevolutions(c)
	    res["entropy"] = entropies
	    res["ts"] = ts
	    elseif data == "diagensemble"
        autocorr = inftimevals(c)
	    res["autocorr"] = autocorr
	else
	    throw("Unknown result type: $data")
    end
    
    res
end

scratchpath = datadir("raw")
if isdir("/scratch/users")
    scratchpath = joinpath("/scratch/users", ENV["USER"], "raw")
end
    
function produceorload(c; subpath=nothing, kwargs...)
    @unpack data, W, L, filling = c
    c["N"] = Int(L * filling)
    delete!(c, "filling")
    if isnothing(subpath)
        subpath = data
    end
    path = joinpath(scratchpath, subpath, "W=$(round(W, sigdigits=3))")
    loadorproduce(produceresult, c, path; kwargs...)
end
    
