### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ bccd3a32-53f9-11ec-1a7e-c12044b6accc
using PlutoTest

# ╔═╡ fa07a6a5-8fa7-41f0-88f8-b233f9d60735
using BitBasis

# ╔═╡ a7dfe927-0b2e-4734-aa69-867ece13c55e
using Combinatorics

# ╔═╡ f77ad6eb-8653-4176-90c1-eb102e1d7c18
using CircularArrays

# ╔═╡ 09e81a27-df5b-4661-bd72-18c38d37593e
using LinearAlgebra

# ╔═╡ 3cca24e9-9d98-4b1b-a017-fe81bbbffdc9
using SparseArrays

# ╔═╡ 7b2dc7c7-8a01-4e4c-ab55-fd30d1a47617
using Random

# ╔═╡ b7f5b3ce-a507-481b-8f6c-0d0e57421105
bittoint(v::CircularArrays.CircularVector) = bittoint(v.data)

# ╔═╡ 49f20b9e-7b35-4414-a15c-e68396513500
bittoint(v::BitArray) = only(v.chunks)::UInt64

# ╔═╡ 8dee19cd-bbd4-4ab6-ac78-7cffad0b197f
basisvector(v::AbstractArray) = basisvector(length(v), bittoint(v))

# ╔═╡ 4e2f0e33-b04c-4b0c-bde8-19b3270790de
function basis(L, N)
	v = vcat(ones(Bool, N), zeros(Bool, L-N))
	sort!(CircularArray.(BitArray.(multiset_permutations(v, length(v)))), by=bittoint)
end

# ╔═╡ cf77d289-3c3f-44e1-8cae-aadc995c2ace
⊗ = kron

# ╔═╡ 24825184-6350-432b-a2d9-a0b3ae3bd416
id(L) = L == 0 ? 1 : Diagonal(ones(2^L))

# ╔═╡ d2319a30-0554-4749-9ffd-8749de29c55f
anhilation(L, i) = id(L - mod1(i, L)) ⊗ sparse([0 0; 1 0]) ⊗ id(mod1(i, L) - 1)

# ╔═╡ 91dbf769-5033-4a45-903b-0a06503d61b7
occupation(L, i) = id(L - mod1(i, L)) ⊗ sparse([0 0; 0 1]) ⊗ id(mod1(i, L) - 1)

# ╔═╡ e5746f41-9ce9-4560-a483-8e0b46e7ef99
@test onehot(3, Int(bit"101"))' * occupation(3, 1) * onehot(3, Int(bit"101")) == 1

# ╔═╡ 1f646620-daae-4344-b49a-d2a3c6f6e950
@test onehot(3, Int(bit"101"))' * occupation(3, 2) * onehot(3, Int(bit"101")) == 0

# ╔═╡ 8552a30d-4fd7-49a8-adc9-0c53a627c1f3
@test onehot(3, Int(bit"001"))' * occupation(3, 1) * onehot(3, Int(bit"001")) == 1

# ╔═╡ 6e182939-51c9-4e69-b89c-4b3b4b67f017
@test onehot(3, Int(bit"111"))' * occupation(3, 2) * onehot(3, Int(bit"111")) == 1

# ╔═╡ 5f1582ef-4ba9-4360-8aef-c4fe15b2d035
function hamiltonian(L, λ)
	H = spzeros(2 ^ L, 2 ^ L)
	for i in 1:L
		σᵢ⁺ = anhilation(L, i)'
		nᵢ = occupation(L, i)
		for j in (i + 1):(i + 2)
			σⱼ = anhilation(L, j)
			nⱼ = occupation(L, j)
			Vⁱʲ = nᵢ * (I - nⱼ) + nⱼ * (I - nᵢ)
			Tⁱʲ = σᵢ⁺ * σⱼ
			Tⁱʲ += Tⁱʲ'
			Cⁱʲ = j == (i + 1) ? 
				I - occupation(L, i - 1) * occupation(L, i + 2) : 
				I - occupation(L, i + 1)
			H -= 0.5 * Cⁱʲ * (λ * Tⁱʲ - (1 - λ) * Vⁱʲ)
		end
	end
	H
end

# ╔═╡ 2a7b8063-0555-4986-a755-cd04d2485178
tobasis(A::AbstractMatrix, basis::AbstractVector{UInt}) = A[basis .+ 1, basis .+ 1]

# ╔═╡ ebdfebd3-3ba3-47ac-b058-1aa439e2bcb5
tobasis(v::AbstractVector, basis::AbstractVector{UInt}) = v[basis .+ 1]

# ╔═╡ e25082e3-595a-440d-9583-b16f2ca513cd
tobasis(A, N::Int) = tobasis(A, basis(Int(log2(size(A, 1))), N))

# ╔═╡ 9de9dd47-3c12-4c33-b4bb-88e306214b3f
tobasis(A, basis) = tobasis(A, bittoint.(basis))

# ╔═╡ 9ce3260c-79c2-458b-999d-0a5fe1976e00
basisvector(L, k) = tobasis(onehot(L, k), count_ones(k))

# ╔═╡ b6cc94f2-4639-476d-849e-a6da803b4665
@test basisvector(3, 1) == [1, 0, 0]

# ╔═╡ a2488d77-f34a-4c7b-b76e-7aaf2adbb99d
@test basisvector(bitarray(1, 3)) == basisvector(3, 1)

# ╔═╡ d8bdf31e-dc62-4e94-bc9d-bdd76f02e914
hamiltonian(L, λ, N) = tobasis(hamiltonian(L, λ), basis(L, N))

# ╔═╡ 04b340b4-92a8-4904-b1a8-334826155001
function applyhopping(i, j, b)
	res = copy(b)
	res[i], res[j] = res[j], res[i]
	return !b[i] && b[j], res
end

# ╔═╡ 0b35ddf6-1e34-4f16-9d4d-64ba93a555af
applyinteraction(i, j, b) = b[i] != b[j]

# ╔═╡ b73e857a-f633-49f4-b16a-6a96530fe67b
function applyconstraint(i, j, b)
	if j == i + 1
		return !b[i-1] || !b[j+1]
	elseif j == i + 2
		return !b[i+1]
	else
		return 0
	end
end

# ╔═╡ 2d7eb900-2eb9-4848-a2ee-6e2d6734d7ce
begin
	mutable struct SparseFactory{T}
		rows::Vector{Int64}
		cols::Vector{Int64}
		vals::Vector{T}
		m
		n
	end
	SparseFactory(m, n, T=Float64) = SparseFactory(Int64[], Int64[], T[], m, n)
end

# ╔═╡ 4b3add98-cfb1-40ba-994a-e5efdcb0de7b
function Base.push!(A::SparseFactory, i, j, v) 
	push!(A.rows, i)
	push!(A.cols, j)
	push!(A.vals, v)
end

# ╔═╡ b745228f-78f9-47eb-a20e-15ef7874d08c
Base.setindex!(A::SparseFactory, v, i, j) = push!(A, i, j, v)

# ╔═╡ e5ba77f2-e388-4d17-96c5-8dcd9b12b949
Base.getindex(A::SparseFactory, keys...) = zero(eltype(A.vals))

# ╔═╡ 10d116fd-e6e9-4ef5-af04-7b5c9a87602d
md"""
# Matrix representation of Hamiltonian
"""

# ╔═╡ 1fe2ea75-b853-4319-845d-82b23a9e3cc4
md"""
The Hamiltonian of the system with $L$ sites is 

$$H = -\frac 12 \sum_{i=1}^L \sum_{j=i+1}^{i+2} C^{ij}(\lambda T^{ij} - (1-\lambda)V^{ij}).$$

The constraint between sites $i$ and $j$ is defined as 

$$C^{ij} = \begin{cases}
1 - n_{i-1}n_{j+2} & j = i+1\\
1 - n_{i+1} & j = i+2\\
0 & \text{else}
\end{cases}$$

where the occupation operator $n_i = \ket{1_i}\bra{1_i}$ for a site $i$ is defined by the basis vector $\ket{1_i}$ of site $i$ occupied.

The hopping term between sites $i$ and $j$ is defined as

$$T^{ij} = \sigma_i^\dagger\sigma_j  + \text{h.c.}$$

with anhilation and creation operators $\sigma_i = \ket{0_i}\bra{1_i}$ and $\sigma_i^\dagger=\ket{1_i}\bra{0_i}$ respectively.

The interaction between sites $i$ and $j$ is defined as

$$V^{ij} = n_i(1-n_j).$$

The hopping frequency parameter $0<\lambda<1$ weights between hopping and interaction strength.
"""

# ╔═╡ d44936cd-98db-4383-b517-801b390d8305
md"""
### Choice of basis
A convenient basis to represent the Hamiltonian is the occupation number basis 

$$\ket{b_1 \dots b_L} = \otimes_{i=1}^L \begin{cases}\ket{0_i} & b_i = 0\\
\ket{1_i} & b_i = 1\end{cases}$$

consisting of all the product states of the per site occupation number basis vectors $\ket{0_i}$ and $\ket{1_i}$ at sites $i$.

We can represent each product state as the binary notation of an integer.

$$\ket{b_1\dots b_L} = \ket{k} \quad \text{where } k = \sum_{i=1}^L 2^{b_i}.$$

This gives a unique ordering of our basis vectors and their vector representation is

$$\ket{k} = (0, \dots, 0, 1, 0, \dots, 0)^T$$

where the one is at the $(k+1)$th position in the vector.
"""

# ╔═╡ c9a29f4d-b01c-4c13-8af3-faece6670b37
md"""
### Conservation of particle count
Since particle count is conserved the Hamiltonian will have block structure for
basisvectors with same particle number.

It therefore is convenient to collect all basis vectors with the same particle number by chosing one basis vector of that particle count, all other basis vectors with the same particle count are unique permutations of bits of that basis vector.
"""

# ╔═╡ 588fff8e-0390-4f1d-96e6-3b57ffc539ba
md"""
## Tensor product notation

Any local operator $A_i'$ of dimension 2 for a site $i$ can be transfered to an oeprator $A_i$ or dimension $2^L$ in Hilbert space of the system by tensor products with unities $I^k$ of dimension $k$ 

$$A_i = I^{2^{i-1}} \otimes A_i' \otimes I^{2^{L-i}}.$$

The local anhilation operator at site $i$ is

$$\sigma_i' = \begin{pmatrix}0 & 0\\ 1&0\end{pmatrix}$$

represented in occupation number basis $\{\ket{0}, \ket{1}\}$.

The creation operator $\sigma_i^\dagger$ is the complex conjugate of the anhilation operator.

Accordingly the matrix representation of the local occupation number is 

$$n_i' = \begin{pmatrix}0 & 0\\ 0&1\end{pmatrix}.$$

This way the operators of the Hamiltonian can be defined as sparse matrices and the Hamiltonian can be directly calculated.
"""

# ╔═╡ 4b9fd419-5322-42e2-8da2-6d6228ae0c2f
md"""
When considering conservation of particles the Hamiltonian can be reduced to entries respecting a certain particle count.
"""

# ╔═╡ 03781a73-44bd-4a19-bbf2-4369287b106c
md"""
## Iterating through basis vectors

When applying single terms $T^{ij}$, $V^{ij}$ or $C^{ij}$ the result is a scaled version of at most one other basis vector. 
For the Interaction- and constraint terms $V^{ij}$ or $C^{ij}$ the result is the same scaled basis vector. 
For the Hopping term $T^{ij}$ it is a different basis vector.

The matrix elements of the Hamiltonian can therefore be cosntructed by summing up the scalings of the terms for the matrix element defined by the source basis vector and the resulting basis vector.

The effect of the Hopping operator

$$T^{ij}\ket{b_1, \dots, b_N} = \begin{cases}
\ket{\dots, b_i+1, \dots, b_j-1, \dots} & b_i = 0 \land b_j = 1\\
0 & \text{else}
\end{cases}$$

to a basis vector $\ket{b_1, \dots, b_N}$ is moving the particle from the $j$th to the $i$th site, if that is possible.


"""

# ╔═╡ 8bdb4cc9-a7ce-47f0-8c8d-1dccf46486f1
md"""
The effect of the Interaction operator

$$V^{ij}\ket{b_1, \dots, b_N} = \begin{cases}
\ket{b_1, \dots, b_N} & b_i = 1 \land b_j = 0\\
0 & \text{else}
\end{cases}$$

is an interaction if the $i$th site is occupied and the $j$th site is free. The resulting basis vector is unchanged.
"""

# ╔═╡ de6b23f6-8185-4baf-870c-baec04aa49d3
md"""
Similar the Contraint operator

$$C^{ij}\ket{b_1, \dots, b_N} = \begin{cases}
\ket{b_1, \dots, b_N} & j=i+1 \wedge b_{i-1} = 1 \lor b_{j+1} = 1\\
\ket{b_1, \dots, b_N} & j=i+2 \wedge b_{i+1} = 1 \\
0 & \text{else}
\end{cases}$$
is allowing interaction and hopping to take effect if the common neighbors of sites $i$ and $j$ are free. The resulting basis vector is unchanged.
"""

# ╔═╡ a6c74b08-5d50-447b-a7e1-7b7c3710389f
md"""
Note that the constraint for a pair of sites $i$ and $j$ is never affecting the same sites as the interaction and hopping term. Therefore it does not matter in which order they are applied to the basis state.
"""

# ╔═╡ 3aa6c31c-a66a-474e-9af3-02a189217dff
realize(A::SparseFactory) = sparse(A.rows, A.cols, A.vals, A.m, A.n)

# ╔═╡ e3fa99e2-37ad-4973-8389-f3edf25fb36b
function hamiltonian(basis::AbstractArray, λ)
	L = length(basis[1])
	M = length(basis)
	H = SparseFactory(M, M)
	for (k, b) in enumerate(basis), i in 1:L, j in i+1:i+2
		if applyconstraint(i, j, b)
			v = applyinteraction(i, j, b)
			t, b2 = applyhopping(i, j, b)
			v && push!(H, k, k, (1 - λ)*v)
			if t
				l = searchsortedfirst(basis, b2, by=bittoint)
				push!(H, k, l, -λ*t)
				push!(H, l, k, -λ*t)
			end
		end
	end
	0.5 * realize(H)
end

# ╔═╡ ca34fbb5-4001-48d5-9358-aa93aa94ea15
function allbasis(L)
	# generate all basis vectors for system with L lattice sites
	CircularArray.(bitarray.(0:2^L-1, L))
end

# ╔═╡ 1c8455d7-5ea6-4dbc-96ff-0f15e3b462c3
md"""
## Disorder and independent hopping/interaction hamiltonian
"""

# ╔═╡ 29cc7f2d-1105-4499-b777-6b8d0d2cc57a
begin
function hamiltoniandisorder(basis, V, W, seed=0; J=1, unconstrained=false)
	L = length(basis[1])
	M = length(basis)
	H = SparseFactory(M, M)
	for (k, b) in enumerate(basis), i in 1:L, j in i+1:i+2
		if unconstrained || applyconstraint(i, j, b)
			v = applyinteraction(i, j, b)
			t, b2 = applyhopping(i, j, b)
			v && push!(H, k, k, V*v)
			if t
				l = searchsortedfirst(basis, b2, by=bittoint)
				push!(H, k, l, -J*t)
				push!(H, l, k, -J*t)
			end
		end
	end
	H₀ = 0.5 * realize(H)
	H₀ + Diagonal(rand(MersenneTwister(seed), M) * 2W .- W)
end
hamiltoniandisorderλ(basis, λ, W=0, seed=0) = hamiltoniandisorder(basis, 1-λ, W, seed; J=λ)
end

# ╔═╡ 0cd7bf45-4b62-4491-a855-6a9fef016405
function hamiltoniandisordernew(basis, V, W, seed=0; J=1, unconstrained=false)
	L = length(basis[1])
	ϵᵢ = rand(MersenneTwister(seed), L) * 2W .- W	
	hamiltoniandisorderfrompot(basis, V, ϵᵢ; J=J, unconstrained=unconstrained)
end

function hamiltoniandisorderfrompot(basis, V, pot; J=1, unconstrained=false)
	L = length(basis[1])
	M = length(basis)
	H = SparseFactory(M, M)
	for (k, b) in enumerate(basis), i in 1:L, j in i+1:i+2
		if unconstrained || applyconstraint(i, j, b)
			v = applyinteraction(i, j, b)
			t, b2 = applyhopping(i, j, b)
			v && push!(H, k, k, V*v)
			if t
				l = searchsortedfirst(basis, b2, by=bittoint)
				push!(H, k, l, -J*t)
				push!(H, l, k, -J*t)
			end
		end
	end
	H₀ = realize(H)
	# H₀ = 0.5 * realize(H)
	ϵᵢ = pot
	H₀ + sum([tobasis(occupation(L, i), basis) * x for (i, x) in enumerate(ϵᵢ)])
end

# ╔═╡ dd56efac-37a3-44d6-a3eb-b10c0a47857a
md"""
#### Representative
To construct the momentum basis vectors we can iterate through all occupation basis vectors and calculate their momentum state representation for all compatible momenta. This would give us $R_a$ different momentum basis vectors for each $R_a$-periodic occupation basis vector $\ket a$. However all vectors resulting from $r\in\{1\to R_a\}$ translations of that vector $\ket a$ would yield momentum basis vectors parallel to the ones of $\ket a$.

Thus it is sufficient to chose one representative from the $R_a$ shifted vectors. And only construct the basis from those. We can chose the vector with the smallest integer representation as representative.
"""

# ╔═╡ dd4eb9fa-e2bc-4007-8c4b-1ccb9cd6325a
function representator(v::BitArray)
	L = length(v)
	states = [circshift(v, s) for s in 1:L]
	nums = [only(s.chunks) for s in states]
	return states[argmin(nums)]	
end

# ╔═╡ 0455bec5-421b-41dd-bd6f-7f80cd3e5029
isrepresentator(v::BitArray) = isrepresentator(v, similar(v))

# ╔═╡ fb903b76-8e50-4954-9dd5-bf03c7b2b103
function isrepresentator(v::BitArray, v2)
	L = length(v)
	i = only(v.chunks)
	for s in 1:L
		i2 = only(circshift!(v2, v, s).chunks)
		if i2 < i
			return false
		elseif i2 == i
			return true
		end
	end
	return true
end

# ╔═╡ 0bb504f7-4b1e-4790-bd82-0f6114d16ea1
function isrepresentator(v::Vector)
	v2 = similar(first(v).data)
	[isrepresentator(x.data, v2) for x in v]
end

# ╔═╡ fc7c52b8-ade1-458e-8159-a39223e698fb
representator(v::CircularVector{Bool, BitVector}) = CircularVector(representator(v.data))

# ╔═╡ 3e158d8b-96e2-42eb-ab4a-b49ded9d719c
isrepresentator(v::CircularVector{Bool, BitVector}) = isrepresentator(v.data)

# ╔═╡ a02ca109-4900-4d9e-bbc4-46d819864411
md"""
# Tests
"""

# ╔═╡ d41dce8a-5bf3-4386-8cc3-c34f2fc3997d
md"""
## identity matrix creation `id`
"""

# ╔═╡ d3d8c9c8-c59e-4e20-8f0c-ba82fb0781f3
@test id(0) == 1

# ╔═╡ 1b0ea891-14c4-41d2-93db-a75521dbcbfe
@test id(1) == [1 0; 0 1]

# ╔═╡ 64213987-41a6-4ef4-8573-8d92fa6f0a06
md"""
## anhilation operator creation `anhilation`
"""

# ╔═╡ 461dec51-777a-4437-aff2-db94a18431a8
@test anhilation(1, 1) == [0 0; 1 0]

# ╔═╡ e728da5f-a226-4906-a7e9-224e80a25a79
@test anhilation(2, 3) == anhilation(2, 1)

# ╔═╡ fc754518-ea71-4b67-b6f4-d2a1b46ca593
md"""
## Occupation operator creation `occupation`
"""

# ╔═╡ 356d6b42-19ec-4f27-9ff8-4ba3d0b0c038
@test occupation(2, 0) == occupation(2, 2)

# ╔═╡ 3f594112-2e84-4d38-b506-c6e008a88e1b
md"""
## Tensor product Hamiltonian `hamiltonian(L)`
"""

# ╔═╡ c3beb708-4b3c-4e83-88c4-2bbbf32e5c1a
@test issymmetric(hamiltonian(4, 0.5))

# ╔═╡ a4e942e9-d41f-44a4-9b9d-4a7063f444b1
@test hamiltonian(6, 0.5)[Int(bit"101101")+1, :] == zeros(2^6) # two isolated holes

# ╔═╡ b989637d-3f35-430a-9674-3ca207648812
@test hamiltonian(6, 0.5)[Int(bit"011011")+1, :] == zeros(2^6) # two isolated holes

# ╔═╡ 5da9a9fb-b914-4b50-bf20-ffc868d33e44
md"""
## Basis vectors from integer `basisvector`
"""

# ╔═╡ 2f4d7a14-d72f-4ad8-81d7-0fcca6c8183c
@test basisvector(1, 0) == [1]

# ╔═╡ 757d3b67-7f29-4568-b32d-09f466ac74bd
@test basisvector(1, 1) == [1]

# ╔═╡ f05b2133-75bb-4d0b-84c9-362d4342b824
@test basisvector(3, 0) == [1]

# ╔═╡ dc86534b-4c10-4303-9c90-32b74596cdb1
@test basisvector(3, 7) == [1]

# ╔═╡ 67c51a46-95a8-42d1-8a1f-635d980ee4b4
md"""
## basis for certain particle number `basis`
"""

# ╔═╡ 8bbff18f-2b73-4511-bd09-8a56d99167ad
@test only(basis(4, 4)) == ones(Bool, 4)

# ╔═╡ 88cdcf03-08dd-4a35-afe2-aae1e8d65efc
@test only(basis(4, 0)) == zeros(Bool, 4)

# ╔═╡ 8a7f59c9-d7d4-4dd2-82a1-7a67aba9f633
@test basis(2, 1) == [[1, 0], [0, 1]]

# ╔═╡ 030c63c0-c106-4122-949a-50394530b7f2
md"""
## Hamiltonian for certain particle count `hamiltonian(L, N)`
"""

# ╔═╡ dd961753-00b3-4811-bfa1-3f3ad4c3c8c7
@test hamiltonian(4, 0.5, 3) == zeros(4, 4) # a single free site allows no interactions

# ╔═╡ ff0085b9-faa7-4277-9775-c0faa7c01c15
@test hamiltonian(8, 0.5, 7) == zeros(8, 8) # a single free site allows no interactions

# ╔═╡ b58af3bd-9e1f-43ed-aab7-095d848145c7
@test hamiltonian(6, 0.5, 3)[2, 2] == 1.25 # specific case that checks a previous bug

# ╔═╡ 42adb6e2-4bbe-4ed7-880b-0085aa9a03c0
md"""
## Apply hopping operator to basis vector `applyhopping`
"""

# ╔═╡ fe813e92-02d0-4500-a336-5d0e735a3b72
@test applyhopping(1, 2, BitArray([0, 1, 1])) == (1, [1, 0, 1])

# ╔═╡ 38beff73-1c65-4f6c-a6fe-f3e286c96fb1
@test applyhopping(2, 3, BitArray([0, 1, 1]))[1] == 0

# ╔═╡ 5198905b-1bc5-49f7-9347-66d84d054b03
@test applyhopping(3, 4, CircularArray(BitArray([0, 1, 1])))[1] == 0

# ╔═╡ c4d9b1ff-db23-4224-b06c-336d57e223f0
md"""
## `applyinteraction`
"""

# ╔═╡ dd3d6d18-5238-4310-bab4-a2c07c2a9074
@test applyinteraction(1, 2, BitArray([1, 0]))

# ╔═╡ 0deb2d76-aa3c-4bde-bffb-13192c6561a7
@test !applyinteraction(1, 3, CircularArray(BitArray([1, 0])))

# ╔═╡ 1554f951-e241-4486-a887-d77a15a80971
md"""
## `applyconstraint`
"""

# ╔═╡ cb9c8eda-4731-4d7f-bbab-dbd4b8ea2d1b
@test !applyconstraint(2, 3, BitArray([1, 0, 1, 1, 1]))

# ╔═╡ 963fad51-f2a6-4e23-9eed-0b4ce3da436f
@test applyconstraint(2, 3, BitArray([0, 1, 1, 1]))

# ╔═╡ c1f2efd7-cd4b-4f16-b5ff-caf54ade10c3
@test !applyconstraint(1, 3, BitArray([0, 1, 1, 1]))

# ╔═╡ 1d556577-4c83-42ee-9ff0-361ce141a2e1
@test applyconstraint(4, 6, CircularArray(BitArray([0, 1, 1, 1])))

# ╔═╡ ee2d6f9a-52a2-42e2-8129-c75187cba079
md"""
## create Hamiltonian for set of basis vectors `hamiltonian(basis)`
"""

# ╔═╡ 1f705aaf-be50-4ce7-ac87-5e2c8640721d
@test hamiltonian(basis(2,1), 0.5) == hamiltonian(2, 0.5, 1)

# ╔═╡ e726a652-3898-48c0-bdb6-91f22eb5ad7e
@test hamiltonian(basis(4,2), 0.5) == hamiltonian(4, 0.5, 2)

# ╔═╡ d3d21cd9-85ea-4546-844d-6ae794567a2e
@test hamiltonian(basis(4,3), 0.5) == hamiltonian(4, 0.5, 3)

# ╔═╡ 0e0f74fd-e3df-4b0d-96cf-a928751cb43d
@test hamiltonian(basis(4,4), 0.5) == hamiltonian(4, 0.5, 4)

# ╔═╡ 0b4e6503-7a1f-458f-b292-b50e37f9dffa
@test hamiltonian(basis(6,3), 0.5) ≈ hamiltonian(6, 0.5, 3)

# ╔═╡ f99f8d50-3473-4997-9696-6d1e7b4be18b
@test hamiltonian(allbasis(10), 0.5) == hamiltonian(10, 0.5)

# ╔═╡ 62d63a09-99f9-437d-a21a-97d1e7ba95d4
md"""
## disorder hamiltonian
"""

# ╔═╡ c002c274-4c5f-4ff3-843d-c6fe7dbfd88c
@test hamiltoniandisorderλ(basis(2,1), 0.5) == hamiltonian(basis(2,1), 0.5)

# ╔═╡ b499da64-f4cb-4b30-a466-81779894e948
@test hamiltoniandisorderλ(basis(4,2), 0.5) == hamiltonian(basis(4,2), 0.5)

# ╔═╡ 4159b8c0-d902-4f05-ab42-95187e9cbb96
@test hamiltoniandisorderλ(basis(4,3), 0.5) == hamiltonian(basis(4,3), 0.5)

# ╔═╡ b66e3e83-ed63-4ebf-92e5-6a289195c2e6
@test hamiltoniandisorderλ(basis(4,4), 0.5) == hamiltonian(basis(4,4), 0.5)

# ╔═╡ d0d11d65-0815-4fed-b6d1-c4baf57f99d8
@test hamiltoniandisorderλ(basis(6,3), 0.5) ≈ hamiltonian(basis(6,3), 0.5)

# ╔═╡ d46b5197-056c-4f26-abf9-5ae7cf8720b1
@test hamiltoniandisorderλ(basis(6,3), 0.2) ≈ hamiltonian(basis(6,3), 0.2)

# ╔═╡ 5d76fd04-9d1b-4efd-86c4-c257fb4add1c
@test diag(hamiltoniandisorderλ(basis(4,2), 1)) ≈ zeros(6)

# ╔═╡ ebaf9d60-5b95-4ad8-b073-c81c95370e96
@test !(diag(hamiltoniandisorderλ(basis(4,2), 1, 1, 0)) ≈ zeros(6))

# ╔═╡ 780fa8fe-d24c-443d-ad2a-928003d77069
@test hamiltoniandisorderλ(allbasis(10), 0.5) == hamiltonian(allbasis(10), 0.5)

# ╔═╡ dae7b533-6d3a-4fef-834e-29330c0c3b8d
md"""
## Unconstrained hamiltonian
"""

# ╔═╡ 36298439-0568-42b8-a60c-5f2057189c36
@test hamiltoniandisorder(basis(4,1), 0.5, 0, 0; J=1) == hamiltoniandisorder(basis(4,1), 0.5, 0, 0; J=1, unconstrained=true)

# ╔═╡ c6ed5970-676a-44c8-91a1-391f2f00e1b2
@test hamiltoniandisorder(basis(4,3), 0.5, 0, 0; J=0.5, unconstrained=true) == hamiltoniandisorder(basis(4,1), 0.5, 0, 0; J=0.5)

# ╔═╡ 0755b904-caf3-4a2e-bf9b-54aec0da2fca
md"""
## Representator
"""

# ╔═╡ accf5731-d671-4adc-bd10-0f1a83b468ed
@test isrepresentator(BitArray([1, 0, 1, 0]))

# ╔═╡ 1bce21eb-de21-434d-8173-91891049ee1f
let v = BitArray([1, 0, 1, 0])
	@test only(representator(v).chunks) <= only(v.chunks)
end

# ╔═╡ 292b31e4-e5ad-496c-9892-96175c60c48b
let v = BitArray([0, 1, 0, 1])
	@test only(representator(v).chunks) <= only(v.chunks)
end

# ╔═╡ a7c37145-0c42-46e6-b51b-d457d4783254


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
BitBasis = "50ba71b6-fa0f-514d-ae9a-0916efc90dcf"
CircularArrays = "7a955b69-7140-5f4e-a0ed-f168c5e2e749"
Combinatorics = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[compat]
BitBasis = "~0.7.4"
CircularArrays = "~1.3.0"
Combinatorics = "~1.0.2"
PlutoTest = "~0.2.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.9.1"
manifest_format = "2.0"
project_hash = "6ac185c6691dcffc052ca080bd3ed851d84d800b"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitBasis]]
deps = ["LinearAlgebra", "StaticArrays"]
git-tree-sha1 = "68ce92be119ad7ff44ebbb9ffc0f7a70b1e34c45"
uuid = "50ba71b6-fa0f-514d-ae9a-0916efc90dcf"
version = "0.7.4"

[[deps.CircularArrays]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0598a9ea22c65bfde7f07f21485ebf60deee3302"
uuid = "7a955b69-7140-5f4e-a0ed-f168c5e2e749"
version = "1.3.0"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.2+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.10.11"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "043017e0bdeff61cfbb7afeb558ab29536bbb5ed"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.10.8"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.21+4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.9.0"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.9.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "Pkg", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "5.10.1+6"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "aadb748be58b492045b4f56166b5188aa63ce549"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.7"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╟─10d116fd-e6e9-4ef5-af04-7b5c9a87602d
# ╟─1fe2ea75-b853-4319-845d-82b23a9e3cc4
# ╟─d44936cd-98db-4383-b517-801b390d8305
# ╠═bccd3a32-53f9-11ec-1a7e-c12044b6accc
# ╠═fa07a6a5-8fa7-41f0-88f8-b233f9d60735
# ╠═a7dfe927-0b2e-4734-aa69-867ece13c55e
# ╠═f77ad6eb-8653-4176-90c1-eb102e1d7c18
# ╠═09e81a27-df5b-4661-bd72-18c38d37593e
# ╠═3cca24e9-9d98-4b1b-a017-fe81bbbffdc9
# ╠═7b2dc7c7-8a01-4e4c-ab55-fd30d1a47617
# ╠═b7f5b3ce-a507-481b-8f6c-0d0e57421105
# ╠═49f20b9e-7b35-4414-a15c-e68396513500
# ╠═9ce3260c-79c2-458b-999d-0a5fe1976e00
# ╠═8dee19cd-bbd4-4ab6-ac78-7cffad0b197f
# ╠═b6cc94f2-4639-476d-849e-a6da803b4665
# ╠═a2488d77-f34a-4c7b-b76e-7aaf2adbb99d
# ╟─c9a29f4d-b01c-4c13-8af3-faece6670b37
# ╠═4e2f0e33-b04c-4b0c-bde8-19b3270790de
# ╟─588fff8e-0390-4f1d-96e6-3b57ffc539ba
# ╠═cf77d289-3c3f-44e1-8cae-aadc995c2ace
# ╠═24825184-6350-432b-a2d9-a0b3ae3bd416
# ╠═d2319a30-0554-4749-9ffd-8749de29c55f
# ╠═91dbf769-5033-4a45-903b-0a06503d61b7
# ╠═e5746f41-9ce9-4560-a483-8e0b46e7ef99
# ╠═1f646620-daae-4344-b49a-d2a3c6f6e950
# ╠═8552a30d-4fd7-49a8-adc9-0c53a627c1f3
# ╠═6e182939-51c9-4e69-b89c-4b3b4b67f017
# ╠═5f1582ef-4ba9-4360-8aef-c4fe15b2d035
# ╟─4b9fd419-5322-42e2-8da2-6d6228ae0c2f
# ╠═d8bdf31e-dc62-4e94-bc9d-bdd76f02e914
# ╠═2a7b8063-0555-4986-a755-cd04d2485178
# ╠═ebdfebd3-3ba3-47ac-b058-1aa439e2bcb5
# ╠═e25082e3-595a-440d-9583-b16f2ca513cd
# ╠═9de9dd47-3c12-4c33-b4bb-88e306214b3f
# ╟─03781a73-44bd-4a19-bbf2-4369287b106c
# ╠═04b340b4-92a8-4904-b1a8-334826155001
# ╟─8bdb4cc9-a7ce-47f0-8c8d-1dccf46486f1
# ╠═0b35ddf6-1e34-4f16-9d4d-64ba93a555af
# ╟─de6b23f6-8185-4baf-870c-baec04aa49d3
# ╠═b73e857a-f633-49f4-b16a-6a96530fe67b
# ╟─a6c74b08-5d50-447b-a7e1-7b7c3710389f
# ╠═2d7eb900-2eb9-4848-a2ee-6e2d6734d7ce
# ╠═b745228f-78f9-47eb-a20e-15ef7874d08c
# ╠═4b3add98-cfb1-40ba-994a-e5efdcb0de7b
# ╠═e5ba77f2-e388-4d17-96c5-8dcd9b12b949
# ╠═3aa6c31c-a66a-474e-9af3-02a189217dff
# ╠═e3fa99e2-37ad-4973-8389-f3edf25fb36b
# ╠═ca34fbb5-4001-48d5-9358-aa93aa94ea15
# ╠═1c8455d7-5ea6-4dbc-96ff-0f15e3b462c3
# ╠═29cc7f2d-1105-4499-b777-6b8d0d2cc57a
# ╠═0cd7bf45-4b62-4491-a855-6a9fef016405
# ╠═dd56efac-37a3-44d6-a3eb-b10c0a47857a
# ╠═dd4eb9fa-e2bc-4007-8c4b-1ccb9cd6325a
# ╠═0455bec5-421b-41dd-bd6f-7f80cd3e5029
# ╠═fb903b76-8e50-4954-9dd5-bf03c7b2b103
# ╠═0bb504f7-4b1e-4790-bd82-0f6114d16ea1
# ╠═fc7c52b8-ade1-458e-8159-a39223e698fb
# ╠═3e158d8b-96e2-42eb-ab4a-b49ded9d719c
# ╟─a02ca109-4900-4d9e-bbc4-46d819864411
# ╟─d41dce8a-5bf3-4386-8cc3-c34f2fc3997d
# ╠═d3d8c9c8-c59e-4e20-8f0c-ba82fb0781f3
# ╠═1b0ea891-14c4-41d2-93db-a75521dbcbfe
# ╟─64213987-41a6-4ef4-8573-8d92fa6f0a06
# ╠═461dec51-777a-4437-aff2-db94a18431a8
# ╠═e728da5f-a226-4906-a7e9-224e80a25a79
# ╟─fc754518-ea71-4b67-b6f4-d2a1b46ca593
# ╠═356d6b42-19ec-4f27-9ff8-4ba3d0b0c038
# ╟─3f594112-2e84-4d38-b506-c6e008a88e1b
# ╠═c3beb708-4b3c-4e83-88c4-2bbbf32e5c1a
# ╠═a4e942e9-d41f-44a4-9b9d-4a7063f444b1
# ╠═b989637d-3f35-430a-9674-3ca207648812
# ╟─5da9a9fb-b914-4b50-bf20-ffc868d33e44
# ╠═2f4d7a14-d72f-4ad8-81d7-0fcca6c8183c
# ╠═757d3b67-7f29-4568-b32d-09f466ac74bd
# ╠═f05b2133-75bb-4d0b-84c9-362d4342b824
# ╠═dc86534b-4c10-4303-9c90-32b74596cdb1
# ╟─67c51a46-95a8-42d1-8a1f-635d980ee4b4
# ╠═8bbff18f-2b73-4511-bd09-8a56d99167ad
# ╠═88cdcf03-08dd-4a35-afe2-aae1e8d65efc
# ╠═8a7f59c9-d7d4-4dd2-82a1-7a67aba9f633
# ╟─030c63c0-c106-4122-949a-50394530b7f2
# ╠═dd961753-00b3-4811-bfa1-3f3ad4c3c8c7
# ╠═ff0085b9-faa7-4277-9775-c0faa7c01c15
# ╠═b58af3bd-9e1f-43ed-aab7-095d848145c7
# ╟─42adb6e2-4bbe-4ed7-880b-0085aa9a03c0
# ╠═fe813e92-02d0-4500-a336-5d0e735a3b72
# ╠═38beff73-1c65-4f6c-a6fe-f3e286c96fb1
# ╠═5198905b-1bc5-49f7-9347-66d84d054b03
# ╟─c4d9b1ff-db23-4224-b06c-336d57e223f0
# ╠═dd3d6d18-5238-4310-bab4-a2c07c2a9074
# ╠═0deb2d76-aa3c-4bde-bffb-13192c6561a7
# ╟─1554f951-e241-4486-a887-d77a15a80971
# ╠═cb9c8eda-4731-4d7f-bbab-dbd4b8ea2d1b
# ╠═963fad51-f2a6-4e23-9eed-0b4ce3da436f
# ╠═c1f2efd7-cd4b-4f16-b5ff-caf54ade10c3
# ╠═1d556577-4c83-42ee-9ff0-361ce141a2e1
# ╟─ee2d6f9a-52a2-42e2-8129-c75187cba079
# ╠═1f705aaf-be50-4ce7-ac87-5e2c8640721d
# ╠═e726a652-3898-48c0-bdb6-91f22eb5ad7e
# ╠═d3d21cd9-85ea-4546-844d-6ae794567a2e
# ╠═0e0f74fd-e3df-4b0d-96cf-a928751cb43d
# ╠═0b4e6503-7a1f-458f-b292-b50e37f9dffa
# ╠═f99f8d50-3473-4997-9696-6d1e7b4be18b
# ╟─62d63a09-99f9-437d-a21a-97d1e7ba95d4
# ╠═c002c274-4c5f-4ff3-843d-c6fe7dbfd88c
# ╠═b499da64-f4cb-4b30-a466-81779894e948
# ╠═4159b8c0-d902-4f05-ab42-95187e9cbb96
# ╠═b66e3e83-ed63-4ebf-92e5-6a289195c2e6
# ╠═d0d11d65-0815-4fed-b6d1-c4baf57f99d8
# ╠═d46b5197-056c-4f26-abf9-5ae7cf8720b1
# ╠═5d76fd04-9d1b-4efd-86c4-c257fb4add1c
# ╠═ebaf9d60-5b95-4ad8-b073-c81c95370e96
# ╠═780fa8fe-d24c-443d-ad2a-928003d77069
# ╟─dae7b533-6d3a-4fef-834e-29330c0c3b8d
# ╠═36298439-0568-42b8-a60c-5f2057189c36
# ╠═c6ed5970-676a-44c8-91a1-391f2f00e1b2
# ╟─0755b904-caf3-4a2e-bf9b-54aec0da2fca
# ╠═accf5731-d671-4adc-bd10-0f1a83b468ed
# ╠═1bce21eb-de21-434d-8173-91891049ee1f
# ╠═292b31e4-e5ad-496c-9892-96175c60c48b
# ╠═a7c37145-0c42-46e6-b51b-d457d4783254
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
