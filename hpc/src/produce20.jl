using Distributed
using SlurmClusterManager

try
    addprocs(SlurmManager(), exeflags="--project=$(Base.active_project())")
catch e
    @warn "Could not create SlurmClusterManager, using default cluster manager. Error was: $e"
end

println("Distribute code...")
@everywhere include("produceresults.jl")

always = ("filling" => 3/4, "J" => 1)

println("inftime vals...")
seeds = "seed" => collect(1:5)
configsdefinftime = dict_list(Dict("L" => [20], 
    "V" => 4,
    "W" => collect(10 .^ range(-8, 3, length=100)), seeds,
    "unconstrained" => [false, true],
    "data" => "diagensemble",
    "purpose" => "inftimevals",
    always...))

println("Scaling configs...")
common = ("L" => 20, "V" => 4,
    "data" => "eigendensities",
    "purpose" => "scaling",
    "seed" => collect(1:20),
    always...)

disorder = 
configsscaling = vcat(
    dict_list(Dict(
        "W" => collect(10 .^ range(-2, 2, length=17)) * 2, 
        "unconstrained" => [true, false],
        common...)), 
	dict_list(Dict(
	    "W" => collect(10 .^ (0.5:0.25:1.5)) * 2,
	    "unconstrained" => false,
	    common...)))


println("Paper figure configs...")
common = ("data" => "autocorrelationevolutions", 
    "purpose" => "paperfigures",
    "ts" => "detail", "unconstrained" => [false, true],
    "L" => 20,
    "V" => 4,
    always...)

configsdefpaperfigures = vcat(
    dict_list(Dict("W" => 0 * 2, "seed" => 1, common...)),
    dict_list(Dict("W" => 0.5 * 2, "seed" => collect(1:10), common...)))


println("meta stable configs...")
common = ("purpose" => "metastable",
    "data" => "autocorrelationevolutions",
    "V" => 4,
    "ts" => "meta",
    "L" => 20,
    "seed" => collect(1:5),
    always...)

configsmeta = vcat(
	dict_list(Dict(
	    "W" => collect(10. .^ (-2:0.5:2)) * 2, 
	    "seed" => collect(1:5),
	    "unconstrained" => [false, true],
	    common...)), 
	dict_list(Dict(
	    "W" => collect(10. .^ (-1.5:0.25:0.5)) * 2,
	    "unconstrained" => false,
	    common...)), 
	dict_list(Dict(
	    "W" => collect(10. .^ (0.5:0.25:1.5)) * 2,
	    "unconstrained" => true,
	    common...)))


println("Entropies ...")
common = ("purpose" => "logincrease",
    "data" => "entropyevolutions",
    "V" => 1, "W" => [20, 2, 10, 0.5], 
    "ts" => "rough",
    always...)
configsloginc = vcat(dict_list(Dict("state" => "1110_1110_1110_1110_1110", 
	    "L" => 20, "unconstrained" => [true], 
	    "seed" => collect(1:80),
	    common...)),
        dict_list(Dict("state" => "1101_1011_0110_1011_1111", 
	    "L" => 20, "unconstrained" => [true, false], 
	    "seed" => collect(1:80),
	    common...)))

common = ("purpose" => "eigenentropydistribution",
    "data" => "eigenentropies",
    "V" => 1, "W" => [20, 2, 10, 0.5], 
    "unconstrained" => [true, false],
    always...)
configseigentr = dict_list(Dict(
	    "L" => 20, 
	    "seed" => collect(1:5),
	    common...))


"""
Execute Configs
"""

configsdef = unique(vcat(configsscaling,
    configsdefpaperfigures,
    configsmeta,
    configsloginc,
    configseigentr,
    configsdefinftime
))
	
configsdef = sort(configsdef, by = (x -> (x["seed"] - 1) * binomial(x["L"], Int(x["L"] * x["filling"])) / x["L"]))

@everywhere function calc(c)
    p = c["purpose"]
    delete!(c, "purpose")
    produceorload(c; subpath = joinpath("newdisorder2", p))
    nothing
end

println("Executing $(length(configsdef)) configs:")

pmap(calc, configsdef)




