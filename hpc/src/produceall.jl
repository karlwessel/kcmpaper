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
seeds = "seed" => collect(1:20)
configsdefinftime = dict_list(Dict("L" => [12, 16], 
    "V" => 4,
    "W" => collect(10 .^ range(-8, 3, length=200)), seeds,
    "unconstrained" => [false, true],
    "data" => "diagensemble",
    "purpose" => "inftimevals",
    always...))

println("Scaling configs...")
common = ("V" => 4,
    "data" => "eigendensities",
    "purpose" => "scaling",
    always...)

seeds = "seed" => collect(1:400)
disorder = "W" => collect(10 .^ range(-2, 2, length=17)) * 2
configsdefscaling = dict_list(Dict("L" => 12, 
    disorder, seeds,
    "unconstrained" => false,
    common...))
    
disorder = "W" => unique([collect(10 .^ range(-2, 2, length=9)); collect(10 .^ range(0, 2, length=17))]) * 2
configsdefscaling = vcat(configsdefscaling, dict_list(Dict("L" => 12,      
    disorder, seeds,
    "unconstrained" => true,
    common...)))
    
seeds = "seed" => collect(1:4000)
disorder = "W" => collect(10 .^ range(0.5, 2, length=7)) * 2
configsdefscaling = vcat(configsdefscaling, dict_list(Dict("L" => 12, 
    disorder, seeds,
    "unconstrained" => false,
    common...)))

seeds = "seed" => collect(1:30)
disorder = "W" => collect(10 .^ range(-2, 2, length=17)) * 2
configsdefscaling = vcat(configsdefscaling, dict_list(Dict("L" => 16, 
    disorder, seeds,
    "unconstrained" => true,
    common...)))
    
disorder = "W" => unique([collect(10 .^ range(-2, 2, length=9)); collect(10 .^ range(0, 2, length=17))]) * 2
configsdefscaling = vcat(configsdefscaling, dict_list(Dict("L" => 16,      
    disorder, seeds,
    "unconstrained" => false,
    common...)))
    
seeds = "seed" => collect(1:2000)
disorder = "W" => collect(10 .^ range(0.5, 2, length=7)) * 2
configsdefscaling = vcat(configsdefscaling, dict_list(Dict("L" => 16,
    disorder, seeds,
    "unconstrained" => false,
    common...)))
    
println("WV map configs...")
seeds = "seed" => collect(1:20)
configsdefwvmap = dict_list(Dict("L" => [16], 
    "V" => collect(10 .^ range(-1, 2, length=13)),
    "W" => collect(10 .^ range(-1, 2, length=13)) * 2, seeds,
    "unconstrained" => [false, true],
    "data" => "eigendensities",
    "purpose" => "wvmap",
    always...))

println("Paper figure configs...")
common = ("data" => "autocorrelationevolutions", 
    "purpose" => "paperfigures",
    "ts" => "detail", "unconstrained" => [false, true],
    always...)
configsdefpaperfigures = dict_list(Dict("L" => [12], 
    "V" => [1, 4, 16], "W" => [0, 0.5, 10] * 2, "seed" => collect(1:20), common...))
	
configsdefpaperfigures = vcat(configsdefpaperfigures, dict_list(Dict("L" => [16], "V" => [4], "W" => 0, "seed" => 1, common...)),
    dict_list(Dict("L" => [16], "V" => [4], "W" => [0.5] * 2, "seed" => collect(1:10), common...)))


println("meta stable configs...")
common = ("purpose" => "metastable",
    "data" => "autocorrelationevolutions",
    "V" => 4,
    "ts" => "meta",
    always...)

configsmeta = dict_list(Dict("L" => [12], 
    "W" => unique([collect(10. .^ (-2:0.5:2));
        collect(10. .^ (-1.5:0.25:1))]) * 2, 
    "seed" => collect(1:20), "unconstrained" => [false, true],
    common...))

configsmeta = vcat(configsmeta,
    dict_list(Dict("L" => [16], 
    "W" => unique([collect(10. .^ (-2:0.5:2));
        collect(10. .^ (-1.5:0.25:0.5))]) * 2, 
    "seed" => collect(1:5), 
    "unconstrained" => [false],
    common...)), 
    dict_list(Dict("L" => [16], 
    "W" => unique([collect(10. .^ (-2:0.5:2));
        collect(10. .^ (-0:0.25:1.5))]) * 2, 
    "seed" => collect(1:5), 
    "unconstrained" => [true],
    common...)))


println("Entropies")
common = ("purpose" => "logincrease",
    "data" => "entropyevolutions",
    "V" => 1, "W" => [10, 20, 2, 0.5], 
    "ts" => "rough",
    always...)
configsloginc = vcat(
	dict_list(Dict("state" => "1110_1110_1110", 
	    "L" => 12, "seed" => collect(1:1200), 
        "unconstrained" => [true],
	    common...)),
	dict_list(Dict("state" => "1110_1110_1110_1110", 
	    "L" => 16, 
	    "seed" => collect(1:800), 
        "unconstrained" => [true], 
	    common...)))

configsloginc = vcat(
    dict_list(Dict("state" => "1101_1010_1111", 
        "L" => 12, "seed" => collect(1:1200), 
        "unconstrained" => [true, false],
        common...)),
    dict_list(Dict("state" => "1101_1011_0101_1111", 
        "L" => 16, 
        "seed" => collect(1:800), 
        "unconstrained" => [true, false], 
        common...)))

common = ("purpose" => "logincrease",
    "data" => "entropyevolutions",
    "V" => [0, 0.5, 1, 2, 4], "W" => [10, 20, 2, 0.5], 
    "ts" => "rough",
    always...)
configsloginc = vcat(configsloginc,
	dict_list(Dict("state" => "1110_1110_1110_1110", 
	    "L" => 16, "unconstrained" => [true], 
	    "seed" => collect(1:200), 
	    common...)))
configsloginc = vcat(configsloginc,
    dict_list(Dict("state" => "1101_1011_0101_1111", 
        "L" => 16, "unconstrained" => [true, false], 
        "seed" => collect(1:200), 
        common...)))

common = ("purpose" => "eigenentropydistribution",
    "data" => "eigenentropies",
    "V" => 1, "W" => [10, 20, 2, 0.5], 
    "unconstrained" => [true, false],
    always...)
configseigentr = vcat(
	dict_list(Dict(
	    "L" => 12, "seed" => collect(1:100),
	    common...)),
	dict_list(Dict(
	    "L" => 16, 
	    "seed" => collect(1:30), 
	    common...)))
common = ("purpose" => "eigenentropydistribution",
    "data" => "eigenentropies",
    "V" => [0, 0.5, 1, 2, 4], "W" => [10, 20, 2, 0.5], 
    "unconstrained" => [true, false],
    always...)
configseigentr = vcat(configseigentr,
	dict_list(Dict(
	    "L" => 16, 
	    "seed" => collect(1:20), 
	    common...)))


"""
Execute Configs
"""

configsdef = unique(vcat(configsdefscaling,
    configsdefwvmap,
    configsdefpaperfigures,
    configsmeta,
    configsloginc,
    configseigentr,
    configsdefinftime))

	
configsdef = sort(configsdef, by = (x -> (x["seed"] - 1) * binomial(x["L"], Int(x["L"] * x["filling"])) / x["L"]))

@everywhere function calc(c)
    p = c["purpose"]
    delete!(c, "purpose")
    produceorload(c; subpath = joinpath("newdisorder2", p))
    nothing
end

println("Executing $(length(configsdef)) configs:")

pmap(calc, configsdef)




