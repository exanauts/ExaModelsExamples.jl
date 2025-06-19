module ExaModelsExamples

import JLD2
import Downloads
import ExaModels: ExaModels, NLPModels
import PowerModels: PowerModels, silence
import Random

include("OPF/opf.jl")
include("OPF/scopf.jl")
include("OPF/mpopf.jl")
include("Misc/luksanvlcek.jl")
include("OptimalControl/distillation.jl")
include("OptimalControl/quadrotor.jl")
include("OptimalControl/goddard.jl")
include("COPS/robot.jl")
include("COPS/rocket.jl")
include("COPS/bearing.jl")
include("COPS/camshape.jl")
include("COPS/elec.jl")
include("COPS/steering.jl")
include("COPS/pinene.jl")
include("COPS/marine.jl")
include("COPS/gasoil.jl")
include("COPS/pde_models.jl")

const NAMES = filter(names(ExaModelsExamples; all = true)) do x
    str = string(x)
    endswith(str, "model") && !startswith(str, "#")
end

for name in filter(names(ExaModelsExamples; all = true)) do x
    endswith(string(x), "model")
end
    @eval export $name
end

function __init__()
    if haskey(ENV, "EXA_MODELS_DEPOT")
        global TMPDIR = ENV["EXA_MODELS_DEPOT"]
    else
        global TMPDIR = joinpath(@__DIR__,"..","data")
        mkpath(TMPDIR)
    end
    PowerModels.silence()
end

end # module ExaModelsExamples
