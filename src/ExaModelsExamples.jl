module ExaModelsExamples

import JLD2
import Downloads
import ExaModels: ExaModels, NLPModels
import PowerModels: PowerModels, silence
import Random

include("opf.jl")
include("scopf.jl")
include("mpopf.jl")
include("luksanvlcek.jl")
include("distillation.jl")
include("quadrotor.jl")
include("goddard.jl")
include("robot.jl")
include("rocket.jl")
include("bearing.jl")
include("camshape.jl")
include("elec.jl")
include("steering.jl")
include("pinene.jl")
include("marine.jl")
include("gasoil.jl")
include("pde_models.jl")
include("carmix.jl")
include("chain.jl")
include("channel.jl")
include("glider.jl")
include("minsurf.jl")
include("polygon.jl")
include("tetra.jl")
include("torison.jl")
include("triangle.jl")
include("methanol.jl")

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
