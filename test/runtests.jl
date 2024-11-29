using Test, MadNLP, MadNLPGPU, KernelAbstractions, CUDA, MadNLPHSL
import ExaModelsExamples
const CONFIGS = [
     (Float64, nothing),
    # (Float64, CPU()),
    (Float64, CUDABackend()),
]

# Note: The results may vary depending on the solver used.
# Additionally, by using MadNLP default solver, if the tolerance (tol) is set to less than 1e-4,
# the solver may diverge. Adjust the tolerance accordingly

ExaModels_INSTANCES = [
    (ExaModelsExamples.bearing_model, (50, 50), -1.5482e-1),
    (ExaModelsExamples.chain_model, (400,), 5.0698),
    (ExaModelsExamples.camshape_model, (1000,), 4.27907),
    (ExaModelsExamples.catmix_model, (100,), -4.80556e-2),
    (ExaModelsExamples.channel_model, (200,), 1.0),
    (ExaModelsExamples.elec_model, (50,), 1.0552e3),
    (ExaModelsExamples.gasoil_model, (100,), 5.2366e-3),
    (ExaModelsExamples.glider_model, (100,), 1.25505e3),
    (ExaModelsExamples.goddard_model, (400,), -1.01283),
    (ExaModelsExamples.marine_model, (100,), 1.97462e7),
    (ExaModelsExamples.methanol_model, (100,), 9.02229e-3),
    (ExaModelsExamples.minsurf_model, (50, 50), 2.51488),
    (ExaModelsExamples.minsurf_model, (50, 75), 2.50568),
    (ExaModelsExamples.minsurf_model, (50, 100), 2.50694),
    (ExaModelsExamples.pinene_model, (100,), 1.98721e1),
    (ExaModelsExamples.polygon_model, (100,), -0.674981), # N.B: objective depends on the optimizer used.
    (ExaModelsExamples.robot_model, (200,), 9.14138),
    (ExaModelsExamples.steering_model, (200,), 5.54577e-1),
    (ExaModelsExamples.tetra_duct20_model, (), 4.82685e3),
    (ExaModelsExamples.tetra_duct15_model, (), 1.04951e4),
    (ExaModelsExamples.tetra_foam5_model, (), 6.42560e3),
    (ExaModelsExamples.tetra_gear_model, (), 4.15163e3),
    (ExaModelsExamples.tetra_hook_model, (), 6.05735e3),
    (ExaModelsExamples.torsion_model, (50, 50), -4.18087e-1),
    (ExaModelsExamples.dirichlet_model, (20,), 1.71464e-2),
    (ExaModelsExamples.henon_model, (20,), 6.667736), # N.B: objective depends on the optimizer used.
    (ExaModelsExamples.lane_emden_model, (20,), 9.11000),
    (ExaModelsExamples.triangle_deer_model, (), 2.01174e3),
    (ExaModelsExamples.triangle_pacman_model, (), 1.25045e3),
    (ExaModelsExamples.triangle_turtle_model, (), 4.21523e3),
]

function runtests()
    @testset "ExaModelsExamples test" begin
        for (instance, params, result) in ExaModels_INSTANCES
            for (T, backend) in CONFIGS
                m = eval(instance)(params...; T = T, backend = backend)
                result = madnlp(m;)

                @testset "$instance" begin
                    @test result.status == MadNLP.SOLVE_SUCCEEDED
                end
            end
        end
    end
end

runtests()
