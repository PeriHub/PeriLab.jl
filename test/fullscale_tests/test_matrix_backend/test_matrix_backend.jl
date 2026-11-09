# SPDX-FileCopyrightText: 2026 Christian Willberg
#
# SPDX-License-Identifier: BSD-3-Clause

# Matrix Backend CPU (PD_Matrix, licensed matrix elements) on 1–3 ranks against the
# SparseArrays reference of test_model_reduction_Verlet_matrix.
folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    for cores in 1:3
        run_perilab("verlet_full_matrix_cpu", cores, true, folder_name)
    end
    # blocks with different densities: per-node mass (Reference from the CPU backend)
    run_perilab("two_densities", 1, true, folder_name)
    run_perilab("two_densities", 2, true, folder_name)
    # DCB crack growth (Critical Energy): the matrix follows the bond damage; Reference from
    # 1 rank (damage up to 1 along the crack)
    for cores in 1:2
        run_perilab("dcb_damage", cores, true, folder_name)
    end
    # printing (Simple additive model): status mask, the matrix follows the Active field;
    # inactive nodes keep zero displacement
    PeriLab.run("additive_verlet.yaml"; silent = true)
    exo = ExodusDatabase("additive_verlet.e", "r")
    for step in eachindex(read_times(exo))
        inactive = read_values(exo, NodalVariable, step, 1, "Active") .< 0.5
        for name in ("Displacementsx", "Displacementsy")
            @test all(iszero, read_values(exo, NodalVariable, step, 1, name)[inactive])
        end
    end
    close(exo)
    rm("additive_verlet.e")
    for cores in 1:2
        run_perilab("additive_verlet", cores, true, folder_name)
    end
    # Linear Static and Newmark with BiCGStab against the LU references
    for cores in 1:2
        run_perilab("linear_static_bicgstab", cores, true, folder_name * "_implicit")
        run_perilab("newmark_bicgstab", cores, true, folder_name * "_implicit")
    end
    # one iteration does not converge: the first step fails (StepFailed, no fallback), the
    # run is cancelled before any step is written
    PeriLab.run("newmark_bicgstab_itmax.yaml"; silent = true)
    exo = ExodusDatabase("newmark_bicgstab_itmax.e", "r")
    @test length(read_times(exo)) <= 1
    close(exo)
    rm.(filter(f -> startswith(f, "newmark_bicgstab_itmax") &&
                    (endswith(f, ".e") || endswith(f, ".log")), readdir()))
    # GPU and Mixed (rank 0 on the GPU, rank 1 on the CPU): only with CUDA.jl installed
    # and a functional GPU (PeriLab loads CUDA.jl itself)
    PD = PeriLab.Solver_Manager.PD_Matrix
    PD.Backend.load_cuda!()
    hook = PD.GPU_HOOK[]
    if hook !== nothing && Base.invokelatest(hook.n_devices) > 0
        run_perilab("verlet_full_matrix_gpu", 1, true, folder_name)
        run_perilab("verlet_full_matrix_mixed", 2, true, folder_name)
    else
        @test_skip "Matrix Backend GPU / Mixed (needs CUDA.jl and a functional GPU)"
    end
end
