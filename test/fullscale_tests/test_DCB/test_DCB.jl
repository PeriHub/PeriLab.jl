# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    run_perilab("DCBmodel_unified_bb", 1, true, folder_name)
    run_perilab("DCBmodel_correspondence", 1, true, folder_name)
    run_perilab("DCBmodel_UMAT", 1, true, folder_name)
    run_perilab("DCBmodel_PD_solid", 1, true, folder_name)
    run_perilab("DCBmodel_PD_solid_fallback", 1, true, folder_name)
    run_perilab("DCBmodel_PD_solid_local_damping", 1, true, folder_name)
    # run_perilab("DCBmodel_PD_solid_switch", 1, true, folder_name)
    run_perilab("DCBmodel_PD_solid_temp_depended", 1, true, folder_name)
    run_perilab("DCBmodel_correspondence_bond_associated", 1, true, folder_name)
    run_perilab("DCBmodel_correspondence_surface_correction", 1, true, folder_name)
end

# pre-cracked DCB (bond filter with contact: damaged from t = 0), arms released with an
# initial velocity, no loads: local damping must take kinetic energy out
function dcb_kinetic_energy(coefficient)
    deck = "DCBmodel_PD_solid_local_damping_precrack"
    file = "$(deck)_c$(coefficient).yaml"
    write(file,
          replace(read("$deck.yaml", String),
                  "Damping coefficient: 0.45" => "Damping coefficient: $coefficient",
                  "Output Filename: $deck" => "Output Filename: $(deck)_c$(coefficient)"))
    PeriLab.run(file; silent = true)
    rm(file)
    DM = PeriLab.Data_Manager
    v = DM.get_field("Velocity", "N")
    return sum(DM.get_field("Density") .* DM.get_field("Volume") .* sum(v .^ 2; dims = 2)) /
           2
end

cd("fullscale_tests/" * folder_name) do
    undamped = dcb_kinetic_energy(0.0)
    damped = dcb_kinetic_energy(0.45)
    @test isfinite(damped) && isfinite(undamped)
    @test damped < 0.6 * undamped
end
