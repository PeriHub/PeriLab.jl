function get_forces_from_force_density(datamanager)
    force_density = datamanager.get_field("Force Densities", "NP1")
    forces = datamanager.get_field("Forces", "NP1")
    volume = datamanager.get_field("Volume")
    forces[:] = force_density .* volume
    return data
end