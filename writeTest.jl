using Exodus

# Base.Filesystem.rm("./temp_element_variables.e")

exo = ExodusDatabase(
    "./temp_element_variables.e",
    Initialization(3, 4, 2, 1, 0, 0)
)

nodal_names = ["displ_x", "displ_y", "displ_z"]
element_names = ["stress_xx", "stress_yy", "stress_xy"]

coords = [0 0 1 1; 0 1 0 1; 0 0 0 0]
blocks = [Block(1, 4, 1, "point", [1 2 3; 2 4 3])]

write_coordinates(exo, coords)
write_element_blocks(exo, blocks)
write_time(exo, 1, 0.0)
write_number_of_nodal_variables(exo, length(nodal_names))
write_number_of_element_variables(exo, length(element_names))

# for (i, name) in enumerate(nodal_names)
#     write_nodal_variable_name(exo, i, name)
#     write_nodal_variable_values(exo, 1, i, [1, 2, 3, 4])
# end

# for (i, name) in enumerate(element_names)
#     write_element_variable_name(exo, i, name)
#     write_element_variable_values(exo, 1, 1, i, [5, 6, 7, 8])
# end


close(exo)