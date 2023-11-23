# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

file_content = "header: x y block_id volume\n"
file_content_ns = "header: global_id\n"

# Generate x and y values using for loops
num_points = 100  # Adjust the number of points as needed
step_size = 0.1

data = []
ns1 = []
id = 1
for i in range(num_points):
    for j in range(num_points):
        block = 1
        x = i * step_size
        y = j * step_size
        if x >= (num_points * step_size /2) :
            ns1.append(f"{id}")
            block = 2
        data.append(f"{x:.1f} {y:.1f} {block} {step_size*step_size:.2f}")
        id+=1

# Combine header and data
file_content += "\n".join(data)
file_content_ns += "\n".join(ns1)

directory = "examples/Benchmark/"

with open(directory + "mesh" + str(num_points) +".txt", "w") as file:
    file.write(file_content)

with open(directory +"ns"+ str(num_points) + ".txt", "w") as file:
    file.write(file_content_ns)
