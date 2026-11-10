# SPDX-FileCopyrightText: 2023 Christian Willberg <christian.willberg@dlr.de>, Jan-Timo Hesse <jan-timo.hesse@dlr.de>
#
# SPDX-License-Identifier: BSD-3-Clause

# Damage only in the right quarter of the strip: on 2 ranks one rank holds no damage
# block, but all ranks must still exchange the damage model's fields.
folder_name = basename(@__FILE__)[1:(end - 3)]
cd("fullscale_tests/" * folder_name) do
    run_perilab("partial_damage", 2, true, folder_name)
end
