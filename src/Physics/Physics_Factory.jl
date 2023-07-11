include("./Additive/Additive_Factory.jl")
include("./Damage/Damage_Factory.jl")
include("./Material/Material_Factory.jl")
include("./Thermal/Thermal_Factory.jl")
import .Additive
import .Damage
import .Material
import .Thermal
export get_physics

module Physics
function get_physics(params)
    Additive.get_additive(params)
    Damage.get_damage(params)
    Material.get_material(params)
    Thermal.get_thermal(params)
    return
end

end