function thermal_strain(alpha, temperature, strain)
    strain -= alpha * temperature
    return strain
end
