from thermopack.cubic import cubic

# Define the cubic EOS
eos = cubic('CO2,C1', 'SRK')
#eos = cubic('CO2', 'VDW')
x = [0.5,0.5]
#x = [0.05,0.95]
#x = [0.95,0.05]

# Calculate the isotherm
p_iso_T, v_iso_T, s_iso_T, h_iso_T = eos.get_isotherm(
    350, x, minimum_pressure=1e5, maximum_pressure=1.5e7, nmax=100)

# Print header
print(f"{'Volume(m^3/mol)'}{' Pressure(Pa)'} ")
for p, v in zip(p_iso_T, v_iso_T):
    print(f"{v:.10f} {p:.10f}")

