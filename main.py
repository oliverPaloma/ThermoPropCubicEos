from thermopack.cubic import cubic

# Define the cubic EOS

#eos = cubic('CO2,C1', 'pr')         
eos = cubic('C1,C2,C3,C4,NC4,IC5,NC5,NC6,NC7,NC8,NC9,NC10,NC11,NC12,NC13,NC14,NC15,NC16', 'pr')

#eos = cubic('CO2', 'VDW')
#x = [0.5,0.5]
#x = [0.413506,0.040300,0.215300,0.053900,0.054300,0.051500,0.051900,0.1039,0.00347,0.00268,0.00207,0.00159,0.00123,0.00095,0.00073,0.000566,0.000437,0.001671]

x = [0.413506, 0.040300, 0.215300, 0.053900, 0.054300, 0.051500, 
     0.051900, 0.1039, 0.00347, 0.00268, 0.00207, 0.00159, 0.00123, 
     0.00095, 0.00073, 0.000566, 0.000437, 0.001671]



#x = [0.95,0.05]

# Calculate the isotherm
p_iso_T, v_iso_T, s_iso_T, h_iso_T = eos.get_isotherm(
    320, x, minimum_pressure=1e5, maximum_pressure=1.5e7, nmax=100)

# Print header
print(f"{'Volume(m^3/mol)'}{' Pressure(Pa)'} ")
for p, v in zip(p_iso_T, v_iso_T):
    print(f"{v:.10f} {p:.10f}")

