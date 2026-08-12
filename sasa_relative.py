from pymol import cmd

cmd.load("7WUT.cif")

# Calcula SASA relativa por resíduo (valores entre 0 e 1)
sasa_data = cmd.get_sasa_relative("polymer")

# Coleta informações dos resíduos
residue_info = []

cmd.iterate(
    "polymer and name CA",
    'residue_info.append((chain, resn, resi))',
    space={"residue_info": residue_info}
)

output_file = "sasa_relative_7wut.txt"

with open(output_file, "w") as f:
    # Cabeçalho corrigido
    f.write("Chain\tResidue\tNumber\tRelative SASA (%)\n")

    # Assume que a ordem de sasa_data acompanha a ordem dos resíduos selecionados
    for (residue_num, sasa_value), (chain, resn, resi) in zip(sasa_data.items(), residue_info):
        f.write(f"{chain}\t{resn}\t{resi}\t{sasa_value * 100:.2f}\n")

print(f"Resultados salvos em {output_file}")
