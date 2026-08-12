import __main__
__main__.pymol_argv = ['pymol','-qc']

import pymol
from pymol import cmd, stored

pymol.finish_launching()

# Configuração do cálculo de SASA
cmd.set('dot_solvent', 1)
cmd.set('dot_density', 3)

# Carrega a estrutura
cmd.load('7WUT.cif')

# Lista para armazenar cadeia, número e nome do resíduo
stored.residues = []
cmd.iterate('name CA', 'stored.residues.append((chain, resi, resn))')

# Lista para armazenar resultados de SASA
sasa_results = []

# Calcula SASA por resíduo
for chain, residue, resname in stored.residues:
    sasa_value = cmd.get_area(f'chain {chain} and resi {residue}')
    sasa_results.append((chain, residue, resname, sasa_value))

# Arquivo de saída
output_file = "sasa_area_7wut.txt"

with open(output_file, "w") as f:
    
    f.write("Chain\tResidue\tNumber\tSASA (Å^2)\n")
    
    for chain, residue, resname, sasa_value in sasa_results:
        f.write(f"{chain}\t{resname}\t{residue}\t{sasa_value:.2f}\n")

    total_sasa = sum(sasa_value for _, _, _, sasa_value in sasa_results)
    f.write(f"\nTotal SASA\t{total_sasa:.2f} Å^2\n")

print(f"Resultados salvos em '{output_file}'")
