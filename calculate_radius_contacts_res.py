import os
import argparse
import csv
from Bio.PDB import PDBParser
import numpy as np
from Bio import PDB
import matplotlib.pyplot as plt
import pandas as pd

def calculate_radius_of_gyration_per_residue(structure):
    atom_coordinates = []
    residue_count = 0
    for model in structure:
        for chain in model:
            for residue in chain:
                residue_count += 1
                for atom in residue:
                    atom_coordinates.append(atom.get_coord())
    coordinates_array = np.array(atom_coordinates)
    center_of_mass = np.mean(coordinates_array, axis=0)
    squared_distances = np.sum((coordinates_array - center_of_mass) ** 2, axis=1)
    radius_of_gyration = np.sqrt(np.mean(squared_distances))
    radius_of_gyration_per_residue = radius_of_gyration / residue_count
    return radius_of_gyration_per_residue, residue_count

def contar_residuos_aminoacidos(arquivo_pdb):
    parser = PDB.PDBParser()
    estrutura = parser.get_structure('estrutura', arquivo_pdb)
    contador_residuos = 0
    for modelo in estrutura:
        for cadeia in modelo:
            for residuo in cadeia:
                if PDB.is_aa(residuo):
                    contador_residuos += 1
    return contador_residuos

def calculate_contacts(pdb_file, cutoff_distance=6.0):
    parser = PDBParser()
    structure = parser.get_structure('protein', pdb_file)
    ca_atoms = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.has_id('CA'):
                    ca_atoms.append(residue['CA'])
    ca_array = np.array([atom.get_coord() for atom in ca_atoms])
    distance_matrix = np.linalg.norm(ca_array[:, np.newaxis] - ca_array, axis=-1)
    num_contacts = np.sum(distance_matrix < cutoff_distance) - len(ca_atoms)
    return num_contacts

def process_directory(input_dir, output_file):
    results = []
    with open(output_file, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Filename", "Number of Residues",
                            "Radius of Gyration per Residue (80-84)", "Radius of Gyration per Residue (85-89)",
                            "Radius of Gyration per Residue (90-94)", "Radius of Gyration per Residue (95-100)",
                            "Contacts per Residue"])
        for filename in os.listdir(input_dir):
            if filename.endswith(".pdb"):
                pdb_file = os.path.join(input_dir, filename)
                try:
                    num_residues = contar_residuos_aminoacidos(pdb_file)
                    parser = PDBParser(QUIET=True)
                    structure = parser.get_structure("protein", pdb_file)
                    rg_per_residue, residue_count = calculate_radius_of_gyration_per_residue(structure)
                    num_contacts = calculate_contacts(pdb_file)
                    contacts_per_residue = num_contacts / num_residues if num_residues else 0
                    rg_80_84 = rg_85_89 = rg_90_94 = rg_95_100 = None
                    if 80 <= num_residues <= 84:
                        rg_80_84 = rg_per_residue
                    elif 85 <= num_residues <= 89:
                        rg_85_89 = rg_per_residue
                    elif 90 <= num_residues <= 94:
                        rg_90_94 = rg_per_residue
                    elif 95 <= num_residues <= 100:
                        rg_95_100 = rg_per_residue
                    csvwriter.writerow([filename, num_residues, rg_80_84, rg_85_89, rg_90_94, rg_95_100, contacts_per_residue])
                    print(f'Processado {filename}')
                except Exception as e:
                    print(f'Erro ao processar {filename}: {str(e)}')
    print(f'Os resultados foram salvos em {output_file}')

def plot_scatter_charts(csv_file):
    df = pd.read_csv(csv_file)
    intervals = {
        "80-84": "green",
        "85-89": "yellow",
        "90-94": "orange",
        "95-100": "red"
    }

    for interval, color in intervals.items():
        plt.figure(figsize=(10, 6))
        x = df["Contacts per Residue"]
        if interval == "80-84":
            y = df["Radius of Gyration per Residue (80-84)"]
        elif interval == "85-89":
            y = df["Radius of Gyration per Residue (85-89)"]
        elif interval == "90-94":
            y = df["Radius of Gyration per Residue (90-94)"]
        elif interval == "95-100":
            y = df["Radius of Gyration per Residue (95-100)"]
        mask = ~y.isna()
        x = x[mask]
        y = y[mask]
        count = len(y)
        plt.scatter(x, y, color=color, label=f'{interval} ({count})')
        plt.xlabel("Contacts per Residue")
        plt.ylabel("Radius of Gyration per Residue")
        plt.legend(title="Residue Interval")
        plt.title(f"Scatter Plot of Radius of Gyration per Residue vs Contacts per Residue ({interval})")
        plt.grid(True)
        plt.savefig(f'plot_{interval}_aa_ns1_denv1.png', format='png')
        plt.show()

    # Plot combined chart
    plt.figure(figsize=(10, 6))
    for interval, color in intervals.items():
        x = df["Contacts per Residue"]
        if interval == "80-84":
            y = df["Radius of Gyration per Residue (80-84)"]
        elif interval == "85-89":
            y = df["Radius of Gyration per Residue (85-89)"]
        elif interval == "90-94":
            y = df["Radius of Gyration per Residue (90-94)"]
        elif interval == "95-100":
            y = df["Radius of Gyration per Residue (95-100)"]
        mask = ~y.isna()
        x = x[mask]
        y = y[mask]
        count = len(y)
        plt.scatter(x, y, color=color, label=f'{interval} ({count})')
    plt.xlabel("Contacts per Residue")
    plt.ylabel("Radius of Gyration per Residue")
    plt.legend(title="Residue Interval")
    plt.title("Scatter Plot of Radius of Gyration per Residue vs Contacts per Residue (All Intervals)")
    plt.grid(True)
    plt.savefig('plot_all_intervals_ns1_denv1.png', format='png')
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Análise de arquivos PDB e plotagem de gráficos de dispersão.')
    parser.add_argument('diretorio_pdb', help='Caminho para o diretório contendo arquivos PDB')
    parser.add_argument('csv_file', help='Caminho para o arquivo CSV para plotar os gráficos de dispersão')

    args = parser.parse_args()

    input_directory = args.diretorio_pdb
    output_file = 'radius_contacts_res_ns1_denv1.csv'

    process_directory(input_directory, output_file)
    plot_scatter_charts(output_file)

