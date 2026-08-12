import os
import json
from Bio import PDB

def extract_residues_info(peptide, pdb_directory, output_file):
    # Dicionário para armazenar as informações de saída
    result = {}

    # Parser de estrutura PDB
    parser = PDB.PDBParser(QUIET=True)

    # Lista de resíduos do peptídeo para busca
    peptide_residues = list(peptide)

    # Dicionário de resíduos padrão para checar
    residue_dict = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS',
                    'E': 'GLU', 'Q': 'GLN', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE',
                    'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE', 'P': 'PRO',
                    'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

    # Converter a lista de resíduos do peptídeo para seus códigos de três letras
    peptide_residues = [residue_dict[res] for res in peptide_residues]

    # Percorre todos os arquivos no diretório pdb_directory
    for pdb_file in os.listdir(pdb_directory):
        if pdb_file.endswith(".pdb"):
            print(f"Processando arquivo: {pdb_file}")
            structure = parser.get_structure(pdb_file, os.path.join(pdb_directory, pdb_file))

            # Nome do arquivo sem extensão
            file_name = os.path.splitext(pdb_file)[0]
            # Dicionário para armazenar informações de cadeias e resíduos encontrados para o arquivo atual
            file_info = {}

            for model in structure:
                for chain in model:
                    chain_id = chain.get_id()
                    chain_residues = []

                    for residue in chain:
                        res_name = residue.get_resname()

                        # Checa se o nome do resíduo está na lista de resíduos do peptídeo
                        if res_name in peptide_residues:
                            residue_number = residue.get_id()[1]  # Obtém o número do resíduo
                            chain_residues.append((res_name, residue_number))

                    # Verifica se a sequência de resíduos na cadeia corresponde ao peptídeo
                    for i in range(len(chain_residues) - len(peptide_residues) + 1):
                        if [res[0] for res in chain_residues[i:i+len(peptide_residues)]] == peptide_residues:
                            positions = [res[1] for res in chain_residues[i:i+len(peptide_residues)]]
                            file_info[chain_id] = positions
                            break

            if file_info:
                result[file_name] = file_info

    # Escreve o resultado no arquivo JSON sem espaços desnecessários
    with open(output_file, 'w') as outfile:
        json.dump(result, outfile, separators=(',', ':'))

# Parâmetros de entrada
peptide = "XXXXXXXXXXXXXXXX"  # Exemplo de peptídeo (todos os 20 aminoácidos padrões)
pdb_directory = "../filter_radius_contacts/ns1_denv1/selected_ns1_denv1_rg0.03_ct0.3"  # Diretório contendo os arquivos PDB
output_file = "outputs/ns1_denv1/fixed_positions_ns1_denv1_pdbs.jsonl"  # Nome do arquivo de saída

# Chama a função para extrair as informações
extract_residues_info(peptide, pdb_directory, output_file)
