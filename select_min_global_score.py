import os

def parse_fasta_file(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        current_sequence = None
        for line in file:
            if line.startswith('>'):
                if current_sequence:
                    sequences.append(current_sequence)
                current_sequence = {'header': line.strip(), 'sequence': ''}
            else:
                current_sequence['sequence'] += line.strip()
        if current_sequence:
            sequences.append(current_sequence)
    return sequences

def write_fasta_file(output_file, header, sequence, input_file_name):
    with open(output_file, 'w') as file:
        modified_header = f"{header} | source_file={input_file_name}"
        file.write(modified_header + '\n')
        file.write(sequence + '\n')

def main(input_directory_path, output_directory_path):
    # Verifica se o diretório de saída existe, caso contrário, cria
    if not os.path.exists(output_directory_path):
        os.makedirs(output_directory_path)

    for file_name in os.listdir(input_directory_path):
        if file_name.endswith('.fa'):
            file_path = os.path.join(input_directory_path, file_name)
            sequences = parse_fasta_file(file_path)
            if sequences:
                sequences = sequences[1:]  # Ignorando a primeira sequência
                min_global_score_seq = min(sequences, key=lambda x: float(x['header'].split(", global_score=")[1].split(",")[0]))                                                                                                                                                               
                output_file = os.path.join(output_directory_path, os.path.splitext(file_name)[0] + '_selected.fasta')
                write_fasta_file(output_file, min_global_score_seq['header'], min_global_score_seq['sequence'], file_name)

# Exemplo de uso:
input_directory_path = '../protein_mpnn/outputs/ns1_denv3/seqs'
output_directory_path = 'ns1_denv3'
main(input_directory_path, output_directory_path)

