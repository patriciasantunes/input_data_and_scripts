import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os
import shutil

# Função para ler o CSV, filtrar os dados, plotar o gráfico de dispersão e copiar arquivos
def plot_scatter_chart(csv_file, input_dir):
    # Ler o arquivo CSV
    df = pd.read_csv(csv_file)
    
    # Definir as cores para os intervalos de resíduos
    colors = {
        "80-84": "green",
        "85-89": "yellow",
        "90-94": "orange",
        "95-100": "red"
    }

    # Inicializar uma lista para armazenar todos os valores mínimos de raio de giração
    min_radius_values = []

    # Identificar o menor raio de giração considerando todos os intervalos
    for interval in colors.keys():
        col_name = f"Radius of Gyration per Residue ({interval})"
        min_radius_values.append(df[col_name].min())

    # Menor valor de raio de giração entre todos os intervalos
    global_min_radius_value = min(min_radius_values)

    # Identificar o maior valor de contatos por resíduo no DataFrame original
    global_max_contacts_value = df["Contacts per Residue"].max()

    # Inicializar DataFrame vazio para armazenar os resultados filtrados
    filtered_dfs = []

    # Filtrar os dados para cada grupo de raio de giração
    for interval in colors.keys():
        col_name = f"Radius of Gyration per Residue ({interval})"
        
        # Filtrar os dados de acordo com os critérios especificados
        filtered_df = df[
            (abs(df[col_name] - global_min_radius_value) <= 0.03)
        ]

        # Adicionar os dados filtrados ao DataFrame de resultados
        filtered_dfs.append(filtered_df)

    # Combinar todos os DataFrames filtrados
    filtered_df = pd.concat(filtered_dfs).drop_duplicates().reset_index(drop=True)

    # Filtrar os dados com base no valor global de contatos por resíduo
    final_filtered_df = filtered_df[
        abs(filtered_df["Contacts per Residue"] - global_max_contacts_value) <= 0.3
    ]

    # Definir o nome do arquivo de saída
    output_csv_file = "selected_ns1_denv1_rg0.03_ct0.3.csv"

    # Salvar os dados filtrados em um novo arquivo CSV
    final_filtered_df.to_csv(output_csv_file, index=False)

    # Copiar arquivos da coluna Filename para o diretório de saída
    output_dir = "selected_ns1_denv1_rg0.03_ct0.3"
    os.makedirs(output_dir, exist_ok=True)

    for filename in final_filtered_df["Filename"]:
        src_path = os.path.join(input_dir, filename)
        dst_path = os.path.join(output_dir, filename)
        if os.path.exists(src_path):
            shutil.copy(src_path, dst_path)
        else:
            print(f"File {filename} not found in the input directory.")

    # Criar o gráfico de dispersão
    plt.figure(figsize=(10, 6))

    # Plotar os dados para cada intervalo de resíduos
    for interval, color in colors.items():
        x = final_filtered_df["Contacts per Residue"]
        y = final_filtered_df[f"Radius of Gyration per Residue ({interval})"]

        # Remover valores nulos para não causar erros no gráfico
        mask = ~y.isna()
        x = x[mask]
        y = y[mask]

        # Plotar os dados com a cor correspondente
        plt.scatter(x, y, color=color, label=interval)

    # Configurações do gráfico
    plt.xlabel("Contacts per Residue")
    plt.ylabel("Radius of Gyration per Residue")
    plt.legend(title="Residue Interval")
    plt.title("Scatter Plot of Radius of Gyration per Residue vs Contacts per Residue")
    plt.grid(True)

    # Salvar o gráfico em um formato compatível com o Google Sheets
    plt.savefig('plot_selected_ns1_denv1_rg0.03_ct0.3.png', format='png')

    # Exibir o gráfico
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot a scatter chart from a CSV file and copy selected files.')
    parser.add_argument('csv_file', help='Path to the CSV file')
    parser.add_argument('input_dir', help='Path to the input directory containing PDB files')

    args = parser.parse_args()

    plot_scatter_chart(args.csv_file, args.input_dir)

