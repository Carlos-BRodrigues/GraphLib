#include "graph_tools2.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <chrono>
#include <random>
#include <bits/stdc++.h>
using namespace std;

void analisar_grafo(const std::string& nome_arquivo) {
    std::cout << "\n==============================================" << std::endl;
    std::cout << "Analisando o grafo: " << nome_arquivo << std::endl;
    std::cout << "==============================================" << std::endl;

    graph_tools_lib::Graph graph(nome_arquivo, graph_tools_lib::RepresentationType::ADJACENCY_LIST);
    int num_vertices = graph.getVertexCount();
    std::string output_filename = "analise_" + nome_arquivo;

    std::cout << "Resultados detalhados serão salvos em: " << output_filename << std::endl;

    // Apaga o conteúdo antigo do arquivo de relatório, se houver
    std::ofstream(output_filename, std::ofstream::trunc).close();

    graph.writeResults(output_filename);

    // Gera e salva a árvore BFS a partir dos nós 1 e 10
    //if (num_vertices >= 1) graph.generateBfsReport(1, output_filename);
    
    // Gera e salva a árvore DFS a partir do nó 1
    //if (num_vertices >= 1) graph.generateDfsReport(1, output_filename);
    
    // Gera e salva o relatório de componentes conexas
    graph.generateConnectedComponentsReport(output_filename);
    
    // Cálculo de Distância e Caminho mínimo (Questão 1)

    for (int target : {20, 30, 40, 50, 60}) {
    auto path = graph.getPath(target, 10);

    std::cout << "Caminho " << 10 << " -> " << target << ": ";
    for (int v : path){
        std::cout << v + 1 << " ";
    }
    std::cout << std::endl;
    std::cout << "Distância entre o vértice 1 e o vértice "<< target <<": "<< graph.getDistance(10, target) << std::endl;
    }

    // Análise de Tempo: Dijkstra - Vetor (Questão 2)

    std::cout << "\n 2. Análise de Tempo: Dikjstra - Vetor (100 execuções)" << std::endl;
    {
        std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> dist(0, num_vertices - 1);
        
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 100; ++i) {
            graph.dijkstra(dist(rng), graph_tools_lib::DijkstraImplType::VECTOR);
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << "Vetor: Tempo médio por Dijkstra: " << duration.count() / 100.0 << " millissegundos" << std::endl;
    }

    // Análise de Tempo: Dijkstra - Heap (Questão 2)
    std::cout << "\n 3. Análise de Tempo: Dijkstra - Heap (100 execuções)" << std::endl;
    {
        std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> dist(0, num_vertices - 1);

        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 100; ++i) {
            graph.dijkstra(dist(rng));
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        std::cout << "Heap: Tempo médio por Dijkstra: " << duration.count() / 100.0 << " millissegundos" << std::endl;
    }



    // Rede de colaboração entre pesquisadores da area de Computação (Questão 3)

    ifstream file("rede_colaboracao_vertices.txt");
    if (!file.is_open()) {
        cerr << "Erro ao abrir nomes.txt\n";
        cout << 1;
    }

    unordered_map<string, int> name_to_id; // nome → índice
    unordered_map<int, string> id_to_name; // índice → nome

    string line;
    while (getline(file, line)) {
        if (line.empty()) continue;
        size_t pos = line.find(',');
        if (pos == string::npos) continue; // linha mal formatada
        int id = stoi(line.substr(0, pos));
        string name = line.substr(pos + 1);
        id_to_name[id] = name;
        name_to_id[name] = id;
    }
    file.close();
    
    auto inicio = name_to_id["Edsger W. Dijkstra"];
    for (int target : {name_to_id["Alan M. Turing"], name_to_id ["J. B. Kruskal"], name_to_id["Jon M. Kleinberg"] , name_to_id["Éva Tardos"], name_to_id["Daniel R. Figueiredo"]}) {
    auto path = graph.getPath(target, inicio);

    std::cout << "Caminho " << "Edsger W. Dijkstra" << " -> " << id_to_name[target] << ": ";
    for (size_t i = 0; i < path.size(); ++i) {
    std::cout << id_to_name[path[i] + 1];
    if (i + 1 < path.size()) std::cout << " -> ";
    }
    std::cout << std::endl;
    std::cout << "Distância entre Edsger W. Dijkstra e "<< id_to_name[target] <<": "<< graph.getDistance(inicio, target) << std::endl;
    std::cout << std::endl;
    }

    
    std::cout << "Análise de '" << nome_arquivo << "' concluída." << std::endl;
}

int main() {
    std::vector<std::string> arquivos_de_grafos = {"rede_colaboracao.txt"};

    for (const auto& arquivo : arquivos_de_grafos) {
        try {
            analisar_grafo(arquivo);
        } catch (const std::exception& e) {
            std::cerr << "ERRO ao processar " << arquivo << ": " << e.what() << std::endl;
        }
    }
    return 0;
}



