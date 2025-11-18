#include "graph_tools3.hpp"
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <queue>
#include <stack>
#include <limits>
#include <sstream>
#include <chrono> // Para medição de tempo
#include <random> // Para geração de números aleatórios

// Certifique-se de incluir os namespaces necessários no seu ambiente
using namespace std;
using namespace graph_tools_lib;

// Adicionando a função Bellman_Ford à classe Graph (apenas para referência, você precisa 
// implementá-la em Graph.cpp, mas a chamada será feita aqui)
// SearchResult Graph::Bellman_Ford(int start_node) const; // Supondo que você use start_node, não end_node.

void analisar_grafo(const std::string& nome_arquivo) {
    std::cout << "\n==============================================" << std::endl;
    std::cout << "Analisando o grafo: " << nome_arquivo << std::endl;
    std::cout << "==============================================" << std::endl;

    // Usando ADJACENCY_LIST para todos os testes
    // Mude NO_Direction para Direction para testar grafos direcionados
    graph_tools_lib::Graph graph(nome_arquivo, graph_tools_lib::RepresentationType::ADJACENCY_LIST, graph_tools_lib::Direction_Graph::Direction);
    int num_vertices = graph.getVertexCount();
    std::string output_filename = "analise_" + nome_arquivo;

    std::cout << "Resultados detalhados serão salvos em: " << output_filename << std::endl;

    // Apaga o conteúdo antigo do arquivo de relatório, se houver
    std::ofstream(output_filename, std::ofstream::trunc).close();

    graph.writeResults(output_filename);

    auto grafo_reverso = graph.reverseEdges();

    
    // Estudo de caso 1
    std::cout << "\n--- 1. Teste de Bellman-Ford ---" << std::endl;

    auto res = graph.Bellman_Ford(100-1);
    for (int start : {10, 20, 30}) {
        if (num_vertices > start) {
            std::cout << "Distância entre o vértice " << start << " e o vértice 100: "<< res.distance[start-1] << std::endl;
        }
    }
    
    // Estudo de caso 2 ---

    std::cout << "\n 2. Análise de Tempo: Bellman-Ford (10 execuções)" << std::endl;
    {
        std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> dist(0, num_vertices - 1);
        
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 10; ++i) {
            graph.Bellman_Ford(dist(rng));
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        std::cout << "Tempo médio por Bellman_Ford: " << duration.count() / 10.0 << "segundos" << std::endl;
    }

    // Estudo de caso 3

    std::cout << "\n--- 3. Teste Dijkstra ---" << std::endl;
    if (!graph.hasNegativeWeights()) {
    
    auto res2 = grafo_reverso.dijkstra(100-1);
    for (int start : {10, 20, 30}) {
        if (num_vertices > start) {
            std::cout << "Distância entre o vértice " << start << " e o vértice 100: "<< res2.distance[start-1] << std::endl;
        }
    }

    std::cout << "\n  Análise de Tempo: Dijkstra (10 execuções)" << std::endl;
    {
        std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> dist(0, num_vertices - 1);
        
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 10; ++i) {
            grafo_reverso.dijkstra(dist(rng));
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        std::cout << "Tempo médio por dijkstra: " << duration.count() / 10.0 << "segundos" << std::endl;
    }
    }
    else{
        std::cout << "O Grafo possui pesos negativos" << std::endl;
    }
}

int main() {
    std::vector<std::string> arquivos_de_grafos = {"grafo_W_1.txt", "grafo_W_2.txt", "grafo_W_3.txt", "grafo_W_4.txt", "grafo_W_5.txt"};
    
    for (const auto& arquivo : arquivos_de_grafos) {
        try {
            analisar_grafo(arquivo);
        } catch (const std::exception& e) {
            std::cerr << "ERRO ao processar " << arquivo << ": " << e.what() << std::endl;
        }
    }
    return 0;
}
