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
    graph_tools_lib::Graph graph(nome_arquivo, graph_tools_lib::RepresentationType::ADJACENCY_MATRIX, graph_tools_lib::Direction_Graph::Direction);
    int num_vertices = graph.getVertexCount();
    std::string output_filename = "analise_" + nome_arquivo;

    std::cout << "Resultados detalhados serão salvos em: " << output_filename << std::endl;

    // Apaga o conteúdo antigo do arquivo de relatório, se houver
    std::ofstream(output_filename, std::ofstream::trunc).close();

    graph.writeResults(output_filename);
    
    // --- 1. Teste de BFS e DFS ---
    std::cout << "\n--- 1. Teste de Travessias (BFS e DFS) ---" << std::endl;
    
    // Gera e salva a árvore BFS a partir do nó 1 (vértice 0)
    if (num_vertices >= 1) {
        std::cout << "Executando BFS a partir do nó 1..." << std::endl;
        graph.generateBfsReport(1, output_filename);
    }
    
    // Gera e salva a árvore DFS a partir do nó 1
    if (num_vertices >= 1) {
        std::cout << "Executando DFS a partir do nó 1..." << std::endl;
        graph.generateDfsReport(1, output_filename);
    }
    
    // Gera e salva o relatório de componentes conexas
    graph.generateConnectedComponentsReport(output_filename);
    
    graph.print();

    
    // Gera e salva o relatório de Dijkstra a partir do nó 10
    /** if (num_vertices >= start_node_test) {
        std::cout << "Executando Dijkstra a partir do nó " << start_node_test << "..." << std::endl;
        graph.generateDijkstraReport(start_node_test, output_filename);
    }
    
    // Exemplo de cálculo de caminho e distância
    for (int target : {20, 30, 40, 50, 60}) {
        if (num_vertices > target) {
            auto path = graph.getPath(target, start_node_test); // Assumindo que getPath calcula com base no último Dijkstra/Bellman-Ford
            
            std::cout << "Caminho " << start_node_test << " -> " << target << ": ";
            for (int v : path){
                std::cout << v + 1 << " ";
            }
            std::cout << std::endl;
            std::cout << "Distância entre o vértice " << start_node_test << " e o vértice "<< target <<": "<< graph.getDistance(start_node_test, target) << std::endl;
        }
    }
    **/
    // --- 3. Teste de Bellman-Ford ---
    std::cout << "\n--- 3. Teste de Bellman-Ford ---" << std::endl;

    if (graph.hasNegativeWeights()) {
        std::cout << "Pesos negativos detectados. Executando Bellman-Ford a partir do nó 3..." << std::endl;
        // Chamada direta, assumindo que Bellman_Ford retorna SearchResult
        try {
            graph.Bellman_Ford(3);
            std::cout << "Bellman-Ford concluído com sucesso (sem ciclos negativos)." << std::endl;
            graph.generateBellman_FordReport(3, output_filename);
        } catch (const std::exception& e) {
            std::cerr << "Bellman-Ford ERRO: " << e.what() << std::endl;
        }
    } else {
        std::cout << "Sem pesos negativos. Bellman-Ford não é necessário/executado para teste." << std::endl;
    }
    
}

int main() {
    std::vector<std::string> arquivos_de_grafos = {"test2.txt"};
    
    for (const auto& arquivo : arquivos_de_grafos) {
        try {
            analisar_grafo(arquivo);
        } catch (const std::exception& e) {
            std::cerr << "ERRO ao processar " << arquivo << ": " << e.what() << std::endl;
        }
    }
    return 0;
}
