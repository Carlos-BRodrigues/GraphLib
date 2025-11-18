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

    auto result = graph.Run_Search(1);
    
}

int main() {
    std::vector<std::string> arquivos_de_grafos = {"test.txt","test2.txt", "test3.txt"};
    
    for (const auto& arquivo : arquivos_de_grafos) {
        try {
            analisar_grafo(arquivo);
        } catch (const std::exception& e) {
            std::cerr << "ERRO ao processar " << arquivo << ": " << e.what() << std::endl;
        }
    }
    return 0;
}
