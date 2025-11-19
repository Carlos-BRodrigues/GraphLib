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
#include <chrono> 
#include <random> 

using namespace std;
using namespace graph_tools_lib;

void analisar_grafo(const std::string& nome_arquivo) {
    std::cout << "\n" << std::endl;
    std::cout << "Analisando o grafo: " << nome_arquivo << std::endl;
    std::cout << " " << std::endl;

    graph_tools_lib::Graph graph(nome_arquivo, graph_tools_lib::RepresentationType::ADJACENCY_LIST, graph_tools_lib::Direction_Graph::Direction);
    int num_vertices = graph.getVertexCount();
    std::string output_filename = "analise_" + nome_arquivo;

    std::cout << "Resultados detalhados serão salvos em: " << output_filename << std::endl;

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

