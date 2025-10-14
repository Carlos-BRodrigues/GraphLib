#include "graph_tools.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <fstream>

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
    
    std::cout<<graph.getDistance(1 , 5)<<std::endl;
    
    std::cout<<graph.getApproximateDiameter()<<std::endl;

    graph.generateDijkstraReport(1 , output_filename, graph_tools_lib::DijkstraImplType::VECTOR);


    std::cout << "Análise de '" << nome_arquivo << "' concluída." << std::endl;
}

int main() {
    std::vector<std::string> arquivos_de_grafos = {"test.txt"};

    for (const auto& arquivo : arquivos_de_grafos) {
        try {
            analisar_grafo(arquivo);
        } catch (const std::exception& e) {
            std::cerr << "ERRO ao processar " << arquivo << ": " << e.what() << std::endl;
        }
    }
    return 0;
}

