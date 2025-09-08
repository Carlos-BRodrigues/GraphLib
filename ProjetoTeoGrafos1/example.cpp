#include "graph_tools.hpp"
#include <iostream>

int main() {
    try {
        std::cout << "--- Testando a representacao com LISTA DE ADJACENCIA ---" << std::endl;
        // O usuário escolhe a representação aqui
        graph_tools_lib::Graph graph_list("test2.txt", graph_tools_lib::RepresentationType::ADJACENCY_LIST);
        
        // E usa o objeto Graph normalmente, sem se preocupar com os detalhes
        graph_list.print();
        graph_list.writeResults("saida2.txt");
        //std::vector<double> stats_list = graph_list.getDegreeStats();
        //...

        std::cout << "\n--- Testando a representacao com MATRIZ DE ADJACENCIA ---" << std::endl;
        // O usuário pode criar outro grafo com a outra representação
        graph_tools_lib::Graph graph_matrix("test.txt", graph_tools_lib::RepresentationType::ADJACENCY_MATRIX);
        graph_matrix.print();
        graph_matrix.writeResults("saida.txt");

    } catch (const std::exception& e) {
        std::cerr << "Ocorreu um erro: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}