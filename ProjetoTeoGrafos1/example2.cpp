#include "graph_tools.hpp"
#include <iostream>

int main() {
    try {
        
        //Cria o grafo como a lista de adjacência
        graph_tools_lib::Graph graph_list("test.txt", graph_tools_lib::RepresentationType::ADJACENCY_LIST);
        //Cria o grafo como a matriz de adjacência


        //Imprime e escreve no arquivo de saída
        graph_list.print();
        graph_list.writeResults("saida2.txt");

        int start_node_bfs = 1;
        graph_list.BFS(start_node_bfs, "saida2.txt");
        
        //Busca em Profundidade (DFS)
        int start_node_dfs = 1; // Exemplo de vértice inicial
        graph_list.DFS(start_node_dfs, "saida2.txt");
        
        //Componentes Conexas
        graph_list.writeConnectedComponents("saida2.txt");
        
        //Diâmetro do Grafo
        int diameter_list = graph_list.getDiameter();
        std::cout << "Diametro do grafo (Lista): " << diameter_list << std::endl;

        int apdiameter_list = graph_list.getDiameterapprox();
        std::cout << "Diametro aproximado do grafo (Lista): " << apdiameter_list << std::endl;

        //Distância entre dois vértices (exemplo: 1 e 5)
        int distance_list = graph_list.getDistance(1, 5);
        std::cout << "Distancia entre 1 e 5 (Lista): " << distance_list << std::endl;
        
        std::cout << "\n--- Testando a representacao com MATRIZ DE ADJACENCIA ---" << std::endl;
        
        //Cria o grafo com a matriz de adjacência
        graph_tools_lib::Graph graph_matrix("test.txt", graph_tools_lib::RepresentationType::ADJACENCY_MATRIX);
        
        //Imprime o grafo e as estatísticas de grau no arquivo de saída
        graph_matrix.print();
        graph_matrix.writeResults("saida.txt");
        
        
        graph_matrix.BFS(1, "saida.txt");
        
        graph_matrix.DFS(1, "saida.txt");
        
        graph_matrix.writeConnectedComponents("saida.txt");
        
        int diameter_matrix = graph_matrix.getDiameter();
        std::cout << "Diametro do grafo (Matriz): " << diameter_matrix << std::endl;

        int apdiameter_matrix = graph_matrix.getDiameterapprox();
        std::cout << "Diametro aproximado do grafo (Matriz): " << apdiameter_matrix << std::endl;

        int distance_matrix = graph_matrix.getDistance(1, 5);
        std::cout << "Distancia entre 1 e 5 (Matriz): " << distance_matrix << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Ocorreu um erro: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
