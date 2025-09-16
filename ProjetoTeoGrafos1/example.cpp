#include "graph_tools.hpp"
#include <iostream>

int main() {
    try {
        std::cout << "--- Testando a representacao com LISTA DE ADJACENCIA ---" << std::endl;
        
        // Crie o grafo com a lista de adjacência
        graph_tools_lib::Graph graph_list("test2.txt", graph_tools_lib::RepresentationType::ADJACENCY_LIST);
        
        // Imprime o grafo e as estatísticas de grau no arquivo de saída
        graph_list.print();
        graph_list.writeResults("saida2.txt");
        
        // --- Testando as novas funcionalidades com Lista de Adjacência ---
        
        // 1. Busca em Largura (BFS)
        int start_node_bfs = 1; // Exemplo de vértice inicial
        graph_list.BFS(start_node_bfs, "saida2.txt");
        
        // 2. Busca em Profundidade (DFS)
        int start_node_dfs = 1; // Exemplo de vértice inicial
        graph_list.DFS(start_node_dfs, "saida2.txt");
        
        // 3. Componentes Conexas
        graph_list.writeConnectedComponents("saida2.txt");
        
        // 4. Diâmetro do Grafo
        int diameter_list = graph_list.getDiameter();
        std::cout << "Diametro do grafo (Lista): " << diameter_list << std::endl;

        // 5. Distância entre dois vértices (exemplo: 1 e 5)
        int distance_list = graph_list.getDistance(1, 4);
        std::cout << "Distancia entre 1 e 4 (Lista): " << distance_list << std::endl;
        
        std::cout << "\n--- Testando a representacao com MATRIZ DE ADJACENCIA ---" << std::endl;
        
        // Crie o grafo com a matriz de adjacência
        graph_tools_lib::Graph graph_matrix("test.txt", graph_tools_lib::RepresentationType::ADJACENCY_MATRIX);
        
        // Imprime o grafo e as estatísticas de grau no arquivo de saída
        graph_matrix.print();
        graph_matrix.writeResults("saida.txt");
        
        // --- Testando as novas funcionalidades com Matriz de Adjacência ---
        
        // 1. Busca em Largura (BFS)
        graph_matrix.BFS(1, "saida.txt");
        
        // 2. Busca em Profundidade (DFS)
        graph_matrix.DFS(1, "saida.txt");
        
        // 3. Componentes Conexas
        graph_matrix.writeConnectedComponents("saida.txt");
        
        // 4. Diâmetro do Grafo
        int diameter_matrix = graph_matrix.getDiameter();
        std::cout << "Diametro do grafo (Matriz): " << diameter_matrix << std::endl;

        // 5. Distância entre dois vértices (exemplo: 1 e 5)
        int distance_matrix = graph_matrix.getDistance(1, 5);
        std::cout << "Distancia entre 1 e 5 (Matriz): " << distance_matrix << std::endl;

    } catch (const std::exception& e) {
        std::cerr << "Ocorreu um erro: " << e.what() << std::endl;
        return 1;
    }
    return 0;
}
