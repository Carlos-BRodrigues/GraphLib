#include "graph_tools.hpp"
#include <iostream>
#include <vector>
#include <string>
#include <chrono>   // Para medir o tempo
#include <random>   // Para gerar vértices aleatórios
#include <numeric>  // Para std::accumulate
#include <algorithm>// Para std::min_element, std::max_element

void analisar_grafo(const std::string& nome_arquivo) {
    std::cout << "\n==================================================" << std::endl;
    std::cout << "Analisando o grafo: " << nome_arquivo << std::endl;
    std::cout << "==================================================" << std::endl;

    // --- Análise de Memória (Questão 1) ---
    std::cout << "\n--- 1. Análise de Memória ---" << std::endl;
    {
        graph_tools_lib::Graph graph_list(nome_arquivo, graph_tools_lib::RepresentationType::ADJACENCY_LIST);
        std::cout << "LISTA: Grafo com " << graph_list.getVertexCount() << " vertices carregado. Verifique a memória e pressione Enter." << std::endl;
        std::cin.get();
    }
    {
        graph_tools_lib::Graph graph_matrix(nome_arquivo, graph_tools_lib::RepresentationType::ADJACENCY_MATRIX);
        std::cout << "MATRIZ: Grafo com " << graph_matrix.getVertexCount() << " vertices carregado. Verifique a memória e pressione Enter." << std::endl;
        std::cin.get();
    }

    // Criar os grafos para os testes de tempo e algoritmos
    graph_tools_lib::Graph graph_list(nome_arquivo, graph_tools_lib::RepresentationType::ADJACENCY_LIST);
    graph_tools_lib::Graph graph_matrix(nome_arquivo, graph_tools_lib::RepresentationType::ADJACENCY_MATRIX);
    
    int num_vertices = graph_list.getVertexCount();
    if (num_vertices == 0) {
        std::cout << "Grafo vazio, pulando análises." << std::endl;
        return;
    }

    // --- Análise de Tempo: BFS (Questão 2) ---
    std::cout << "\n--- 2. Análise de Tempo: BFS (100 execuções) ---" << std::endl;
    {
        std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> dist(0, num_vertices - 1);
        
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 100; ++i) {
            graph_list.bfs(dist(rng));
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "LISTA: Tempo médio por BFS: " << duration.count() / 100.0 << " microssegundos" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 100; ++i) {
            graph_matrix.bfs(dist(rng));
        }
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "MATRIZ: Tempo médio por BFS: " << duration.count() / 100.0 << " microssegundos" << std::endl;
    }
    
    // --- Análise de Tempo: DFS (Questão 3) ---
    std::cout << "\n--- 3. Análise de Tempo: DFS (100 execuções) ---" << std::endl;
    {
        std::mt19937 rng(std::chrono::steady_clock::now().time_since_epoch().count());
        std::uniform_int_distribution<int> dist(0, num_vertices - 1);

        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 100; ++i) {
            graph_list.dfs(dist(rng));
        }
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "LISTA: Tempo médio por DFS: " << duration.count() / 100.0 << " microssegundos" << std::endl;

        start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < 100; ++i) {
            graph_matrix.dfs(dist(rng));
        }
        stop = std::chrono::high_resolution_clock::now();
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        std::cout << "MATRIZ: Tempo médio por DFS: " << duration.count() / 100.0 << " microssegundos" << std::endl;
    }
    
    // --- Pais na Árvore Geradora (Questão 4) ---
    std::cout << "\n--- 4. Pais na Árvore Geradora ---" << std::endl;
    if (num_vertices >= 30) {
        std::vector<int> start_nodes = {1, 2, 3};
        std::vector<int> target_nodes = {10, 20, 30};
        for (int start_node : start_nodes) {
            auto bfs_res = graph_list.bfs(start_node - 1);
            auto dfs_res = graph_list.dfs(start_node - 1);
            for (int target_node : target_nodes) {
                int bfs_parent = bfs_res.parent[target_node - 1];
                int dfs_parent = dfs_res.parent[target_node - 1];
                std::cout << "Iniciando em " << start_node << " -> Pai de " << target_node << " | BFS: " << (bfs_parent == -1 ? "Nenhum" : std::to_string(bfs_parent + 1))
                          << ", DFS: " << (dfs_parent == -1 ? "Nenhum" : std::to_string(dfs_parent + 1)) << std::endl;
            }
        }
    } else {
        std::cout << "Grafo com menos de 30 vértices, pulando teste." << std::endl;
    }

    // --- Distância entre Vértices (Questão 5) ---
    std::cout << "\n--- 5. Distância entre Vértices ---" << std::endl;
    if (num_vertices >= 30) {
        std::cout << "Distância (10, 20): " << graph_list.getDistance(10 - 1, 20 - 1) << std::endl;
        std::cout << "Distância (10, 30): " << graph_list.getDistance(10 - 1, 30 - 1) << std::endl;
        std::cout << "Distância (20, 30): " << graph_list.getDistance(20 - 1, 30 - 1) << std::endl;
    } else {
        std::cout << "Grafo com menos de 30 vértices, pulando teste." << std::endl;
    }

    // --- Componentes Conexas (Questão 6) ---
    std::cout << "\n--- 6. Componentes Conexas ---" << std::endl;
    auto componentes = graph_list.getConnectedComponents();
    std::cout << "Número de componentes conexas: " << componentes.size() << std::endl;
    if (!componentes.empty()) {
        auto min_max_it = std::minmax_element(componentes.begin(), componentes.end(), 
            [](const auto& a, const auto& b) {
                return a.size() < b.size();
            });
        std::cout << "Tamanho da menor componente: " << min_max_it.first->size() << std::endl;
        std::cout << "Tamanho da maior componente: " << min_max_it.second->size() << std::endl;
    }

    // --- Diâmetro do Grafo (Questão 7) ---
    std::cout << "\n--- 7. Diâmetro do Grafo ---" << std::endl;
    std::cout << "Diâmetro (aproximado): " << graph_list.getDiameter() << std::endl;
}

int main() {
    // Coloque os nomes corretos dos seus 6 arquivos de grafo aqui
    /*std::vector<std::string> arquivos_de_grafos = {
        "grafo_1.txt", "grafo_2.txt", "grafo_3.txt", 
        "grafo_4.txt", "grafo_5.txt", "grafo_6.txt"
    };*/

    std::vector<std::string> arquivos_de_grafos = {
        "grafo_1.txt", "grafo_2.txt"
    };

    for (const auto& arquivo : arquivos_de_grafos) {
        try {
            analisar_grafo(arquivo);
        } catch (const std::exception& e) {
            std::cerr << "Erro ao processar " << arquivo << ": " << e.what() << std::endl;
        }
    }

    return 0;
}