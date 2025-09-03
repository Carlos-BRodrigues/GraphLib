#include <iostream>
#include <vector>
#include <fstream> // !! Importante para arquivos
#include <string>
#include <tuple>
#include <algorithm>
#include <numeric>
#include <stdexcept> // For throwing errors
//#include "graph_tools.hpp" // Link to our own interface file

namespace graph_tools_lib {

bool Graph::ReadFile(const std::string& filename){

    std::ifstream arquivo_grafo(filename);

    if (!arquivo_grafo.is_open()) {
        std::cerr << "Erro ao abrir o arquivo para leitura!" << std::endl;
        return 1;
    }
    std::vector<std::string> lines;
    std::string line;

    while(std::getline(arquivo_grafo, line)){
        lines.push_back(line);
    }

    arquivo_grafo.close();

    int num_vertices; int num_arestas;

    num_vertices = lines[0][0] - '0';
    num_arestas = lines.size()-1;
    
    std::cout << num_vertices << std::endl;
    std::cout << num_arestas << std::endl;

    for(double j : grau(lines,num_vertices)){
        std::cout << j << std::endl;
    }

    std::ofstream arquivo_saida("saida.txt");

    if (!arquivo_saida.is_open()) {
        std::cerr << "Erro ao abrir o arquivo para escrita!" << std::endl;
        return 1;
    }
    
    arquivo_saida << num_vertices << std::endl;
    arquivo_saida << num_arestas << std::endl;
    for(double j : Graph::grau(lines,num_vertices)){
        arquivo_saida << j << std::endl;
    }

    arquivo_saida.close();

    return 0;
}

std::vector<double> grau(const std::vector<std::string>& graph, int n) const {

    std::vector<double> degrees(n, 0.0);

    for (int i = 1; i < graph.size(); ++i) {
        int node1 = graph[i][0] - '0';
        int node2 = graph[i][2] - '0';
        
        degrees[node1 - 1] += 1;
        degrees[node2 - 1] += 1;
    }

    double min_degree = *std::min_element(degrees.begin(), degrees.end());
    double max_degree = *std::max_element(degrees.begin(), degrees.end());
    double sum = std::accumulate(degrees.begin(), degrees.end(), 0.0);
    double mean_degree =  sum / degrees.size();

    std::sort(degrees.begin(), degrees.end());

    double median_degree;
    if (n % 2 == 1) {
        median_degree = degrees[n / 2];
    } else {
        median_degree = (degrees[n / 2 - 1] + degrees[n / 2]) / 2.0;
    }
    
    return {min_degree, max_degree, mean_degree, median_degree};
}

void print(const std::vector<std::string>& graph) const {
    for(std::string line : graph){
        std::cout << line << std::endl;
    }
}
}