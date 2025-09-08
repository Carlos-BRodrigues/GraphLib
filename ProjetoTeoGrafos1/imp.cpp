#include "graph_tools.hpp"
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <stdexcept>

namespace graph_tools_lib {

// --- MUDANÇA 1: Implementação da AdjacencyList ---
AdjacencyList::AdjacencyList(int num_vertices) : num_vertices_(num_vertices), num_edges_(0) {
    adj_list_.assign(num_vertices, std::vector<int>());
}

void AdjacencyList::addEdge(int u, int v) {
    adj_list_[u].push_back(v+1);
    adj_list_[v].push_back(u+1);
    num_edges_++;
}

int AdjacencyList::getVertexCount() const { return num_vertices_; }

int AdjacencyList::getEdgeCount() const { return num_edges_; }

int AdjacencyList::getDegree(int v) const { return adj_list_[v].size(); }

void AdjacencyList::print() const {

    for (int i = 0; i<this->getVertexCount(); i++){
        std::cout<<i+1<<" => ";
        for(int j=0; j<this->getDegree(i); j++){
            std::cout<<adj_list_[i][j]<<' ';
        }
        std::cout<<std::endl;
    }

 }


// --- MUDANÇA 2: Implementação da AdjacencyMatrix ---
AdjacencyMatrix::AdjacencyMatrix(int num_vertices) : num_vertices_(num_vertices), num_edges_(0) {
    adj_matrix_.assign(num_vertices, std::vector<int>(num_vertices, 0));
}

void AdjacencyMatrix::addEdge(int u, int v) {
    if (adj_matrix_[u][v] == 0) { // Evita contar arestas duplicadas
        num_edges_++;
    }
    adj_matrix_[u][v] = 1;
    adj_matrix_[v][u] = 1;
}

int AdjacencyMatrix::getVertexCount() const { return num_vertices_; }

int AdjacencyMatrix::getEdgeCount() const { return num_edges_; }

int AdjacencyMatrix::getDegree(int v) const {
    return std::accumulate(adj_matrix_[v].begin(), adj_matrix_[v].end(), 0);
}

void AdjacencyMatrix::print() const { 

    for (int i = 0; i<this->getVertexCount(); i++){
        std::cout<<i+1<<" => ";
        for(int j=0; j<this->getVertexCount(); j++){
            if(adj_matrix_[i][j]==1){
                std::cout<<j+1<<' ';
            }
        }
        std::cout<<std::endl;
    }
    
 }


// --- MUDANÇA 3: Refatoração do Construtor e Métodos da Graph ---

// O construtor agora é uma "fábrica" que constrói a representação correta.
Graph::Graph(const std::string& filename, RepresentationType type) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: Could not open file " + filename);
    }
    
    int num_vertices;
    file >> num_vertices;

    if (type == RepresentationType::ADJACENCY_LIST) {
        representation_ = std::make_unique<AdjacencyList>(num_vertices);
    } else { // ADJACENCY_MATRIX
        representation_ = std::make_unique<AdjacencyMatrix>(num_vertices);
    }

    int node1, node2;
    while (file >> node1 >> node2) {
        representation_->addEdge(node1 - 1, node2 - 1); // Delega a adição da aresta
    }
    file.close();
}

// Os métodos da Graph agora simplesmente delegam as chamadas.
std::vector<double> Graph::getDegreeStats() const {
    int n = representation_->getVertexCount();
    std::vector<double> degrees;
    degrees.reserve(n);
    for (int i = 0; i < n; ++i) {
        degrees.push_back(representation_->getDegree(i));
    }
    // ... o resto da lógica (min, max, media, etc.) permanece o mesmo ...

     if (degrees.empty()) {
        return {0, 0, 0, 0}; // Return zero for all stats if graph is empty.
    }
    
    double min_degree = *std::min_element(degrees.begin(), degrees.end());
    double max_degree = *std::max_element(degrees.begin(), degrees.end());
    double sum = std::accumulate(degrees.begin(), degrees.end(), 0.0);
    double mean_degree = sum / degrees.size();

    std::sort(degrees.begin(), degrees.end());
    
    double median_degree;
    auto n2 = degrees.size();
    if (n2 % 2 == 1) {
        median_degree = degrees[n2 / 2];
    } else {
        median_degree = (degrees[n2 / 2 - 1] + degrees[n2 / 2]) / 2.0;
    }
    
    return {min_degree, max_degree, mean_degree, median_degree};
}


void Graph::print() const {
    // A classe Graph não precisa saber como imprimir, ela apenas pede!
    representation_->print();
}

bool Graph::writeResults(const std::string& output_filename) const {
    std::ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        std::cerr << "Error opening file for writing: " << output_filename << std::endl;
        return false;
    }
    
    output_file << "Number of vertices: " << representation_->getVertexCount() << std::endl;
    output_file << "Number of edges: " << representation_->getEdgeCount() << std::endl;
    
    std::vector<double> stats = this->getDegreeStats();
    output_file << "Min degree: " << stats[0] << std::endl;
    output_file << "Max degree: " << stats[1] << std::endl;
    output_file << "Mean degree: " << stats[2] << std::endl;
    output_file << "Median degree: " << stats[3] << std::endl;
    
    output_file.close();
    return true;
}

} // namespace