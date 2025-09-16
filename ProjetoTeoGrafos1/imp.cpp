#include "graph_tools.hpp"
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <queue>
#include <stack>

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

void AdjacencyList::BFS(int start_node, std::vector<int>& parent, std::vector<int>& level) const {
    parent.assign(num_vertices_, -1);
    level.assign(num_vertices_, -1);
    std::vector<bool> visited(num_vertices_, false);
    std::queue<int> q;

    parent[start_node] = -1;
    level[start_node] = 0;
    visited[start_node] = true;
    q.push(start_node);

    while (!q.empty()) {
        int u = q.front();
        q.pop();

        for (int v : adj_list_[u]) {
            int v_idx = v - 1;
            if (!visited[v_idx]) {
                visited[v_idx] = true;
                parent[v_idx] = u;
                level[v_idx] = level[u] + 1;
                q.push(v_idx);
            }
        }
    }
}

void AdjacencyList::DFS(int start_node, std::vector<int>& parent, std::vector<int>& level) const {
    parent.assign(num_vertices_, -1);
    level.assign(num_vertices_, -1);
    std::vector<bool> visited(num_vertices_, false);
    
    std::stack<int> s;
    s.push(start_node);
    visited[start_node] = true;
    parent[start_node] = -1; // O vértice inicial não tem pai
    level[start_node] = 0;
    
    while (!s.empty()) {
        int u = s.top();
        s.pop();

        // Percorre os vizinhos em ordem inversa para que a ordem de visitação
        // seja a mesma da implementação recursiva (depende da ordem de inserção)
        for (auto it = adj_list_[u].rbegin(); it != adj_list_[u].rend(); ++it) {
            int v_idx = (*it) - 1;
            if (!visited[v_idx]) {
                visited[v_idx] = true;
                parent[v_idx] = u;
                level[v_idx] = level[u] + 1;
                s.push(v_idx);
            }
        }
    }
}

int AdjacencyList::getDistance(int u, int v) const {
    std::vector<int> parent, level;
    BFS(u, parent, level);
    return level[v];
}

int AdjacencyList::getDiameter() const {
    int max_distance = 0;
    for (int i = 0; i < num_vertices_; ++i) {
        std::vector<int> parent, level;
        BFS(i, parent, level);
        for (int j = 0; j < num_vertices_; ++j) {
            if (level[j] > max_distance) {
                max_distance = level[j];
            }
        }
    }
    return max_distance;
}

std::vector<std::vector<int>> AdjacencyList::getConnectedComponents() const {
    std::vector<std::vector<int>> components;
    std::vector<bool> visited(num_vertices_, false);
    
    for (int i = 0; i < num_vertices_; ++i) {
        if (!visited[i]) {
            std::vector<int> current_component;
            std::queue<int> q;
            q.push(i);
            visited[i] = true;
            
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                current_component.push_back(u + 1);
                
                for (int v : adj_list_[u]) {
                    int v_idx = v - 1;
                    if (!visited[v_idx]) {
                        visited[v_idx] = true;
                        q.push(v_idx);
                    }
                }
            }
            components.push_back(current_component);
        }
    }
    
    std::sort(components.begin(), components.end(), [](const auto& a, const auto& b) {
        return a.size() > b.size();
    });
    
    return components;
}

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

void AdjacencyMatrix::BFS(int start_node, std::vector<int>& parent, std::vector<int>& level) const {
    parent.assign(num_vertices_, -1);
    level.assign(num_vertices_, -1);
    std::vector<bool> visited(num_vertices_, false);
    std::queue<int> q;

    parent[start_node] = -1;
    level[start_node] = 0;
    visited[start_node] = true;
    q.push(start_node);

    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (int v = 0; v < num_vertices_; ++v) {
            if (adj_matrix_[u][v] == 1 && !visited[v]) {
                visited[v] = true;
                parent[v] = u;
                level[v] = level[u] + 1;
                q.push(v);
            }
        }
    }
}

void AdjacencyMatrix::DFS(int start_node, std::vector<int>& parent, std::vector<int>& level) const {
    parent.assign(num_vertices_, -1);
    level.assign(num_vertices_, -1);
    std::vector<bool> visited(num_vertices_, false);

    std::stack<int> s;
    s.push(start_node);
    visited[start_node] = true;
    parent[start_node] = -1;
    level[start_node] = 0;

    while (!s.empty()) {
        int u = s.top();
        s.pop();

        for (int v = 0; v < num_vertices_; ++v) {
            if (adj_matrix_[u][v] == 1 && !visited[v]) {
                visited[v] = true;
                parent[v] = u;
                level[v] = level[u] + 1;
                s.push(v);
            }
        }
    }
}

int AdjacencyMatrix::getDistance(int u, int v) const {
    std::vector<int> parent, level;
    BFS(u, parent, level);
    return level[v];
}

int AdjacencyMatrix::getDiameter() const {
    int max_distance = 0;
    for (int i = 0; i < num_vertices_; ++i) {
        std::vector<int> parent, level;
        BFS(i, parent, level);
        for (int j = 0; j < num_vertices_; ++j) {
            if (level[j] > max_distance) {
                max_distance = level[j];
            }
        }
    }
    return max_distance;
}

std::vector<std::vector<int>> AdjacencyMatrix::getConnectedComponents() const {
    std::vector<std::vector<int>> components;
    std::vector<bool> visited(num_vertices_, false);
    
    for (int i = 0; i < num_vertices_; ++i) {
        if (!visited[i]) {
            std::vector<int> current_component;
            std::queue<int> q;
            q.push(i);
            visited[i] = true;
            
            while (!q.empty()) {
                int u = q.front();
                q.pop();
                current_component.push_back(u + 1);
                
                for (int v = 0; v < num_vertices_; ++v) {
                    if (adj_matrix_[u][v] == 1 && !visited[v]) {
                        visited[v] = true;
                        q.push(v);
                    }
                }
            }
            components.push_back(current_component);
        }
    }
    
    std::sort(components.begin(), components.end(), [](const auto& a, const auto& b) {
        return a.size() > b.size();
    });
    
    return components;
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
void Graph::BFS(int start_node, const std::string& output_file) const {
    std::vector<int> parent, level;
    representation_->BFS(start_node - 1, parent, level);
    writeSearchTree(output_file, "Busca em Largura (BFS)", start_node, parent, level);
}

void Graph::DFS(int start_node, const std::string& output_file) const {
    std::vector<int> parent, level;
    representation_->DFS(start_node - 1, parent, level);
    writeSearchTree(output_file, "Busca em Profundidade (DFS)", start_node, parent, level);
}

void Graph::writeSearchTree(const std::string& output_file, const std::string& algorithm,
                            int start_node, const std::vector<int>& parent,
                            const std::vector<int>& level) const {
    std::ofstream out(output_file, std::ios_base::app);
    if (!out.is_open()) {
        throw std::runtime_error("Erro ao abrir o arquivo para escrita.");
    }
    
    out << "\n--- Resultados da " << algorithm << " a partir do vertice " << start_node << " ---\n";
    for (size_t i = 0; i < parent.size(); ++i) {
        out << "Vertice: " << i + 1 << ", ";
        if (parent[i] != -1) {
            out << "Pai: " << parent[i] + 1 << ", ";
        } else {
            out << "Pai: Nenhum, ";
        }
        out << "Nivel: " << level[i] << std::endl;
    }
    out.close();
}

int Graph::getDistance(int u, int v) const {
    return representation_->getDistance(u - 1, v - 1);
}

int Graph::getDiameter() const {
    return representation_->getDiameter();
}

std::vector<std::vector<int>> Graph::getConnectedComponents() const {
    return representation_->getConnectedComponents();
}

bool Graph::writeConnectedComponents(const std::string& output_filename) const {
    std::ofstream out(output_filename, std::ios_base::app);
    if (!out.is_open()) {
        std::cerr << "Erro ao abrir o arquivo para escrita: " << output_filename << std::endl;
        return false;
    }

    auto components = getConnectedComponents();
    out << "\n--- Componentes Conexas ---\n";
    out << "Numero de componentes: " << components.size() << std::endl;

    for (const auto& component : components) {
        out << "Tamanho: " << component.size() << " | Vertices: ";
        for (int v : component) {
            out << v << " ";
        }
        out << std::endl;
    }
    out.close();
    return true;
}

} // namespace
