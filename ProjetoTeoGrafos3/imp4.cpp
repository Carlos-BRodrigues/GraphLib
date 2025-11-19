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

namespace graph_tools_lib {


AdjacencyList::AdjacencyList(int num_vertices, Direction_Graph direction) : num_vertices_(num_vertices), num_edges_(0), direction_type_(direction) {
    adj_list_.assign(num_vertices, std::vector<Edge>());
}
void AdjacencyList::addEdge(int u, int v, double weight) {
    adj_list_[u].push_back({v, weight});
    if (direction_type_ == Direction_Graph::NO_Direction){ // adiciona o inverso da aresta se não for direcionado
        adj_list_[v].push_back({u, weight});
    }
    num_edges_++;
}
int AdjacencyList::getVertexCount() const { return num_vertices_; }
int AdjacencyList::getEdgeCount() const { return num_edges_; }
std::vector<Edge> AdjacencyList::getNeighbors(int v) const {
    return adj_list_[v];
}
void AdjacencyList::print() const {
    for (int i = 0; i < num_vertices_; ++i) {
        std::cout << i + 1 << " => ";
        for (const auto& edge : adj_list_[i]) {
            std::cout << "[" << edge.target + 1 << ", w:" << edge.weight << "] ";
        }
        std::cout << std::endl;
    }
}
AdjacencyMatrix::AdjacencyMatrix(int num_vertices, Direction_Graph direction) : num_vertices_(num_vertices), num_edges_(0), direction_type_(direction) {
    adj_matrix_.assign(num_vertices, std::vector<double>(num_vertices, std::numeric_limits<double>::infinity()));
}
void AdjacencyMatrix::addEdge(int u, int v, double weight) {
    if (adj_matrix_[u][v] == std::numeric_limits<double>::infinity()) num_edges_++;
    adj_matrix_[u][v] = weight;
    if (direction_type_ == Direction_Graph::NO_Direction){ // adiciona o inverso da aresta se não for direcionado
        adj_matrix_[v][u] = weight;
    }
}
int AdjacencyMatrix::getVertexCount() const { return num_vertices_; }
int AdjacencyMatrix::getEdgeCount() const { return num_edges_; }
std::vector<Edge> AdjacencyMatrix::getNeighbors(int v) const {
    std::vector<Edge> neighbors;
    for (int i = 0; i < num_vertices_; ++i) {
        if (adj_matrix_[v][i] != 0.0) {
            neighbors.push_back({i, adj_matrix_[v][i]});
        }
    }
    return neighbors;
}
void AdjacencyMatrix::print() const {
    for (int i = 0; i < num_vertices_; ++i) {
        std::cout << i + 1 << " => ";
        for (int j = 0; j < num_vertices_; ++j) {
            if (adj_matrix_[i][j] != 0.0) {
                std::cout << "[" << j + 1 << ", w:" << adj_matrix_[i][j] << "] ";
            }
        }
        std::cout << std::endl;
    }
}




Graph::Graph(const std::string& filename, RepresentationType type, Direction_Graph direction) : has_negative_weights_(false), has_weights(false) {
    std::ifstream file(filename);
    if (!file.is_open()) throw std::runtime_error("Erro: não pôde abrir o arquivo " + filename);
    int num_vertices;
    file >> num_vertices;
    if (type == RepresentationType::ADJACENCY_LIST) {
        representation_ = std::make_unique<AdjacencyList>(num_vertices, direction);
    } else {
        representation_ = std::make_unique<AdjacencyMatrix>(num_vertices, direction);
    }
    std::string line;
    std::getline(file, line);
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        int node1, node2;
        double weight = 1.0;
        ss >> node1 >> node2;
        if (ss >> weight) {
            if (weight != 1) has_weights = true;
            if (weight < 0) has_negative_weights_ = true;
        }
        if (node1 > 0 && node2 > 0) {
            representation_->addEdge(node1 - 1, node2 - 1, weight);
        }
    }
    file.close();
}

int Graph::getVertexCount() const { return representation_->getVertexCount(); }
int Graph::getEdgeCount() const { return representation_->getEdgeCount(); }
void Graph::print() const { representation_->print(); }
bool Graph::hasNegativeWeights() const { return has_negative_weights_; }
bool Graph::hasWeights() const { return has_weights; }

SearchResult Graph::bfs(int start_node) const {
    int n = getVertexCount();
    SearchResult result;
    result.parent.assign(n, -1);
    result.distance.assign(n, std::numeric_limits<double>::infinity());
    std::queue<int> q;
    result.distance[start_node] = 0.0;
    q.push(start_node);
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        for (const auto& edge : representation_->getNeighbors(u)) {
            int v = edge.target;
            if (result.distance[v] < 0 || result.distance[v] == std::numeric_limits<double>::infinity()) {
                result.distance[v] = result.distance[u] + 1;
                result.parent[v] = u;
                q.push(v);
            }
        }
    }
    return result;
}

SearchResult Graph::dfs(int start_node) const {
    int n = getVertexCount();
    SearchResult result;
    result.parent.assign(n, -1);
    result.distance.assign(n, std::numeric_limits<double>::infinity());
    std::stack<int> s;
    std::vector<bool> visited(n, false);

    s.push(start_node);
    visited[start_node] = true;
    result.distance[start_node] = 0.0;

    while (!s.empty()) {
        int u = s.top();
        s.pop();
        for (const auto& edge : representation_->getNeighbors(u)) {
            int v = edge.target;
            if (!visited[v]) {
                visited[v] = true;
                result.parent[v] = u;
                result.distance[v] = result.distance[u] + 1;
                s.push(v);
            }
        }
    }
    return result;
}

/**
 * @brief Implementação unificada de Dijkstra.
 * O usuário pode escolher entre a versão com vetor (lenta, O(V^2)) ou
 * com heap (rápida, O(E log V)). A versão com heap é a padrão.
 */
SearchResult Graph::dijkstra(int start_node, DijkstraImplType impl_type) const {
    if (has_negative_weights_) {
        throw std::runtime_error("O algoritmo de Dijkstra não suporta pesos negativos.");
    }
    
    int n = getVertexCount();
    SearchResult result;
    result.distance.assign(n, std::numeric_limits<double>::infinity());
    result.parent.assign(n, -1);
    result.distance[start_node] = 0.0;

    //Heap

    if (impl_type == DijkstraImplType::HEAP) {
        using PDI = std::pair<double, int>;
        std::priority_queue<PDI, std::vector<PDI>, std::greater<PDI>> pq;
        pq.push({0.0, start_node});

        while (!pq.empty()) {
            double d = pq.top().first;
            int u = pq.top().second;
            pq.pop();
            if (d > result.distance[u]) continue;

            for (const auto& edge : representation_->getNeighbors(u)) {
                int v = edge.target;
                double weight = edge.weight;
                if (result.distance[u] + weight < result.distance[v]) {
                    result.distance[v] = result.distance[u] + weight;
                    result.parent[v] = u;
                    pq.push({result.distance[v], v});
                }
            }
        }
    } 
    
    //Vetor

    else {
        std::vector<bool> visited(n, false);
        for (int i = 0; i < n; ++i) {
            int u = -1;
            double min_dist = std::numeric_limits<double>::infinity();
            for (int j = 0; j < n; ++j) {
                if (!visited[j] && result.distance[j] < min_dist) {
                    min_dist = result.distance[j];
                    u = j;
                }
            }
            if (u == -1) break;
            visited[u] = true;

            for (const auto& edge : representation_->getNeighbors(u)) {
                int v = edge.target;
                double weight = edge.weight;
                if (result.distance[u] + weight < result.distance[v]) {
                    result.distance[v] = result.distance[u] + weight;
                    result.parent[v] = u;
                }
            }
        }
    }
    return result;
}

// Algoritmo de Bellman_Ford já otimizado
// From all to sink
SearchResult Graph::Bellman_Ford(int end_node) const{
    int n = getVertexCount();
    SearchResult result;
    result.distance.assign(n, std::numeric_limits<double>::infinity());
    result.parent.assign(n, -1); // Tem a termologia parent mas ele guarda o sucessor
    result.distance[end_node] = 0.0;
    for (int i = 0; i < n - 1; i++){
        bool k = false;
        for(int v = 0; v < n; v++){
            for (const auto& edge : representation_->getNeighbors(v)) {
                int w = edge.target;
                double weight = edge.weight;
                if (result.distance[w] + weight < result.distance[v]) {
                    result.distance[v] = result.distance[w] + weight;
                    result.parent[v] = w;
                    k = true;
                }
            }
        }
        if (k == false){
            break;
        }
    }
    for(int v = 0; v < n; v++){
            for (const auto& edge : representation_->getNeighbors(v)) {
                int w = edge.target;
                double weight = edge.weight;
                if (result.distance[w] + weight < result.distance[v]) {
                   throw std::runtime_error("O grafo possui ciclos negativos");
                }
            }
        }

    return result;
}
    // Outra versão do Bellman_Ford que funciona igual dijkstra ele calcula a distância de um vértice para os outros
SearchResult Graph::Bellman_Ford_reversed(int start_node) const{
    int n = getVertexCount();
    SearchResult result;
    result.distance.assign(n, std::numeric_limits<double>::infinity());
    result.parent.assign(n, -1); //parent
    result.distance[start_node] = 0.0;
    for (int i = 0; i < n - 1; i++){
        bool k = false;
        for(int j = 0; j < n; j++){
            for (const auto& edge : representation_->getNeighbors(j)) {
                int v = edge.target;
                double weight = edge.weight;
                if (result.distance[j] + weight < result.distance[v]) {
                    result.distance[v] = result.distance[j] + weight;
                    result.parent[v] = j;
                    k = true;
                }
            }
        }
        if (k == false){
            break;
        }
    }
    for(int j = 0; j < n; j++){
            for (const auto& edge : representation_->getNeighbors(j)) {
                int v = edge.target;
                double weight = edge.weight;
                if (result.distance[j] + weight < result.distance[v]) {
                   throw std::runtime_error("O grafo possui ciclos negativos");
                }
            }
        }

    return result;
}

SearchResult Graph::Floyd_Warshall() const {
    int n = getVertexCount();
    const double INF = std::numeric_limits<double>::infinity();
    SearchResult result;
    result.distMatrix.assign(n, std::vector<double>(n, INF));
    result.predMatrix.assign(n, std::vector<int>(n, -1));

    // 2. Distância zero para cada vértice para si mesmo
     for (int i = 0; i < n; i++){
        result.distMatrix[i][i] = 0.0;
        result.predMatrix[i][i] = i;
     }

    // Pesos das arestas
    for (int u = 0; u < n; u++) {
        for (const auto& edge : representation_->getNeighbors(u)) {
            int v = edge.target;
            result.distMatrix[u][v] = edge.weight;
            result.predMatrix[u][v] = u;
        }
    }

    // Floyd–Warshall
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (result.distMatrix[i][k] < INF && result.distMatrix[k][j] < INF) {
                    result.distMatrix[i][j] = std::min(result.distMatrix[i][j], result.distMatrix[i][k] + result.distMatrix[k][j]);
                    result.predMatrix[i][j] = result.predMatrix[k][j];
                }
            }
        }
    }

    for (int i = 0; i < n; i++){
        if (result.distMatrix[i][i] < 0){
            throw std::runtime_error("O grafo possui ciclos negativos");
        }
    }

    return result;
}


// Função para reconstruir o caminha por meio dos pais
std::vector<int> Graph::getPath(int target, int u) { //1 - based
    SearchResult result;
    result = this->Run_Search(u - 1);
    std::vector<int> path;
    int current = target - 1;

    if (result.distance[current] == std::numeric_limits<double>::infinity() || result.distance[current] < 0) {
        // Unreachable or error state (e.g., disconnected or negative cycle issue)
        return path; 
    }

    while (current != -1) {
        path.push_back(current);
        if (current == u - 1) { 
            break;
        }
        current = result.parent[current];
    }
    if (path.back() != u - 1) {
        return {};
    }

    std::reverse(path.begin(), path.end());
    
    // Path contains 0-based vertices. The function returns 0-based indices.
    return path;
}

SearchResult Graph::Run_Search(int start_node) const{
    SearchResult result;
    int s = start_node;
    if(!hasWeights()){
        result = this->bfs(s);
    }
    if (hasNegativeWeights()){
        result = this->Bellman_Ford_reversed(s);
    }
    if (!hasNegativeWeights()){
        result = this->dijkstra(s);
    }
    return result;

}


double Graph::getDistance(int s, int t) const {
    // s = start_node (1-based), t = end_node (1-based)
    
    // Run the appropriate search from the source node (u-1)
    SearchResult result = this->Run_Search(s - 1);
    
    // Return the distance to the target node (v-1)
    return result.distance[t - 1];
}


double Graph::getDiameter() const {
    int n = getVertexCount();
    double max_distance = 0.0;
    
    for (int i = 0; i < n; ++i) {
        try {
            SearchResult result = this->Run_Search(i); 

            for (double dist : result.distance) {
                // If any distance is infinity, the graph is disconnected (diameter is undefined/infinite)
                if (dist == std::numeric_limits<double>::infinity()) {
                    return -1.0;
                }
                if (dist > max_distance) {
                    max_distance = dist;
                }
            }
        } catch (const std::runtime_error& e) {
            // Catches negative cycle error from Bellman-Ford
            std::cerr << "Erro ao calcular diâmetro: " << e.what() << std::endl;
            return std::numeric_limits<double>::infinity(); // Saída de Erro 
        }
    }
    return max_distance;
}

double Graph::getApproximateDiameter() const {
    int n = getVertexCount();
    if (n < 2) return 0.0;

    // Search from an arbitrary starting node (0)
    SearchResult res1 = this->Run_Search(0);
    
    int u = 0;
    double max_dist = 0.0;
    for (int i = 0; i < n; ++i) {
        if (res1.distance[i] == std::numeric_limits<double>::infinity()) return -1.0; // Disconnected
        if (res1.distance[i] > max_dist) {
            max_dist = res1.distance[i];
            u = i; // Node furthest from the initial start
        }
    }

    // Search from the furthest node (u)
    SearchResult res2 = this->Run_Search(u);
    
    // Return the max distance found in the second search
    return *std::max_element(res2.distance.begin(), res2.distance.end());
}


// Funções Auxiliares 
std::vector<double> Graph::getDegreeStats() const {
    int n = getVertexCount();
    std::vector<double> degrees;
    degrees.reserve(n);
    for (int i = 0; i < n; ++i) {
        // Para calcular o grau
        degrees.push_back(representation_->getNeighbors(i).size());
    }
    
    if (degrees.empty()) return {0, 0, 0, 0};
    
    double min_degree = *std::min_element(degrees.begin(), degrees.end());
    double max_degree = *std::max_element(degrees.begin(), degrees.end());
    double sum = std::accumulate(degrees.begin(), degrees.end(), 0.0);
    double mean_degree = sum / degrees.size();
    std::sort(degrees.begin(), degrees.end());
    double median_degree = (degrees.size() % 2 == 1) ? 
                           degrees[degrees.size() / 2] : 
                           (degrees[degrees.size() / 2 - 1] + degrees[degrees.size() / 2]) / 2.0;
    
    return {min_degree, max_degree, mean_degree, median_degree};
}

bool Graph::writeResults(const std::string& output_filename) const {
    std::ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        std::cerr << "Erro ao abrir arquivo para escrita " << output_filename << std::endl;
        return false;
    }
    
    output_file << "Número de vértices: " << representation_->getVertexCount() << std::endl;
    output_file << "Número de arestas: " << representation_->getEdgeCount() << std::endl;
    
    std::vector<double> stats = this->getDegreeStats();
    output_file << "Grau mínimo: " << stats[0] << std::endl;
    output_file << "Grau máximo: " << stats[1] << std::endl;
    output_file << "Grau médio: " << stats[2] << std::endl;
    output_file << "Grau mediano: " << stats[3] << std::endl;
    
    output_file.close();
    return true;
}

void Graph::generateReport(int start_node, const std::string& output_filename, AlgorithmType algo_type) const {
    std::cout << "Gerando relatório para o nó " << start_node << "..." << std::endl;
    
    // Determine which algorithm to run
    SearchResult result;
    std::string algo_name;
    
    // Convert 1-based start_node to 0-based for internal functions
    int start_node_0based = start_node - 1;

    // Use a switch statement to select the algorithm
    switch (algo_type) {
        case AlgorithmType::BFS:
            result = this->bfs(start_node_0based);
            algo_name = "Busca em Largura (BFS)";
            break;
        case AlgorithmType::DFS:
            result = this->dfs(start_node_0based);
            algo_name = "Busca em Profundidade (DFS)";
            break;
        case AlgorithmType::DIJKSTRA_HEAP:
            result = this->dijkstra(start_node_0based, DijkstraImplType::HEAP);
            algo_name = "Dijkstra (Heap)";
            break;
        case AlgorithmType::DIJKSTRA_VECTOR:
            result = this->dijkstra(start_node_0based, DijkstraImplType::VECTOR);
            algo_name = "Dijkstra (Vector)";
            break;
        case AlgorithmType::BELLMAN_FORD_REVERSED:
            result = this->Bellman_Ford_reversed(start_node_0based);
            algo_name = "Bellman-Ford (Source-to-All)";
            break;
        case AlgorithmType::AUTO_SEARCH:
        default:
            result = this->Run_Search(start_node_0based);
            std::string algorithm_name_part = (hasWeights() ? (hasNegativeWeights() ? "Bellman-Ford" : "Dijkstra") : "BFS");
            algo_name = std::string("Busca Automática (") + algorithm_name_part + ")";
            break;
    }

    // Write the results to the file
    std::ofstream out(output_filename, std::ios_base::app);
    if (!out.is_open()) {
        throw std::runtime_error("Erro: não foi possível abrir o arquivo para escrita: " + output_filename);
    }

    out << "\n--- Árvore de " << algo_name << " a partir do Vértice " << start_node << " ---\n";
    
    const std::string parent_label = "Pai: ";
    
    for (size_t i = 0; i < result.parent.size(); ++i) {
        out << "Vértice: " << i + 1 << ", "; // 1-based vertex
        if (result.parent[i] != -1) {
            out << parent_label << result.parent[i] + 1 << ", "; // 1-based parent
        } else {
            out << parent_label << "Nenhum, ";
        }
        
        if (result.distance[i] == std::numeric_limits<double>::infinity()) {
            out << "Distância: INFINITO" << std::endl;
        } else {
            out << "Distância: " << result.distance[i] << std::endl;
        }
    }

    out.close();
    std::cout << "Relatório salvo em '" << output_filename << "'." << std::endl;
}

void Graph::generateBellman_FordReport(int end_node, const std::string& output_filename) const {
    std::cout << "Gerando relatório Bellman_Ford (De todos para o nó final) para o nó" << end_node << "..." << std::endl;

    SearchResult result = this->Bellman_Ford(end_node - 1);

    std::ofstream out(output_filename, std::ios_base::app);
    if (!out.is_open()) {
        throw std::runtime_error("Erro: não foi possível abrir o arquivo para escrita: " + output_filename);
    }

    out << "\n--- Árvore do Algoritmo de Bellman_Ford (Sink Tree) a partir do Vértice " << end_node << " ---\n";

    for (size_t i = 0; i < result.parent.size(); ++i) {
        out << "Vértice: " << i + 1 << ", "; // Exibe o vértice como 1-based
        if (result.parent[i] != -1) {
            out << "Sucessor: " << result.parent[i] + 1 << ", "; // Exibe o pai como 1-based
        } else {
            out << "Sucessor: Nenhum, ";
        }
        out << "Distância: " << result.distance[i] << std::endl;
    }

    out.close();
    std::cout << "Relatório salvo em '" << output_filename << "'." << std::endl;
}

std::vector<std::vector<int>> Graph::getConnectedComponents() const {
    int n = getVertexCount();
    std::vector<std::vector<int>> components;

    // Caso NÃO DIRECIONADO → componentes conexos normais
    if (direction_type_ == Direction_Graph::NO_Direction) {
        std::vector<bool> visited(n, false);

        for (int i = 0; i < n; ++i) {
            if (!visited[i]) {
                SearchResult r = this->bfs(i);
                std::vector<int> comp;

                for (int j = 0; j < n; ++j) {
                    if (r.distance[j] >= 0 && r.distance[j] != std::numeric_limits<double>::infinity() && !visited[j]) {
                        visited[j] = true;
                        comp.push_back(j);
                    }
                }
                components.push_back(comp);
            }
        }

        return components;
    }

    // Caso DIRECIONADO → SCC usando apenas BFS (slow but correct)

    Graph rev = reverseEdges();
    std::vector<bool> assigned(n, false);  // para não repetir SCCs

    for (int start = 0; start < n; ++start) {

        if (assigned[start]) continue;

        // BFS no grafo normal
        SearchResult fromStart = this->bfs(start);

        // BFS no grafo reverso
        SearchResult toStart = rev.bfs(start);

        std::vector<int> comp;

        // Interseção: vértices alcançados nos DOIS sentidos
        for (int v = 0; v < n; ++v) {
            if (!assigned[v] &&
                fromStart.distance[v] >= 0 && 
                toStart.distance[v] >= 0 &&
                fromStart.distance[v] != std::numeric_limits<double>::infinity() && 
                toStart.distance[v] != std::numeric_limits<double>::infinity()) {

                assigned[v] = true;
                comp.push_back(v);
            }
        }

        components.push_back(comp);
    }

    return components;
}

// E para as componentes conexas
void Graph::generateConnectedComponentsReport(const std::string& output_filename) const {
    auto components = getConnectedComponents();
    
    std::ofstream out(output_filename, std::ios_base::app);
    if (!out.is_open()) {
        throw std::runtime_error("Erro: não foi possível abrir o arquivo para escrita: " + output_filename);
    }

    out << "\n--- Componentes Conexas ---\n";
    out << "Número de componentes: " << components.size() << std::endl;

    std::sort(components.begin(), components.end(), [](const auto& a, const auto& b) {
        return a.size() > b.size();
    });

    if (!components.empty()){
        out << "Tamanho da maior componente: " << components.front().size() << std::endl;
        out << "Tamanho da menor componente: " << components.back().size() << std::endl;
    }

    for (const auto& component : components) {
        out << "\nComponentes:\nTamanho: " << component.size() << " | Vertices: ";
        for (int v : component) {
            out << v + 1 << " ";
        }
        out << std::endl;
    }

    out.close();
}

Graph Graph::reverseEdges() const {
    Graph reversed; 
    //Inverte o grafo em lista (por enquanto só em lista)
    reversed.representation_ = std::make_unique<AdjacencyList>(
        this->getVertexCount(),
        Direction_Graph::Direction // Mantém a direção no grafo reverso
    );
    
    reversed.has_weights = this->hasWeights(); 
    reversed.has_negative_weights_ = this->hasNegativeWeights(); 

    // Para cada aresta u -> v, insere v -> u
    for (int u = 0; u < getVertexCount(); u++) {
        for (const Edge &e : representation_->getNeighbors(u)) {
            // Insere a aresta invertida (e.target -> u) no novo grafo
            reversed.representation_->addEdge(e.target, u, e.weight);
        }
    }

    return reversed;
}

}




