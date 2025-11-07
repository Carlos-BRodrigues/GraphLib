#ifndef GRAPH_TOOLS_HPP
#define GRAPH_TOOLS_HPP

#include <string>
#include <vector>
#include <memory>

namespace graph_tools_lib {

enum class RepresentationType {
    ADJACENCY_LIST,
    ADJACENCY_MATRIX
};

// Enum para escolher a implementação de Dijkstra.
enum class DijkstraImplType {
    VECTOR, // Implementação O(V^2)
    HEAP    // Implementação O(E log V)
};

// Paramêtro que define se um grafo é direcionado
enum class Direction_Graph{
    Direction,
    NO_Direction
};

struct Edge {
    int target;
    double weight;
};

struct SearchResult {
    std::vector<int> parent;
    std::vector<double> distance;
};

class IGraphRepresentation {
public:
    virtual ~IGraphRepresentation() = default;
    virtual void addEdge(int u, int v, double weight) = 0;
    virtual int getVertexCount() const = 0;
    virtual int getEdgeCount() const = 0;
    virtual std::vector<Edge> getNeighbors(int v) const = 0;
    virtual void print() const = 0;
};

class AdjacencyList : public IGraphRepresentation {
private:
    int num_vertices_;
    int num_edges_;
    std::vector<std::vector<Edge>> adj_list_;
    Direction_Graph direction_type_;
public:
    AdjacencyList(int num_vertices, Direction_Graph direction);
    void addEdge(int u, int v, double weight) override;
    int getVertexCount() const override;
    int getEdgeCount() const override;
    std::vector<Edge> getNeighbors(int v) const override;
    void print() const override;
};

class AdjacencyMatrix : public IGraphRepresentation {
private:
    int num_vertices_;
    int num_edges_;
    std::vector<std::vector<double>> adj_matrix_;
    Direction_Graph direction_type_;
public:
    AdjacencyMatrix(int num_vertices, Direction_Graph direction);
    void addEdge(int u, int v, double weight) override;
    int getVertexCount() const override;
    int getEdgeCount() const override;
    std::vector<Edge> getNeighbors(int v) const override;
    void print() const override;
};

class Graph {
public:
    Graph(const std::string& filename, RepresentationType type, Direction_Graph direction);

    int getVertexCount() const;
    int getEdgeCount() const;
    std::vector<double> getDegreeStats() const;
    bool writeResults(const std::string& output_filename) const;
    void print() const;
    bool hasNegativeWeights() const;

    // Algoritmos
    SearchResult bfs(int start_node) const;
    SearchResult dfs(int start_node) const;
    SearchResult dijkstra(int start_node, DijkstraImplType impl_type = DijkstraImplType::HEAP) const;
    SearchResult Bellman_Ford(int start_node) const;

    std::vector<std::vector<int>> getConnectedComponents() const;

    std::vector<int> getPath(int target, int v);
    
    double getDistance(int u, int v) const;
    double getDiameter() const;
    double getApproximateDiameter() const;

    // Funções de relatório
    void generateBfsReport(int start_node, const std::string& output_filename) const;
    void generateDfsReport(int start_node, const std::string& output_filename) const;
    void generateDijkstraReport(int start_node, const std::string& output_filename, DijkstraImplType impl_type = DijkstraImplType::HEAP) const;
    void generateBellman_FordReport(int start_node, const std::string& output_filename) const;
    void generateConnectedComponentsReport(const std::string& output_filename) const;

private:
    std::unique_ptr<IGraphRepresentation> representation_;
    bool has_negative_weights_;
};

} // namespace graph_tools_lib

#endif

