#ifndef GRAPH_TOOLS_HPP
#define GRAPH_TOOLS_HPP

#include <string>
#include <vector>
#include <memory>

namespace graph_tools_lib {

//Usuário escolhe representação, mais elegante eu imagino do que uma classe só com tudo dentro
enum class RepresentationType {
    ADJACENCY_LIST,
    ADJACENCY_MATRIX
};

struct SearchResult {
    std::vector<int> parent;
    std::vector<int> distance; // Renomeado de 'level' para 'distance' para clareza
};

//Interface abstrata
class IGraphRepresentation {
public:
    virtual ~IGraphRepresentation() = default;
    virtual void addEdge(int u, int v) = 0;
    virtual int getVertexCount() const = 0;
    virtual int getEdgeCount() const = 0;
    virtual int getDegree(int v) const = 0;
    
    //Métodos de busca
    virtual void BFS(int start_node, std::vector<int>& parent, std::vector<int>& level) const = 0;
    virtual void DFS(int start_node, std::vector<int>& parent, std::vector<int>& level) const = 0;
    
    virtual SearchResult bfs(int start_node) const = 0;
    virtual SearchResult dfs(int start_node) const = 0;

    //Diâmetro e componentes
    virtual int getDistance(int u, int v) const = 0;
    virtual int getDiameter() const = 0;
    virtual std::vector<std::vector<int>> getConnectedComponents() const = 0;
    virtual void print() const = 0;
};

//Classes concretas
class AdjacencyList : public IGraphRepresentation {
private:
    int num_vertices_;
    int num_edges_;
    std::vector<std::vector<int>> adj_list_;
public:
    AdjacencyList(int num_vertices);
    void addEdge(int u, int v) override;
    int getVertexCount() const override;
    int getEdgeCount() const override;
    int getDegree(int v) const override;

    void BFS(int start_node, std::vector<int>& parent, std::vector<int>& level) const override;
    void DFS(int start_node, std::vector<int>& parent, std::vector<int>& level) const override;
    SearchResult bfs(int start_node) const override;
    SearchResult dfs(int start_node) const override;
    
    int getDistance(int u, int v) const override;
    int getDiameter() const override;
    std::vector<std::vector<int>> getConnectedComponents() const override;
    void print() const override;
};

class AdjacencyMatrix : public IGraphRepresentation {
private:
    int num_vertices_;
    int num_edges_;
    std::vector<std::vector<int>> adj_matrix_;
public:
    AdjacencyMatrix(int num_vertices);
    void addEdge(int u, int v) override;
    int getVertexCount() const override;
    int getEdgeCount() const override;
    int getDegree(int v) const override;
    
    void BFS(int start_node, std::vector<int>& parent, std::vector<int>& level) const override;
    void DFS(int start_node, std::vector<int>& parent, std::vector<int>& level) const override;

    SearchResult bfs(int start_node) const override;
    SearchResult dfs(int start_node) const override;
    
    int getDistance(int u, int v) const override;
    int getDiameter() const override;
    std::vector<std::vector<int>> getConnectedComponents() const override;
    void print() const override;
};


//Classe Graph
class Graph {
public:
    //Construtor recebe arquivo e tipo de representação
    Graph(const std::string& filename, RepresentationType type);

    //Funções públicas
    std::vector<double> getDegreeStats() const;
    bool writeResults(const std::string& output_filename) const;
    int getVertexCount() const;
    int getEdgeCount() const;
    void BFS(int start_node, const std::string& output_file) const;
    void DFS(int start_node, const std::string& output_file) const;
    SearchResult bfs(int start_node) const;
    SearchResult dfs(int start_node) const;
    void writeSearchTree(const std::string& output_file, const std::string& algorithm, int start_node, const std::vector<int>& parent, const std::vector<int>& level) const;
    int getDistance(int u, int v) const;
    int getDiameter() const;
    std::vector<std::vector<int>> getConnectedComponents() const;
    bool writeConnectedComponents(const std::string& output_filename) const;
    void print() const;
private:
    //Ponteiro inteligente para a representação
    std::unique_ptr<IGraphRepresentation> representation_;
};

}

#endif
