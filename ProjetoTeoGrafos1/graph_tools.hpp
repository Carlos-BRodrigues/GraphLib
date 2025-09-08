#ifndef GRAPH_TOOLS_HPP
#define GRAPH_TOOLS_HPP

#include <string>
#include <vector>
#include <memory> // ESSENCIAL: para std::unique_ptr

namespace graph_tools_lib {

// MUDANÇA 1: Adicionar um enum para a escolha do usuário
enum class RepresentationType {
    ADJACENCY_LIST,
    ADJACENCY_MATRIX
};

// MUDANÇA 2: Definir a interface abstrata (o blueprint)
class IGraphRepresentation {
public:
    virtual ~IGraphRepresentation() = default;
    virtual void addEdge(int u, int v) = 0;
    virtual int getVertexCount() const = 0;
    virtual int getEdgeCount() const = 0;
    virtual int getDegree(int v) const = 0;
    virtual void print() const = 0;
    // Adicione outras funções essenciais aqui se precisar (ex: getNeighbors)
};

// MUDANÇA 3: Declarar as classes concretas que implementam a interface
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
    void print() const override;
};


// MUDANÇA 4: Simplificar drasticamente a classe Graph
class Graph {
public:
    // O construtor agora aceita a escolha do tipo
    Graph(const std::string& filename, RepresentationType type);

    // As funções públicas não mudam para o usuário!
    std::vector<double> getDegreeStats() const;
    bool writeResults(const std::string& output_filename) const;
    void print() const;

private:
    // REMOVER os membros antigos:
    // int num_vertices_;
    // int num_edges_;
    // std::vector<std::vector<int>> adj_list_;

    // ADICIONAR um único ponteiro inteligente para a representação
    std::unique_ptr<IGraphRepresentation> representation_;
};

} // namespace

#endif