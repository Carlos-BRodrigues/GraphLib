#ifndef GRAPH_TOOLS_HPP
#define GRAPH_TOOLS_HPP

#include <string>
#include <vector>

namespace graph_tools_lib {

class Graph {
public:
    // Constructor: Creates and initializes a Graph object directly from a file.
    explicit Graph(const std::string& filename);

    // Member function to get degree statistics.
    // 'const' promises this function won't change the graph.
    std::vector<double> getDegreeStats() const;

    // Member function to write the results to a file.
    bool writeResults(const std::string& output_filename) const;

    void print() const;

private:
    // Private member variables to STORE the graph's data.
    // The trailing underscore is a common convention for private members.
    int num_vertices_;
    int num_edges_;
    std::vector<std::vector<int>> adj_list_;

    // A private helper function to do the actual file parsing.
    void parseFile(const std::string& filename);
};

} // namespace graph_tools_lib

#endif // GRAPH_TOOLS_HPP