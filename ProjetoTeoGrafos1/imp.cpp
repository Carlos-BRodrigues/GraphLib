#include "graph_tools.hpp"
#include <fstream>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <stdexcept>

namespace graph_tools_lib {

// The constructor's job is to call the parser to set up the object.
Graph::Graph(const std::string& filename) {
    parseFile(filename);
}

// Private helper to populate the object's member variables.
void Graph::parseFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Error: Could not open file " + filename);
    }

    // 1. Read ONLY the number of vertices from the first line.
    file >> this->num_vertices_;
    this->num_edges_ = 0; // Initialize edge count to zero.

    // 2. Initialize the adjacency list to the correct size.
    this->adj_list_.assign(num_vertices_, std::vector<int>());

    int node1, node2;
    // 3. Loop until the end of the file, reading each edge.
    while (file >> node1 >> node2) {
        // Store as 0-based indices.
        // Also, add checks to prevent crashing if the file contains invalid vertex numbers.
        if (node1 > 0 && node1 <= num_vertices_ && node2 > 0 && node2 <= num_vertices_) {
            adj_list_[node1 - 1].push_back(node2 - 1);
            adj_list_[node2 - 1].push_back(node1 - 1);
            this->num_edges_++; // Increment the edge count for each valid edge read.
        }
    }
    file.close();
}

// 'getDegreeStats' is now a member function. It's 'const' because it only reads member data.
std::vector<double> Graph::getDegreeStats() const {
    std::vector<double> degrees;
    degrees.reserve(this->num_vertices_);
    for (const auto& neighbors : this->adj_list_) {
        degrees.push_back(neighbors.size());
    }

    if (degrees.empty()) {
        return {0, 0, 0, 0}; // Return zero for all stats if graph is empty.
    }
    
    double min_degree = *std::min_element(degrees.begin(), degrees.end());
    double max_degree = *std::max_element(degrees.begin(), degrees.end());
    double sum = std::accumulate(degrees.begin(), degrees.end(), 0.0);
    double mean_degree = sum / degrees.size();

    std::sort(degrees.begin(), degrees.end());
    
    double median_degree;
    auto n = degrees.size();
    if (n % 2 == 1) {
        median_degree = degrees[n / 2];
    } else {
        median_degree = (degrees[n / 2 - 1] + degrees[n / 2]) / 2.0;
    }
    
    return {min_degree, max_degree, mean_degree, median_degree};
}

// This function now has a single responsibility: writing output.
bool Graph::writeResults(const std::string& output_filename) const {
    std::ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        std::cerr << "Error opening file for writing: " << output_filename << std::endl;
        return false;
    }
    
    output_file << "Number of vertices: " << this->num_vertices_ << std::endl;
    output_file << "Number of edges: " << this->num_edges_ << std::endl;
    
    std::vector<double> stats = this->getDegreeStats();
    output_file << "Min degree: " << stats[0] << std::endl;
    output_file << "Max degree: " << stats[1] << std::endl;
    output_file << "Mean degree: " << stats[2] << std::endl;
    output_file << "Median degree: " << stats[3] << std::endl;
    
    output_file.close();
    return true;
}

    void Graph::print() const {
    std::cout << "Graph Adjacency List:" << std::endl;
    std::cout << "Vertices: " << this->num_vertices_ << std::endl;
    std::cout << "Edges: " << this->num_edges_ << std::endl;
    std::cout << "--------------------" << std::endl;

    // Loop through each vertex from 0 to num_vertices_ - 1
    for (int i = 0; i < this->num_vertices_; ++i) {
        // We print the vertex number (i + 1 for 1-based display)
        std::cout << i + 1 << " => ";

        // Directly access the vector of neighbors for the current vertex 'i'
        for (int neighbor : this->adj_list_[i]) {
            // We print each neighbor (neighbor + 1 for 1-based display)
            std::cout << neighbor + 1 << " ";
        }
        std::cout << std::endl;
    }
}

} // namespace graph_tools_lib