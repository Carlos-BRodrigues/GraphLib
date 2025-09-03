#include <iostream>
#include <vector>
#include "graph_tools.hpp"

int main() {
    std::cout << "--- Program Started ---" << std::endl; // <-- ADD THIS

    try {
        graph_tools_lib::Graph my_graph("test.txt");
        my_graph.print();
        my_graph.writeResults("saida.txt");

    } catch (const std::exception& e) {
        std::cerr << "An error occurred: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "--- Program Finished ---" << std::endl; // <-- AND THIS
    return 0;
}