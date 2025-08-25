#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>

int main() {
    // Input and output file paths
    const std::string input_file = "wp_pwgevents_14058.lhe";
    const std::string output_file = "testOutput.out";
    const std::string weight_file = "weights.txt";
    const std::string event_line_file = "event_lines.txt";

    // Read weights
    std::vector<std::string> weights;
    {
        std::ifstream weight_stream(weight_file);
        std::string weight;
        while (std::getline(weight_stream, weight)) {
            weights.push_back(weight);
        }
    }

    // Read event lines
    std::vector<int> event_lines;
    {
        std::ifstream event_line_stream(event_line_file);
        int line_number;
        while (event_line_stream >> line_number) {
            event_lines.push_back(line_number);
        }
    }

    // Process input file
    std::ifstream input_stream(input_file);
    std::ofstream output_stream(output_file);
    std::string line;
    int current_line_number = 0;

    while (std::getline(input_stream, line)) {
        current_line_number++;

        // Check if this is a line to modify
        auto it = std::find(event_lines.begin(), event_lines.end(), current_line_number - 1);
        if (it != event_lines.end()) {
            // Calculate index for weight replacement
            size_t index = std::distance(event_lines.begin(), it);
            
            // Split the line
            std::istringstream iss(line);
            std::vector<std::string> tokens;
            std::string token;
            
            while (iss >> token) {
                tokens.push_back(token);
            }

            // Replace third token (index 2)
            if (tokens.size() >= 3) {
                tokens[2] = weights[index];
            }

            // Reconstruct line
            std::ostringstream oss;
            for (const auto& t : tokens) {
                oss << t << " ";
            }
            line = oss.str();
        }

        output_stream << line << std::endl;
    }

    return 0;
}

