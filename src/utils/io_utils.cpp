#include "../../include/utils/io_utils.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <filesystem>
#include <stdexcept>

namespace lbm {
namespace utils {

bool IOUtils::createDirectory(const std::string& directory_path) {
    try {
        std::filesystem::create_directories(directory_path);
        return std::filesystem::exists(directory_path);
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating directory " << directory_path << ": " << e.what() << std::endl;
        return false;
    }
}

std::string IOUtils::generateTimestepFilename(const std::string& base_name, 
                                             int timestep, 
                                             const std::string& extension) {
    std::ostringstream oss;
    oss << base_name << "_t" << timestep << extension;
    return oss.str();
}

void IOUtils::writeCSVHeader(std::ofstream& file, const std::vector<std::string>& headers) {
    if (!file.is_open()) {
        throw std::runtime_error("File is not open for writing");
    }
    
    for (size_t i = 0; i < headers.size(); ++i) {
        file << headers[i];
        if (i < headers.size() - 1) {
            file << ",";
        }
    }
    file << "\n";
}

void IOUtils::writeCSVRow(std::ofstream& file, const std::vector<double>& data) {
    if (!file.is_open()) {
        throw std::runtime_error("File is not open for writing");
    }
    
    for (size_t i = 0; i < data.size(); ++i) {
        file << formatDouble(data[i]);
        if (i < data.size() - 1) {
            file << ",";
        }
    }
    file << "\n";
}

void IOUtils::writeFieldToCSV(const std::string& filename, 
                              double** field, 
                              int Nx, int Ny,
                              const std::string& field_name,
                              bool include_coordinates) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    // Write header
    std::vector<std::string> headers;
    if (include_coordinates) {
        headers.push_back("x");
        headers.push_back("y");
    }
    headers.push_back(field_name);
    writeCSVHeader(file, headers);
    
    // Write data
    for (int y = 0; y < Ny; ++y) {
        for (int x = 0; x < Nx; ++x) {
            std::vector<double> row_data;
            if (include_coordinates) {
                row_data.push_back(static_cast<double>(x));
                row_data.push_back(static_cast<double>(y));
            }
            row_data.push_back(field[x][y]);
            writeCSVRow(file, row_data);
        }
    }
    
    file.close();
}

void IOUtils::writeVelocityProfile(const std::string& filename,
                                  const std::vector<int>& y_positions,
                                  const std::vector<double>& velocities,
                                  const std::vector<double>& analytical_velocities) {
    if (y_positions.size() != velocities.size()) {
        throw std::invalid_argument("Y positions and velocities must have same size");
    }
    
    if (!analytical_velocities.empty() && analytical_velocities.size() != velocities.size()) {
        throw std::invalid_argument("Analytical velocities size mismatch");
    }
    
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    // Write header
    std::vector<std::string> headers = {"y", "velocity"};
    if (!analytical_velocities.empty()) {
        headers.push_back("analytical_velocity");
    }
    writeCSVHeader(file, headers);
    
    // Write data
    for (size_t i = 0; i < y_positions.size(); ++i) {
        std::vector<double> row_data = {
            static_cast<double>(y_positions[i]), 
            velocities[i]
        };
        if (!analytical_velocities.empty()) {
            row_data.push_back(analytical_velocities[i]);
        }
        writeCSVRow(file, row_data);
    }
    
    file.close();
}

void IOUtils::writeMetadata(const std::string& filename,
                           const std::vector<std::pair<std::string, std::string>>& metadata) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for writing: " + filename);
    }
    
    file << "# Simulation Metadata\n";
    for (const auto& pair : metadata) {
        file << pair.first << "=" << pair.second << "\n";
    }
    
    file.close();
}

std::vector<std::vector<double>> IOUtils::readCSV(const std::string& filename, bool skip_header) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file for reading: " + filename);
    }
    
    std::vector<std::vector<double>> data;
    std::string line;
    bool first_line = true;
    
    while (std::getline(file, line)) {
        if (first_line && skip_header) {
            first_line = false;
            continue;
        }
        
        std::vector<std::string> tokens = split(line, ',');
        std::vector<double> row;
        
        for (const std::string& token : tokens) {
            try {
                row.push_back(stringToDouble(token));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Warning: Could not parse value '" << token << "' in file " << filename << std::endl;
                row.push_back(0.0);  // Use default value
            }
        }
        
        data.push_back(row);
        first_line = false;
    }
    
    file.close();
    return data;
}

bool IOUtils::fileExists(const std::string& filename) {
    return std::filesystem::exists(filename);
}

long IOUtils::getFileSize(const std::string& filename) {
    try {
        return std::filesystem::file_size(filename);
    } catch (const std::filesystem::filesystem_error&) {
        return -1;
    }
}

bool IOUtils::createBackup(const std::string& filename) {
    if (!fileExists(filename)) {
        return false;
    }
    
    std::string backup_filename = filename + ".backup";
    try {
        std::filesystem::copy_file(filename, backup_filename, 
                                  std::filesystem::copy_options::overwrite_existing);
        return true;
    } catch (const std::filesystem::filesystem_error& e) {
        std::cerr << "Error creating backup: " << e.what() << std::endl;
        return false;
    }
}

bool IOUtils::validateOutputDirectory(const std::string& output_dir) {
    if (!createDirectory(output_dir)) {
        return false;
    }
    
    // Test write permissions by creating a temporary file
    std::string test_file = output_dir + "/test_write_permissions.tmp";
    std::ofstream test(test_file);
    if (!test.is_open()) {
        return false;
    }
    test.close();
    
    // Clean up test file
    try {
        std::filesystem::remove(test_file);
    } catch (const std::filesystem::filesystem_error&) {
        // Ignore cleanup errors
    }
    
    return true;
}

std::string IOUtils::formatDouble(double value, int precision) {
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(precision) << value;
    return oss.str();
}

void IOUtils::printProgress(int current_step, int total_steps, const std::string& additional_info) {
    double percentage = 100.0 * current_step / total_steps;
    
    std::cout << "\rProgress: " << std::fixed << std::setprecision(1) << percentage << "% "
              << "(" << current_step << "/" << total_steps << ")";
    
    if (!additional_info.empty()) {
        std::cout << " - " << additional_info;
    }
    
    std::cout << std::flush;
    
    if (current_step == total_steps) {
        std::cout << std::endl;  // New line when complete
    }
}

std::vector<std::string> IOUtils::split(const std::string& str, char delimiter) {
    std::vector<std::string> tokens;
    std::stringstream ss(str);
    std::string token;
    
    while (std::getline(ss, token, delimiter)) {
        tokens.push_back(token);
    }
    
    return tokens;
}

double IOUtils::stringToDouble(const std::string& str) {
    try {
        return std::stod(str);
    } catch (const std::invalid_argument& e) {
        throw std::invalid_argument("Cannot convert '" + str + "' to double");
    } catch (const std::out_of_range& e) {
        throw std::out_of_range("Value '" + str + "' is out of range for double");
    }
}

} // namespace utils
} // namespace lbm