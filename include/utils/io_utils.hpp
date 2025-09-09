#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <memory>

namespace lbm {
namespace utils {

/**
 * @brief Input/Output utility functions for LBM simulations
 */
class IOUtils {
public:
    /**
     * @brief Create output directory if it doesn't exist
     * @param directory_path Path to directory to create
     * @return true if directory exists or was created successfully
     */
    static bool createDirectory(const std::string& directory_path);

    /**
     * @brief Generate timestamped filename
     * @param base_name Base filename without extension
     * @param timestep Current timestep
     * @param extension File extension (with dot)
     * @return Formatted filename
     */
    static std::string generateTimestepFilename(const std::string& base_name, 
                                               int timestep, 
                                               const std::string& extension = ".csv");

    /**
     * @brief Write CSV header to file
     * @param file Output file stream
     * @param headers Vector of header names
     */
    static void writeCSVHeader(std::ofstream& file, const std::vector<std::string>& headers);

    /**
     * @brief Write CSV data row to file
     * @param file Output file stream
     * @param data Vector of data values
     */
    static void writeCSVRow(std::ofstream& file, const std::vector<double>& data);

    /**
     * @brief Write 2D field data to CSV file
     * @param filename Output filename
     * @param field 2D field data [Nx][Ny]
     * @param Nx Grid size in x-direction
     * @param Ny Grid size in y-direction
     * @param field_name Name of the field for header
     * @param include_coordinates Whether to include x,y coordinates
     */
    static void writeFieldToCSV(const std::string& filename, 
                                double** field, 
                                int Nx, int Ny,
                                const std::string& field_name,
                                bool include_coordinates = true);

    /**
     * @brief Write velocity profile data to CSV file
     * @param filename Output filename
     * @param y_positions Y-coordinate positions
     * @param velocities Velocity values
     * @param analytical_velocities Analytical velocity values (optional)
     */
    static void writeVelocityProfile(const std::string& filename,
                                    const std::vector<int>& y_positions,
                                    const std::vector<double>& velocities,
                                    const std::vector<double>& analytical_velocities = {});

    /**
     * @brief Write simulation metadata to file
     * @param filename Output filename
     * @param metadata Key-value pairs of simulation parameters
     */
    static void writeMetadata(const std::string& filename,
                             const std::vector<std::pair<std::string, std::string>>& metadata);

    /**
     * @brief Read CSV file into vector of vectors
     * @param filename Input filename
     * @param skip_header Whether to skip the first row
     * @return Vector of rows, each row is a vector of values
     */
    static std::vector<std::vector<double>> readCSV(const std::string& filename, 
                                                   bool skip_header = true);

    /**
     * @brief Check if file exists
     * @param filename File path to check
     * @return true if file exists and is readable
     */
    static bool fileExists(const std::string& filename);

    /**
     * @brief Get file size in bytes
     * @param filename File path
     * @return File size in bytes, -1 if file doesn't exist
     */
    static long getFileSize(const std::string& filename);

    /**
     * @brief Create backup of existing file
     * @param filename Original filename
     * @return true if backup was created successfully
     */
    static bool createBackup(const std::string& filename);

    /**
     * @brief Validate output directory and create if necessary
     * @param output_dir Output directory path
     * @return true if directory is valid and writable
     */
    static bool validateOutputDirectory(const std::string& output_dir);

    /**
     * @brief Format double value for CSV output with specified precision
     * @param value Value to format
     * @param precision Number of decimal places
     * @return Formatted string
     */
    static std::string formatDouble(double value, int precision = 6);

    /**
     * @brief Write simulation progress to console
     * @param current_step Current timestep
     * @param total_steps Total number of timesteps
     * @param additional_info Additional information to display
     */
    static void printProgress(int current_step, int total_steps, 
                             const std::string& additional_info = "");

private:
    /**
     * @brief Split string by delimiter
     * @param str String to split
     * @param delimiter Delimiter character
     * @return Vector of split strings
     */
    static std::vector<std::string> split(const std::string& str, char delimiter);

    /**
     * @brief Convert string to double with error checking
     * @param str String to convert
     * @return Converted double value
     * @throws std::invalid_argument if conversion fails
     */
    static double stringToDouble(const std::string& str);
};

} // namespace utils
} // namespace lbm