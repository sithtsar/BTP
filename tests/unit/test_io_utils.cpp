#include <gtest/gtest.h>
#include <fstream>
#include <filesystem>
#include <vector>
#include "../../include/utils/io_utils.hpp"

using namespace lbm::utils;

class IOUtilsTest : public ::testing::Test {
protected:
    void SetUp() override {
        test_dir = "test_output";
        test_file = test_dir + "/test_file.csv";
        
        // Clean up any existing test files
        if (std::filesystem::exists(test_dir)) {
            std::filesystem::remove_all(test_dir);
        }
    }
    
    void TearDown() override {
        // Clean up test files
        if (std::filesystem::exists(test_dir)) {
            std::filesystem::remove_all(test_dir);
        }
    }
    
    std::string test_dir;
    std::string test_file;
};

TEST_F(IOUtilsTest, CreateDirectory) {
    // Test creating new directory
    EXPECT_TRUE(IOUtils::createDirectory(test_dir));
    EXPECT_TRUE(std::filesystem::exists(test_dir));
    
    // Test creating existing directory (should still return true)
    EXPECT_TRUE(IOUtils::createDirectory(test_dir));
    
    // Test creating nested directory
    std::string nested_dir = test_dir + "/nested/deep";
    EXPECT_TRUE(IOUtils::createDirectory(nested_dir));
    EXPECT_TRUE(std::filesystem::exists(nested_dir));
}

TEST_F(IOUtilsTest, GenerateTimestepFilename) {
    std::string filename = IOUtils::generateTimestepFilename("velocity", 1000, ".csv");
    EXPECT_EQ(filename, "velocity_t1000.csv");
    
    std::string filename_no_ext = IOUtils::generateTimestepFilename("data", 500);
    EXPECT_EQ(filename_no_ext, "data_t500.csv");
}

TEST_F(IOUtilsTest, WriteAndReadCSV) {
    // Create test directory
    ASSERT_TRUE(IOUtils::createDirectory(test_dir));
    
    // Test writing CSV with headers and data
    std::ofstream file(test_file);
    ASSERT_TRUE(file.is_open());
    
    std::vector<std::string> headers = {"x", "y", "value"};
    IOUtils::writeCSVHeader(file, headers);
    
    std::vector<double> row1 = {1.0, 2.0, 3.5};
    std::vector<double> row2 = {4.0, 5.0, 6.7};
    IOUtils::writeCSVRow(file, row1);
    IOUtils::writeCSVRow(file, row2);
    
    file.close();
    
    // Test reading CSV
    auto data = IOUtils::readCSV(test_file, true);  // Skip header
    ASSERT_EQ(data.size(), 2);
    ASSERT_EQ(data[0].size(), 3);
    ASSERT_EQ(data[1].size(), 3);
    
    EXPECT_DOUBLE_EQ(data[0][0], 1.0);
    EXPECT_DOUBLE_EQ(data[0][1], 2.0);
    EXPECT_DOUBLE_EQ(data[0][2], 3.5);
    EXPECT_DOUBLE_EQ(data[1][0], 4.0);
    EXPECT_DOUBLE_EQ(data[1][1], 5.0);
    EXPECT_DOUBLE_EQ(data[1][2], 6.7);
}

TEST_F(IOUtilsTest, WriteFieldToCSV) {
    // Create test directory
    ASSERT_TRUE(IOUtils::createDirectory(test_dir));
    
    // Create test field
    int Nx = 3, Ny = 3;
    double** field = new double*[Nx];
    for (int x = 0; x < Nx; x++) {
        field[x] = new double[Ny];
        for (int y = 0; y < Ny; y++) {
            field[x][y] = x + y * 0.1;
        }
    }
    
    // Write field to CSV
    std::string field_file = test_dir + "/field.csv";
    IOUtils::writeFieldToCSV(field_file, field, Nx, Ny, "test_field", true);
    
    // Verify file exists and has correct content
    EXPECT_TRUE(IOUtils::fileExists(field_file));
    
    auto data = IOUtils::readCSV(field_file, true);  // Skip header
    EXPECT_EQ(data.size(), Nx * Ny);
    EXPECT_EQ(data[0].size(), 3);  // x, y, field_value
    
    // Clean up
    for (int x = 0; x < Nx; x++) {
        delete[] field[x];
    }
    delete[] field;
}

TEST_F(IOUtilsTest, WriteVelocityProfile) {
    // Create test directory
    ASSERT_TRUE(IOUtils::createDirectory(test_dir));
    
    std::vector<int> y_positions = {0, 1, 2, 3, 4};
    std::vector<double> velocities = {0.0, 0.1, 0.2, 0.1, 0.0};
    std::vector<double> analytical = {0.0, 0.12, 0.18, 0.12, 0.0};
    
    std::string profile_file = test_dir + "/velocity_profile.csv";
    IOUtils::writeVelocityProfile(profile_file, y_positions, velocities, analytical);
    
    // Verify file exists and has correct content
    EXPECT_TRUE(IOUtils::fileExists(profile_file));
    
    auto data = IOUtils::readCSV(profile_file, true);  // Skip header
    EXPECT_EQ(data.size(), y_positions.size());
    EXPECT_EQ(data[0].size(), 3);  // y, velocity, analytical_velocity
    
    // Test without analytical velocities
    std::string profile_file2 = test_dir + "/velocity_profile2.csv";
    IOUtils::writeVelocityProfile(profile_file2, y_positions, velocities);
    
    auto data2 = IOUtils::readCSV(profile_file2, true);
    EXPECT_EQ(data2.size(), y_positions.size());
    EXPECT_EQ(data2[0].size(), 2);  // y, velocity only
}

TEST_F(IOUtilsTest, WriteMetadata) {
    // Create test directory
    ASSERT_TRUE(IOUtils::createDirectory(test_dir));
    
    std::vector<std::pair<std::string, std::string>> metadata = {
        {"grid_size", "256x64"},
        {"tau", "1.0"},
        {"gravity", "1e-6"},
        {"max_timesteps", "10000"}
    };
    
    std::string metadata_file = test_dir + "/metadata.txt";
    IOUtils::writeMetadata(metadata_file, metadata);
    
    // Verify file exists
    EXPECT_TRUE(IOUtils::fileExists(metadata_file));
    
    // Read and verify content
    std::ifstream file(metadata_file);
    std::string line;
    std::getline(file, line);
    EXPECT_EQ(line, "# Simulation Metadata");
    
    std::getline(file, line);
    EXPECT_EQ(line, "grid_size=256x64");
    
    file.close();
}

TEST_F(IOUtilsTest, FileOperations) {
    // Test file existence check
    EXPECT_FALSE(IOUtils::fileExists("nonexistent_file.txt"));
    
    // Create test file
    ASSERT_TRUE(IOUtils::createDirectory(test_dir));
    std::ofstream test_file_stream(test_file);
    test_file_stream << "test content\n";
    test_file_stream.close();
    
    // Test file existence
    EXPECT_TRUE(IOUtils::fileExists(test_file));
    
    // Test file size
    long size = IOUtils::getFileSize(test_file);
    EXPECT_GT(size, 0);
    
    // Test backup creation
    EXPECT_TRUE(IOUtils::createBackup(test_file));
    EXPECT_TRUE(IOUtils::fileExists(test_file + ".backup"));
    
    // Test backup of nonexistent file
    EXPECT_FALSE(IOUtils::createBackup("nonexistent_file.txt"));
}

TEST_F(IOUtilsTest, ValidateOutputDirectory) {
    // Test creating and validating new directory
    EXPECT_TRUE(IOUtils::validateOutputDirectory(test_dir));
    EXPECT_TRUE(std::filesystem::exists(test_dir));
    
    // Test validating existing directory
    EXPECT_TRUE(IOUtils::validateOutputDirectory(test_dir));
}

TEST_F(IOUtilsTest, FormatDouble) {
    EXPECT_EQ(IOUtils::formatDouble(3.14159, 2), "3.14");
    EXPECT_EQ(IOUtils::formatDouble(3.14159, 4), "3.1416");
    EXPECT_EQ(IOUtils::formatDouble(1.0, 0), "1");
    EXPECT_EQ(IOUtils::formatDouble(0.0, 3), "0.000");
}

TEST_F(IOUtilsTest, ErrorHandling) {
    // Test writing to closed file
    std::ofstream closed_file;
    std::vector<std::string> headers = {"test"};
    EXPECT_THROW(IOUtils::writeCSVHeader(closed_file, headers), std::runtime_error);
    
    std::vector<double> data = {1.0};
    EXPECT_THROW(IOUtils::writeCSVRow(closed_file, data), std::runtime_error);
    
    // Test reading nonexistent file
    EXPECT_THROW(IOUtils::readCSV("nonexistent_file.csv"), std::runtime_error);
    
    // Test writing field with invalid filename (assuming we can't write to root)
    int Nx = 2, Ny = 2;
    double** field = new double*[Nx];
    for (int x = 0; x < Nx; x++) {
        field[x] = new double[Ny];
        for (int y = 0; y < Ny; y++) {
            field[x][y] = 1.0;
        }
    }
    
    // This might not throw on all systems, so we'll just test it doesn't crash
    try {
        IOUtils::writeFieldToCSV("/invalid/path/file.csv", field, Nx, Ny, "test", true);
    } catch (const std::runtime_error&) {
        // Expected behavior
    }
    
    // Clean up
    for (int x = 0; x < Nx; x++) {
        delete[] field[x];
    }
    delete[] field;
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}