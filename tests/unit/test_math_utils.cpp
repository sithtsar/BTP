#include <gtest/gtest.h>
#include <cmath>
#include <vector>
#include "../../include/utils/math_utils.hpp"

using namespace lbm::utils;

class MathUtilsTest : public ::testing::Test {
protected:
    void SetUp() override {
        // Set up test data
        test_distributions.resize(9);
        for (int i = 0; i < 9; i++) {
            test_distributions[i] = 0.1 + 0.01 * i;  // Valid distribution values
        }
        
        // Create test field
        Nx = 5;
        Ny = 5;
        test_field = new double*[Nx];
        for (int x = 0; x < Nx; x++) {
            test_field[x] = new double[Ny];
            for (int y = 0; y < Ny; y++) {
                test_field[x][y] = x + y * 0.1;  // Simple test pattern
            }
        }
    }
    
    void TearDown() override {
        // Clean up test field
        for (int x = 0; x < Nx; x++) {
            delete[] test_field[x];
        }
        delete[] test_field;
    }
    
    std::vector<double> test_distributions;
    double** test_field;
    int Nx, Ny;
};

TEST_F(MathUtilsTest, CalculateEquilibrium) {
    // Test valid equilibrium calculation
    double density = 1.0;
    double velocity_x = 0.1;
    double velocity_y = 0.05;
    
    // Test direction 0 (rest particle)
    double eq0 = MathUtils::calculateEquilibrium(0, density, velocity_x, velocity_y);
    EXPECT_GT(eq0, 0.0);
    EXPECT_TRUE(std::isfinite(eq0));
    
    // Test direction 1 (east)
    double eq1 = MathUtils::calculateEquilibrium(1, density, velocity_x, velocity_y);
    EXPECT_GT(eq1, 0.0);
    EXPECT_TRUE(std::isfinite(eq1));
    
    // Test invalid direction
    EXPECT_THROW(MathUtils::calculateEquilibrium(-1, density, velocity_x, velocity_y), 
                 std::invalid_argument);
    EXPECT_THROW(MathUtils::calculateEquilibrium(9, density, velocity_x, velocity_y), 
                 std::invalid_argument);
    
    // Test non-finite inputs
    EXPECT_THROW(MathUtils::calculateEquilibrium(0, NAN, velocity_x, velocity_y), 
                 std::invalid_argument);
    EXPECT_THROW(MathUtils::calculateEquilibrium(0, density, INFINITY, velocity_y), 
                 std::invalid_argument);
}

TEST_F(MathUtilsTest, CalculateDensity) {
    // Test valid density calculation
    double density = MathUtils::calculateDensity(test_distributions.data());
    EXPECT_GT(density, 0.0);
    EXPECT_TRUE(std::isfinite(density));
    
    // Test null pointer
    EXPECT_THROW(MathUtils::calculateDensity(nullptr), std::invalid_argument);
    
    // Test with NaN values
    std::vector<double> bad_distributions = test_distributions;
    bad_distributions[0] = NAN;
    EXPECT_THROW(MathUtils::calculateDensity(bad_distributions.data()), std::runtime_error);
}

TEST_F(MathUtilsTest, CalculateVelocity) {
    double density = 1.0;
    
    // Test valid velocity calculation
    auto velocity = MathUtils::calculateVelocity(test_distributions.data(), density);
    EXPECT_TRUE(std::isfinite(velocity.first));
    EXPECT_TRUE(std::isfinite(velocity.second));
    
    // Test null pointer
    EXPECT_THROW(MathUtils::calculateVelocity(nullptr, density), std::invalid_argument);
    
    // Test zero density
    auto zero_vel = MathUtils::calculateVelocity(test_distributions.data(), 0.0);
    EXPECT_EQ(zero_vel.first, 0.0);
    EXPECT_EQ(zero_vel.second, 0.0);
    
    // Test very small density
    auto small_vel = MathUtils::calculateVelocity(test_distributions.data(), 1e-15);
    EXPECT_EQ(small_vel.first, 0.0);
    EXPECT_EQ(small_vel.second, 0.0);
}

TEST_F(MathUtilsTest, CalculateGradient) {
    // Test valid gradient calculation
    auto gradient = MathUtils::calculateGradient(test_field, Nx, Ny, 2, 2);
    EXPECT_TRUE(std::isfinite(gradient.first));
    EXPECT_TRUE(std::isfinite(gradient.second));
    
    // Test null pointer
    EXPECT_THROW(MathUtils::calculateGradient(nullptr, Nx, Ny, 2, 2), 
                 std::invalid_argument);
    
    // Test out of bounds
    EXPECT_THROW(MathUtils::calculateGradient(test_field, Nx, Ny, -1, 2), 
                 std::out_of_range);
    EXPECT_THROW(MathUtils::calculateGradient(test_field, Nx, Ny, Nx, 2), 
                 std::out_of_range);
    EXPECT_THROW(MathUtils::calculateGradient(test_field, Nx, Ny, 2, -1), 
                 std::out_of_range);
    EXPECT_THROW(MathUtils::calculateGradient(test_field, Nx, Ny, 2, Ny), 
                 std::out_of_range);
}

TEST_F(MathUtilsTest, CalculateLaplacian) {
    // Test valid Laplacian calculation
    double laplacian = MathUtils::calculateLaplacian(test_field, Nx, Ny, 2, 2);
    EXPECT_TRUE(std::isfinite(laplacian));
    
    // Test null pointer
    EXPECT_THROW(MathUtils::calculateLaplacian(nullptr, Nx, Ny, 2, 2), 
                 std::invalid_argument);
    
    // Test out of bounds
    EXPECT_THROW(MathUtils::calculateLaplacian(test_field, Nx, Ny, -1, 2), 
                 std::out_of_range);
    EXPECT_THROW(MathUtils::calculateLaplacian(test_field, Nx, Ny, Nx, 2), 
                 std::out_of_range);
}

TEST_F(MathUtilsTest, Clamp) {
    // Test clamping
    EXPECT_EQ(MathUtils::clamp(5.0, 0.0, 10.0), 5.0);
    EXPECT_EQ(MathUtils::clamp(-1.0, 0.0, 10.0), 0.0);
    EXPECT_EQ(MathUtils::clamp(15.0, 0.0, 10.0), 10.0);
    
    // Test with integers
    EXPECT_EQ(MathUtils::clamp(5, 0, 10), 5);
    EXPECT_EQ(MathUtils::clamp(-1, 0, 10), 0);
    EXPECT_EQ(MathUtils::clamp(15, 0, 10), 10);
}

TEST_F(MathUtilsTest, IsFinite) {
    EXPECT_TRUE(MathUtils::isFinite(1.0));
    EXPECT_TRUE(MathUtils::isFinite(0.0));
    EXPECT_TRUE(MathUtils::isFinite(-1.0));
    EXPECT_FALSE(MathUtils::isFinite(NAN));
    EXPECT_FALSE(MathUtils::isFinite(INFINITY));
    EXPECT_FALSE(MathUtils::isFinite(-INFINITY));
}

TEST_F(MathUtilsTest, SafeDivision) {
    // Test normal division
    EXPECT_DOUBLE_EQ(MathUtils::safeDivision(10.0, 2.0), 5.0);
    
    // Test division by zero
    EXPECT_DOUBLE_EQ(MathUtils::safeDivision(10.0, 0.0, -1.0), -1.0);
    
    // Test division by very small number
    EXPECT_DOUBLE_EQ(MathUtils::safeDivision(10.0, 1e-15, -1.0), -1.0);
    
    // Test non-finite inputs
    EXPECT_DOUBLE_EQ(MathUtils::safeDivision(NAN, 2.0, -1.0), -1.0);
    EXPECT_DOUBLE_EQ(MathUtils::safeDivision(10.0, NAN, -1.0), -1.0);
}

TEST_F(MathUtilsTest, AnalyticalPoiseuilleVelocity) {
    double channel_height = 10.0;
    double body_force = 1e-6;
    double viscosity = 1.0/6.0;
    
    // Test center of channel (should be maximum velocity)
    double vel_center = MathUtils::analyticalPoiseuilleVelocity(5, channel_height, body_force, viscosity);
    EXPECT_GT(vel_center, 0.0);
    
    // Test at walls (should be zero)
    double vel_wall1 = MathUtils::analyticalPoiseuilleVelocity(0, channel_height, body_force, viscosity);
    double vel_wall2 = MathUtils::analyticalPoiseuilleVelocity(10, channel_height, body_force, viscosity);
    EXPECT_NEAR(vel_wall1, 0.0, 1e-10);
    EXPECT_NEAR(vel_wall2, 0.0, 1e-10);
    
    // Test invalid parameters
    EXPECT_THROW(MathUtils::analyticalPoiseuilleVelocity(5, 0.0, body_force, viscosity), 
                 std::invalid_argument);
    EXPECT_THROW(MathUtils::analyticalPoiseuilleVelocity(5, channel_height, body_force, 0.0), 
                 std::invalid_argument);
}

TEST_F(MathUtilsTest, CalculateL2Norm) {
    std::vector<double> field1 = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> field2 = {1.1, 2.1, 3.1, 4.1, 5.1};
    
    double norm = MathUtils::calculateL2Norm(field1.data(), field2.data(), field1.size());
    EXPECT_NEAR(norm, 0.1, 1e-10);
    
    // Test null pointers
    EXPECT_THROW(MathUtils::calculateL2Norm(nullptr, field2.data(), field1.size()), 
                 std::invalid_argument);
    EXPECT_THROW(MathUtils::calculateL2Norm(field1.data(), nullptr, field1.size()), 
                 std::invalid_argument);
    
    // Test invalid size
    EXPECT_THROW(MathUtils::calculateL2Norm(field1.data(), field2.data(), 0), 
                 std::invalid_argument);
}

TEST_F(MathUtilsTest, CalculateMaxDifference) {
    std::vector<double> field1 = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> field2 = {1.1, 2.05, 3.2, 4.1, 5.0};
    
    double max_diff = MathUtils::calculateMaxDifference(field1.data(), field2.data(), field1.size());
    EXPECT_NEAR(max_diff, 0.2, 1e-10);
    
    // Test null pointers
    EXPECT_THROW(MathUtils::calculateMaxDifference(nullptr, field2.data(), field1.size()), 
                 std::invalid_argument);
    EXPECT_THROW(MathUtils::calculateMaxDifference(field1.data(), nullptr, field1.size()), 
                 std::invalid_argument);
    
    // Test invalid size
    EXPECT_THROW(MathUtils::calculateMaxDifference(field1.data(), field2.data(), 0), 
                 std::invalid_argument);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}