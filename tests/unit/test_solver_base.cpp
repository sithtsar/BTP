#include <gtest/gtest.h>
#include "../../include/lbm/solver_base.hpp"

using namespace lbm;

// Concrete implementation of SolverBase for testing
class TestSolver : public SolverBase {
public:
    TestSolver(int Nx, int Ny) : SolverBase(Nx, Ny), initialized_(false), stable_(true) {}
    
    void initialize() override {
        initialized_ = true;
    }
    
    void step() override {
        if (!initialized_) {
            throw std::runtime_error("Solver not initialized");
        }
        updateTimestep(getCurrentTimestep() + 1);
    }
    
    bool isStable() const override {
        return stable_;
    }
    
    void writeOutput(int timestep) override {
        // Mock implementation
        last_output_timestep_ = timestep;
    }
    
    // Test helpers
    void setStable(bool stable) { stable_ = stable; }
    bool isInitialized() const { return initialized_; }
    int getLastOutputTimestep() const { return last_output_timestep_; }
    
private:
    bool initialized_;
    bool stable_;
    int last_output_timestep_ = -1;
};

class SolverBaseTest : public ::testing::Test {
protected:
    void SetUp() override {
        solver = std::make_unique<TestSolver>(10, 5);
    }
    
    void TearDown() override {
        solver.reset();
    }
    
    std::unique_ptr<TestSolver> solver;
};

TEST_F(SolverBaseTest, Construction) {
    // Test valid construction
    EXPECT_EQ(solver->getGridDimensions().first, 10);
    EXPECT_EQ(solver->getGridDimensions().second, 5);
    EXPECT_EQ(static_cast<const TestSolver*>(solver.get())->getCurrentTimestep(), 0);
    
    // Test invalid construction
    EXPECT_THROW(TestSolver(-1, 5), std::invalid_argument);
    EXPECT_THROW(TestSolver(10, 0), std::invalid_argument);
    EXPECT_THROW(TestSolver(0, 0), std::invalid_argument);
}

TEST_F(SolverBaseTest, Initialization) {
    // Test initialization
    EXPECT_FALSE(solver->isInitialized());
    solver->initialize();
    EXPECT_TRUE(solver->isInitialized());
}

TEST_F(SolverBaseTest, Timestepping) {
    solver->initialize();
    
    // Test initial timestep
    EXPECT_EQ(static_cast<const TestSolver*>(solver.get())->getCurrentTimestep(), 0);
    
    // Test stepping
    solver->step();
    EXPECT_EQ(static_cast<const TestSolver*>(solver.get())->getCurrentTimestep(), 1);
    
    solver->step();
    EXPECT_EQ(static_cast<const TestSolver*>(solver.get())->getCurrentTimestep(), 2);
    
    // Test stepping without initialization
    TestSolver uninitialized_solver(5, 5);
    EXPECT_THROW(uninitialized_solver.step(), std::runtime_error);
}

TEST_F(SolverBaseTest, StabilityMonitoring) {
    // Test default stability
    EXPECT_TRUE(solver->isStable());
    
    // Test stability monitoring enabled by default
    EXPECT_TRUE(solver->testCheckStability());
    
    // Test setting instability
    solver->setStable(false);
    EXPECT_FALSE(solver->isStable());
    
    // Test disabling stability monitoring
    solver->setStabilityMonitoring(false);
    EXPECT_TRUE(solver->testCheckStability());  // Should return true when disabled
    
    // Test re-enabling stability monitoring
    solver->setStabilityMonitoring(true);
    // Note: checkAndHandleStability() behavior depends on actual stability data
    // which is not available in this mock implementation
}

TEST_F(SolverBaseTest, OutputHandling) {
    // Test output directory setting
    solver->setOutputDirectory("test_output");
    
    // Test write output
    solver->writeOutput(100);
    EXPECT_EQ(solver->getLastOutputTimestep(), 100);
    
    solver->writeOutput(200);
    EXPECT_EQ(solver->getLastOutputTimestep(), 200);
}

TEST_F(SolverBaseTest, StabilityMonitorAccess) {
    // Test accessing stability monitor
    const auto& monitor = solver->getStabilityMonitor();
    
    // Test that we can access the monitor (basic functionality test)
    // More detailed testing would require actual stability data
    EXPECT_FALSE(monitor.hasStabilityIssues());  // Should be false initially
}

TEST_F(SolverBaseTest, DataValidation) {
    // Test data validation (mock implementation)
    EXPECT_TRUE(solver->testValidateData("test data"));
}

TEST_F(SolverBaseTest, LBMConstants) {
    // Test that LBM constants are accessible and correct
    // These are protected members, so we test them indirectly through derived class
    
    // Test grid dimensions are stored correctly
    auto dims = solver->getGridDimensions();
    EXPECT_EQ(dims.first, 10);
    EXPECT_EQ(dims.second, 5);
    
    // The actual LBM constants (Q, cs2, ex, ey, w) are tested implicitly
    // through the math utilities tests
}

TEST_F(SolverBaseTest, PolymorphicBehavior) {
    // Test polymorphic behavior through base class pointer
    std::unique_ptr<SolverBase> base_ptr = std::make_unique<TestSolver>(8, 4);
    
    // Test virtual function calls
    base_ptr->initialize();
    EXPECT_EQ(static_cast<const SolverBase*>(base_ptr.get())->getCurrentTimestep(), 0);
    
    base_ptr->step();
    EXPECT_EQ(static_cast<const SolverBase*>(base_ptr.get())->getCurrentTimestep(), 1);
    
    EXPECT_TRUE(base_ptr->isStable());
    
    base_ptr->writeOutput(50);
    // Can't test the result directly since writeOutput is pure virtual,
    // but we can verify it doesn't crash
    
    EXPECT_EQ(base_ptr->getGridDimensions().first, 8);
    EXPECT_EQ(base_ptr->getGridDimensions().second, 4);
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}