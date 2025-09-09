#pragma once

#include "solver_base.hpp"
#include "../utils/math_utils.hpp"
#include "../utils/io_utils.hpp"
#include "../analysis/h_theorem.hpp"
#include <memory>
#include <string>
#include <vector>

namespace lbm {

/**
 * @brief Single-phase LBM solver for Poiseuille flow with BGK collision operators
 * 
 * This class implements a single-phase Lattice Boltzmann Method solver
 * with support for both standard and entropic BGK collision operators.
 * It includes comprehensive stability monitoring and error handling.
 */
class SinglePhaseSolver : public SolverBase {
public:
    /**
     * @brief Configuration parameters for single-phase simulation
     */
    struct SimulationConfig {
        // Grid parameters
        int Nx = 20;                    ///< Grid size in x-direction
        int Ny = 21;                    ///< Grid size in y-direction
        
        // Physical parameters
        double tau = 1.0;               ///< Relaxation time
        double rho0 = 1.0;              ///< Reference density
        double gravity = 1e-6;          ///< Body force (gravity)
        
        // Simulation parameters
        int max_timesteps = 10000;      ///< Maximum number of timesteps
        int output_interval = 1000;     ///< Output data every N timesteps
        bool use_entropic_bgk = false;  ///< Use entropic BGK collision operator
        
        // Stability parameters
        double max_velocity_limit = 0.1;    ///< Maximum allowed velocity
        double min_density_limit = 1e-6;    ///< Minimum allowed density
        int stability_check_interval = 100; ///< Check stability every N timesteps
        
        // Output parameters
        std::string output_prefix = "single_phase";  ///< Prefix for output files
        bool write_analytical_comparison = true;     ///< Include analytical solution in output
        
        // H-theorem analysis parameters
        bool enable_h_theorem_analysis = true;       ///< Enable H-theorem monitoring
        int h_analysis_interval = 1;                 ///< Analyze H-theorem every N timesteps
    };

    /**
     * @brief Constructor
     * @param config Simulation configuration parameters
     */
    explicit SinglePhaseSolver(const SimulationConfig& config);

    /**
     * @brief Destructor - cleans up allocated memory
     */
    ~SinglePhaseSolver();

    // Inherited from SolverBase
    void initialize() override;
    void step() override;
    bool isStable() const override;
    void writeOutput(int timestep) override;

    /**
     * @brief Run complete simulation
     * @return true if simulation completed successfully, false if terminated due to instability
     */
    bool runSimulation();

    /**
     * @brief Get current simulation configuration
     * @return Reference to configuration
     */
    const SimulationConfig& getConfig() const { return config_; }

    /**
     * @brief Get macroscopic density at grid point
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     * @return Density value
     */
    double getDensity(int x, int y) const;

    /**
     * @brief Get macroscopic velocity at grid point
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     * @return Pair of (velocity_x, velocity_y)
     */
    std::pair<double, double> getVelocity(int x, int y) const;

    /**
     * @brief Calculate analytical Poiseuille velocity profile
     * @param y Grid position in y-direction
     * @return Analytical velocity
     */
    double getAnalyticalVelocity(int y) const;

    /**
     * @brief Get H-function value for entropy analysis
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     * @return H-function value
     */
    double getHFunction(int x, int y) const;

    /**
     * @brief Get H-theorem analyzer
     * @return Reference to H-theorem analyzer
     */
    const analysis::HTheoremAnalyzer& getHTheoremAnalyzer() const { return *h_analyzer_; }

    /**
     * @brief Get latest H-theorem metrics
     * @return Latest H-theorem metrics
     */
    analysis::HTheoremAnalyzer::HTheoremMetrics getLatestHMetrics() const;

private:
    /**
     * @brief Allocate memory for distribution functions
     */
    void allocateMemory();

    /**
     * @brief Deallocate memory for distribution functions
     */
    void deallocateMemory();

    /**
     * @brief Initialize distribution functions to equilibrium
     */
    void initializeDistributions();

    /**
     * @brief Perform collision step
     */
    void collision();

    /**
     * @brief Perform standard BGK collision
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     */
    void standardBGKCollision(int x, int y);

    /**
     * @brief Perform entropic BGK collision
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     */
    void entropicBGKCollision(int x, int y);

    /**
     * @brief Perform streaming step
     */
    void streaming();

    /**
     * @brief Apply boundary conditions
     */
    void applyBoundaryConditions();

    /**
     * @brief Update macroscopic variables from distribution functions
     */
    void updateMacroscopicVariables();

    /**
     * @brief Check numerical stability of current state
     * @return true if stable, false otherwise
     */
    bool checkNumericalStability();

    /**
     * @brief Calculate force term for collision operator
     * @param direction Velocity direction index
     * @param ux Velocity in x-direction
     * @param uy Velocity in y-direction
     * @return Force term value
     */
    double calculateForceTerm(int direction, double ux, double uy) const;

    /**
     * @brief Find optimal alpha parameter for entropic BGK
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     * @param feq Equilibrium distribution functions
     * @return Optimal alpha value
     */
    double findOptimalAlpha(int x, int y, const double* feq) const;

    /**
     * @brief Calculate H-function (entropy) at a grid point
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     * @return H-function value
     */
    double calculateHFunction(int x, int y) const;

    /**
     * @brief Write velocity profile data
     * @param timestep Current timestep
     */
    void writeVelocityProfile(int timestep);

    /**
     * @brief Write H-function evolution data
     * @param timestep Current timestep
     */
    void writeHFunctionData(int timestep);

    /**
     * @brief Perform H-theorem analysis
     * @param timestep Current timestep  
     */
    void performHTheoremAnalysis(int timestep);

    /**
     * @brief Write H-theorem analysis results
     * @param timestep Current timestep
     */
    void writeHTheoremResults(int timestep);

    // Configuration
    SimulationConfig config_;

    // Distribution functions
    double*** f_;        ///< Distribution functions [Q][Nx][Ny]
    double*** f_temp_;   ///< Temporary distribution functions for collision

    // Macroscopic variables
    double** rho_;       ///< Density field [Nx][Ny]
    double** ux_;        ///< Velocity x-component [Nx][Ny]
    double** uy_;        ///< Velocity y-component [Nx][Ny]

    // Derived parameters
    double omega_;       ///< Collision frequency (1/tau)
    double nu_;          ///< Kinematic viscosity
    double cs2_;         ///< Speed of sound squared

    // Stability tracking
    bool is_stable_;
    int last_stability_check_;

    // H-theorem analysis
    std::unique_ptr<analysis::HTheoremAnalyzer> h_analyzer_;

    // Memory management flag
    bool memory_allocated_;
};

} // namespace lbm