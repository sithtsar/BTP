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
 * @brief Multiphase LBM solver with phase-field method for two-phase flows
 * 
 * This class implements a multiphase Lattice Boltzmann Method solver using
 * the phase-field approach for simulating two-phase flows with surface tension.
 * It includes comprehensive stability monitoring and interface tracking validation.
 */
class MultiphaseSolver : public SolverBase {
public:
    /**
     * @brief Configuration parameters for multiphase simulation
     */
    struct SimulationConfig {
        // Grid parameters
        int Nx = 256;                   ///< Grid size in x-direction
        int Ny = 64;                    ///< Grid size in y-direction
        
        // Phase field parameters
        double phi_L = 0.0;             ///< Light fluid phase field value
        double phi_H = 1.0;             ///< Heavy fluid phase field value
        double phi_0 = 0.5;             ///< Interface location
        
        // Fluid properties
        double rho_L = 1.0;             ///< Light fluid density
        double rho_H = 1000.0;          ///< Heavy fluid density
        double mu_L = 0.01;             ///< Light fluid viscosity
        double mu_H = 1.0;              ///< Heavy fluid viscosity
        
        // Surface tension parameters
        double xi = 4.0;                ///< Interface thickness
        double sigma = 0.01;            ///< Surface tension coefficient
        double beta = 0.0;              ///< Calculated from sigma and xi
        double kappa = 0.0;             ///< Calculated from sigma and xi
        
        // Simulation parameters
        double tau_phi = 0.7;           ///< Phase field relaxation time
        double gravity = 1e-6;          ///< Gravitational acceleration
        int max_timesteps = 20000;      ///< Maximum number of timesteps
        int output_interval = 1000;     ///< Output data every N timesteps
        
        // Stability parameters
        double max_velocity_limit = 0.1;    ///< Maximum allowed velocity
        double min_density_limit = 1e-6;    ///< Minimum allowed density
        double max_density_ratio = 1e6;     ///< Maximum density ratio
        double min_tau_limit = 0.51;        ///< Minimum relaxation time
        double max_tau_limit = 5.0;         ///< Maximum relaxation time
        int stability_check_interval = 100; ///< Check stability every N timesteps
        
        // Interface tracking parameters
        double interface_thickness_tolerance = 0.1;  ///< Allowed interface thickness variation
        double phase_conservation_tolerance = 1e-6;  ///< Phase conservation tolerance
        double interface_position_tolerance = 0.5;   ///< Interface position tracking tolerance
        
        // Output parameters
        std::string output_prefix = "multiphase";    ///< Prefix for output files
        bool write_full_fields = true;               ///< Write complete field data
        bool write_averaged_profiles = true;         ///< Write y-averaged profiles
        bool write_analytical_comparison = true;     ///< Include analytical solution
        
        // H-theorem analysis parameters
        bool enable_h_theorem_analysis = true;       ///< Enable H-theorem monitoring
        int h_analysis_interval = 10;                ///< Analyze H-theorem every N timesteps
    };

    /**
     * @brief Constructor
     * @param config Simulation configuration parameters
     */
    explicit MultiphaseSolver(const SimulationConfig& config);

    /**
     * @brief Destructor - cleans up allocated memory
     */
    ~MultiphaseSolver();

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
     * @brief Get phase field value at grid point
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     * @return Phase field value
     */
    double getPhaseField(int x, int y) const;

    /**
     * @brief Get density at grid point
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     * @return Density value
     */
    double getDensity(int x, int y) const;

    /**
     * @brief Get velocity at grid point
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     * @return Pair of (velocity_x, velocity_y)
     */
    std::pair<double, double> getVelocity(int x, int y) const;

    /**
     * @brief Get pressure at grid point
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     * @return Pressure value
     */
    double getPressure(int x, int y) const;

    /**
     * @brief Calculate analytical two-layer Poiseuille velocity
     * @param y Grid position in y-direction
     * @return Analytical velocity
     */
    double getAnalyticalVelocity(int y) const;

    /**
     * @brief Get interface position at given x-coordinate
     * @param x Grid position in x-direction
     * @return Interface y-position
     */
    double getInterfacePosition(int x) const;

    /**
     * @brief Calculate total mass of each phase
     * @return Pair of (light_phase_mass, heavy_phase_mass)
     */
    std::pair<double, double> calculatePhaseMasses() const;

private:
    /**
     * @brief Allocate memory for all field arrays
     */
    void allocateMemory();

    /**
     * @brief Deallocate memory for all field arrays
     */
    void deallocateMemory();

    /**
     * @brief Initialize phase field with layered configuration
     */
    void initializePhaseField();

    /**
     * @brief Initialize hydrodynamic fields
     */
    void initializeHydrodynamics();

    /**
     * @brief Update fluid properties from phase field
     */
    void updateFluidProperties();

    /**
     * @brief Calculate gradients using isotropic finite differences
     */
    void calculateGradients();

    /**
     * @brief Calculate chemical potential for phase field
     * @param x Grid position in x-direction
     * @param y Grid position in y-direction
     * @return Chemical potential value
     */
    double calculateChemicalPotential(int x, int y) const;

    /**
     * @brief Perform phase field collision step
     */
    void phaseFieldCollision();

    /**
     * @brief Perform hydrodynamic collision step
     */
    void hydrodynamicCollision();

    /**
     * @brief Perform streaming step for both distribution functions
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
     * @brief Check interface tracking stability
     * @return true if interface is stable, false otherwise
     */
    bool checkInterfaceStability();

    /**
     * @brief Check phase conservation
     * @return true if phases are conserved within tolerance
     */
    bool checkPhaseConservation();

    /**
     * @brief Write complete field data to files
     * @param timestep Current timestep
     */
    void writeFullFieldData(int timestep);

    /**
     * @brief Write y-averaged profile data
     * @param timestep Current timestep
     */
    void writeAveragedProfiles(int timestep);

    // Configuration
    SimulationConfig config_;

    // Distribution functions
    double*** h_;        ///< Phase field distribution functions [Q][Nx][Ny]
    double*** h_temp_;   ///< Temporary phase field distributions
    double*** g_;        ///< Velocity distribution functions [Q][Nx][Ny]
    double*** g_temp_;   ///< Temporary velocity distributions

    // Macroscopic variables
    double** phi_;       ///< Phase field [Nx][Ny]
    double** rho_;       ///< Density field [Nx][Ny]
    double** mu_;        ///< Dynamic viscosity field [Nx][Ny]
    double** tau_;       ///< Relaxation time field [Nx][Ny]
    double** ux_;        ///< Velocity x-component [Nx][Ny]
    double** uy_;        ///< Velocity y-component [Nx][Ny]
    double** p_;         ///< Pressure field [Nx][Ny]

    // Gradient fields
    double** grad_phi_x_; ///< Phase field gradient x-component [Nx][Ny]
    double** grad_phi_y_; ///< Phase field gradient y-component [Nx][Ny]
    double** lapl_phi_;   ///< Phase field Laplacian [Nx][Ny]

    // Derived parameters
    double M_;           ///< Mobility parameter
    double cs2_;         ///< Speed of sound squared

    // Stability tracking
    bool is_stable_;
    bool interface_stable_;
    bool phase_conserved_;
    int last_stability_check_;
    std::pair<double, double> initial_phase_masses_;

    // Memory management flag
    bool memory_allocated_;
};

} // namespace lbm