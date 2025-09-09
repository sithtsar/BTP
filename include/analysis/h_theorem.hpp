#pragma once

#include <vector>
#include <string>
#include <memory>
#include <map>

namespace lbm {
namespace analysis {

/**
 * @brief Comprehensive H-theorem analysis for LBM simulations
 * 
 * This class provides detailed H-theorem verification and entropy analysis
 * for both standard and entropic BGK collision operators. It tracks H-function
 * evolution, verifies monotonicity, and calculates entropy production rates.
 */
class HTheoremAnalyzer {
public:
    /**
     * @brief Metrics for H-theorem analysis at a single timestep
     */
    struct HTheoremMetrics {
        double h_function_value = 0.0;      ///< Current H-function value
        double h_function_change = 0.0;     ///< Change from previous timestep
        bool monotonic_decrease = true;     ///< Whether H decreases monotonically
        double entropy_production = 0.0;    ///< Rate of entropy production
        double total_mass = 0.0;           ///< Total system mass
        double kinetic_energy = 0.0;       ///< Total kinetic energy
        int timestep = 0;                  ///< Associated timestep
        
        // Statistical measures
        double h_function_variance = 0.0;   ///< Spatial variance of H-function
        double h_function_min = 0.0;        ///< Minimum H-function value
        double h_function_max = 0.0;        ///< Maximum H-function value
        
        // Stability indicators
        bool has_nan_values = false;        ///< Whether NaN values detected
        bool has_inf_values = false;        ///< Whether infinite values detected
        double max_local_change = 0.0;      ///< Maximum local H-function change
    };

    /**
     * @brief Evolution data for H-theorem analysis over time
     */
    struct HEvolutionData {
        std::vector<int> timesteps;                    ///< Timestep values
        std::vector<double> h_values;                  ///< H-function values
        std::vector<double> h_changes;                 ///< H-function changes
        std::vector<double> entropy_production_rates;  ///< Entropy production rates
        std::vector<bool> monotonicity_violations;     ///< Monotonicity violation flags
        
        // Statistical evolution
        std::vector<double> h_variances;               ///< H-function spatial variances
        std::vector<double> kinetic_energies;          ///< Kinetic energy evolution
        std::vector<double> total_masses;              ///< Mass conservation tracking
        
        // Derived metrics
        double average_entropy_production = 0.0;       ///< Time-averaged entropy production
        int monotonicity_violation_count = 0;          ///< Number of violations
        double h_function_decay_rate = 0.0;            ///< Exponential decay rate
        bool overall_monotonic = true;                 ///< Overall monotonicity assessment
    };

    /**
     * @brief Configuration for H-theorem analysis
     */
    struct AnalysisConfig {
        bool enable_spatial_analysis = true;           ///< Enable spatial H-function analysis
        bool enable_monotonicity_check = true;         ///< Enable monotonicity verification
        bool enable_entropy_production = true;         ///< Enable entropy production calculation
        bool enable_statistical_analysis = true;       ///< Enable statistical measures
        
        double monotonicity_tolerance = 1e-12;         ///< Tolerance for monotonicity violations
        int analysis_interval = 1;                     ///< Analyze every N timesteps
        bool log_violations = true;                    ///< Log monotonicity violations
        
        // Output configuration
        bool write_evolution_data = true;              ///< Write H-evolution to file
        std::string output_prefix = "h_theorem";       ///< Output file prefix
    };

    /**
     * @brief Constructor
     * @param Nx Grid size in x-direction
     * @param Ny Grid size in y-direction
     * @param config Analysis configuration
     */
    HTheoremAnalyzer(int Nx, int Ny, const AnalysisConfig& config);

    /**
     * @brief Constructor with default configuration
     * @param Nx Grid size in x-direction
     * @param Ny Grid size in y-direction
     */
    HTheoremAnalyzer(int Nx, int Ny);

    /**
     * @brief Destructor
     */
    ~HTheoremAnalyzer() = default;

    /**
     * @brief Calculate H-function for standard BGK collision operator
     * @param f Distribution functions [Q][Nx][Ny]
     * @param rho Density field [Nx][Ny]
     * @param ux Velocity x-component [Nx][Ny]
     * @param uy Velocity y-component [Nx][Ny]
     * @param timestep Current timestep
     * @return H-theorem metrics
     */
    HTheoremMetrics calculateHFunctionStandardBGK(
        double*** f, double** rho, double** ux, double** uy, int timestep);

    /**
     * @brief Calculate H-function for entropic BGK collision operator
     * @param f Distribution functions [Q][Nx][Ny]
     * @param rho Density field [Nx][Ny]
     * @param ux Velocity x-component [Nx][Ny]
     * @param uy Velocity y-component [Nx][Ny]
     * @param timestep Current timestep
     * @return H-theorem metrics
     */
    HTheoremMetrics calculateHFunctionEntropicBGK(
        double*** f, double** rho, double** ux, double** uy, int timestep);

    /**
     * @brief Record H-theorem metrics for evolution tracking
     * @param metrics H-theorem metrics to record
     */
    void recordHEvolution(const HTheoremMetrics& metrics);

    /**
     * @brief Verify monotonic decrease of H-function
     * @param tolerance Tolerance for numerical errors
     * @return true if H-function decreases monotonically
     */
    bool verifyMonotonicDecrease(double tolerance = 1e-12) const;

    /**
     * @brief Get complete evolution data
     * @return Reference to evolution data
     */
    const HEvolutionData& getEvolutionData() const { return evolution_data_; }

    /**
     * @brief Get latest H-theorem metrics
     * @return Latest metrics, or default if no data recorded
     */
    HTheoremMetrics getLatestMetrics() const;

    /**
     * @brief Calculate entropy production rate
     * @param current_h Current H-function value
     * @param previous_h Previous H-function value
     * @param dt Time step size
     * @return Entropy production rate
     */
    double calculateEntropyProductionRate(double current_h, double previous_h, double dt) const;

    /**
     * @brief Analyze H-function spatial distribution
     * @param h_field H-function field [Nx][Ny]
     * @return Statistical measures of spatial distribution
     */
    struct SpatialStatistics {
        double mean = 0.0;
        double variance = 0.0;
        double min_value = 0.0;
        double max_value = 0.0;
        double standard_deviation = 0.0;
    };
    SpatialStatistics analyzeSpatialDistribution(double** h_field) const;

    /**
     * @brief Write H-theorem evolution data to file
     * @param filename Output filename
     * @return true if successful
     */
    bool writeEvolutionData(const std::string& filename) const;

    /**
     * @brief Generate diagnostic report
     * @return String containing detailed analysis report
     */
    std::string generateDiagnosticReport() const;

    /**
     * @brief Reset analysis data
     */
    void reset();

    /**
     * @brief Set analysis configuration
     * @param config New configuration
     */
    void setConfig(const AnalysisConfig& config) { config_ = config; }

    /**
     * @brief Get current configuration
     * @return Reference to current configuration
     */
    const AnalysisConfig& getConfig() const { return config_; }

private:
    /**
     * @brief Calculate equilibrium distribution function
     * @param direction Velocity direction index
     * @param rho Density
     * @param ux Velocity x-component
     * @param uy Velocity y-component
     * @return Equilibrium distribution value
     */
    double calculateEquilibrium(int direction, double rho, double ux, double uy) const;

    /**
     * @brief Calculate H-function at a single grid point
     * @param f Distribution functions at grid point [Q]
     * @param feq Equilibrium distribution functions at grid point [Q]
     * @return Local H-function value
     */
    double calculateLocalHFunction(const double* f, const double* feq) const;

    /**
     * @brief Calculate H-function using logarithmic form
     * @param f Distribution functions at grid point [Q]
     * @return Local H-function value (logarithmic form)
     */
    double calculateLocalHFunctionLog(const double* f) const;

    /**
     * @brief Check for numerical anomalies in H-function field
     * @param h_field H-function field [Nx][Ny]
     * @return true if anomalies detected
     */
    bool checkNumericalAnomalies(double** h_field) const;

    /**
     * @brief Update evolution statistics
     */
    void updateEvolutionStatistics();

    /**
     * @brief Log monotonicity violation
     * @param timestep Timestep where violation occurred
     * @param h_current Current H-function value
     * @param h_previous Previous H-function value
     */
    void logMonotonicityViolation(int timestep, double h_current, double h_previous) const;

    // Grid parameters
    int Nx_, Ny_;
    static constexpr int Q = 9;  // D2Q9 lattice

    // LBM constants
    static constexpr double cs2 = 1.0 / 3.0;  // Speed of sound squared
    static constexpr int ex[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
    static constexpr int ey[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
    static constexpr double w[Q] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 
                                   1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};

    // Configuration
    AnalysisConfig config_;

    // Evolution tracking
    HEvolutionData evolution_data_;
    std::vector<HTheoremMetrics> metrics_history_;

    // Working arrays for H-function calculation
    std::unique_ptr<double[]> h_field_flat_;
    std::unique_ptr<double[]> feq_temp_;

    // Previous timestep data for change calculation
    double previous_h_value_;
    bool has_previous_data_;
    int last_analysis_timestep_;
};

} // namespace analysis
} // namespace lbm