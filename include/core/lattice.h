#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <array>
#include <cmath>

namespace elbm {

// D2Q9 Lattice structure
class D2Q9 {
public:
    static constexpr int D = 2;  // Dimensions
    static constexpr int Q = 9;  // Number of velocities

    // Discrete velocities
    static constexpr std::array<std::array<int, 2>, 9> c = {{
        {0, 0},   // 0: rest
        {1, 0},   // 1: east
        {0, 1},   // 2: north
        {-1, 0},  // 3: west
        {0, -1},  // 4: south
        {1, 1},   // 5: northeast
        {-1, 1},  // 6: northwest
        {-1, -1}, // 7: southwest
        {1, -1}   // 8: southeast
    }};

    // Weights for equilibrium distribution
    static constexpr std::array<double, 9> w = {
        4.0/9.0,   // 0: rest
        1.0/9.0,   // 1-4: cardinal directions
        1.0/9.0,
        1.0/9.0,
        1.0/9.0,
        1.0/36.0,  // 5-8: diagonal directions
        1.0/36.0,
        1.0/36.0,
        1.0/36.0
    };

    // Speed of sound squared
    static constexpr double cs2 = 1.0/3.0;
    static constexpr double cs = 0.57735026918962576451; // 1/sqrt(3)

    // Opposite direction indices (for bounce-back)
    static constexpr std::array<int, 9> opposite = {0, 3, 4, 1, 2, 7, 8, 5, 6};

    // Get velocity component
    static inline int cx(int i) { return c[i][0]; }
    static inline int cy(int i) { return c[i][1]; }

    // Get weight
    static inline double weight(int i) { return w[i]; }
};

// D3Q19 Lattice structure
class D3Q19 {
public:
    static constexpr int D = 3;  // Dimensions
    static constexpr int Q = 19; // Number of velocities

    // Discrete velocities
    static constexpr std::array<std::array<int, 3>, 19> c = {{
        {0, 0, 0},    // 0: rest
        {1, 0, 0},    // 1-6: face directions
        {-1, 0, 0},
        {0, 1, 0},
        {0, -1, 0},
        {0, 0, 1},
        {0, 0, -1},
        {1, 1, 0},    // 7-18: edge directions
        {-1, -1, 0},
        {1, -1, 0},
        {-1, 1, 0},
        {1, 0, 1},
        {-1, 0, -1},
        {1, 0, -1},
        {-1, 0, 1},
        {0, 1, 1},
        {0, -1, -1},
        {0, 1, -1},
        {0, -1, 1}
    }};

    // Weights for equilibrium distribution
    static constexpr std::array<double, 19> w = {
        1.0/3.0,   // 0: rest
        1.0/18.0,  // 1-6: face directions
        1.0/18.0,
        1.0/18.0,
        1.0/18.0,
        1.0/18.0,
        1.0/18.0,
        1.0/36.0,  // 7-18: edge directions
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
        1.0/36.0
    };

    // Speed of sound squared
    static constexpr double cs2 = 1.0/3.0;
    static constexpr double cs = 0.57735026918962576451; // 1/sqrt(3)

    // Opposite direction indices
    static constexpr std::array<int, 19> opposite = {
        0, 2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11, 14, 13, 16, 15, 18, 17
    };

    // Get velocity components
    static inline int cx(int i) { return c[i][0]; }
    static inline int cy(int i) { return c[i][1]; }
    static inline int cz(int i) { return c[i][2]; }

    // Get weight
    static inline double weight(int i) { return w[i]; }
};

} // namespace elbm

#endif // LATTICE_H
