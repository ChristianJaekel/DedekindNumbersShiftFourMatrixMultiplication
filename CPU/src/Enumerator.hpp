#pragma once

#include <IVLCalc.hpp>
#include <Uint128.hpp>
#include <vector>

using V8T       = std::vector<uint8_t>;
using V16T      = std::vector<uint16_t>;
using V32T      = std::vector<uint32_t>;
using V64T      = std::vector<uint64_t>;
using IntervalT = std::vector<V64T>;

/// @brief Class that holds the algorithm to enumerate the fee distributive lattice, as described in
/// the preprint "A computation of the ninth Dedekind Number" by Christian Jäkel, 2023.
class Enumerator {
public:
    /// C'tor that initializes all data necessary for the enumeration.
    /// That is generating FDL(n-4), an index map for the elements of FDL(n-4) and the interval size
    /// calculator.
    /// @param n FDL(n)
    Enumerator(unsigned int n);

    /// @brief Performs the computation using the CPU.
    /// @param verbose if true, prints progress and details of the computation
    uint128T doEnumerationCPU(bool verbose = false);

private:
    /// @brief Clusters intervals [x,y] by their bottom value x.
    /// E.g. [0,1],[0,3],[0,15],[1,3],[1,5],[3,15] --> [[0,1],[0,3],[0,15]],[[1,3],[1,5]],[[3,15]]
    std::vector<IntervalT> clusterIntervals(const IntervalT& intervals);

    /// @brief Returns the index of the element x.
    /// index computation, instead of lookup, can be implemented here.
    inline unsigned int index(uint64_t x) const { return mIdx[x]; };

    /// @brief Reads intervals from file.
    /// @details Intervals represent an equivalence class. Hence, the classes cardinality is stored
    /// in a third value: bottom, top, cardinality
    /// @param filename name/path of the file
    IntervalT readIntervals(const std::string& filename) const;

    /// @brief Reads ab-values from file.
    /// @details Values represent an equivalence class. Hence, the classes cardinality is stored
    /// in a third value: a, b, cardinality.
    /// Values are order w.r.t. a and delta compressed w.r.t. b.
    /// @param filename name/path of the file
    IntervalT readABValues(const std::string& filename) const;

    /// @brief Generates the interval [x,y] as a vector of elements.
    /// @param x interval bottom
    /// @param y interval top
    V16T generateInterval(uint64_t x, uint64_t y) const;

    // elements of FDL(n-4)
    V64T mFdl{};

    // index map for the elements of FDL(n-4)
    V8T mIdx{};

    // number of generators of FDL(n-4)
    uint64_t mNumOfGens{};

    // interval size calculator
    IVLCalc mIvlCalc{};
};
