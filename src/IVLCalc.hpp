#pragma once
#include <array>
#include <stdint.h>
#include <vector>

/// @brief Class that enables interval length computation in FDL(n), via intervals of FDL(n-1),
/// FDL(n-2) or FDL(n-3).
class IVLCalc {
public:
    IVLCalc() = default;

    /// C'tor that initializes all data necessary for the interval length calculation.
    /// @param n number of generators
    /// @param difference use FDL(n-difference) to calculate interval length
    IVLCalc(unsigned n, unsigned difference);

    /// @brief Calculates the length of the interval v = [x,y] in FDL(n).
    /// @param vector containing upper bound and lower bound
    uint64_t intervalLength(const std::vector<uint64_t>& v) const;

private:
    /// @brief calculates length of interval [x,y] via FDL(n-1)
    uint64_t intervalLength1(uint64_t a, uint64_t b, uint64_t c, uint64_t d) const;

    /// @brief calculates length of interval [x,y] via FDL(n-2)

    uint64_t intervalLength2(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e, uint64_t f,
                             uint64_t g, uint64_t h) const;

    /// @brief calculates length of interval [x,y] via FDL(n-3)
    uint64_t intervalLength3(uint64_t a, uint64_t b, uint64_t c, uint64_t d, uint64_t e, uint64_t f,
                             uint64_t g, uint64_t h, uint64_t i, uint64_t j, uint64_t k, uint64_t l,
                             uint64_t m, uint64_t n, uint64_t o, uint64_t p) const;

    /// @brief initialization of mIntRed (Intervals of FDL(n-1), FDL(n-2) or FDL(n-3))
    void initIntervals();

    /// @brief splits x in two equal parts w.r.t. mBitsize
    std::array<uint64_t, 2> splitIn2(uint64_t x) const;

    /// @brief splits x in four equal parts w.r.t. mBitsize
    std::array<uint64_t, 4> splitIn4(uint64_t x) const;

    /// @brief splits x in eight equal parts w.r.t. mBitsize
    std::array<uint64_t, 8> splitIn8(uint64_t x) const;

    /// @brief generates the interval [x,y] w.r.t. mFdlRed
    std::vector<uint64_t> generateInterval(uint64_t x, uint64_t y) const;

    /// @brief generates the interval interval [x,mIndex[y]] w.r.t. mFdlRed
    std::vector<uint64_t> generateIntervalMixed(uint64_t x, unsigned int ind) const;

    /// @brief initialization of mIndex
    void initIndices();

    // nr of bits needed to represent elements of FDL(n)
    unsigned int mNrOfBits;

    // difference value for FDL(n-difference), to switch interval length calculation method
    unsigned int mDifference;

    // nr of generators of FDL(n-difference)
    unsigned int mGeneratorsRed;

    // elements of FDL(n-difference)
    std::vector<uint64_t> mFdlRed;
    // intervals of FDL(n-difference)
    std::vector<uint64_t> mIntRed;
    // index map FDL(n-difference)
    std::vector<unsigned int> mIndex;
};
