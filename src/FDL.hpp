#pragma once
#include <cstdint>
#include <vector>

/// @brie Class that represents a free distributive lattice.
class FreeDistLat {
public:
    /// @brief c'tor that generates all elements of the free distributive lattice up to n<=6 in
    /// ascending order
    /// @param n number of generators
    FreeDistLat(const uint64_t n);

    /// @brief writes elements as uints to cout
    void printElements() const;

    /// @brief writes elements in binary to cout
    void printElementsBinary() const;

    /// @brief writes number of elements to cout
    void printNumberOfElements() const;

    /// @brief writes detailed neighborhood analysis of elements to cout
    void printElementsDetailed() const;

    /// @brief returns the number of elements
    uint64_t getNumberOfElements() const;

    /// @brief returns the number of generators
    uint64_t getNumberOfGenerators() const;

    /// @brief returns the number of bits needed to represent the elements
    uint64_t getNumberOfBits() const;

    /// @brief returns a vector containing all elements of FDL(n)
    std::vector<uint64_t> getElements() const;

private:
    std::vector<uint64_t> mElements;
    uint64_t              mNumberOfGenerators;
    uint64_t              mNumberOfBits;
};
