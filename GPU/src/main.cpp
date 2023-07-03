#include "Enumerator.hpp"
#include <FDL.hpp>
#include <Timer.hpp>
#include <iostream>

/// @brief main function
/// @details If a command-line argument is provided, the program will print
/// additional information about the enumeration process.
int main(int argc, char* argv[]) {
    bool verbose = false;

    // check if there are command-line arguments
    if (argc > 1) {
        verbose = true;
    }

    // enumeration of FDL(6), FDL(7) and FDL(8)
    for (unsigned i : {6, 7, 8}) {
        Enumerator e(i);

        Timer t("FDL(" + std::to_string(i) + ")");

        const auto nrOfElements = e.doEnumerationGPU(verbose);
        t.stop();

        std::cout << "Nr of Elements: " << nrOfElements << std::endl;
        std::cout << std::endl;
    }

    return 0;
}
