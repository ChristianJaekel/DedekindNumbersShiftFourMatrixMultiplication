#pragma once

#include <iosfwd>

using uint128T = unsigned __int128;

// enables to write unsigned __int128_t to std::cout
std::ostream& operator<<(std::ostream& os, uint128T value);