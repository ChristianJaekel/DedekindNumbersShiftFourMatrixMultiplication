#include "Uint128.hpp"
#include <iostream>

std::ostream& operator<<(std::ostream& os, uint128T value) {
    // convert the value to a string
    std::string strValue{};
    do {
        char digit = static_cast<char>(value % 10);
        strValue.insert(strValue.begin(), '0' + digit);
        value /= 10;
    } while (value != 0);

    // print the string representation
    os << strValue;

    return os;
}