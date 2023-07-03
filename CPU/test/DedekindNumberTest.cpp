#include <Enumerator.hpp>
#include <Uint128.hpp>
#include <gtest/gtest.h>

static uint128T getDedekindNumber(unsigned int n) {
    switch (n) {
        case 0:
            return 2;
        case 1:
            return 3;
        case 2:
            return 6;
        case 3:
            return 20;
        case 4:
            return 168;
        case 5:
            return 7581;
        case 6:
            return 7828354;
        case 7:
            return 2414682040998;
        case 8: {
            const auto value =
                (static_cast<uint128T>(3042) << 64) |
                static_cast<uint128T>(15441756463101891916); // = 56130437228687557907788;

            return value;
        }
        default:
            return 0;
    }
}

class DedekindNumberComputation : public testing::TestWithParam<unsigned int> {};

TEST_P(DedekindNumberComputation, DedekindNumber) {
    auto e = Enumerator(GetParam());

    const auto result = e.doEnumerationCPU();

    ASSERT_TRUE(result == getDedekindNumber(GetParam()));
}

INSTANTIATE_TEST_SUITE_P(ParamValues, DedekindNumberComputation, testing::Values(6, 7, 8));
