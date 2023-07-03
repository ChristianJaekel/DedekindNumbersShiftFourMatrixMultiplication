#include <FDL.hpp>
#include <IVLCalc.hpp>
#include <cstdint>
#include <gtest/gtest.h>

/// @brief class that computes interval length via naive counting
class IntervalLengthForTesting {
public:
    IntervalLengthForTesting(unsigned int n) : mFdl(FreeDistLat(n)) {}

    uint64_t intervalLength(uint64_t bottom, uint64_t top) const {
        uint64_t length = 0;

        for (auto&& x : mFdl.getElements())
            if ((bottom == (x & bottom)) && (top == (x | top)))
                ++length;

        return length;
    }

private:
    FreeDistLat mFdl;
};

class IntervalLengthTest : public ::testing::TestWithParam<std::tuple<unsigned int, unsigned int>> {
};

TEST_P(IntervalLengthTest, IntervalLength) {
    const auto params = GetParam();

    const auto ivl = IVLCalc(std::get<0>(params), std::get<1>(params));

    const auto ivlTest = IntervalLengthForTesting(std::get<0>(params));

    const auto fdl = FreeDistLat(std::get<0>(params));

    for (auto&& bottom : fdl.getElements())
        for (auto&& top : fdl.getElements())
            ASSERT_EQ(ivl.intervalLength({bottom, top}), ivlTest.intervalLength(bottom, top));
}

const std::vector<unsigned int> generators = {3, 4};
const std::vector<unsigned int> difference = {1, 2, 3};

INSTANTIATE_TEST_SUITE_P(ParamValues, IntervalLengthTest,
                         ::testing::Combine(::testing::ValuesIn(generators),
                                            ::testing::ValuesIn(difference)));
