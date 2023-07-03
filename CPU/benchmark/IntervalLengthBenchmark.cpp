#include <FDL.hpp>
#include <IVLCalc.hpp>
#include <benchmark/benchmark.h>

static void ComputeIntervalLength(benchmark::State& state) {
    const auto ivl = IVLCalc(state.range(0), state.range(1));

    const auto fdl = FreeDistLat(state.range(0)).getElements();

    for (auto _ : state) {
        ivl.intervalLength({0, fdl.back()});
    }
}

BENCHMARK(ComputeIntervalLength)->Unit(benchmark::kNanosecond)->ArgsProduct({{3, 4, 5, 6}, {2, 3}});