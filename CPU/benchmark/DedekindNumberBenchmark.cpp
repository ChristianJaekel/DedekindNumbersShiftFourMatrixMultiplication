#include <Enumerator.hpp>
#include <benchmark/benchmark.h>
#include <tbb/global_control.h>

static void ComputeDedekindNumber(benchmark::State& state) {
    Enumerator e(state.range(0));

    for (auto _ : state) {
        e.doEnumerationCPU();
    }
}

BENCHMARK(ComputeDedekindNumber)
    ->Unit(benchmark::kMillisecond)
    ->Arg(6)
    ->Arg(7)
    ->Arg(8)
    ->MinTime(10);

static void ComputeDedekindNumberOneThread(benchmark::State& state) {
    tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, 1);

    Enumerator e(state.range(0));

    for (auto _ : state) {
        e.doEnumerationCPU();
    }
}

BENCHMARK(ComputeDedekindNumberOneThread)->Unit(benchmark::kMillisecond)->Arg(6)->Arg(7)->Arg(8);