#include <benchmark/benchmark.h>
#include <cmath>
#include <cstring>
#include <random>
#include "math/hilbert.h"

static void BM_Hilbert2D(benchmark::State& state) {
  std::random_device                     r;
  std::default_random_engine             e1(r());
  std::uniform_real_distribution<double> uniform_dist(0, 1);

  for(auto _ : state) {
    //    state.PauseTiming(); // takes more time than the random num gen...
    VectorD<2> quadrant1 = {uniform_dist(e1), uniform_dist(e1)}; // only relevant at hilbertLevel = 1
    //    state.ResumeTiming(); // takes more time than the random num gen...
    hilbert::index<2>(quadrant1, state.range(0)); //linear behaviour -> level * 1.2 ns
  }
}
// Register the function as a benchmark
BENCHMARK(BM_Hilbert2D)->DenseRange(1, 10, 1);

static void BM_Hilbert3D(benchmark::State& state) {
  std::random_device                     r;
  std::default_random_engine             e1(r());
  std::uniform_real_distribution<double> uniform_dist(0, 1);

  for(auto _ : state) {
    //    state.PauseTiming(); // takes more time than the random num gen...
    VectorD<3> quadrant1 = {uniform_dist(e1), uniform_dist(e1), uniform_dist(e1)}; // only relevant at hilbertLevel = 1
    //    state.ResumeTiming(); // takes more time than the random num gen...
    hilbert::index<3>(quadrant1, state.range(0)); //linear behaviour -> level * 1.4 ns
  }
}
// Register the function as a benchmark
BENCHMARK(BM_Hilbert3D)->DenseRange(1, 10, 1);

// example
static void BM_memcpy(benchmark::State& state) {
  char* src = new char[state.range(0)];
  char* dst = new char[state.range(0)];
  memset(src, 'x', state.range(0));
  for(auto _ : state)
    memcpy(dst, src, state.range(0));
  state.SetBytesProcessed(int64_t(state.iterations()) * int64_t(state.range(0)));
  delete[] src;
  delete[] dst;
}
BENCHMARK(BM_memcpy)->Arg(8)->Arg(64)->Arg(512)->Arg(1 << 10)->Arg(8 << 10);

BENCHMARK_MAIN();