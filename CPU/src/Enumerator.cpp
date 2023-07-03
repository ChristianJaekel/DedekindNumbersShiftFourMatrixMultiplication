#include "Enumerator.hpp"

#include "DataDir.h"
#include <FDL.hpp>
#include <Timer.hpp>
#include <fstream>
#include <iostream>
#include <tbb/combinable.h>
#include <tbb/parallel_for.h>

#include <Eigen/Dense>

using namespace Eigen;

using MatU128T = Matrix<uint128T, Dynamic, Dynamic>;
using MatDT    = Matrix<double, Dynamic, Dynamic>;

Enumerator::Enumerator(unsigned int n) {
    n = n - 4;

    mNumOfGens = n;

    FreeDistLat fdln(n);
    mFdl = fdln.getElements();

    mIdx.resize(mFdl.back() + 1);

    for (uint8_t j = 0; j < mFdl.size(); ++j)
        mIdx[mFdl[j]] = j;

    mIvlCalc = IVLCalc(n, 2);
}

IntervalT Enumerator::readIntervals(const std::string& fileName) const {
    IntervalT intervals{};

    std::ifstream t(fileName.c_str());
    if (t.fail())
        throw std::runtime_error("intervals file not found");

    std::string line{};

    while (getline(t, line)) {
        intervals.push_back({});
        std::istringstream ss(line);
        uint64_t           number = 0;
        while (ss >> number)
            intervals.back().push_back(number);
    }

    return intervals;
}

IntervalT Enumerator::readABValues(const std::string& fileName) const {
    IntervalT intervals{};

    std::ifstream t(fileName.c_str());

    if (t.fail())
        throw std::runtime_error("abValues file not found");

    std::string line{};

    uint64_t front = 0;
    uint64_t back  = 0;

    while (getline(t, line)) {
        V64T               v = {};
        std::istringstream ss(line);
        uint64_t           number = 0;

        while (ss >> number)
            v.push_back(number);

        if (v.size() > 1) {
            intervals.push_back({front, v[0], v[1]});
            if (intervals.size() > 1 && intervals[intervals.size() - 2].front() == front) {
                intervals.back()[1] += intervals[intervals.size() - 2][1];
            }
        } else {
            front = v[0];
        }
    }

    return intervals;
}

V16T Enumerator::generateInterval(uint64_t x, uint64_t y) const {
    V16T interval{};
    for (unsigned i = index(x); i <= index(y); ++i)
        if ((x == (x & mFdl[i])) && (y == (y | mFdl[i])))
            interval.push_back(static_cast<uint16_t>(mFdl[i]));

    return interval;
}

std::vector<IntervalT> Enumerator::clusterIntervals(const IntervalT& intervals) {
    std::vector<IntervalT> separatedIntervals;

    uint64_t topElement = std::numeric_limits<uint64_t>::max();
    for (const auto& i : intervals) {
        if (topElement != i.front()) {
            if (separatedIntervals.size() > 0)
                std::sort(separatedIntervals.back().begin(), separatedIntervals.back().end(),
                          [](const V64T& a, const V64T& b) { return a[1] < b[1]; });

            separatedIntervals.push_back(IntervalT{});
            topElement = i.front();
        }

        separatedIntervals.back().push_back(i);
    }

    std::sort(separatedIntervals.back().begin(), separatedIntervals.back().end(),
              [](const V64T& a, const V64T& b) { return a[1] < b[1]; });

    return separatedIntervals;
}

uint128T Enumerator::doEnumerationCPU(bool verbose /* = false */) {
    uint128T sum = 0;

    // load intervals
    const auto allIntervals =
        readIntervals(DATA_DIR + std::to_string(mNumOfGens) + "/intervals.txt");
    const auto intervals = clusterIntervals(allIntervals);

    size_t totalCounter = 1;

    // backwards iterate over all interval clusters --> smaller intervals are treated firstly
    for (int i = intervals.size() - 1; i >= 0; --i) {
        const auto bottom = intervals[i].front().front();

        // values for bottom operator
        V64T distBottom(mFdl.size());
        tbb::parallel_for(tbb::blocked_range<unsigned int>(index(bottom), mFdl.size()),
                          [this, &distBottom, bottom](const tbb::blocked_range<unsigned int>& r) {
                              for (unsigned int j = r.begin(); j < r.end(); j++)
                                  distBottom[j] = mIvlCalc.intervalLength({bottom, mFdl[j]});
                          });

        uint128T clusterInterval = 0;
        for (auto&& interval : intervals[i]) {
            Timer t{};

            const auto top = interval[1];

            // values for top operator
            V64T distTop(mFdl.size());
            tbb::parallel_for(tbb::blocked_range<unsigned int>(index(bottom), index(top) + 1),
                              [this, &distTop, top](const tbb::blocked_range<unsigned int>& r) {
                                  for (unsigned int j = r.begin(); j < r.end(); j++)
                                      distTop[j] = mIvlCalc.intervalLength({mFdl[j], top});
                              });

            const auto elements = generateInterval(bottom, top);

            const auto s = elements.size();

            const auto abValues =
                readABValues(DATA_DIR + std::to_string(mNumOfGens) + "/Values_" +
                             std::to_string(bottom) + "_" + std::to_string(interval[1]) + ".txt");

            if (verbose) {
                std::cout << "Total Counter:  " << totalCounter << std::endl;
                std::cout << "Interval:       " << bottom << " " << top << std::endl;
                std::cout << "Nr of Matrices: " << abValues.size() << std::endl;
            }

            tbb::combinable<uint128T> batchSum{};
            tbb::parallel_for(
                tbb::blocked_range<unsigned int>(0, abValues.size()),
                [this, &batchSum, &abValues, &elements, &distBottom, &distTop,
                 s](const tbb::blocked_range<unsigned int>& r) {
                    MatDT    A(s, s);
                    MatDT    B(s, s);
                    MatU128T C(s, s);

                    for (unsigned int j = r.begin(); j < r.end(); j++) {
                        const auto a = abValues[j][0];
                        const auto b = abValues[j][1];

                        // fill matrices
                        for (unsigned int ind1 = 0; ind1 < s; ++ind1) {
                            const auto c   = elements[ind1];
                            const auto cMa = c & a;
                            const auto cJa = c | a;
                            const auto cMb = c & b;
                            const auto cJb = c | b;

                            for (unsigned int ind2 = ind1; ind2 < s; ++ind2) {
                                const auto d     = elements[ind2];
                                const auto aMcMd = d & cMa;
                                const auto bJcJd = d | cJb;

                                const auto bMcMd = d & cMb;
                                const auto aJcJd = d | cJa;

                                const auto alpha = distBottom[index(aMcMd)] * distTop[index(bJcJd)];
                                const auto beta  = distBottom[index(bMcMd)] * distTop[index(aJcJd)];

                                A.col(ind1)[ind2] = alpha;
                                B.col(ind1)[ind2] = beta;
                            }
                        }

                        // compute matrix product
                        A = A.selfadjointView<Lower>();
                        B = B.selfadjointView<Lower>();

                        C = (A * B).cast<uint128T>();

                        // compute trace of C^2
                        uint128T trace = 0;
                        for (unsigned idx = 0; idx < s; ++idx) {
                            trace += (static_cast<uint128T>(2) *
                                      (C.row(idx)
                                           .tail(s - idx - 1)
                                           .dot(C.col(idx).tail(s - idx - 1).transpose())));
                            trace += (C.col(idx)[idx] * C.row(idx)[idx]);
                        }

                        batchSum.local() += trace * static_cast<uint128T>(abValues[j][2]);
                    }
                },
                tbb::static_partitioner());

            const uint128T batchSumWeighted =
                batchSum.combine(std::plus<uint128T>()) * static_cast<uint128T>(interval[2]);
            clusterInterval += batchSumWeighted;

            if (verbose) {
                std::cout << "Interval Batch Sum: " << batchSumWeighted << std::endl;
                std::cout << "Cluster Sum:        " << clusterInterval << std::endl;
                std::cout << "Accumulated Sum:    " << sum + clusterInterval << std::endl;
                t.stop();
                std::cout << std::endl;
            }

            totalCounter++;
        }
        sum += clusterInterval;
    }
    return sum;
}
