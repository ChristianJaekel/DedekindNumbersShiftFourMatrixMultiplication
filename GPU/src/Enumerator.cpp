#include "Enumerator.hpp"
#include "Timer.hpp"
#include <FDL.hpp>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <tbb/parallel_for.h>

Enumerator::Enumerator(unsigned int n) {
    mNumOfGens = n;

    FreeDistLat fdln(n);
    mFdl = fdln.getElements();

    mIdx.resize(mFdl.back() + 1);

    for (uint8_t j = 0; j < mFdl.size(); ++j)
        mIdx[mFdl[j]] = j;

    mIvlCalc = IVLCalc(n, 2);

    mDE = std::make_unique<DeviceEnumerator>(mIdx);
}

IntervalT Enumerator::readIntervals(const std::string& fileName) const {
    IntervalT intervals;

    std::ifstream t(fileName.c_str());
    std::string   line;

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
    IntervalT intervals;

    std::ifstream t(fileName.c_str());
    std::string   line;

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

uint128T Enumerator::doEnumerationGPU(bool verbose /* = false */) {
    uint128T sum = 0;

    // load intervals
    const auto allIntervals =
        readIntervals(DATA_DIR + std::to_string(mNumOfGens) + "/intervals.txt");
    const auto intervals = clusterIntervals(allIntervals);

    const auto memory = mDE->getAvailableGPUMemory();

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
        mDE->setBottom(distBottom);

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
            mDE->setTop(distTop);

            const auto elements = generateInterval(bottom, top);
            mDE->setElements(elements);

            const auto s  = elements.size();
            const auto ss = s * s;

            // nr of batches that gpu memory can hold (memory for indices, elements, bottom operator
            // and top operator has to be deducted)
            const auto possibleNumberOfBatches =
                memory / (3 * ss * sizeof(double)) - (mFdl.back() + 1) * sizeof(uint8_t) -
                mFdl.size() * sizeof(uint16_t) - 2 * mFdl.size() * sizeof(uint64_t);

            const auto abValues =
                readABValues(DATA_DIR + std::to_string(mNumOfGens) + "/Values_" +
                             std::to_string(bottom) + "_" + std::to_string(interval[1]) + ".txt");

            const auto abValuesSize = abValues.size();

            mDE->initMatrices(s, ss * std::min(possibleNumberOfBatches, abValuesSize));

            uint128T batchSum{};

            if (verbose) {
                std::cout << "Total Counter:  " << totalCounter << std::endl;
                std::cout << "Interval:       " << bottom << " " << top << std::endl;
                std::cout << "Nr of Matrices: " << abValuesSize << std::endl;
            }

            if (possibleNumberOfBatches < abValuesSize) {
                for (size_t j = 0; j < abValuesSize; j += possibleNumberOfBatches)
                    batchSum += mDE->doEnumerationGPU(
                        j, std::min(j + possibleNumberOfBatches, abValuesSize), abValues, s);

            } else {
                batchSum += mDE->doEnumerationGPU(0, abValuesSize, abValues, s);
            }

            const uint128T batchSumWeighted = batchSum * static_cast<uint128T>(interval[2]);
            clusterInterval += batchSumWeighted;

            if (verbose) {
                std::cout << "Interval Batch Sum: " << batchSumWeighted << std::endl;
                std::cout << "Cluster Sum:        " << clusterInterval << std::endl;
                std::cout << "Accumulated Sum:    " << sum + clusterInterval << std::endl;
                t.stop();
                std::cout << std::endl;
            }

            mDE->freeMatrices();
            mDE->freeTop();
            mDE->freeElements();
            totalCounter++;
        }
        sum += clusterInterval;
        mDE->freeBottom();
    }
    return sum;
}
