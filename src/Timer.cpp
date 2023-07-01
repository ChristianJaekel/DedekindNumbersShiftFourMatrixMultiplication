#include "Timer.hpp"
#include <fstream>
#include <iostream>

Timer::Timer() : mProcessName("It"), mStartTime(std::chrono::high_resolution_clock::now()) {}

Timer::Timer(const std::string& processName)
    : mProcessName(processName), mStartTime(std::chrono::high_resolution_clock::now()) {}

void Timer::stop() {
    const auto endTime = std::chrono::high_resolution_clock::now();
    const auto timeDelta =
        std::chrono::duration_cast<std::chrono::milliseconds>(endTime - mStartTime).count();

    const auto t = getTime(timeDelta);
    std::cout << mProcessName + " took " << t.first << t.second;
}

void Timer::stop(std::ofstream& os) {
    const auto endTime = std::chrono::high_resolution_clock::now();
    const auto timeDelta =
        std::chrono::duration_cast<std::chrono::milliseconds>(endTime - mStartTime).count();

    const auto t = getTime(timeDelta);
    os << mProcessName + " took " << t.first << t.second;
}

void Timer::start() { mStartTime = std::chrono::high_resolution_clock::now(); }

void Timer::start(const std::string& processName) {
    mProcessName = processName;
    start();
}

std::pair<double, std::string> Timer::getTime(double time) {
    if (time < 1000) {
        return {time, " milliseconds\n"};
    } else if (time >= 1000 && time < 60000) {
        return {time / 1000, " seconds\n"};
    } else if (time >= 60000 && time < 3600000) {
        return {time / (1000 * 60), " minutes\n"};
    } else {
        return {time / (1000 * 60 * 60), " hours\n"};
    }
}