#pragma once
#include <chrono>
#include <iosfwd>
#include <string>

using TimePoint = std::chrono::high_resolution_clock::time_point;

/// @brief Timer class to measure the time of a process and print it to the console or a file.
class Timer {
public:
    /// @brief Default constructor. Starts the timer.
    /// @details Sets the process name to "It".
    Timer();

    /// @brief Constructor with process name. Starts the timer.
    /// @param processName Name of the process.
    Timer(const std::string& processName);

    /// @brief Stops the timer and prints the time to the console.
    void stop();

    /// @brief Stops the timer and prints the time to given ofstream.
    void stop(std::ofstream& os);

    /// @brief Restarts the timer.
    void start();

    /// @brief Restarts the timer and sets the process name.
    /// @param processName Name of the process.
    void start(const std::string& processName);

private:
    /// @brief Converts the time to the appropriate unit.
    std::pair<double, std::string> getTime(double time);

    TimePoint   mStartTime;
    std::string mProcessName;
};
