#ifndef TIMER_H
#define TIMER_H

#include<chrono>

typedef std::chrono::high_resolution_clock hclock;

class Timer {
public:
    static void start() {
        start_time = hclock::now();
    }

    static void stop() {
        end_time = hclock::now();
    }

    // timing in ms
    static double get_timing_in_ms() {
        return std::chrono::nanoseconds(end_time - start_time).count() * 1e-6;
    }

private:
    static std::chrono::time_point<hclock> start_time, end_time;
};

std::chrono::time_point<hclock> Timer::start_time, Timer::end_time;

#endif
