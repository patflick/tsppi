/**
 * @file    cputimer.h
 * @author  Patrick Flick <patrick.flick@gmail.com>
 * @brief   Provides a timing stop watch for benchmarking code.
 *
 * Timing C++ code
 * Copyright 2014 Patrick Flick
 *
 * Licensed under the "THE BEER-WARE LICENSE" (Revision 42):
 * Patrick Flick wrote this file. As long as you retain this notice you
 * can do whatever you want with this stuff. If we meet some day, and you think
 * this stuff is worth it, you can buy me a beer or coffee in return
 */
#ifndef CPU_TIMER_H
#define CPU_TIMER_H

#include <time.h>
#include <stdint.h>
// enable time measurements on MAC OS
#ifdef __MACH__
#include <mach/clock.h>
#include <mach/mach.h>
#endif

class CPUTimer {

    private:
        uint64_t starttime;
        uint64_t measured_time;

        inline long long microsecsonds_ts() {
            struct timespec ts;
#ifdef __MACH__ // OS X does not have clock_gettime, use clock_get_time
            clock_serv_t cclock;
            mach_timespec_t mts;
            host_get_clock_service(mach_host_self(), SYSTEM_CLOCK, &cclock);
            clock_get_time(cclock, &mts);
            mach_port_deallocate(mach_task_self(), cclock);
            ts.tv_sec = mts.tv_sec;
            ts.tv_nsec = mts.tv_nsec;
#else
            clock_gettime(CLOCK_MONOTONIC, &ts);
#endif
            return (uint64_t)ts.tv_sec * 1000000LL + (uint64_t)ts.tv_nsec / 1000LL;
        }
    public:
        inline void start() {
            starttime = microsecsonds_ts();
        }

        inline void stop() {
            measured_time = microsecsonds_ts() - starttime;
        }

        // Returns time in microseconds
        uint64_t getTimeMicro() {
            return measured_time;
        }

        // returns time in seconds as double
        double getTime() {
            return (double)measured_time/1000000.0;
        }

};

#endif
