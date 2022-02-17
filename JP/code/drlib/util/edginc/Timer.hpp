//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Timer.hpp
//
//   Small class that computes elapsed time since its creation
//
//   Author      : Regis S Guichard
//
//   Date        : 11 Jan 2005
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_TIMER_HPP
#define QLIB_TIMER_HPP

#include "edginc/config.hpp"
#include <time.h>

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL Timer{
public:
    Timer();

    double calcTime() const;

private:
#ifdef UNIX
    time_t  startSec;
#else
    clock_t startTime;
#endif
};

// another version of performance counter
// apparently there is a name conflict between windows.h
// and qlib math so i am hiding the data representation
class UTIL_DLL PerfTimer {
    struct Impl;
    Impl*  data;
public:
    PerfTimer();
    ~PerfTimer();

    void start();
    void stop();
    double getTimeLapse();
};

DRLIB_END_NAMESPACE

#endif
