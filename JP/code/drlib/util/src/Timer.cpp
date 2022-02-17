//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Timer.cpp
//
//   Small class that computes elapsed time since its creation
//
//   Author      : Regis S Guichard
//
//   Date        : 11 Jan 2005
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Timer.hpp"

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#include <unistd.h>
#include <time.h>
#endif

DRLIB_BEGIN_NAMESPACE

#ifdef UNIX
Timer::Timer():
startSec(time(0)){}

double Timer::calcTime() const{
    time_t endSec = time(0);
    // UNIX overflows after 2147 seconds
    return difftime(endSec, startSec);
}
#else
Timer::Timer():
startTime(clock()){}

double Timer::calcTime() const{
    clock_t endTime = clock();
    return static_cast<double>(endTime - startTime) / CLOCKS_PER_SEC;
}
#endif

// another version of performance timer
#ifdef _WIN32

struct PerfTimer::Impl {
    LARGE_INTEGER _tstart;
    LARGE_INTEGER _tend;
    LARGE_INTEGER _freq;
};

void PerfTimer::start()
{
    QueryPerformanceFrequency(&data->_freq);
    QueryPerformanceCounter(&data->_tstart);    
}

void PerfTimer::stop()
{
    QueryPerformanceCounter(&data->_tend);
}

double PerfTimer::getTimeLapse()
{
    return ((double)data->_tend.QuadPart - 
        (double)data->_tstart.QuadPart)/((double)data->_freq.QuadPart);
}

#else

#if defined( _POSIX_TIMERS) && (_POSIX_TIMERS > 0 ) && defined (_POSIX_CPUTIME)
/** POSIX hi-res timer: see man clock_gettime().
    Measures user time of the process, as it's more meaningful than the wall clock time in most cases.
*/

struct PerfTimer::Impl {
  struct timespec _tstart;
  struct timespec _tend;
};

void PerfTimer::start() {
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, & data->_tstart);
}

void PerfTimer::stop()
{
  clock_gettime(CLOCK_PROCESS_CPUTIME_ID, & data->_tend);
}

double PerfTimer::getTimeLapse()
{
  return (data->_tend.tv_sec - data->_tstart.tv_sec) + (double) (data->_tend.tv_nsec - data->_tstart.tv_nsec)/1e+9;
}

#else // !WIN32 and  no POSIX timers available 
struct PerfTimer::Impl {
    struct timeval _tstart;
    struct timeval _tend;
    struct timezone _tz;
};

void PerfTimer::start()
{
    gettimeofday(&data->_tstart, &data->_tz);   
}

void PerfTimer::stop()
{
    gettimeofday(&data->_tend, &data->_tz);
}

double PerfTimer::getTimeLapse()
{
    double t1, t2;
    t1 =  (double)data->_tstart.tv_sec + (double)data->_tstart.tv_usec/(1000*1000);
    t2 =  (double)data->_tend.tv_sec + (double)data->_tend.tv_usec/(1000*1000);
    return t2-t1;
}
#endif // _POSIX_TIMERS
#endif // Win32

PerfTimer::PerfTimer()
{
    data = new Impl();
}

PerfTimer::~PerfTimer()
{
    delete data;
}

DRLIB_END_NAMESPACE
