#include <esl_timer.h>

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

double eslNow()
{
    static LARGE_INTEGER count, freq;
    static int initflag;

    if (!initflag)
    {
        QueryPerformanceFrequency(&freq);
        initflag++;
    }

    QueryPerformanceCounter(&count);
    return (double)count.QuadPart / (double)freq.QuadPart;
}

#else
#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>

double eslNow()
{
    struct timeval tv;
    gettimeofday(&tv, 0);
    return (double)tv.tv_sec + (double)tv.tv_usec / 1000000;
}

#endif
