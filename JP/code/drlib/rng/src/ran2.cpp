/***

    C++ wrapping over SC_ran2 algorithm
    (Produces Uniform RNG)

*/
#include "edginc/ran2.h"
#include "edginc/SCException.h"

#include <cstdio>
#include <cassert>
#include <fstream>
#include <cmath>

CORE_BEGIN_NAMESPACE

#ifndef CHECKRNG
#define CHECKRNG 0
#endif

using namespace std;


#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)



#if defined(LINUX)
#include <signal.h>
#include <execinfo.h>

/** show_stackframe is a Linux specific way to print the current call stack. Used mainly for debugging purposes.
 */

static void show_stackframe()
{
    void *trace[16];
    char **messages = (char **)NULL;
    int i, trace_size = 0;

    trace_size = backtrace(trace, 16);
    messages = backtrace_symbols(trace, trace_size);
    printf("[bt] Execution path:\n");
    for (i=0; i<trace_size; ++i)
        printf("[bt] %s\n", messages[i]);
}
#else
static void show_stackframe()
{}
#endif

class RNG_DLL BlackBox
{
public:
    virtual void observe(long seed, double rnd) = 0;
    virtual ~BlackBox()
    {}
}
;

/** BlackBoxReader is a helper class to match numbers produced by different refactorings of the library. It reads RNGs from a file named "seed_rnd.txt" in the current directory, compares it to the one it is about to return and crashes if there is a discrepancy.
 */

class RNG_DLL BlackBoxReader : public BlackBox
{
    std::ifstream in;
public:
    void observe(long seed, double rnd)
    {
        double r;
        long   s;
        long line;
        string s1, s2, s3;
        in >> line >> s >> r;
        in >> s1 >> s2 >> s3;

        if (seed != s || ! (fabs(rnd-r) < 1e-5))
            fprintf(stdout, "[%ld] Expected: (%ld %f) Got: (%ld %f)\n", line, s, r, seed, rnd);

        assert(seed == s);
        assert(fabs(rnd-r) < 1e-5);
        //assert(line != 53101);
    }
    BlackBoxReader(const char * fn = "seed_rnd.txt") : in(fn)
    {
        if (!in)
            assert("Could not open BlackBoxReader file!" == NULL);
    }
};

/** \return uniformly distributed number in [0,1] */
double Ran2Gen::fetch()
{
#if defined(CHECKRNG) && CHECKRNG > 0
    static BlackBoxReader bb;
#endif

    int j;
    long k;
    double temp;
    long seed_in = seed;
    long * idum = & seed;
    double result;




    if (*idum <= 0) {
        if (-(*idum) < 1)
            *idum=1;
        else
            *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0)
                *idum += IM1;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0)
        *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0)
        idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1)
        iy += IMM1;
    if ((temp=AM*iy) > RNMX)
        result = RNMX;
    else
        result = temp;
#if defined (CHECKRNG) && CHECKRNG > 0

    bb.observe(seed_in, result);
#endif

    return result;
}

/////////////////////////// Ran2StaticGen ///////////////////////////////
/** Ran2StaticGen shares inner state between all instances of this class */

long Ran2StaticGen::idum2 = 123456789L;
long Ran2StaticGen::iy = 0L;
long Ran2StaticGen::iv[NTAB];

double Ran2StaticGen::fetch()
{
#if CHECKRNG
    static BlackBoxReader bb;
#endif

    int j;
    long k;
    double temp;
    long seed_in = seed;
    long * idum = & seed;
    double result;

    if (*idum <= 0) {
        if (-(*idum) < 1)
            *idum=1;
        else
            *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0)
                *idum += IM1;
            if (j < NTAB)
                iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0)
        *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0)
        idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1)
        iy += IMM1;
    if ((temp=AM*iy) > RNMX)
        result = RNMX;
    else
        result = temp;
#if CHECKRNG

    bb.observe(seed_in, result);
#endif

    return result;
}



void SC_ran2Gen::superCubeAdjustment()
{

    if( seed > 0 ) {
        seed *= -1;
    }
    if(seed == 0)
        throw SCException(__FILE__,__LINE__,
                          "zero seed passed in. Will lead to problems in pseudo-random number generator");
}

// Dark magic from the original SuperCube.cpp
// the client should clone the generator at the end of the path and then seek to the given path
void SC_ran2Gen::seekToPath(int iPath)
{
    const long seedStride = 1000;
    const long seedMod = 16061968;
    long seed_in = seed;
    seed = (seed +  seedStride* (long)(1+iPath) )% seedMod + iPath;
    if(seed > 0)
        seed *= -1;
    // std::fflush(NULL);
    //    std::fprintf(stdout, "seekToPath: iPath= %d seed= %ld => seed= %ld\n", iPath, seed_in, seed);
}

void SC_ran2StaticGen::superCubeAdjustment()
{

    if( seed > 0 ) {
        seed *= -1;
    }
    if(seed == 0)
        throw SCException(__FILE__,__LINE__,
                          "zero seed passed in. Will lead to problems in pseudo-random number generator");
}

// Dark magic from the original SuperCube.cpp
// the client should clone the generator at the end of the path and then seek to the given path
void SC_ran2StaticGen::seekToPath(int iPath)
{
    const long seedStride = 1000;
    const long seedMod = 16061968;
    long seed_in = seed;
    seed = (seed +  seedStride* (long)(1+iPath) )% seedMod + iPath;
    if(seed > 0)
        seed *= -1;
    // std::fflush(NULL);
    //    std::fprintf(stdout, "seekToPath: iPath= %d seed= %ld => seed= %ld\n", iPath, seed_in, seed);
}

#undef IM1
#undef IM2
#undef AM
#undef IMM1
#undef IA1
#undef IA2
#undef IQ1
#undef IQ2
#undef IR1
#undef IR2
#undef NTAB
#undef NDIV
#undef EPS
#undef RNMX
/* (C) Copr. 1986-92 Numerical Recipes Software *0-12'=. */
CORE_END_NAMESPACE
