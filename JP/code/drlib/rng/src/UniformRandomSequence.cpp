/*************************************************************************
 UniformRandomSequence.cpp
 Author: M.Huq, CPG DR
 
 Based on Numerical Recipies ran2
 *************************************************************************/

#include "edginc/UniformRandomSequence.h"
#include "edginc/SCException.h"

CORE_BEGIN_NAMESPACE

// Internal constants
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

static const double INVQ1 = 1./IQ1;
static const double INVQ2 = 1./IQ2;
static const long   AQR1 = (IA1*IQ1 + IR1);
static const long   AQR2 = (IA2*IQ2 + IR2);
static const double RNMXbyAM = RNMX / AM ;

void UniformRandomSequence::initializeVariables()
{
    long j,k;
    m_idum2 = 123456789;
    m_iy    = 0;
    // Initialize random number generator
    // Handle negative seeds
    if (-seed < 1)
        seed=1;
    else
        seed = -seed;
    m_idum2=seed;
    // Initialize iv table
    for (j=NTAB+7;j>=0;j--) {
        k=seed/IQ1;
        seed=IA1*(seed-k*IQ1)-k*IR1;
        if (seed < 0)
            seed += IM1;
        if (j < NTAB)
            m_iv[j] = seed;
    }
    m_iy=m_iv[0];
}// UniformRandomSequence::initializeVariables()

void UniformRandomSequence::initializeObject(long theSeed, int numDimensions)
{
    seed = theSeed;
    m_initSeed = seed;

    nDimensions = numDimensions;
    vector.assign(numDimensions, 0);

    m_a0 = 1;
    m_a1 = 2147483647;
    m_a2 = 0;
    m_a3 = 1;

    if(seed == 0) {
        throw SCException(__FILE__, __LINE__, "Bad value for seed (=0)");
    }
    if(nDimensions <=0) {
        throw SCException(__FILE__, __LINE__, "Bad value for numDimensions (<=0)");
    }
    /*  if(!vector)
        vector = new double [nDimensions];*/
    //
    initializeVariables();
}// UniformRandomSequence::initializeObject(long theSeed, int numDimensions)

// Method for reinitializing sequence
void UniformRandomSequence::reinitialize(const long _seed)
{
    seed = _seed;
    m_initSeed = seed;

    if(seed == 0) {
        throw SCException(__FILE__, __LINE__, "Bad value for seed (=0)");
    }

    if(vector.empty())
        throw SCException(__FILE__, __LINE__,
                          "UniformRandomSequence::reinitialize : vector not initialized. Use initializeObject instead.");

    initializeVariables();

}// UniformRandomSequence::reinitialize(const long _seed)


double UniformRandomSequence::fetchDeviate()
{

    long j,k1, k2, s1;
    //   double temp;

    k1= seed /IQ1;
    seed=IA1*seed - k1*AQR1;
    s1 = seed + IM1;
    if (seed < 0)
        seed = s1;

    k2=  m_idum2 /IQ2;
    m_idum2=IA2*m_idum2-k2*AQR2;
    if ( m_idum2 < 0)
        m_idum2 += IM2;

    j=m_iy/NDIV;
    m_iy=m_iv[j]-m_idum2;
    m_iv[j] = seed;
    if (m_iy < 1)
        m_iy += IMM1;

    if (m_iy < RNMXbyAM)
        return AM*m_iy;
    else
        return RNMX;

}// UniformRandomSequence::fetchDeviate()
//
// No reference to path here.
void UniformRandomSequence::populateVector()
{
    for(int iFactor = 0; iFactor < nDimensions; iFactor++) {
        vector[iFactor] = fetchDeviate();
    }//iFactor
}// void UniformRandomSequence::populateVector()

// Method for advancing sequence up by numberElements.
// This makes sense only in the context of using populateVector()
void UniformRandomSequence::advanceSequence(const int numberElements)
{
    if(nDimensions == 0)
        throw SCException (__FILE__, __LINE__,
                           "UniformRandomSequence not initialized : nDimensions =0");
    double *tmpx;

    if (numberElements<=0)
        return;

    for(int iSeq = 0; iSeq < numberElements-1; iSeq++) {
        skipVector();
    }
    // make one regular populate just in case the next call starts reading the vector[]
    populateVector();
    tmpx = getVector();

}// void UniformRandomSequence::advanceSequence(const int numberElements)

//===========================================================
// populateVector(const int iPath)
// Populates vector of random deviates based on path index.
// Seed is deterministically determined for each path
void UniformRandomSequence::populateVector(int iPath)
{

    seed = (abs(m_initSeed) +  m_a0* (long)(1+iPath) )% m_a1 + m_a2*(long)iPath + m_a3;
    seed *= -1;

    initializeVariables();

    for(int iFactor = 0; iFactor < nDimensions; iFactor++) {
        vector[iFactor] = fetchDeviate();
    }//iFactor
}// void UniformRandomSequence(const int iPath)

UniformRandomSequence::~UniformRandomSequence()
{
    //   if(vector){
    //     delete [] vector;
    //     vector = 0;
    //   }
}

void UniformRandomSequence::setSeedParameters(long _a0,
        long _a1,
        long _a2,
        long _a3)
{
    m_a0 = _a0;
    m_a1 = _a1;
    m_a2 = _a2;
    m_a3 = _a3;
}


void UniformRandomSequence::skipDeviate()
{
    long j,k1, k2, s1;
    //     double temp;

    k1= seed /IQ1;
    seed=IA1*seed - k1*AQR1;
    s1 = seed + IM1;
    if (seed < 0)
        seed = s1;

    k2=  m_idum2 /IQ2;
    m_idum2=IA2*m_idum2-k2*AQR2;
    if ( m_idum2 < 0)
        m_idum2 += IM2;

    j=m_iy/NDIV;
    m_iy=m_iv[j]-m_idum2;
    m_iv[j] = seed;
    if (m_iy < 1)
        m_iy += IMM1;

}


void UniformRandomSequence::skipVector()
{
    for(int iFactor = 0; iFactor < nDimensions; iFactor++) {
        skipDeviate();
    }
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

CORE_END_NAMESPACE
