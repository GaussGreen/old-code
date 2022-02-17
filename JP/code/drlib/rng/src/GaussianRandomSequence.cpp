/************************************************************************
 GaussianRandomSequence.cpp
 Author: M.Huq, Derivatives Research
 
 C++ version of gasdev2
 ************************************************************************/

#include "edginc/GaussianRandomSequence.h"
#include "edginc/UniformRandomSequence.h"

CORE_BEGIN_NAMESPACE

GaussianRandomSequence::GaussianRandomSequence(const long _seed,
        const int  _nDimensions)
{
    initializeObject(_seed, _nDimensions);
}/*
GaussianRandomSequence::GaussianRandomSequence(const long _seed,
            const int  _nDimensions);
  */

void GaussianRandomSequence::initializeObject(const long _seed,
        const int _nDimensions)
{
    seed = _seed;
    m_initSeed = seed;
    nDimensions = _nDimensions;
    vector.assign(_nDimensions, 0.0);
    m_a0 = 1;
    m_a1 = 2147483647;
    m_a2 = 0;
    m_a3 = 1;

    if(seed == 0) {
        throw SCException(__FILE__, __LINE__, "Bad value for _seed (=0)");
    }
    if(nDimensions <=0) {
        throw SCException(__FILE__, __LINE__, "Bad value for _nDimensions (<=0)");
    }
    /*  if(!vector)
        vector = new double [nDimensions];*/
    //

    m_cachedAvail = 0;
    // Initialize the uniform random number sequence as a 2-d sequence.
    m_uniformSeq.initializeObject(seed, 2);

}// UniformRandomSequence::initializeObject(long theSeed, int numDimensions)

// Method to fetch a single gasdev deviate
double GaussianRandomSequence::fetchDeviate()
{
    double fac,rsq,v1,v2;

    if  (m_cachedAvail == 0) {
        do {
            m_uniformSeq.populateVector();
            seed = m_uniformSeq.getSeed();
            double *uniVect = m_uniformSeq.getVector();
            v1=2.0*uniVect[0]-1.0;
            v2=2.0*uniVect[1]-1.0;
            rsq=v1*v1+v2*v2;
        } while (rsq >= 1.0 || rsq == 0.0);
        fac=sqrt(-2.0*log(rsq)/rsq);
        m_cachedValue=v1*fac;
        m_cachedAvail=1;
        return v2*fac;
    } else {
        m_cachedAvail=0;
        return m_cachedValue;
    }
}// double GaussianRandomSequence::fetchDeviate()

// Method for populating vector
void GaussianRandomSequence::populateVector()
{
    for(int iFactor = 0; iFactor < nDimensions; iFactor++) {
        vector[iFactor] = fetchDeviate();
    }//iFactor
}// void GaussianRandomSequence::populateVector()


//===========================================================
// populateVector(const int iPath)
// Populates vector of random deviates based on path index.
// Seed is deterministically determined for each path
void GaussianRandomSequence::populateVector(int iPath)
{

    seed = (abs(m_initSeed) +  m_a0* (long)(1+iPath) )% m_a1 + m_a2*(long)iPath + m_a3;
    seed *= -1;

    // reinitialize the uniform random number generator
    m_uniformSeq.reinitialize(seed);
    m_cachedAvail = 0;

    for(int iFactor = 0; iFactor < nDimensions; iFactor++) {
        vector[iFactor] = fetchDeviate();
    }//iFactor
}// void GaussianRandomSequence(const int iPath)
// Method for advancing sequence up by numberElements.
// This makes sense only in the context of using populateVector()
void GaussianRandomSequence::advanceSequence(const int numberElements)
{
    if(nDimensions == 0)
        throw SCException (__FILE__, __LINE__,
                           "GaussianRandomSequence not initialized : nDimensions =0");
    double *tmpx;

    if (numberElements<=0)
        return;

    for(int iSeq = 0; iSeq < numberElements-1; iSeq++) {
        skipVector();
    }//iSeq
    populateVector();
    tmpx = getVector();

}// void GaussianRandomSequence::advanceSequence(const int numberElements)


// Method for reinitializing sequence
void GaussianRandomSequence::reinitialize(const long _seed)
{
    seed = _seed;
    m_initSeed = seed;

    if(seed == 0) {
        throw SCException(__FILE__, __LINE__, "Bad value for seed (=0)");
    }

    if(vector.empty())
        throw SCException(__FILE__, __LINE__,
                          "UniformRandomSequence::reinitialize : vector not initialized. Use initializeObject instead.");

    m_uniformSeq.reinitialize(seed);
    m_cachedAvail = 0;

}// GaussianRandomSequence::reinitialize(const long _seed )


GaussianRandomSequence::~GaussianRandomSequence()
{
    /*  if(vector){
        delete [] vector;
        vector = 0;
      }*/
}


void GaussianRandomSequence::setSeedParameters(long _a0, long _a1, long _a2, long _a3)
{
    m_a0 = _a0;
    m_a1 = _a1;
    m_a2 = _a2;
    m_a3 = _a3;
}



// This version skips pairs of deviates!
void GaussianRandomSequence::skipDeviate()
{
    double /*fac,*/rsq,v1,v2;

    do {
        m_uniformSeq.populateVector();
        seed = m_uniformSeq.getSeed();
        double *uniVect = m_uniformSeq.getVector();
        v1=2.0*uniVect[0]-1.0;
        v2=2.0*uniVect[1]-1.0;
        rsq=v1*v1+v2*v2;
    } while (rsq >= 1.0 || rsq == 0.0);

    m_cachedAvail=0;

}

void GaussianRandomSequence::skipVector()
{

    int iFactor = 0;
    if (nDimensions>0 && m_cachedAvail==1) {
        fetchDeviate();
        iFactor++;
    }
    for(; iFactor < nDimensions-2; iFactor += 2) {
        //now skip in pairs
        skipDeviate();
    }
    for(; iFactor < nDimensions; iFactor++) {
        //now skip as usual
        fetchDeviate();
    }

}// void GaussianRandomSequence::populateVector()

CORE_END_NAMESPACE
