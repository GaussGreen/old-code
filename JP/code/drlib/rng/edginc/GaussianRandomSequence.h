/************************************************************************
 GaussianRandomSequence.h
 Author: M.Huq, Derivatives Research
 
 C++ version of gasdev2
 ************************************************************************/
#ifndef _GAUSSIANRANDOMSEQUENCE__H_
#define _GAUSSIANRANDOMSEQUENCE__H_


#include "edginc/Sequence.h"
#include "edginc/UniformRandomSequence.h"
#include "edginc/SCException.h"
#include "edginc/DECLARESP.h"
//#include <math.h>

CORE_BEGIN_NAMESPACE

class  RNG_DLL GaussianRandomSequence : public Sequence
{
private:
    long   m_initSeed;
    long   m_cachedAvail;
    double m_cachedValue;
    long   m_a0;
    long   m_a1;
    long   m_a2;
    long   m_a3;
    // Internal uniform random sequence
    UniformRandomSequence m_uniformSeq;
    double fetchDeviate();

    long seed;

public:
    GaussianRandomSequence()
    {
        //     vector = 0;
        //     nDimensions = 0;
        seed = 0;
        m_initSeed = 0;
        m_cachedAvail = 0;
        m_cachedValue = 0.0;
        m_a0 = 1;
        m_a1 = 2147483647;
        m_a2 = 0;
        m_a3 = 1;
    };
    GaussianRandomSequence(const long _seed,
                           const int numDimensions);
    virtual void initializeObject(const long _seed,
                                  const int numDimensions);
    virtual ~GaussianRandomSequence();

    // Usual methods overloaded from Sequence
    virtual void populateVector();
    virtual void populateVector(int iPath);
    virtual void advanceSequence(const int numberElements);
    virtual void reinitialize(const long _seed);
    virtual void setSeedParameters(long _a0, long _a1, long _a2, long _a3);

    long getSeed() const
    {
        return seed;
    }
    void setSeed(long thisSeed)
    {
        seed =thisSeed;
    }
private:

    void skipDeviate();
    void skipVector();


};
DECLARESP(GaussianRandomSequence);
CORE_END_NAMESPACE

#endif// _GAUSSIANRANDOMSEQUENCE__H_
