/*************************************************************************
 UniformRandomSequence.h
 Author: M.Huq, CPG DR
 
 Definition for uniform pseudo-random sequence class - modified version
 of ran2 which does not require static variables. Eliminates some of the
 problems associated with uninitialized seeds etc.
 *************************************************************************/
#ifndef _UNIFORMRANDOMSEQUENCE__H_
#define _UNIFORMRANDOMSEQUENCE__H_

#include "edginc/Sequence.h"

CORE_BEGIN_NAMESPACE

class RNG_DLL UniformRandomSequence : public Sequence
{
protected:
    // variables
    long m_idum2;
    long m_iy;
    long m_iv[32];  // This table needs to be initialized on first use.
    long m_initSeed;
    long m_a0; // Parameters for random seed when using populateVector(int)
    long m_a1;
    long m_a2;
    long m_a3;
    // methods unique to UniformRandomSequence
    virtual void initializeVariables();
    virtual double fetchDeviate();


public:
    // Constructors
    UniformRandomSequence()
    {
        seed=-1;
        /*    nDimensions=0;
            vector=0; */
        m_a0= 1;
        m_a1= 2147483647;
        m_a2= 0;
        m_a3= 1;
    };
    UniformRandomSequence(long theSeed, int numDimensions)
    {
        initializeObject(theSeed, numDimensions);
    };
    // Initialization method assocated with above constructors
    virtual void initializeObject(long theSeed, int numDimensions);
    // Destructor
    virtual ~UniformRandomSequence();
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

private:

    void skipDeviate();
    void skipVector();

    void setSeed( long thisSeed)
    {
        seed =thisSeed;
    }
    long seed;


};

CORE_END_NAMESPACE

#endif// _UNIFORMRANDOMSEQUENCE__H_

