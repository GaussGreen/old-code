

//----------------------------------------------------------------------------
//
//   Group       : Cross Asset
//
//   Filename    : MCRandomLite.hpp
//
//   Description : lightweight random number generators
//
//   Author      : Henrik Rasmussen
//
//   Date        : 10 April 2006
//
//
//----------------------------------------------------------------------------

#ifndef QLIB_MCRANLITE_H
#define QLIB_MCRANLITE_H

#include "edginc/config.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

const double SMALL_RADIUS = 1.e-10;

// default random number generator seed
const unsigned long DEFAULT_RNG_SEED = 123456789;

// different random number generators
typedef enum {
	no_rng,
	gaussianGfsr4,
	gaussianKnuth,
	gaussianZiggurat,
	gaussianNR,
	gaussianSobol,
	gaussianNiederr,
	randomBit
} RNGName;

// interface class for random number generators
class RanNumClass : public virtual VirtualDestructorBase
{

public:

	// initialize generator to obtain a new sequence
	// of random number
	virtual void setSeed(unsigned long seed) = 0;

	// reset the generator to obtain the same sequence
	// of random numbers
	virtual void reset() = 0;

	virtual unsigned long getSeed() const = 0;

	// get one random number from the generator
	virtual double getRandNumber() =0 ;

	// get one random vector, or dim independent random
	// numbers from the generator
	virtual void getRandVector(double *vec, long dim) =0 ;

	// routines to use log file for writing intermidiate
	// results/testing/debugging
	// void setFileStream(std::ofstream *outFile);
	// virtual void printObjectState() const;

};


// Very long period (10^2917) uniform  random number generator
// It is a generalized feedback shift-register (gfsr) generator
// i.e. the generator's state is updated by taking the XOR-sum of
// past lagged values. An implementation which performs well in tests
// uses four past lagged values i.e. the n-th random number in the sequence
// is obtained from:
// state[n] = state[n-A]^state[n-B]^state[n-C]^state[n-D]
// The state is an array of 2^{14} longs
// For more details see GNU Scientific Library documentation and source code,
// the implementation is based on the G.S.L implementation
// http://gnuwin32.sourceforge.net/packages/gsl.htm

class RngUniformGfsr4Class : public RanNumClass {

public:

	explicit RngUniformGfsr4Class(unsigned long seed=DEFAULT_RNG_SEED);
	~RngUniformGfsr4Class();

	// copy constructor
	RngUniformGfsr4Class(const RngUniformGfsr4Class &otherRng);
	// assignment
	RngUniformGfsr4Class &operator=(const RngUniformGfsr4Class &otherRng);

	// class is copyable, compiler generated assignment operator
	// and copy ctor work fine

	void setSeed(unsigned long seed);
	void reset();
	unsigned long getSeed() const;

	double getRandNumber();
	void getRandVector(double *vec, long dim);

private:

	// backward shifts used in the state update relation 
	unsigned long m_A;
    unsigned long m_B;
    unsigned long m_C;
    unsigned long m_D;

	// size of state array = smallest power of 2 larger than
	// the largest shift D
    unsigned long m_M;

	unsigned long m_seed,
		          m_index; // current index of the state array

	// internal state of the generator: array of 2^14 longs
	unsigned long *m_state;

	// inverse of largest long, we devide the generator output to
	// obtain a real number in [0,1]
	double m_invMax;

	// function used to initialize of state array
	// from the seed
	unsigned long Lcg(unsigned long n) 
	{
		return ( (69069*n) & 0xffffffffUL );
	}

	// initializes the state array
	void initializeState(unsigned long seed);

	// updates the state to produce the next random number in the sequence
	unsigned long getNumber()
	{
		// index must always be within the state array limits
		// we use bitwise AND with the array size to ensure
		// m_index is in 0,...,M
		m_index= (m_index+1)&m_M;

		if ( m_index > m_M )
		{
			throw("Index error in RngUniformGfsr4Class");
		}

        m_state[m_index] = m_state[(m_index+(m_M+1-m_A))&m_M]^
                           m_state[(m_index+(m_M+1-m_B))&m_M]^
                           m_state[(m_index+(m_M+1-m_C))&m_M]^
                           m_state[(m_index+(m_M+1-m_D))&m_M];

		return m_state[m_index];
	}
};


// Gaussian random number generator using the Box-Muller method
// with the gnu scientific library long period Gfsr4 generator
// as the source of uniform deviates
// For more details on Box-Muller method see: Numerical Recipes
//   in C++, chapter 7
class RngGaussianGfsr4Class : public RanNumClass {

public:

	explicit RngGaussianGfsr4Class(unsigned long seed=DEFAULT_RNG_SEED);
	~RngGaussianGfsr4Class();

	// copy constructor
	RngGaussianGfsr4Class(const RngGaussianGfsr4Class &otherRng);
	// assignment
	RngGaussianGfsr4Class &operator=(const RngGaussianGfsr4Class &otherRng);

	// class is copyable, compiler generated assignment operator
	// and copy ctor work fine

	void setSeed(unsigned long seed);
	void reset();
	unsigned long getSeed() const;

	double getRandNumber();
	void getRandVector(double *vec, long dim);

private:

	// backward shifts used in the state update relation 
	unsigned long m_A;
    unsigned long m_B;
    unsigned long m_C;
    unsigned long m_D;

	// size of state array = smallest power of 2 larger than
	// the largest shift D
    unsigned long m_M;

	unsigned long m_seed,
		          m_index; // current index of the state array

	// internal state of the generator: array of 2^14 longs
	unsigned long *m_state;

	// inverse of largest long. We devide to
	// obtain a uniform deviate [0,1]
	double m_invMax;

	// flag which denotes whether we draw a new
	// point for Gaussian deviates or not in the 
	// Box-Muller algorithm
	int m_draw;

	// function used to initialize of state array
	// from the seed
	unsigned long Lcg(unsigned long n) 
	{
		return ( (69069*n) & 0xffffffffUL );
	}

	// initializes the state array
	void initializeState(unsigned long seed);

	// updates the state to produce the next random number in the sequence
	unsigned long getNumber()
	{
		// index must always be within the state array limits
		// we use bitwise AND with the array size to ensure
		// m_index is in 0,...,M
		m_index= (m_index+1)&m_M;

		if ( m_index > m_M )
		{
			throw("Index error in RngGaussianGfsr4Class");
		}

        m_state[m_index] = m_state[(m_index+(m_M+1-m_A))&m_M]^
                           m_state[(m_index+(m_M+1-m_B))&m_M]^
                           m_state[(m_index+(m_M+1-m_C))&m_M]^
                           m_state[(m_index+(m_M+1-m_D))&m_M];
		return m_state[m_index];
	}
};

// Impementation of Uniform  generator using
// Knuth's subtractive method (see numerical recipes in C++) instead
// of linear congruential methods. The idea is to construct each new
// random number in a sequence by taking the difference of of two numbers
// stored in an array and then store it in the array
// see: Numerical Recipes in C++, chapter 7
class RngUniformKnuthClass : public RanNumClass {

public:

	explicit RngUniformKnuthClass(unsigned long seed=DEFAULT_RNG_SEED);
	~RngUniformKnuthClass();

	// copy constructor
	RngUniformKnuthClass(const RngUniformKnuthClass &otherRng);
	// assignment operator
	RngUniformKnuthClass &operator=(const RngUniformKnuthClass &otherRng);

	// class is copyable, compiler generated assignment operator
	// and copy ctor work fine

	void setSeed(unsigned long seed);
	void reset();
	unsigned long getSeed() const;

	double getRandNumber();
	void getRandVector(double *vec, long dim);

private:
	long m_seed,
		 m_rndNum;
	int m_index1, // indicies of random number array
		m_index2;

	// array which holds the last 56 random number in a sequence
	// The value 56 is algorithm specific and should not be changed
	long m_buffer[56];

	// necessary constants used in the algorithm
	long mBig,
		 mSeed;

	// inverse of largest long. We devide to
	// obtain a uniform deviate in [0,1]
	double m_invMax;

	// initialize array of random number
	void initializeState(long seed);

	// get a new random number and update array
	unsigned long getNumber()
	{
		// increase indicies, make sure the are in range
		m_index1++;
	    if ( m_index1==56 )
		{
			m_index1 = 1;
		}

	    m_index2++;
	    if ( m_index2==56 )
		{
			m_index2 = 1;
		}

		// genrate a new random number subtractively
	    m_rndNum = m_buffer[m_index1]-m_buffer[m_index2];
	    if ( m_rndNum < 0 )
		{
			// make sure it is in range
			m_rndNum += mBig;
		}

		// store it in array
	    m_buffer[m_index1]=m_rndNum;

	    return m_rndNum;
	}
};


// Gaussian random number generator using the Box-Muller method
// and Knuth's substractive generator as source of uniform deviates
// see: Numerical Recipes in C++, chapter 7
class RngGaussianKnuthClass : public RanNumClass {

public:

	explicit RngGaussianKnuthClass(unsigned long seed=DEFAULT_RNG_SEED);
	~RngGaussianKnuthClass();

	// copy constructor
	RngGaussianKnuthClass(const RngGaussianKnuthClass &otherRng);
	// assingment operator
	RngGaussianKnuthClass &operator=(const RngGaussianKnuthClass &otherRng);

	// class is copyable, compiler generated assignment operator
	// and copy ctor work fine

	void setSeed(unsigned long seed);
	void reset();
	unsigned long getSeed() const;

	double getRandNumber();
	void getRandVector(double *vec, long dim);

private:
	long m_seed,
		 m_rndNum;
	int m_index1,
		m_index2; // indicies of random number array

	// array which holds the last 56 random number in a sequence
	// The value 56 is algorithm specific and should not be changed
	long m_buffer[56];

	// necessary constants used in the algorithm
	long mBig,
		 mSeed;

	// inverse of largest long. We devide to
	// obtain a uniform deviate in [0,1]
	double m_invMax;

	// flag which denotes whether we draw a new
	// point for Gaussian deviates or not in the 
	// Box-Muller algorithm
	int m_draw;

	// initialize array of random number
	void initializeState(long seed);

	// get a new random number and update array
	unsigned long getNumber()
	{
		// increase indicies, make sure the are in range
		m_index1++;
	    if ( m_index1==56 )
		{
			m_index1 = 1;
		}

	    m_index2++;
	    if ( m_index2==56 )
		{
			m_index2 = 1;
		}

		// genrate a new random number subtractively
	    m_rndNum = m_buffer[m_index1]-m_buffer[m_index2];
	    if ( m_rndNum < 0 )
		{
			// make sure it is in range
			m_rndNum += mBig;
		}

		// store it in array
	    m_buffer[m_index1]=m_rndNum;

	    return m_rndNum;
	}
};


// Gaussian random number generator using the Ziggurat algorithm
// supposed to be much faster than the Box-Muller based generators
// For more details see preprints:
// http://www.jstatsoft.org/v05/i08/ziggurat.pdf
// http://www.doornik.com/research/ziggurat.pdf
class RngGaussianZigguratClass : public RanNumClass {

public:

	explicit RngGaussianZigguratClass(unsigned long seed=DEFAULT_RNG_SEED);
	~RngGaussianZigguratClass();

	// copy constructor
    RngGaussianZigguratClass(const RngGaussianZigguratClass &otherRng);
	// assignment
	RngGaussianZigguratClass &operator=(const RngGaussianZigguratClass &otherRng);

	// class is copyable, compiler generated assignment operator
	// and copy ctor work fine

	void setSeed(unsigned long seed);
	void reset();
	unsigned long getSeed() const;
	double getRandNumber();
	void getRandVector(double *vec, long dim);

private:

	unsigned long m_seed,
		          m_rngULongState,
				  m_rngULongStatePrev;

	// array of 32-bit integers used into the step which forms the x's that
	// are to be returned
	unsigned long m_int32bitGaussian[256];

	// array of doubles which is used to form the random gaussian numbers
	// that are returned
	double m_wGaussian[256];

	// array with values of the normal density (modulo the normalization)
	// at the rectangle points x_i
	double m_gaussianDen[256];

	unsigned long getULong()
    {
	    m_rngULongStatePrev= m_rngULongState;
	    m_rngULongState ^= ( m_rngULongState<<13 );
	    m_rngULongState ^= ( m_rngULongState>>17 );
	    m_rngULongState ^= ( m_rngULongState<<5 );

       return m_rngULongState+m_rngULongStatePrev;
	}
   
    double getUniform()
	{
		return ( getULong()*0.465661e-9 );
	}
};



// Another Box-Muller gaussian random number generator
// using the long period ran2 from numerical recipies in C++
// as source of uniform variates
// see: Numerical Recipes in C++, chapter 7
class RngGaussianNRClass : public RanNumClass {

public:

	explicit RngGaussianNRClass(unsigned long seed=DEFAULT_RNG_SEED);
	~RngGaussianNRClass();

	// copy constructor
	RngGaussianNRClass(const RngGaussianNRClass &otherRng);
	// assignment
	RngGaussianNRClass &operator=(const RngGaussianNRClass &otherRng);

	// class is copyable, compiler generated assignment operator
	// and copy ctor work fine

	void setSeed(unsigned long seed);
	void reset();
	unsigned long getSeed() const;
	double getRandNumber();
	void getRandVector(double *vec, long dim);

private:

	// constants of the two linear congruential generators
	// in getNumber()
	long m_im1, m_im2, m_ia1, m_ia2, m_iq1, m_iq2,
		          m_ir1,m_ir2,m_imm1,m_ndiv,m_iy;

	double m_eps, m_rnmx, m_am;

	// shuffle table
	long m_iv[32];

	long m_seed,    // seed of combined generator
		 m_state1,  // state of linear congruential generator 1
		 m_state2;  // state of linear congruential generator 2

	// flag which denotes whether we draw a new
	// point for Gaussian deviates or not in the 
	// Box-Muller algorithm
	int m_draw;

	void initializeState(long seed);

	// combines two linear congruential generators to incease the
	// period and shuffles the output to reduce serial correlations
	double getNumber()
	{
		long k,j;
		double temp;

		// compute state = ( ia * state) % im 
		// with out overflow (see numerical Recipies in C++)
		k = m_state1 / m_iq1;
		m_state1 = m_ia1 * (m_state1 - k*m_iq1) - k*m_ir1;
		if ( m_state1 < 0 )
		{
			m_state1 += m_im1;
		}

		k = m_state2 / m_iq2;
		m_state2 = m_ia2 * (m_state2 - k*m_iq2) - k*m_ir2;
		if ( m_state2 < 0 )
		{
			m_state2 += m_im2;
		}

		j=m_iy/m_ndiv;

		// state1 is shuffled and state1 and state2 are combined to
		// generate the output
		m_iy = m_iv[j] - m_state2;
		m_iv[j] = m_state1;

		if ( m_iy < 1 )
		{
			m_iy += m_imm1;
		}

		if ( (temp=m_am*m_iy) > m_rnmx )
		{
			return m_rnmx;
		}
		else
		{
			return temp;
		}
	}

};

// binary generator, returns -1.0, 1.0 with probability 1/2
// has small period
// see: Numerical Recipes in C++
class RngRandomBitClass : public RanNumClass {

public:

	explicit RngRandomBitClass(unsigned long seed = DEFAULT_RNG_SEED);

	~RngRandomBitClass();

	// class is copyable, compiler generated assignment operator
	// and copy ctor work fine

	void setSeed(unsigned long seed);
	void reset();
	unsigned long getSeed() const;
	double getRandNumber();
	void getRandVector(double *vec, long dim);

private:

	unsigned long m_seed,
		          m_state,
				  m_output;
};


//declare smartPtr versions	
DECLARE(RanNumClass);

DRLIB_END_NAMESPACE

#endif
