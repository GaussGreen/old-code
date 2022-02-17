

#include "edginc/config.hpp"
#include "edginc/MCRandomLite.hpp"

DRLIB_BEGIN_NAMESPACE

// TODO :
// 1) add references
// 2) check and fix, if needed, assignment and copy constructors
//    i.e., all memeber variables should be copied over
//    (NRCLass already fixed)

// ********************* functions of RngUniformGfsr4Class - start ******************
RngUniformGfsr4Class::RngUniformGfsr4Class(unsigned long seed)
{
	// set constsnts
	m_A=471;
    m_B=1586;
    m_C=6988;
    m_D=9689;
	m_M=16383; // = 2^14
	m_invMax= 1.0/4294967296.0;

	// allocate memory for state array
	m_state = new unsigned long[m_M+1];

	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
	initializeState(m_seed);

}

RngUniformGfsr4Class::RngUniformGfsr4Class(const RngUniformGfsr4Class &otherRng)
{
	long i = -1;

	if ( this != &otherRng )
	{
		m_seed    = otherRng.m_seed;
		m_index   = otherRng.m_index;
		m_A = otherRng.m_A;
		m_B = otherRng.m_B;
		m_C = otherRng.m_C;
		m_D = otherRng.m_D;
		m_M = otherRng.m_M;
		m_invMax  = otherRng.m_invMax;

		m_state = new unsigned long[m_M+1];

		for (i=0; i<=(long)m_M; i++)
		{
			m_state[i] = otherRng.m_state[i];
		}
	}
}

RngUniformGfsr4Class &RngUniformGfsr4Class::operator=(const RngUniformGfsr4Class &otherRng)
{
	long i = -1;
	if ( this != &otherRng )
	{
		m_seed    = otherRng.m_seed;
		m_index   = otherRng.m_index;
		m_A = otherRng.m_A;
		m_B = otherRng.m_B;
		m_C = otherRng.m_C;
		m_D = otherRng.m_D;
		m_M = otherRng.m_M;
		m_invMax  = otherRng.m_invMax;

		delete[] m_state;
		m_state = new unsigned long[m_M+1];

		for (i=0; i<=(long)m_M; i++)
		{
			m_state[i] = otherRng.m_state[i];
		}
	}

	return *this;
}

RngUniformGfsr4Class::~RngUniformGfsr4Class()
{
	delete[] m_state;
}

void RngUniformGfsr4Class::setSeed(unsigned long seed)
{
	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
	initializeState(m_seed);
}

void RngUniformGfsr4Class::reset()
{
	initializeState(m_seed);
}

unsigned long RngUniformGfsr4Class::getSeed() const
{
	return m_seed;
}

double RngUniformGfsr4Class::getRandNumber()
{
	return getNumber()*m_invMax;
}

void RngUniformGfsr4Class::getRandVector(double *vec, long dim)
{
	for (long i=0; i<dim; i++)
	{
		vec[i]= getNumber()*m_invMax;
	}
}

void RngUniformGfsr4Class::initializeState(unsigned long seed)
{
	// Masks for turning on the diagonal bit and turning off the
    // leftmost bits
    unsigned long msb = 0x80000000UL;
    unsigned long mask = 0xffffffffUL;
	int i;

	if (seed == 0)
	{
		seed = DEFAULT_RNG_SEED;  
	}

    // We use the congruence s_{n+1} = (69069*s_n) mod 2^32
	// i.e. the private function Lcg(n) to initialize the state from the seed. 
	// This works because ANSI-C unsigned long
    // integer arithmetic is automatically modulo 2^32 (or a higher
    // power of two), so we can safely ignore overflow.

    // state array initialization loop using the linear congruence
	// Lcg(n) and the masks to avoid low bit correlations
    for (i = 0; i <= (int)m_M; i++)
    {
		unsigned long t = 0UL;
        unsigned long bit = msb ;

        for (int j = 0; j < 32; j++)
        {
			seed = Lcg(seed) ;
            if (seed & msb) 
			{
				t |= bit ;
			}
			bit >>= 1 ;
		}
		m_state[i] = t ;
	}

    // Perform the "orthogonalization" of the matrix
    for (i=0; i<32; ++i)
	{
		int k=7+i*3;
        m_state[k] &= mask;     // Turn off bits left of the diagonal 
        m_state[k] |= msb;      // Turn on the diagonal bit           
        mask >>= 1;
        msb >>= 1;
	}
	
	// initialize the state array index
	m_index = i;
}

// ******************** functions of RngUniformGfsr4Class - end *********************



// ******************* functions of RngGaussianGfsr4Class - start *******************
RngGaussianGfsr4Class::RngGaussianGfsr4Class(unsigned long seed)
{
	// set constsnts
	m_A=471;
    m_B=1586;
    m_C=6988;
    m_D=9689;
	m_M=16383; // = 2^14
	m_invMax= 1.0/4294967296.0;

	// allocate memory for state array
	m_state = new unsigned long[m_M+1];

	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
	initializeState(m_seed);
}

RngGaussianGfsr4Class::RngGaussianGfsr4Class(const RngGaussianGfsr4Class &otherRng)
{
	long i = -1;

	if ( this != &otherRng )
	{
		m_seed    = otherRng.m_seed;
		m_index   = otherRng.m_index;
		m_A = otherRng.m_A;
		m_B = otherRng.m_B;
		m_C = otherRng.m_C;
		m_D = otherRng.m_D;
		m_M = otherRng.m_M;
		m_invMax  = otherRng.m_invMax;
		m_draw = otherRng.m_draw;

		/*
		if ( m_state != 0 )
		{
			delete[] m_state;
		}
		*/

		m_state = new unsigned long[m_M+1];

		for (i=0; i<=(long)m_M; i++)
		{
			m_state[i] = otherRng.m_state[i];
		}
	}
}

RngGaussianGfsr4Class &RngGaussianGfsr4Class::operator=(const RngGaussianGfsr4Class &otherRng)
{
	long i = -1;

	if ( this != &otherRng )
	{
		m_seed    = otherRng.m_seed;
		m_index   = otherRng.m_index;
		m_A = otherRng.m_A;
		m_B = otherRng.m_B;
		m_C = otherRng.m_C;
		m_D = otherRng.m_D;
		m_M = otherRng.m_M;
		m_invMax  = otherRng.m_invMax;

		if ( m_state )
		{
			delete[] m_state;
		}

		m_state = new unsigned long[m_M+1];

		for (i=0; i<=(long)m_M; i++)
		{
			m_state[i] = otherRng.m_state[i];
		}
	}
	
	return *this;
}

RngGaussianGfsr4Class::~RngGaussianGfsr4Class()
{
	delete[] m_state;
}

void RngGaussianGfsr4Class::setSeed(unsigned long seed)
{
	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
	initializeState(m_seed);
}

void RngGaussianGfsr4Class::reset()
{
	initializeState(m_seed);
}

unsigned long RngGaussianGfsr4Class::getSeed() const
{
	return m_seed;
}

double RngGaussianGfsr4Class::getRandNumber()
{
	// uses Box-Muller alrorithm
	// numerical recipes in C++ implementation
	static double gset;
	double fac,radius,v1,v2;

	if ( m_draw==0 )
	{
		do
		{
			// pick to random numbers in [-1,1]^2
			v1=(2.0*getNumber()*m_invMax)-1.0;
			v2=(2.0*getNumber()*m_invMax)-1.0;

			// check they are in the unit circle
			radius= v1*v1+v2*v2;
		} while ( (radius>=1.0) || (radius < SMALL_RADIUS) );

		// perform Box-Muller transformation to get two
		// gaussian deviates. return one, save the other for
		// next call and set the flag to one
		fac= sqrt(-2.0*log(radius)/radius);
		gset=v1*fac;
		m_draw=1;

		return v2*fac;
	}
	else
	{
		// m_draw flag == 1, we don't need to pick a new point
		// in the unit circle, just return the gaussian
		// number from the previous call
		m_draw=0;
		return gset;
	}
}

void RngGaussianGfsr4Class::getRandVector(double *vec, long dim)
{
	static double gset;
	double fac,radius,v1,v2;

	for (long i=0; i<dim; i++)
	{
		// for every component of random vector use Box-Muller
		if ( m_draw==0 )
	    {
			do
		    {
				v1=(2.0*getNumber()*m_invMax)-1.0;
				v2=(2.0*getNumber()*m_invMax)-1.0;
			    radius= v1*v1+v2*v2;
			} while ( (radius>=1.0) || (radius < SMALL_RADIUS) );
		    fac= sqrt(-2.0*log(radius)/radius);

		    gset=v1*fac;
		    m_draw=1;
		    vec[i]= v2*fac;
		}
	    else
		{
			m_draw=0;
			vec[i]= gset;
		}
	}
}

void RngGaussianGfsr4Class::initializeState(unsigned long seed)
{
	// Masks for turning on the diagonal bit and turning off the
    // leftmost bits
    unsigned long msb = 0x80000000UL;
    unsigned long mask = 0xffffffffUL;
	int i;

	if (seed == 0)
	{
		seed = DEFAULT_RNG_SEED;   // the default seed is 4357
	}

    // We use the congruence s_{n+1} = (69069*s_n) mod 2^32
	// i.e. the private function Lcg(n) to initialize the state from the seed. 
	// This works because ANSI-C unsigned long
    // integer arithmetic is automatically modulo 2^32 (or a higher
    // power of two), so we can safely ignore overflow.

    // state array initialization loop using the linear congruence
	// Lcg(n) and the masks to avoid low bit correlations
    for (i = 0; i <= (int)m_M; i++)
    {
		unsigned long t = 0UL;
        unsigned long bit = msb ;

        for (int j = 0; j < 32; j++)
        {
			seed = Lcg(seed) ;
            if (seed & msb) 
			{
				t |= bit ;
			}
			bit >>= 1 ;
		}
		m_state[i] = t ;
	}

    // Perform the "orthogonalization" of the matrix
    for (i=0; i<32; ++i)
	{
		int k=7+i*3;
        m_state[k] &= mask;     // Turn off bits left of the diagonal 
        m_state[k] |= msb;      // Turn on the diagonal bit           
        mask >>= 1;
        msb >>= 1;
	}
	
	// initialize the state array index
	m_index = i;

	// initialize the flag for Box-Muller
	m_draw = 0;
}

// ************************* functions of RngGaussianGfsr4Class - end ***************



// ************************* functions of RngUniformKnuthClass - start **************
RngUniformKnuthClass::RngUniformKnuthClass(unsigned long seed)
{
	// set constants
	mBig= 1000000000;
	mSeed= 161803398;
	m_invMax= 0.000000001;

	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
	initializeState(m_seed);
}

RngUniformKnuthClass::RngUniformKnuthClass(const RngUniformKnuthClass &otherRng)
{
	long i = -1;

	if ( this != &otherRng )
	{
		m_seed = otherRng.m_seed;
		m_rndNum = otherRng.m_rndNum;
		m_index1 = otherRng.m_index1;
		m_index2 = otherRng.m_index2;
		mBig = otherRng.mBig;
		mSeed = otherRng.mSeed;
		m_invMax = otherRng.m_invMax;

		for (i=0; i<56; i++)
		{
			m_buffer[i] = otherRng.m_buffer[i];
		}
	}
}

RngUniformKnuthClass &RngUniformKnuthClass::operator=(const RngUniformKnuthClass &otherRng)
{
	long i = -1;

	if ( this != &otherRng )
	{
		m_seed = otherRng.m_seed;
		m_rndNum = otherRng.m_rndNum;
		m_index1 = otherRng.m_index1;
		m_index2 = otherRng.m_index2;
		mBig = otherRng.mBig;
		mSeed = otherRng.mSeed;
		m_invMax = otherRng.m_invMax;

		for (i=0; i<56; i++)
		{
			m_buffer[i] = otherRng.m_buffer[i];
		}
	}

	return *this;
}

RngUniformKnuthClass::~RngUniformKnuthClass()
{
}

void RngUniformKnuthClass::setSeed(unsigned long seed)
{
	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
	initializeState(m_seed);
}

void RngUniformKnuthClass::reset()
{
	initializeState(m_seed);
}

unsigned long RngUniformKnuthClass::getSeed() const
{
	return m_seed;
}

double RngUniformKnuthClass::getRandNumber()
{
	return getNumber()*m_invMax;
}

void RngUniformKnuthClass::getRandVector(double *vec, long dim)
{
	for (long i=0; i<dim; i++)
		vec[i]= getNumber()*m_invMax;
}

void RngUniformKnuthClass::initializeState(long seed)
{
	long j= labs((mSeed-labs(seed))) % mBig;

	// initialize first and last entries of array
	m_buffer[0]=0;
	m_buffer[55]=j;

	// initialize the rest in a slightly random order,
	// with number that are not especially random
	long k=1;
	for (int i=1; i<55; i++)
	{
		int n = (21*i) % 55;

		m_buffer[n] = k;
		k = j-k;
		if ( k<0 )
		{
			k += mBig;
		}

		j= m_buffer[n];
	}

	// we "randomize" the array
	for (int i1=0; i1<4; i1++)
	{
		for (int i=1; i<56; i++)
		{
			long n = m_buffer[i] - m_buffer[1+(i+30)%55];
			if ( n<0 )
			{
				n += mBig;
			}
			m_buffer[i] = n;
		}
	}

    // initialize array indicies for the first generated number.
	// value 31 is special and should not be changed
	m_index1 = 0;
	m_index2 = 31;
}

// ************************** functions of RngUniformKnuthClass - end ***************



// ************************* functions of RngGaussianKnuthClass - start *************
RngGaussianKnuthClass::RngGaussianKnuthClass(unsigned long seed)
{
	// set constants
	mBig= 1000000000;
	mSeed= 161803398;
	m_invMax= 0.000000001;

	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
	initializeState(m_seed);
}

RngGaussianKnuthClass::RngGaussianKnuthClass(const RngGaussianKnuthClass &otherRng)
{
	long i = -1;

	if ( this != &otherRng )
	{
		m_seed = otherRng.m_seed;
		m_rndNum = otherRng.m_rndNum;
		m_index1 = otherRng.m_index1;
		m_index2 = otherRng.m_index2;
		mBig = otherRng.mBig;
		mSeed = otherRng.mSeed;
		m_invMax = otherRng.m_invMax;
		m_draw = otherRng.m_draw;

		for (i=0; i<56; i++)
		{
			m_buffer[i] = otherRng.m_buffer[i];
		}
	}
}

RngGaussianKnuthClass &RngGaussianKnuthClass::operator=(const RngGaussianKnuthClass &otherRng)
{
	long i = -1;

	if ( this != &otherRng )
	{
		m_seed = otherRng.m_seed;
		m_rndNum = otherRng.m_rndNum;
		m_index1 = otherRng.m_index1;
		m_index2 = otherRng.m_index2;
		mBig = otherRng.mBig;
		mSeed = otherRng.mSeed;
		m_invMax = otherRng.m_invMax;
		m_draw = otherRng.m_draw;

		for (i=0; i<56; i++)
		{
			m_buffer[i] = otherRng.m_buffer[i];
		}
	}

	return *this;
}

RngGaussianKnuthClass::~RngGaussianKnuthClass()
{
}

void RngGaussianKnuthClass::setSeed(unsigned long seed)
{
	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
	initializeState(m_seed);
}

void RngGaussianKnuthClass::reset()
{
	initializeState(m_seed);
}

unsigned long RngGaussianKnuthClass::getSeed() const
{
	return m_seed;
}

double RngGaussianKnuthClass::getRandNumber()
{
	// uses Box-Muller algorithm
	static double gset;
	double fac,radius,v1,v2;

	if ( m_draw==0 )
	{
		do
		{
			// pick to random numbers in [-1,1]^2
			v1=(2.0*getNumber()*m_invMax)-1.0;
			v2=(2.0*getNumber()*m_invMax)-1.0;

			// check they are in the unit circle
			radius= v1*v1+v2*v2;
		} while ( (radius>=1.0) || (radius < SMALL_RADIUS) );

		// perform Box-Muller transformation to get two
		// gaussian deviates. return one, save the other for
		// next call and set the flag to one
		fac= sqrt(-2.0*log(radius)/radius);
		gset=v1*fac;
		m_draw=1;

		return v2*fac;
	}
	else
	{
		// draw flag == 1, we don't need to pick a new point
		// in the unit circle, just return the gaussian
		// number from the previous call
		m_draw=0;
		return gset;
	}
}

void RngGaussianKnuthClass::getRandVector(double *vec, long dim)
{
	static double gset;
	double fac,radius,v1,v2;

	for (long i=0; i<dim; i++)
	{
		// do Box-Muller for each component of the random vector
		if ( m_draw==0 )
	    {
			do
		    {
				v1=(2.0*getNumber()*m_invMax)-1.0;
				v2=(2.0*getNumber()*m_invMax)-1.0;
			    radius= v1*v1+v2*v2;
			} while ( (radius>=1.0) || (radius < SMALL_RADIUS) );
		    fac= sqrt(-2.0*log(radius)/radius);

		    gset=v1*fac;
		    m_draw=1;
		    vec[i]= v2*fac;
		}
	    else
		{
			m_draw=0;
			vec[i]= gset;
		}
	}
}

void RngGaussianKnuthClass::initializeState(long seed)
{
	long j= labs((mSeed-labs(seed))) % mBig;

	// initialize first and last entries of array
	m_buffer[0]=0;
	m_buffer[55]=j;

	// initialize the rest in a slightly random order,
	// with number that are not especially random
	long k=1;
	for (int i=1; i<55; i++)
	{
		int n = (21*i) % 55;

		m_buffer[n] = k;
		k = j-k;
		if ( k<0 )
		{
			k += mBig;
		}

		j= m_buffer[n];
	}

	// we "randomize" the array
	for (int i1=0; i1<4; i1++)
	{
		for (int i=1; i<56; i++)
		{
			long n = m_buffer[i] - m_buffer[1+(i+30)%55];
			if ( n<0 )
			{
				n += mBig;
			}
			m_buffer[i] = n;
		}
	}

    // initialize array indicies for the first generated number.
	// value 31 is special and should not be changed
	m_index1 = 0;
	m_index2 = 31;

	// initialize the flag for Box-Muller
	m_draw = 0;
}

// ********************* functions of RngGaussianKnuthClass - end *******************



// ********************* functions of RngGaussianZigguratClass - start **************
RngGaussianZigguratClass::RngGaussianZigguratClass(unsigned long seed) 
{
	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
	m_rngULongState = m_seed;

	// set up arrays
	const double rectAreaGaussian   = 0.00492867323399;
	const double m= 2147483648.0; // =2^31
	double rightMostXgaussian = 3.6541528853610088; // =x_255
	double q=rectAreaGaussian*exp(0.5*rightMostXgaussian*rightMostXgaussian);

	m_int32bitGaussian[0]= static_cast<unsigned long> ((rightMostXgaussian/q)*m);
	m_int32bitGaussian[1]=0;
	
	m_wGaussian[0]=q/m;
	m_wGaussian[255]=rightMostXgaussian/m;

	m_gaussianDen[0]=1.0;
	m_gaussianDen[255]=exp(-0.5*rightMostXgaussian*rightMostXgaussian);

	double previousX=rightMostXgaussian;
	// fill up arrays starting from right most point and
	// moving towards x_0=0
	for (int i=254; i>=1; i--)
	{
		// find x_i from x_{i+1} and make this current rightMost point
		rightMostXgaussian=sqrt(-2.0*log( (rectAreaGaussian/rightMostXgaussian)+
			                             exp(-0.5*rightMostXgaussian*rightMostXgaussian)));
		m_int32bitGaussian[i+1]=static_cast<unsigned long> ( (rightMostXgaussian/previousX)*m);
		previousX=rightMostXgaussian;
		m_gaussianDen[i]=exp(-0.5*rightMostXgaussian*rightMostXgaussian);
		m_wGaussian[i]=rightMostXgaussian/m;
	}

}

RngGaussianZigguratClass::RngGaussianZigguratClass(const RngGaussianZigguratClass &otherRng)
{
	if ( this != &otherRng )
	{
		m_seed = otherRng.m_seed;
		m_rngULongState = otherRng.m_rngULongState;

		// set up arrays
		const double rectAreaGaussian   = 0.00492867323399;
		const double m= 2147483648.0; // =2^31
		double rightMostXgaussian = 3.6541528853610088; // =x_255
		double q=rectAreaGaussian*exp(0.5*rightMostXgaussian*rightMostXgaussian);

		m_int32bitGaussian[0]= static_cast<unsigned long> ((rightMostXgaussian/q)*m);
		m_int32bitGaussian[1]=0;

		m_wGaussian[0]=q/m;
		m_wGaussian[255]=rightMostXgaussian/m;

		m_gaussianDen[0]=1.0;
		m_gaussianDen[255]=exp(-0.5*rightMostXgaussian*rightMostXgaussian);

		double previousX=rightMostXgaussian;
		// fill up arrays starting from right most point and
		// moving towards x_0=0
		for (int i=254; i>=1; i--)
		{
			// find x_i from x_{i+1} and make this current rightMost point
			rightMostXgaussian=sqrt(-2.0*log( (rectAreaGaussian/rightMostXgaussian)+
				exp(-0.5*rightMostXgaussian*rightMostXgaussian)));
			m_int32bitGaussian[i+1]=static_cast<unsigned long> ( (rightMostXgaussian/previousX)*m);
			previousX=rightMostXgaussian;
			m_gaussianDen[i]=exp(-0.5*rightMostXgaussian*rightMostXgaussian);
			m_wGaussian[i]=rightMostXgaussian/m;
		}
	}
}

RngGaussianZigguratClass &RngGaussianZigguratClass::operator=(const RngGaussianZigguratClass &otherRng)
{
	if ( this != &otherRng )
	{
		m_seed = otherRng.m_seed;
		m_rngULongState = otherRng.m_rngULongState;

		// set up arrays
		const double rectAreaGaussian   = 0.00492867323399;
		const double m= 2147483648.0; // =2^31
		double rightMostXgaussian = 3.6541528853610088; // =x_255
		double q=rectAreaGaussian*exp(0.5*rightMostXgaussian*rightMostXgaussian);

		m_int32bitGaussian[0]= static_cast<unsigned long> ((rightMostXgaussian/q)*m);
		m_int32bitGaussian[1]=0;

		m_wGaussian[0]=q/m;
		m_wGaussian[255]=rightMostXgaussian/m;

		m_gaussianDen[0]=1.0;
		m_gaussianDen[255]=exp(-0.5*rightMostXgaussian*rightMostXgaussian);

		double previousX=rightMostXgaussian;
		// fill up arrays starting from right most point and
		// moving towards x_0=0
		for (int i=254; i>=1; i--)
		{
			// find x_i from x_{i+1} and make this current rightMost point
			rightMostXgaussian=sqrt(-2.0*log( (rectAreaGaussian/rightMostXgaussian)+
				exp(-0.5*rightMostXgaussian*rightMostXgaussian)));
			m_int32bitGaussian[i+1]=static_cast<unsigned long> ( (rightMostXgaussian/previousX)*m);
			previousX=rightMostXgaussian;
			m_gaussianDen[i]=exp(-0.5*rightMostXgaussian*rightMostXgaussian);
			m_wGaussian[i]=rightMostXgaussian/m;
		}
	}

	return *this;
}

RngGaussianZigguratClass::~RngGaussianZigguratClass()
{
}

void RngGaussianZigguratClass::setSeed(unsigned long seed)
{
	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
	m_rngULongState=m_seed;
}

void RngGaussianZigguratClass::reset()
{
	m_rngULongState=m_seed;
}

unsigned long RngGaussianZigguratClass::getSeed() const
{
	return m_seed;
}

double RngGaussianZigguratClass::getRandNumber()
{
	long j;
	int i;
	static double x,y;
	const double rightMost=3.6541528853610088f;
	const double invRightMost=0.273661237329759f;

	// generate a random 32 bit integer j
	j=getULong();
	// take the rightmost 8 bits to get a random index i
	// ( STEP 1 of the algorithm )
	i= (j&255);

	// ( STEP 2 of the algorithm )
	if ( labs(j)<(long)m_int32bitGaussian[i] )
		return j*m_wGaussian[i];
	else
	{
		for (;;)
		{
			x=j*m_wGaussian[i];
			// if rectangle index is zero return from tail
			// ( STEP 3 of algorithm )
			if ( i==0 )
			{
				do
				{
					x= -log(getUniform())*invRightMost;
					y= -log(getUniform());
				 } while ( y+y < x*x );
				 return ( j>0 ) ? rightMost+x : -rightMost-x;
			}
			// ( STEP 4 of the algorithm )
			if ( (m_gaussianDen[i]+getUniform()*
				 (m_gaussianDen[i-1]-m_gaussianDen[i])) < exp(-0.5*x*x) )
				return x;
			// Try again if not successfull yet
			// ( STEP 5 of the algorithm)
			j=getULong();
			i=(j&255);
			if ( labs(j)<(long)m_int32bitGaussian[i] )
				return j*m_wGaussian[i];
		}
	}
}

void RngGaussianZigguratClass::getRandVector(double *vec, long dim)
/*{
	for (long i=0; i<dim; i++)
		vec[i]=GetRandNumber();
}*/
{
	long k,j;
	int i;
	static double x,y;
	const double rightMost=3.6541528853610088f;
	const double invRightMost=0.273661237329759f;

	for (k=0; k<dim; k++)
	{
		// generate a random 32 bit integer j
	    j=getULong();
	    // take the rightmost 8 bits to get a random index i
	    // ( STEP 1 of the algorithm )
	    i= (j&255);

	    // ( STEP 2 of the algorithm )
	    if ( labs(j)<(long)m_int32bitGaussian[i] )
		   vec[k] = j*m_wGaussian[i];
	    else
	    {
			for (;;)
		    {
				x=j*m_wGaussian[i];
			    // if rectangle index is zero return from tail
			    // ( STEP 3 of algorithm )
			    if ( i==0 )
			    {
					do
				    {
						x= -log(getUniform())*invRightMost;
					    y= -log(getUniform());
					} while ( y+y < x*x );
				    vec[k] = ( j>0 ) ? rightMost+x : -rightMost-x;
					break;
				}
			    // ( STEP 4 of the algorithm )
			    if ( (m_gaussianDen[i]+getUniform()*
				     (m_gaussianDen[i-1]-m_gaussianDen[i])) < exp(-0.5*x*x) )
				{
					 vec[k]= x;
					 break;
				}
			    // Try again if not successfull yet
			    // ( STEP 5 of the algorithm)
			    j=getULong();
			    i=(j&255);
			    if ( labs(j)<(long)m_int32bitGaussian[i] )
				{
					vec[k]= j*m_wGaussian[i];
					break;
				}
			} // end forever for loop
		}
	} // for k=0 to dim loop
}
// ******************* functions of RngGaussianZigguratClass - end ******************



// ****************** functions of RngGaussianNRClass - start ***********************
RngGaussianNRClass::RngGaussianNRClass(unsigned long seed)
{
	// set algorithm's constants
	m_im1 = 2147483563;
	m_im2 = 2147483399;
	m_ia1 = 40014;
	m_ia2 = 40692;
	m_iq1 = 53668;
	m_iq2 = 52774;
	m_ir1 = 12211;
	m_ir2 = 3791;
	m_imm1 = m_im1-1;
	m_ndiv = 1+ m_imm1/32;
	m_eps = 3.0e-16;
	m_rnmx = 1.0-m_eps;
	m_am = 1.0/m_im1;

	if ( seed != 0)
	{
		m_seed = seed;
		initializeState(m_seed);
	}
	else
	{
		m_seed = DEFAULT_RNG_SEED;
		initializeState(m_seed);
	}
	
}

RngGaussianNRClass::RngGaussianNRClass(const RngGaussianNRClass &otherRng)
{
	if ( this != &otherRng )
	{
		// set algorithm's constants
		m_im1  = otherRng.m_im1;
		m_im2  = otherRng.m_im2;
		m_ia1  = otherRng.m_ia1;
		m_ia2  = otherRng.m_ia2;
		m_iq1  = otherRng.m_iq1;
		m_iq2  = otherRng.m_iq2;
		m_ir1  = otherRng.m_ir1;
		m_ir2  = otherRng.m_ir2;
		m_imm1 = otherRng.m_imm1;
		m_ndiv = otherRng.m_ndiv;
		m_iy   = otherRng.m_iy;
		m_eps  = otherRng.m_eps;
		m_rnmx = otherRng.m_rnmx;
		m_am   = otherRng.m_am;
		m_seed = otherRng.m_seed;
		m_state1 = otherRng.m_state1;
		m_state2 = otherRng.m_state2;
		m_draw = otherRng.m_draw;
 
	    for(long i = 0; i < 32; i++)
		{
	       m_iv[i] = otherRng.m_iv[i];
		}
	}
}

RngGaussianNRClass &RngGaussianNRClass::operator=(const RngGaussianNRClass &otherRng)
{
	if ( this != &otherRng )
	{
		// set algorithm's constants
		// set algorithm's constants
		m_im1  = otherRng.m_im1;
		m_im2  = otherRng.m_im2;
		m_ia1  = otherRng.m_ia1;
		m_ia2  = otherRng.m_ia2;
		m_iq1  = otherRng.m_iq1;
		m_iq2  = otherRng.m_iq2;
		m_ir1  = otherRng.m_ir1;
		m_ir2  = otherRng.m_ir2;
		m_imm1 = otherRng.m_imm1;
		m_ndiv = otherRng.m_ndiv;
		m_iy   = otherRng.m_iy;
		m_eps  = otherRng.m_eps;
		m_rnmx = otherRng.m_rnmx;
		m_am   = otherRng.m_am;
		m_seed = otherRng.m_seed;
		m_state1 = otherRng.m_state1;
		m_state2 = otherRng.m_state2;
		m_draw = otherRng.m_draw;
 
	    for(long i = 0; i < 32; i++)
		{
	       m_iv[i] = otherRng.m_iv[i];
		}
	}

	return *this;
}

RngGaussianNRClass::~RngGaussianNRClass()
{
}

void RngGaussianNRClass::setSeed(unsigned long seed)
{
	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
	initializeState(m_seed);
}

void RngGaussianNRClass::reset()
{
	initializeState(m_seed);
}

unsigned long RngGaussianNRClass::getSeed() const
{
	return m_seed;
}

double RngGaussianNRClass::getRandNumber()
{
	// uses Box-Muller algorithm
	static double gset;
	double fac,radius,v1,v2;

	if ( m_draw==0 )
	{
		do
		{
			v1=(2.0*getNumber())-1.0;
			v2=(2.0*getNumber())-1.0;

			radius= v1*v1+v2*v2;
		} while ( (radius>=1.0) || (radius < SMALL_RADIUS) );

		fac= sqrt(-2.0*log(radius)/radius);
		gset=v1*fac;
		m_draw=1;

		return v2*fac;
	}
	else
	{
		m_draw=0;
		return gset;
	}
}

void RngGaussianNRClass::getRandVector(double *vec, long dim)
{
	static double gset;
	double fac,radius,v1,v2;

	for (long i=0; i<dim; i++)
	{
		// do Box-Muller for each component
		if ( m_draw==0 )
	    {
			do
		    {
				v1=(2.0*getNumber())-1.0;
				v2=(2.0*getNumber())-1.0;
			    radius= v1*v1+v2*v2;
			} while ( (radius>=1.0) || (radius < SMALL_RADIUS) );
		    fac= sqrt(-2.0*log(radius)/radius);

		    gset=v1*fac;
		    m_draw=1;
		    vec[i]= v2*fac;
		}
	    else
		{
			m_draw=0;
			vec[i]= gset;
		}
	}
}

void RngGaussianNRClass::initializeState(long seed)
{
	long k;

	m_state1 = seed;
	m_state2 = m_state1;

	// load the shuffle table after 8 warm ups
	// using the first linear congruential generator
	for (int j=32+7; j>=0; j--)
	{
		k = m_state1/m_iq1;
		m_state1 = m_ia1*(m_state1-k*m_iq1)-k*m_ir1;
		if ( m_state1 < 0 )
		{
			m_state1 += m_im1;
		}
		if ( j < 32 )
		{
			m_iv[j]= m_state1;
		}
	}
	
	m_iy = m_iv[0];

	// initialize the flag for Box-Muller
	m_draw = 0;
}
// ******************** functions of RngGaussianNRClass - end ***********************



// ******************** functions of RngRandomBitClass - start **********************
RngRandomBitClass::RngRandomBitClass(unsigned long seed)
{
	if ( seed == 0 )
	{
		m_seed = DEFAULT_RNG_SEED;
	}
	else
	{
		m_seed = seed;
	}

	m_state = m_seed;

}

RngRandomBitClass::~RngRandomBitClass()
{
}

void RngRandomBitClass::setSeed(unsigned long seed)
{
	m_seed = ( seed == 0 ) ? DEFAULT_RNG_SEED : seed ;
}

void RngRandomBitClass::reset()
{
	m_state = m_seed;
}

unsigned long RngRandomBitClass::getSeed(void) const
{
	return m_seed;
}

double RngRandomBitClass::getRandNumber()
{
	m_output = ( (m_state >> 17)& 1 )^
		       ( (m_state >> 4) & 1 )^
			   ( (m_state >> 1) & 1)^
			   ( m_state & 1);
	m_state  = (m_state << 1) | m_output;

	if ( m_output )
	{
		return 1.0;
	}
	else
	{
		return -1.0;
	}
}

void RngRandomBitClass::getRandVector(double *vec, long dim)
{
	for (long i=0; i<dim; i++)
	{
		m_output = ( (m_state >> 17) & 1 )^
			       ( (m_state >> 4)  & 1 )^
				   ( (m_state >> 1)  & 1 )^
				   ( m_state & 1);
	    m_state  = (m_state << 1) | m_output;

	    if ( m_output )
	    {
		    vec[i] = 1.0;
		}
	    else
		{
			vec[i] = -1.0;
		}
	}
}

// ********************** functions of RngRandomBitClass - end **********************


DRLIB_END_NAMESPACE
