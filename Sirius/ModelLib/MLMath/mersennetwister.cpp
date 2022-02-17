
#include "stdafx.h"
#include "mersennetwister.h"
#include "mleqnormal.h"
#include "MlEqMaths.h"

// Period parameters 
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   // constant vector a 
#define UPPER_MASK 0x80000000UL // most significant w-r bits 
#define LOWER_MASK 0x7fffffffUL // least significant r bits 

static unsigned long mt[N]; /* the array for the state vector  */
static int mti = N+1; /* mti==N+1 means mt[N] is not initialized */
/*
const long nRandNumbers = 20;//0000000;	// 20 million numbers so far...be careful 
bool mt_initialized = false;
unsigned long mt_nbCalls = 0;
double *pRandomNumber = CMersenneTwister::setup_random_numbers();
double *random_numbers[nRandNumbers];


double *CMersenneTwister::setup_random_numbers()
{
	if( !mt_initialized )
	{
//		CMersenneTwister mt(0);
//		mt.generate_numbers();
		mt_initialized = true;
	}

	return random_numbers;
}
*/
CMersenneTwister::CMersenneTwister( int dim )
{
	initialize( dim );

	
	m_nSteps = m_nFactors = m_nAssets =	0;
	m_correl.resize(0);
}

void CMersenneTwister::initialize(int dimension)
{
	m_numberDimensions = dimension;

	init_genrand(19650218UL);
//	unsigned long init[4] = {0x123, 0x234, 0x345, 0x456};	
//	init_by_array(init, 4);

//	pRandomNumber = random_numbers;
}

void CMersenneTwister::initialize(int dimension, int numberScenariosToBeStored,int numberRandomsPerScenario)
{
	initialize( dimension );

	long dummy = -123;
	randomGenerator::initialize(dummy,numberScenariosToBeStored,numberRandomsPerScenario);
}

/*
void CMersenneTwister::generate_numbers()
{	
	for(int k=0; k<nRandNumbers; k++)
	{
		random_numbers[k] = normal_inv2( genrand_real2() );	
	}
}
*/

/*
inline 
void CMersenneTwister::generateRandoms( CVector& nRandVect, int idim )
{
	for(int i=0; i<m_numberDimensions; ++i){
		nRandVect[i] = *pRandomNumber++ ;
		++mt_nbCalls;
	}
	
}
*/

inline 
void CMersenneTwister::generateRandoms( CVector& nRandVect, int idim )
{
	for(int i=0; i<m_numberDimensions; ++i){
		nRandVect[i] = normal_inv2( genrand_real2() );//*pRandomNumber++ ;
//		++mt_nbCalls;
	}
	
}
/*
double* CMersenneTwister::next()
{
	double* g = pRandomNumber;
	
	pRandomNumber	+= m_numberDimensions;
	mt_nbCalls		+= m_numberDimensions;
	return g;
}
*/

inline 
void CMersenneTwister::generateRandoms( double* nptr )
{
	for(int i=0; i<m_numberDimensions; ++i)
		*nptr++ = normal_inv2( genrand_real2() );
}

inline 
void CMersenneTwister::generateUniformRandoms( CVector& uRandVect, int idim )
{
	for(int i=0; i<m_numberDimensions; ++i)
		uRandVect[i] =  genrand_real2() ;
}



inline 
double CMersenneTwister::generateRandom()
{
	return normal_inv2( genrand_real2() );
}

inline 
double CMersenneTwister::generateUniformRandom()
{
	return  genrand_real2() ;
}


void CMersenneTwister::initCorrelationParameters(GVector < CVector >& randoms,GVector < int > numberOfFactors, Cholesky & cholesky,int ipath,CVector& rndTemp,bool generateGauss)
{
	m_nAssets	=	randoms.getsize();
	m_nFactors	=	numberOfFactors[0] ;
	m_nSteps	=	randoms[0].getsize()/m_nFactors;

//		check for consistency of matrix set up once :

	if ( numberOfFactors.getsize() != m_nAssets ){
		throw("incorrect dimensioning of arrays encountered");
	}

	for (int iasset = 0 ; iasset < m_nAssets; iasset++ )
	{
		int n = randoms[iasset].getsize();
		int m = numberOfFactors[iasset];

		if ( m != m_nFactors || n%m != 0 || m_nSteps != n / m ){
			throw("incorrect dimensioning encountered in generateRandoms");
		}
	}


//	store correlation structure efficiently :

	std::vector<double> correl;		

	int k,l;
	for (int istep = 0 ; istep < m_nSteps; istep++ )
	{
		k = -1;
		for (int iasset = 0 ; iasset < m_nAssets; iasset++ )
		{
			for (int ifactor = 0; ifactor < m_nFactors; ifactor++ )
			{	
				k++;
				
				for ( l = 0; l <= k; l++ )
				{
					double corr = cholesky.getCholesky(istep,k,l) ;
					correl.push_back( corr );
				}					
			}		
		}			
	}

	int dim = correl.size();
	m_correl.resize(dim);

	for(int i=0; i<dim; ++i){		// do this out of laziness...
		m_correl[i] = correl[i];
	}
}



void CMersenneTwister::generateRandoms(GVector < CVector >& randoms,GVector < int > numberOfFactors, Cholesky & cholesky,int ipath,CVector& rndTemp,bool generateGauss)
{

//  randoms[iasset][istep*numberOfFactors[iasset]+jfactor]

	int istep,iasset,ifactor,n,m;

	if ( ipath == 0 )	{
		initCorrelationParameters( randoms, numberOfFactors, cholesky,ipath, rndTemp, generateGauss );
	}

	int ndim = rndTemp.getsize();
	if ( ndim == 0 )
	{		
		ndim = m_nSteps * m_nAssets * m_nFactors ;
		rndTemp.resize(ndim);
	}



	if ( generateGauss ){
		generateRandoms(rndTemp,ipath);
	}
	else{
		generateUniformRandoms(rndTemp,ipath);
	}



	const double* corr = m_correl.getPtr();
	const double* rand = rndTemp.getPtr();

	int factor = m_nAssets * m_nFactors ;
	int k,l;

	for ( istep = 0 ; istep < m_nSteps; istep++ )
	{
		m = istep*factor; 
		k = -1;
		for ( iasset = 0 ; iasset < m_nAssets; iasset++ )
		{
	//		double* corrRand = randoms[iasset].getPtr();

			for (ifactor = 0; ifactor < m_nFactors; ifactor++ )
			{	
				k++;
				n = istep*m_nFactors + ifactor;
			
				double correlatedRandomNumber = 0.0;
				for ( l = 0; l <= k; l++ )
				{	
					correlatedRandomNumber += (*corr++) * ( *(rand+l+m) ); // lazy again...
				}

				randoms[iasset][n] = correlatedRandomNumber ;
		//		*(corrRand + n) = correlatedRandomNumber ;		
			}		
		}			
	}


	
	
}



// original code

/* initializes mt[N] with a seed */
inline 
void CMersenneTwister::init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
inline 
void CMersenneTwister::init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
inline
unsigned long CMersenneTwister::genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
inline
long CMersenneTwister::genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
inline
double CMersenneTwister::genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
inline
double CMersenneTwister::genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
inline
double CMersenneTwister::genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
inline
double CMersenneTwister::genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */


#undef N
#undef M
#undef MATRIX_A
#undef UPPER_MASK
#undef LOWER_MASK