#include	"stdafx.h"


#include "MCInt\INTRng.h"


//----------------------------------------------------------------------------------------------

// Seed for ran2
long int _seed;

// Constants for ran2
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

//----------------------------------------------------------------------------------------------

// ran2 initialisation
void ran2init(long seed)
{
  _seed = seed;
  RAN2();
}

//----------------------------------------------------------------------------------------------

// ran2 random number generator
// From NR in C, p. 282
// Note: Instead of the original parameter long* idum we use a global variable _seed
double RAN2(
	void
)
{
  long *idum = &_seed;

  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];

  int j;
  long k;
  double temp;

  // Initialize
  if (*idum <= 0) 
    {
      // Be sure to prevent idum = 0
      if (-(*idum) < 1) *idum = 1; 
      else *idum = -(*idum);

      idum2 = (*idum);

      // Load the shuffle table (after 8 warm-ups)
      for (j = NTAB + 7; j >= 0; j--) 
	{ 
	  k = (*idum)/IQ1;
	  *idum = IA1*(*idum - k*IQ1) - k*IR1;
	  if (*idum < 0) *idum += IM1;
	  if (j < NTAB) iv[j] = *idum;
	}
      iy = iv[0];
    }

  // Start here when not initializing
  k = (*idum)/IQ1; 

  // Compute idum = (IA1*idum) % IM1 without overflows by Schrage's method
  *idum = IA1*(*idum - k*IQ1) - k*IR1;
  if (*idum < 0) *idum += IM1;

  k = idum2/IQ2;

  // Compute idum2 = (IA2*idum) % IM2 likewise
  idum2 = IA2*(idum2 - k*IQ2) - k*IR2; 
  if (idum2 < 0) idum2 += IM2;

  // Will e in the range 0..NTAB-1
  j = iy/NDIV; 

  // Here idum is shuffled, idum and idum2 are combined to generate output
  iy = iv[j] - idum2; 
  iv[j] = *idum;
  if (iy < 1) iy += IMM1;

  // Because users don't expect endpoint values
  if ((temp=AM*iy) > RNMX) return RNMX; 
  else return temp;
}

//----------------------------------------------------------------------------------------------

