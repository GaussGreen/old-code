/************************************************************************
 * Module:	DRL - RAND
 * Function:	Random Number Generation and Simulation
 *		Sobol Sequence
 * Author:	NRC2 - modified by C. Daher
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include <math.h>
#include <float.h>

#include "drlrand.h"			/* prototype consistency */

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ?\
        (iminarg1) : (iminarg2))

#define MAXBIT 30
#define MAXDIM 6
static unsigned long	inS,
			ixS[MAXDIM+1],
			*iuS[MAXBIT+1];
static unsigned long	mdegS[MAXDIM+1]={0,1,2,3,3,4,4};
static unsigned long	ipS[MAXDIM+1]={0,0,1,1,2,1,4};
static unsigned long	ivS[MAXDIM*MAXBIT+1]={ 0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};



/*f---------------------------------------------------------------------
 * Generates pseudo random Sobol sequence numbers.
 * When $n < 0$, internally initialize a set of direction numbers.
 * When $0<n<\mbox{MAXDIM}=30$, returns as vector <i> x[0..n-1]</i> the
 * next values form $n$ of these sequences
 * ($n$ must not be changed between initializations).
 * Returns 0 iff OK.
 */

DLL_EXPORT(int)
DrlSobolSequence(int *n, double *x)
{
static	char	routine[] = "DrlSobolSequence";
	int		j,k,l;
	unsigned long	i,im,ipp;
static double		fac;
static unsigned long	in,ix[MAXDIM+1],*iu[MAXBIT+1];
static unsigned long	mdeg[MAXDIM+1]={0,1,2,3,3,4,4};
static unsigned long	ip[MAXDIM+1]={0,0,1,1,2,1,4};
static unsigned long	iv[MAXDIM*MAXBIT+1]={ 0,1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};

	/* Because NRC vectors are [1..n] */
	x--;

#ifdef	__DEBUG__
	printf("Sobol: n=%d");
#endif

	if (*n < 0) {
		/*
		 * Reset all static vectors
		 */
#undef	COPYVECT
#define	COPYVECT(x,y,n)	{for(i=0; i<=(n)-1; i++) *((x)+i) = *((y)+i);}
		in = inS;
		COPYVECT(ix, ixS, MAXDIM+1);
		COPYVECT(iu, iuS, MAXBIT+1);
		COPYVECT(mdeg, mdegS, MAXDIM+1);
		COPYVECT(ip, ipS, MAXDIM+1);
		COPYVECT(iv, ivS, MAXDIM*MAXBIT+1);
#undef	COPYVECT


		/*
		 * Initialize
		 */
		for (j=1,k=0;j<=MAXBIT;j++,k+=MAXDIM) {
			iu[j] = &iv[k];
		}
		for (k=1;k<=MAXDIM;k++) {
			for (j=1;j<=(int)mdeg[k];j++) {
				iu[j][k] <<= (MAXBIT-j);
			}
			for (j=mdeg[k]+1;j<=MAXBIT;j++) {
				ipp=ip[k];
				i=iu[j-mdeg[k]][k];
				i ^= (i >> mdeg[k]);
				for (l=mdeg[k]-1;l>=1;l--) {
					if (ipp & 1) i ^= iu[j-l][k];
					ipp >>= 1;
				}
				iu[j][k]=i;
			}
		}
		fac=1.0/(1L << MAXBIT);
		in=0;
	} else {

		if (*n > MAXDIM) {
			GtoErrMsg("%s: maximum dimension %d (got %d).\n",
				routine, MAXDIM, *n);
			return(FAILURE);
		}


		im=in;
		for (j=1;j<=MAXBIT;j++) {
			if (!(im & 1)) break;
			im >>= 1;
		}
		if (j > MAXBIT) {
			GtoErrMsg("%s: MAXBIT too small in sobseq\n",
				routine);
			return(1);
		}
		im=(j-1)*MAXDIM;
		for (k=1;k<=IMIN(*n,MAXDIM);k++) {
			ix[k] ^= iv[im+k];
			x[k]=ix[k]*fac;
		}
		in++;
	}
	return(0);
}
#undef MAXBIT
#undef MAXDIM
#undef NRANSI


