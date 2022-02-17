/*
 * $Log: rand-gen.cpp,v $
 * Revision 1.8  2003/09/23 09:37:46  emezzine
 * Added MRGK5Random()
 *
 * Revision 1.7  2003/07/01 18:30:46  arm
 * abs replaced by labs
 *
 * Revision 1.6  2003/04/28 15:18:33  emezzine
 *  Ajout d'un autre generateur gaussien
 *
 * Revision 1.5  2001/01/24 16:43:03  sgasquet
 * Ajout exception dans la fonction GrayFaureInit
 *
 * Revision 1.4  2000/11/13 14:24:23  sgasquet
 * Cght seed par seed pour eviter conflit avec element class dans seedMT()
 *
 * Revision 1.3  2000/11/13 10:51:50  sgasquet
 * Ajout generateurs MMT et Faure
 *
 * Revision 1.2  1999/10/28 10:26:02  mab
 * comment de : include <iostream.h>
 *
 * Revision 1.1  1998/11/19 11:08:14  nicolasm
 * Initial revision
 *
 */

#include<math.h>
#include<stdlib.h>

/* TMP
#include <iostream.h>
*/


#include "rand-gen.h"
#include "gaussian.h"



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


// Algorithme nrc de generation de variables aleatoires uniformement
// distribuees entre  0 et 1
double ran2( long *idum)
{
    int j;
    long k;
    static long idum2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;

    if (*idum <= 0) 
    {
        if (-(*idum) < 1) *idum=1;
        else *idum = -(*idum);
        idum2=(*idum);
        for (j=NTAB+7;j>=0;j--) 
        {
            k=(*idum)/IQ1;
            *idum=IA1*(*idum-k*IQ1)-k*IR1;
            if (*idum < 0) *idum += IM1;
            if (j < NTAB) iv[j] = *idum;
        }
        iy=iv[0];
    }
    k=(*idum)/IQ1;
    *idum=IA1*(*idum-k*IQ1)-k*IR1;
    if (*idum < 0) *idum += IM1;
    k=idum2/IQ2;
    idum2=IA2*(idum2-k*IQ2)-k*IR2;
    if (idum2 < 0) idum2 += IM2;
    j=iy/NDIV;
    iy=iv[j]-idum2;
    iv[j] = *idum;
    if (iy < 1) iy += IMM1;
    if ((temp=AM*iy) > RNMX) return RNMX;
    else return temp;
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


double MRGK5Random(long FirstUse)
{
  static double M1  = 4294949027.0;
  static double M2  = 4294934327.0;
  static double A12 = 1154721.0;
  static double A14 = 1739991.0;
  static double A15N= 1108499.0;
  static double A21 = 1776413.0;
  static double A23 = 865203.0;
  static double A25N= 1641052.0;
  static double NORM= 2.3283163396834613e-10;
  
  static double x10, x11, x12, x13, x14, x20, x21, x22, x23, x24;
  long k;
  double p1, p2,sample;
  

  
  /* First call to the sequence */
  if (FirstUse) 
    {
      /*Initialization*/
      x10= 231458761.;
      x11= 34125679.;
      x12= 45678213.;
      x13= 7438902.;
      x14= 957345.;

      x20= 57964412.;
      x21= 12365487.;
      x22= 77221456.;
      x23= 816403.;
      x24= 8488912.;
    }
  
  /* For each call to the sequence, computation of a new point */
  /* First generator with Schrage method */
  p1= A12*x13 - A15N*x10;
  
  if(p1> 0.0) 
    p1-= A14*M1;

  p1+= A14*x11;
  k= (long)floor(p1/M1);/*TOCHECK*/
  p1-= k*M1;

  if(p1< 0.0) 
    p1+= M1;
  
  x10= x11;
  x11= x12;
  x12= x13; 
  x13= x14;
  x14= p1;
  
  /* Second generator with Schrage method */
  p2= A21*x24 - A25N*x20;

  if(p2> 0.0)
    p2-= A23*M2;

  p2+= A23*x22;
  k= (long)floor(p2/M2);/*TOCHECK*/
  p2-= k*M2;

  if(p2< 0.0)
    p2+= M2;

  x20= x21;
  x21= x22;
  x22= x23;
  x23= x24;
  x24= p2;

  /*Combination of the two generators */
  if (p1<= p2)
	  sample= (p1- p2+ M1)*NORM;
  else
	  sample=(p1- p2)*NORM;
  return sample;
}


ARM_BaseGenerator::ARM_BaseGenerator(long seed_, long FirstUse)
                   :seed(seed_)
{
  generator_param=-labs(seed);

  itsFirstUse = FirstUse;
}


double ARM_BaseGenerator::Generate()
{
    double rand;

    if(itsFirstUse)
    {
        rand = MRGK5Random(itsFirstUse);
        SetFirstUse(0); 
    }
    else
    {
        rand = MRGK5Random(itsFirstUse);           
    }

  //return ran2(&generator_param);

    return rand;
}

void ARM_BaseGenerator::Restart()
{
  generator_param=-labs(seed);
}





/********************************************************************************/
/*                                                                              */
/*               Monstruous Mersenne Twister                                    */
/*                                                                              */
/********************************************************************************/
ARM_MMTGenerator::ARM_MMTGenerator( uint32 seed_)
{
	seedMT(seed_);
	left = -1;
}


void ARM_MMTGenerator::seedMT(uint32 seed_)
 {

    register uint32 x = (seed_ | 1U) & 0xFFFFFFFFU, *s = state;
    register int    j;

    for(left=0, *s++=x, j=K_N; --j;
        *s++ = (x*=69069U) & 0xFFFFFFFFU);
 }


uint32 ARM_MMTGenerator::reloadMT(void)
 {
    register uint32 *p0=state, *p2=state+2, *pM=state+K_M, s0, s1;
    register int    j;

//    if(left < -1)
//        seedMT(4357U);

    left=K_N-1, next=state+1;

    for(s0=state[0], s1=state[1], j=K_N-K_M+1; --j; s0=s1, s1=*p2++)
        *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K_K : 0U);

    for(pM=state, j=K_M; --j; s0=s1, s1=*p2++)
        *p0++ = *pM++ ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K_K : 0U);

    s1=state[0], *p0 = *pM ^ (mixBits(s0, s1) >> 1) ^ (loBit(s1) ? K_K : 0U);
    s1 ^= (s1 >> 11);
    s1 ^= (s1 <<  7) & 0x9D2C5680U;
    s1 ^= (s1 << 15) & 0xEFC60000U;
    return(s1 ^ (s1 >> 18));
 }


uint32 ARM_MMTGenerator::LGenerate(void)
 {
    uint32 y;

    if(--left < 0)
        return(reloadMT());

    y  = *next++;
    y ^= (y >> 11);
    y ^= (y <<  7) & 0x9D2C5680U;
    y ^= (y << 15) & 0xEFC60000U;
    return(y ^ (y >> 18));
 }

double ARM_MMTGenerator::Generate(void)
{
    uint32 y;

    if(--left < 0)
        return((double)reloadMT() / (unsigned long)0xffffffff );

    y  = *next++;
    y ^= (y >> 11);
    y ^= (y <<  7) & 0x9D2C5680U;
    y ^= (y << 15) & 0xEFC60000U;
	y = y ^ (y >> 18);
    return ( (double)y / (unsigned long)0xffffffff );
}




/********************************************************************************/
/*                                                                              */
/*               Gray Faure multidimentionnal sequence                          */
/*                                                                              */
/********************************************************************************/
int ARM_GFGenerator::GrayFaureInit(int dim,int skip)
{
  if (GFadd)
    FreeMemory();
  if ((dim<2) || (skip<0))
     return 0;

  GFdimension = dim;
  GFbase = Base();
  GFmbit = (int) (log((double) LONG_MAX)/log((double) GFbase));

  if (skip>=GFmbit)
  {
        throw Exception(__LINE__, __FILE__, ERR_INITIAL_VALUE_PB,
                        "Cannot handle this too high dimension pb");
  }

  Tables();
  Pascal();
  FirstPoint(skip);
  return GFbase;
}



void ARM_GFGenerator::FreeMemory(void)
{
  int i;
  int j;
  for (i=0;i<GFbase;i++)
    {
      free(GFadd[i]);
      free(GFsub[i]);
    }
  free(GFadd);
  free(GFsub);
  for (i=0;i<GFmbit;i++)
    {
      free(GFpowbase[i]);
      for (j=0;j<GFdimension;j++)
	free(GFpascal[i][j]);
      free(GFpascal[i]);
    }
  free(GFpowbase);
  free(GFplusone);
  free(GFpascal);
  for (i=0;i<GFdimension;i++)
    free(GFfaure[i]);
  free(GFfaure);
  free(GFb_ary);
  free(GFgray);
  free(GFoutput);
}


int ARM_GFGenerator::Base(void)
{
  int base;
  int found = 0;
  
  if (GFdimension<4)
    return GFdimension;
  if (GFdimension%2==0)
    base = GFdimension+1;
  else
    base = GFdimension;
  while (found==0)
    {
      int quot = 3;
      int qmax = (int) sqrt((double) base);
      found = 1;
      while ((found==1) && (quot<=qmax))
	{
	  if (base%quot==0)
	    {
	      found = 0;
	      base += 2;
	    }
	  quot += 2;
	}
    }
  return base;
}


/* Generation of the tables GFadd, GFsub, GFpowbase and GFplusone  */
/* used by GrayFaureNext to generate quickly the next Faure point. */
/*                                                                 */
/*        GFadd[i][j] = (i+j) (mod GFbase)  i,j=0..GFbase-1        */
/*        GFsub[i][j] = GFbase-1+i-j        i,j=0..GFbase-1        */
/*      GFpowbase[i][j] = (j-GFbase+1)*(GFbase^(GFmbit-i-1))       */
/*                i=0..GFmbit-1  j=0..2*GFbase-2                   */
/*         GFplusone[i] = (i+1) (mod GFbase)  i=0..GFbase-1        */
void ARM_GFGenerator::Tables(void)
{
  int i;
  int j;

  GFadd = (int **) calloc((unsigned) GFbase,sizeof(int *));
  for (i=0;i<GFbase;i++)
    {
      GFadd[i] = (int *) calloc((unsigned) GFbase,sizeof(int));
      for (j=0;j<GFbase;j++)
	GFadd[i][j] = (i+j)%GFbase;
    }
  GFsub = (int **) calloc((unsigned) GFbase,sizeof(int *));
  for (i=0;i<GFbase;i++)
    {
      GFsub[i] = (int *) calloc((unsigned) GFbase,sizeof(int));
      for (j=0;j<GFbase;j++)
	GFsub[i][j] = GFbase-1+i-j;
    }
  GFpowbase = (long **) calloc((unsigned) GFmbit,sizeof(long *));
  for (i=0;i<GFmbit;i++)
    GFpowbase[i] = (long *) calloc((unsigned) 2*GFbase-1,sizeof(long));
  GFpowbase[GFmbit-1][GFbase] = 1;
  for (i=GFmbit-2;i>=0;i--)
    GFpowbase[i][GFbase] = GFpowbase[i+1][GFbase]*GFbase;
  for (i=0;i<GFmbit;i++)
    {
      for (j=GFbase+1;j<2*GFbase-1;j++)
	GFpowbase[i][j] = GFpowbase[i][j-1]+GFpowbase[i][GFbase];
      for (j=GFbase-1;j>=0;j--)
	GFpowbase[i][j] = GFpowbase[i][j+1]-GFpowbase[i][GFbase];	
    }
  GFplusone = (int *) calloc((unsigned) GFbase,sizeof(int));
  for (j=0;j<GFbase;j++)
    GFplusone[j] = (j+1)%GFbase;
}

/* Calculation of GFpascal[a][b][c]    a=0,..,GFmbit-1                  */
/*                                     b=0,..,GFdimension-1             */
/* 				       c=0,..,a                         */
/* GFpascal[a] is a GFdimension(columns)x(a+1)(lines) matrix. Its (b)th */
/* column is equal to the (b)th column of the Pascal matrix raised to   */
/* the (a)th power (b=0..GFdimension-1). Every operation is done modulo */
/* GFbase. Note that this matrix has only (a+1) non-zero lines because  */
/* the Pascal matrix (as its powers) is upper-triangular. Computation   */ 
/* is made by recursion and don't require any factorial calculus or     */
/* exponentiation. The Pascal matrix C contains the binomial coeffs:    */
/*                 C[i][j]=j!/(i!*(j-i)!) for 0<=i<=j                   */
/*                 C[i][j]=0              for 0<=j<i                    */
void ARM_GFGenerator::Pascal(void)
{
  int i;
  int j;
  int k;
  int fact=1;
  int diag;

  GFpascal = (int ***) calloc((unsigned) GFmbit,sizeof(int **));
  for (k=0;k<GFmbit;k++)
    {
      GFpascal[k] = (int **) calloc((unsigned) GFdimension,sizeof(int *));
      for (j=0;j<GFdimension;j++)
	  GFpascal[k][j] = (int *) calloc((unsigned) k+1,sizeof(int));
    }
  for (k=0;k<GFmbit;k++)
    {
      GFpascal[k][0][k] = 1;
      GFpascal[k][1][0] = 1;
      GFpascal[k][1][k] = 1;
    }
  for (k=2;k<GFmbit;k++)
    for (i=1;i<k;i++)
      GFpascal[k][1][i] = GFadd[GFpascal[k-1][1][i-1]][GFpascal[k-1][1][i]];
  for (j=2;j<GFdimension;j++)
    for (k=GFmbit-1;k>=0;k--)
      {
	diag = GFmbit-k-1;
	if (diag==0)
	  fact = 1;
	else
	  fact = (fact*j)%GFbase;
	for (i=0;i<=k;i++)
	  GFpascal[diag+i][j][i] = (fact*GFpascal[diag+i][1][i])%GFbase;
      }
}

/* Initialisation of GFb_ary, a vector which contains the decomposition  */
/* in base GFbase of the index of the first iteration (GFbase^skip-1).   */
/* Initialisation of GFgray, a vector which contains the decomposition   */
/* in base GFbase of the Gray code of the index coded in GFb_ary. It can */
/* be shown that this is a null vector except the (skip-1)th component   */
/* which is (GFbase-1) if skip>0.                                        */
/* Initialisation of the matrix GFfaure[i][j]. The vector GFfaure[i] is  */
/* equal to the product (modulo GFbase) of the (i)th power of the Pascal */
/* matrix by the GFgray vector. It is a representation of the output in  */
/* base GFbase.                                                          */
/* Initialisation of the vector GFoutput[GFdimension] which contains the */
/* first point of the Faure sequence. The last component of this vector  */
/* is the common denominator (GFbase^GFmbit) by which one has to divide  */
/* the numerators contained in the GFdimension first ones. So, the first */
/* point of the Faure sequence is (GFoutput[0]/GFoutput[GFdimension],..  */
/* ..,GFoutput[GFdimension-1]/GFoutput[GFdimension]).                    */
void ARM_GFGenerator::FirstPoint(int skip)
{
  int i;
  int j;

  GFfaure = (int **) calloc((unsigned) GFdimension,sizeof(int *));
  for (i=0;i<GFdimension;i++)
    GFfaure[i] = (int *) calloc((unsigned) GFmbit,sizeof(int));
  GFb_ary = (int *) calloc((unsigned) GFmbit+1,sizeof(int));
  for (i=0;i<skip;i++)
    GFb_ary[i] = GFbase-1;
  GFgray = (int *) calloc((unsigned) GFmbit,sizeof(int));
  if (skip>0)
    {
      GFgray[skip-1] = GFbase-1;
      for (i=0;i<GFdimension;i++)
	for (j=0;j<skip;j++)
	  GFfaure[i][j] = ((GFbase-1)*GFpascal[skip-1][i][j])%GFbase;
    }
  GFoutput = (long *) calloc((unsigned) GFdimension+1,sizeof(long *));
  GFoutput[GFdimension]=GFbase*GFpowbase[0][GFbase];
  for (i=0;i<GFdimension;i++)
    for (j=0;j<skip;j++)
      GFoutput[i] += GFpowbase[j][GFsub[GFfaure[i][j]][0]];
}

/* This routine must be called to generate the next point (GrayFaureInit  */
/* must be called once before the first call to GrayFaureNext).           */
/* It increments GFb_ary, determines which bit of the Gray code has to    */
/* change and increments it. Then it generates (in integer arithmetic and */
/* by table lookup) the coordinates of the next point of the sequence.    */
/* It returns a pointer on the vector GFoutput defined as in FirstPoint,  */
/* or the NULL pointer if the next index can't be represented on GFmbit   */
/* GFbase-bits (i.e. if GFbase*index > long_MAX).                         */
long *ARM_GFGenerator::VLGenerate(void)
{
  int bit = -1;
  int tmp;
  int i;
  int j;

  do
    {
      bit++;
      GFb_ary[bit] = GFplusone[GFb_ary[bit]];
    }
  while (GFb_ary[bit]==0);
  if (bit==GFmbit)
    return NULL;
  GFgray[bit] = GFplusone[GFgray[bit]];

  for (i=0;i<GFdimension;i++)
    for (j=0;j<=bit;j++)
      {
	tmp = GFfaure[i][j];
	GFfaure[i][j] = GFadd[GFpascal[bit][i][j]][tmp];
	GFoutput[i] += GFpowbase[j][GFsub[GFfaure[i][j]][tmp]];
      }
  return GFoutput;
}



double *ARM_GFGenerator::VGenerate(void)
{
	VLGenerate();

	for (int i=0;i<GFdimension;i++)
		DGFoutput[i] = GFoutput[i] / (double) GFoutput[GFdimension];

	return DGFoutput;
}

long ARM_GFGenerator::LGenerate(void)
{
		VLGenerate();

		return GFoutput[0];
}

double ARM_GFGenerator::Generate(void)
{
		VLGenerate();

		return (GFoutput[0] / (double) GFoutput[GFdimension]);
}



ARM_GFGenerator::ARM_GFGenerator(int dim, int skip)
{
	
	DGFoutput = new double[dim+1];
	GFadd = NULL;
	GrayFaureInit(dim,skip);
}

ARM_GFGenerator::~ARM_GFGenerator()
{
	delete DGFoutput;
	FreeMemory();
}







ARM_NormalGenerator::ARM_NormalGenerator()
{
  base_generator=new ARM_BaseGenerator(0);
  own_base_generator=1;
  next_available=0;

  itsGaussianType = K_BOX_MULLER;
}

ARM_NormalGenerator::ARM_NormalGenerator( ARM_BaseGenerator* base_generator_,int gaussianType)
  :base_generator(base_generator_),own_base_generator(0),next_available(0)
{
  itsGaussianType = gaussianType;
}

ARM_NormalGenerator::~ARM_NormalGenerator()
{
  if (own_base_generator)
  {
    delete base_generator;
  }
}

double ARM_NormalGenerator::Generate()
{
    if (next_available)
    {
        next_available=0;
        return next;  
    }
    else
    {
        if(!itsGaussianType)
        {
            double u=base_generator->Generate();
            while (u==0)
            {
                u=base_generator->Generate();
            }
            double v=base_generator->Generate();
            double tmp=sqrt(-2*log(u));
        
            next_available=1;
            next=(tmp*cos(2*PI*v));
        
            return (tmp*sin(2*PI*v));
        }
        else 
        {
            double u     = base_generator->Generate();
            double gauss = INV_PART_FUNC_NOR(u);

            return gauss;
        }
    }
}

ARM_GaussianGenerator::ARM_GaussianGenerator(double mu_, double sigma2_)
  :mu(mu_),sigma(sqrt(sigma2_)){
}


ARM_GaussianGenerator::ARM_GaussianGenerator( ARM_BaseGenerator* base_generator_, 
                                                          double mu_, double sigma2_)
  :mu(mu_),sigma(sqrt(sigma2_)),ng(base_generator_)
{
}

double ARM_GaussianGenerator::Generate()
{
    return mu+sigma*ng.Generate();
}


ARM_GaussianGenerator::~ARM_GaussianGenerator(){
}


