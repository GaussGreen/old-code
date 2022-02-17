
#include "math.h"
#include "srt_h_all.h"
#include "LGMSVPDE.h"
#include "LGMSVMC.h"
#include "LGMSVGrfn.h"
#include "LGMSVUtil.h"
#include "Fx3FCalib.h"
#include "Fx3FUtils.h"
#include "CPDCalib.h"
#include "num_h_bessel.h"
#include "uterror.h"
#include "num_h_random.h"
#include "num_h_bessel.h"
#include "num_h_simpso.h"
#include "num_h_proba.h"
#include "num_h_allhdr.h"
#include "nag.h"
#include "nags.h"
#include "nag_types.h"
#include "num_h_f2c.h"



#define MAX_CPN			600

double gamma_(double *a);

double PolyHermite0( double x)
{
double result=1;
return result;
}


double PolyHermite1( double x)
{
double result;
result =x;
return result;
}


double PolyHermite2( double x)
{
double result= x*x - 1;
return result;
}


double PolyHermite3( double x)
{
double result =x*x*x - 3*x;
return result;
}


double PolyHermite4( double x)
{
double result = x*x*x*x - 6*x*x + 3;
return result;
}


double PolyHermite5( double x)
{
double result = x*x*x*x*x - 10 * x*x*x + 15*x;
return result;
}


double h1( double x)
{
double result=PolyHermite2(x)/6;
return result;
}


double h2( double x)
{
double result=PolyHermite3(x)/24;
return result;
}


double h11( double x)
{
double result=- (2 * PolyHermite3(x) + PolyHermite1(x) )/36;
return result;
}


double h3( double x)
{
double result=PolyHermite4(x)/120;
return result;
}


double h12( double x)
{
double result=-(PolyHermite4(x) + PolyHermite2(x))/24;
return result;
}


double h111( double x)
{
double result=(12 * PolyHermite4(x) + 19 * PolyHermite2(x))/324;
return result;
}


double h4( double x)
{
double result=PolyHermite5(x)/720;
return result;
}

double h22( double x)
{
double result=-( 3 * PolyHermite5(x) + 6 * PolyHermite3(x) + 2 * PolyHermite1(x))/384;
return result;
}

double h13( double x)
{
double result=- (2 * PolyHermite5(x) + 3 * PolyHermite3(x))/180;
return result;
}

double h112( double x)
{
double result=( 14 * PolyHermite5(x) + 37 * PolyHermite3(x) + 8 * PolyHermite1(x))/288;
return result;
}

double h1111( double x)
{
double result=-( 252 * PolyHermite5(x) + 832 * PolyHermite3(x))/7776;
return result;
}


double g1( double x)
{
double result =  pow(x,2) - 1;
return result;
}

double g2( double x)
{
double result =  pow(x,3)  - 3 * x;
return result;
}

double g3( double x)
{
double result =  5 * x - 2 *  pow(x,3) ;
return result;
}

double g4( double x)
{
double result = -6 * pow(x,2) + pow(x,4) +3;
return result;
}

double g5( double x)
{
double result =  30 * pow(x,2)- 6* pow(x,4) -12;
return result;
}

double g6( double x)
{
double result = -106 * pow(x,2) + 24 * pow(x,4) + 34;
return result;
}

double g7( double x)
{
double result = 15 * x - 10 * pow(x,3) + pow(x,5);
return result;
}

double g8( double x)
{
double result = -84 * x + 68 * pow(x,3) - 8 * pow(x,5);
return result;
}

double g9( double x)
{
double result = 321 * x - 309 * pow(x,3) + 42 * pow(x,5);
return result;
}

double g10( double x)
{
double result = -87 * x + 72 * pow(x,3) - 9 * pow(x,5);
return result;
}

double g11( double x)
{
double result = -1511 * x + 1688 - 252 * pow(x,5);
return result;
}

double HermitePolynome(int n, double x)
{

double result;

if(n==0 || n==1)
result=1;
else
result= x * HermitePolynome(n-1,  x) - (n-1) * HermitePolynome(n-2, x);

return result;
}

/*
struct SimulationCIR

{
	BasicParameters basicParameters;
	SimulationParameters simulationParameters;

	SimulationCIR(const BasicParameters & basicParameters0, const SimulationParameters & simulationParameters0)
	
	{		
		basicParameters = basicParameters0;
		simulationParameters = simulationParameters0;
	}
};
  */


// typedef struct always employed
typedef struct _IntegerVector
{
int size;
int *degree;
} IntegerVector;


IntegerVector * CreateIntegerVector( int size0, int *degree )
{
int i;
IntegerVector *res;
res = malloc(sizeof(IntegerVector));

for(i=0; i<size0; i++)
res->degree[i]=degree[i];

return res;
}

void DestroyIntegerVector(IntegerVector ** IntegerVector)
{
	if (*IntegerVector == NULL) return;

	free((*IntegerVector)->degree);
	free(*IntegerVector);
	
	*IntegerVector = NULL;
}

typedef struct _Polynom
{
int degree;
double *coef;
} Polynom;

typedef struct _PolynomVector
{

Polynom *vector;
int size;

}PolynomVector;



typedef struct _PolynomVector
{

Polynom *vector;
int size;

}PolynomVector;



Polynom * CreateZeroPolynom(int degree0)
{
	int iPolynom;	
	Polynom * res;
	res = malloc(sizeof(Polynom));
	res->degree = degree0;
	res->coef = malloc((degree0+1) * sizeof(double));
	
	
	for(iPolynom=0;iPolynom<=degree0;iPolynom++)
	{
		res->coef[iPolynom]=0.0;
	}

	return res;

}


// CreatePolynom is a function like "new" in C++
Polynom * CreatePolynom(int degree0, double *coef0)
{
	int iPolynom;	
	Polynom * res = CreateZeroPolynom(degree0);
	for(iPolynom=0;iPolynom<=degree0;iPolynom++)
	{
		res->coef[iPolynom]=coef0[iPolynom];
	}
	return res;
}

PolynomVector * CreateGenericPolynomVector()
{

PolynomVector *res;
res = malloc(sizeof(PolynomVector));

return res;
}


PolynomVector * CreateZeroPolynomVector(IntegerVector degreeVector)
{
int i;
PolynomVector *res;
res->size = degreeVector.size;
res->vector = malloc(res->size * sizeof(Polynom));
	
	for(i=0; i<res->size; i++)
	{
		res->vector[i] = *CreateZeroPolynom(degreeVector.degree[i]);
	}

return res;
}


// DestroyPolynom is a function like "delete" in C++

void DestroyPolynom(Polynom ** poly)
{
	if (*poly == NULL) return;

	free((*poly)->coef);
	free(*poly);
	
	*poly = NULL;
}


void DestroyPolynomVector(PolynomVector ** polyVector)
{
	int i;
	
	if (*polyVector == NULL) return;
	
	for(i=0; i<=(*polyVector)->size; i++)	
	DestroyPolynom( ( &(*polyVector) -> vector[i]) );
	
	free(*polyVector);
	
	*polyVector = NULL;
}





Polynom * lambdaPolynom( double lambda, Polynom * P)
{
int i;

Polynom *res;
res = CreateZeroPolynom(P->degree);

for(i=0; i<= P->degree; i++)
{
	res->coef[i] = lambda * P->coef[i];
}

return res;

}


Polynom * SumPolynoms( Polynom * P, Polynom * Q)
{
	
	int iPolynom, iNomNullTest=0;	
	int minDegree, maxDegree;
	Polynom * res, *minPQ, *maxPQ;

if(P->degree!=Q->degree)

	{
	if(P->degree<Q->degree)
		{
			minDegree = P->degree;
			maxDegree = Q->degree;
			
			minPQ = CreateZeroPolynom(maxDegree);		
			maxPQ = CreateZeroPolynom(maxDegree);
			
			for(iPolynom = 0;iPolynom <= minDegree; iPolynom++)
			minPQ->coef[iPolynom] = P -> coef[iPolynom];
			maxPQ->coef[iPolynom] = Q -> coef[iPolynom];			
			
			for(iPolynom = minDegree + 1; iPolynom <= maxDegree; iPolynom++)
			minPQ->coef[iPolynom] = 0;	
			maxPQ->coef[iPolynom] = Q -> coef[iPolynom];			
		}

	else
		{
			minDegree = Q->degree;
			maxDegree = P->degree;

			minPQ = CreateZeroPolynom(maxDegree);		
			maxPQ = CreateZeroPolynom(maxDegree);
			
			for(iPolynom = 0;iPolynom <= minDegree; iPolynom++)
			minPQ->coef[iPolynom] = Q -> coef[iPolynom];
			maxPQ->coef[iPolynom] = P -> coef[iPolynom];
			
			for(iPolynom = minDegree + 1; iPolynom <= maxDegree; iPolynom++)
			minPQ->coef[iPolynom] = 0;		
			maxPQ->coef[iPolynom] = P -> coef[iPolynom];		
		}
	}

else
maxDegree = P->degree;
 


	res = CreateZeroPolynom(maxDegree);

	for(iPolynom = 0;iPolynom <= maxDegree; iPolynom++)
	{
		res->coef[iPolynom]=minPQ->coef[iPolynom] + maxPQ->coef[iPolynom];
		iNomNullTest += res->coef[iPolynom]==0 ? 0: 1;
	}

	if(iNomNullTest==0)
	{
		res->degree =0;
		realloc(res->coef,1);
	}

return res;

}



Polynom * ProductPolynoms( Polynom * P, Polynom * Q)
{
	int iPolynom, iConvol;	
	int minDegree, maxDegree, productDegree;
	Polynom * res, *minPQ, *maxPQ;

if( (Q->degree==0 && P->coef[0]==0) || (Q->degree==0 && P->coef[0]==0) )
		{
			res = CreateZeroPolynom(0);
			res->coef[0] = 0;
		}

else
		{
				productDegree = P->degree + Q->degree;
				
				if(P->degree!=Q->degree)

					{
					if(P->degree<Q->degree)
						{
							minDegree = P->degree;
							maxDegree = Q->degree;

							minPQ = CreateZeroPolynom(minDegree);
							maxPQ = CreateZeroPolynom(maxDegree);
							
							for(iPolynom = 0;iPolynom <= minDegree; iPolynom++)
							{
								minPQ->coef[iPolynom] = P -> coef[iPolynom];
								maxPQ->coef[iPolynom] = Q -> coef[iPolynom];
							}

							for(iPolynom = minDegree + 1; iPolynom <= maxDegree; iPolynom++)
							{
								minPQ->coef[iPolynom] = 0;	
								maxPQ->coef[iPolynom] = Q -> coef[iPolynom];
							}

							if(maxPQ->degree < productDegree)
							{
								for(iPolynom = maxDegree + 1; iPolynom <= productDegree; iPolynom++)
								{
									minPQ->coef[iPolynom] = 0;	
									maxPQ->coef[iPolynom] = 0;
								}
							}
							
						
						}


					else
						{
							minDegree = Q->degree;
							maxDegree = P->degree;

							minPQ = CreateZeroPolynom(minDegree);
							maxPQ = CreateZeroPolynom(maxDegree);	
							
							for(iPolynom = 0;iPolynom <= minDegree; iPolynom++)
							{
								minPQ->coef[iPolynom] = Q -> coef[iPolynom];
								maxPQ->coef[iPolynom] = P -> coef[iPolynom];
							}

							for(iPolynom = minDegree + 1; iPolynom <= maxDegree; iPolynom++)
							{
								minPQ->coef[iPolynom] = 0;	
								maxPQ->coef[iPolynom] = P -> coef[iPolynom];
							}

							if(maxPQ->degree < productDegree)
							{
								for(iPolynom = maxDegree + 1; iPolynom <= productDegree; iPolynom++)
								{
									minPQ->coef[iPolynom] = 0;	
									maxPQ->coef[iPolynom] = 0;
								}
							}

						}
					}



			res = CreateZeroPolynom(productDegree);

			for(iPolynom = 0;iPolynom <= maxDegree; iPolynom++)
			{
				res->coef[iPolynom]=0;
					
					for(iConvol = 0;iConvol <= iPolynom; iPolynom++)	
					res->coef[iPolynom]+= minPQ->coef[iConvol]* maxPQ->coef[iPolynom-iConvol];	
			}
		}

return res;

}


double PolynomValueAtx(Polynom *P, double x)
{
double result;
int i;

for(i=0;i<=P->degree;i++) 
{
result += P->coef[i] * pow(x,i);
}

return result;
}


double CornishFisher( double x, int n, double * cumulant)
{
double result;
result =

}
 


/*
Err PolylSum(Polynom P, polynom Q, Polynom &Sum)
{
int  iSum;

Sum.degree = max(P.degree,Q.degree);

for(iSum=0; iSum<=Sum.degree; iSum++)
Sum.coef[k]=P.coef[k]+Q.coef[k];

return 0;
}


Err PolynomialProduct(Polynom &P, polynom &Q, Polynom &Product)
{
int  iSum, iConvol;

if(P.Degree==0;
Product.Degree = P.degree  + Q.degree;

for(iSum=0; iSum<=Product.Degree; iSum++)

	{
		for(iConvol=0; iConvol<=iSum; iConvol++)
		Product.coef[iSum]=P.coef[iConvol] * Q.coef[iSum-iConvol];
	}

return 0;
}


err PolynomialSum(Polynom *polynomArray, int k, Polynom &Sum)
{
int iSum;

for(iSum==1; iSum<=k-1; iSum++)
PolynomialSum(polynomArray[], Polynom


}

double KsikDegree(int k)
{
int result, j;

if(k==1)
result = 1;

else
	{
	result = 1;
		for(j==1; j<=k-1; j++)
		result = max(KsikDegree(k - j) + j + 1 + 1, result);
	}
}


double Ksik( double k, double xAlpha, double *cumulant)
{
int j, degreeKsik =KsikDegree(k) ;

Polynom *p0;
Polynom *p1;
Polynom *p2;
Polynom *p3;

if(k==1)
p0->degree = 2;
p0->coef[2]=cumulant[2];
p0->coef[1]=0;
p0->coef[0]=0;

else
{
p0->degree = degreeKsik;
p0->coef[k+1]=cumulant[k];

	for(j==1; j <= k-1; j++)
	{

	p1->degree = degreeKsik; 
	p2->degree = degreeKsik;
	p3->degree = degreeKsik;

	}

}

}





double CornishFisherExpansionValue(double xAlpha, int n)
{
double result;
int k;
result= xAlpha;

for(k==1; k<=n; k++)
	{
		result += Ksik(xAlpha,k);
	}
return result;
}

*/


double ModifiedBesselFunction( double nu ,double  x)
{


double result;
long nnu = (long)floor(nu);
long sizenu = abs(nnu);


Complex * tabb = malloc((sizenu + 1) * sizeof(Complex));


NagError fail;
memset(&fail,0,sizeof(NagError));
nag_bessel_i_alpha (x, nu - nnu, nnu, tabb, &fail);
result = tabb[sizenu].re;

free(tabb);

return result;


}




Err BesselArrays( double nu ,double  x, double *arraynup1)
{

if (nu>=0)
{
long nnu = (long)floor(nu);
long sizenu = abs(nnu)+2;
int iFillArray;
Complex * tabb = malloc((sizenu) * sizeof(Complex));
NagError fail;
memset(&fail,0,sizeof(NagError));
nag_bessel_i_alpha (x, nu - nnu, nnu+1, tabb, &fail);
for(iFillArray=0; iFillArray < sizenu; iFillArray++)
arraynup1[iFillArray] = tabb[iFillArray].re;
}

else
{
long nnu = (long)floor(nu);
long sizenu = abs(nnu);
int iFillArray;
Complex * tabb = malloc((sizenu + 1) * sizeof(Complex));
NagError fail;
memset(&fail,0,sizeof(NagError));
nag_bessel_i_alpha (x, nu - nnu, nnu, tabb, &fail);
for(iFillArray=0; iFillArray < sizenu + 1; iFillArray++)
arraynup1[iFillArray] = tabb[iFillArray].re;
}


return NULL;

}


Err	 lgmSV_mc_balsam(	
					/*	Time Information  */
					int			iNbTime,
					int			iNbEvent,
					double		*dTime,
					double		*dDate,
					
					int			iNumPaths,					
						
					/*	Model data Information	*/
					double		dLambdaX,
					double		*dSigma,

					double		*dAlpha,
					double		*dLambdaEps,
					double		*dLvlEps,
					double		*dRho,

					/* Parameters for DF(t,T*) reconstruction */
					double		*dff_star,
					double		*gam_star,
					double		*gam2_star,
					
					/* Parameters */
					LGMSVParam	Params,

					/*	Product data */
					void		**func_parm_tab, 
					int			*EvalEvent,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,

					/*	Initialisation function to be called at the beggining of each path or NULL if none */
					void		(*init_func)(),
					
					/*	Payoff function */
					Err (*payoff_func)(	long	path_index,
										double	evt_date,
										double	evt_time,
										void	*func_parm,
										
										double	ft,
										double	psi,
										double	v,
																
										/* Vector of results to be updated */
										int		nprod,
										/* Result	*/
										double	*prod_val,
										int		*stop_path),
					/*	Result */
					int			iNbProduct, 
					double		**res)
{
Err		err = NULL;
clock_t	time1, time2;
long	i, j, k, m;

double	*sum_payoff			= NULL,
		*sum_payoff2		= NULL,
		*path_payoff		= NULL,
		*res_evt			= NULL,				
		**matrix			= NULL,
		*matrixi			= NULL,
				yInterpolated		,

		*muf				= NULL,
		*sigf				= NULL,
		*sigf2				= NULL,
		*sigf3				= NULL,
		*sigf4				= NULL,
		*sigf5				= NULL,
		*sigf6				= NULL,
		*mupsi				= NULL,
		*mupsi2				= NULL,
		*sigpsi				= NULL,
		*muv1				= NULL,
		*muv2				= NULL,
		*sigv				= NULL,
		*sigv2				= NULL,
		*sigv3				= NULL,
		*sigv4				= NULL,
		*sigv5				= NULL,
		*sigv6				= NULL,		
		*sigv7				= NULL,
		*sigv8				= NULL,
		*sigv9				= NULL,			
		*sigv10				= NULL,
		*sigv11				= NULL,
		*sigv12				= NULL,		
		*sigv13				= NULL,					
		*coef1				= NULL,
		*coef2				= NULL,

		*b1					= NULL,
		*b2					= NULL,
		*b3					= NULL,
							
		*c2					= NULL,
		*c3					= NULL,
		*c1					= NULL,
		*c4					= NULL,
		*c5					= NULL,
		*c6					= NULL,
		*c7					= NULL,
		
		*n0					= NULL,
		*n2					= NULL,
		*n4					= NULL,
		*d1					= NULL,
		**gamma				= NULL,
		**cumulant1			= NULL,
		**cumulant2			= NULL,
		**cumulant3			= NULL,
		*nu					= NULL,
		*xAux				= NULL,			
		**coefBesselApprox1	= NULL,	
		**coefBesselApprox2	= NULL,	
		**xValues								= NULL,
		**TabulatedBesselFunction				= NULL,	
		**TabulatedBesselFunctionSecondDerivate	= NULL,
		**ModifiedBesselFunctionValues			= NULL,				
		***save_values		= NULL;




int     iTabCount, NbTabCount=5000;
int		stop_path;
int		evtindex;
int		optim_today;
long     iBesselApprox, besselApproxOrder = 2;
int     NbCumulants = 7;
double	 v0, sqrtV, sumV, difV, x, argGamma1, argGamma2, besselApproxNum, besselApproxDen, SimpleBesselFunction, BigBesselFunction ;

double  xInfBound = 2.0e-308, xSupBound = 1000 + 2.0e-308;
int		sizeNuMain;
double	dt, sqdt ;
double	ft, v, sqv, psi;
double	brow1, brow2;
double	df;
double	log_drift, log_var;

long	seed = -123456789;

	/* --------------------------------------------------------------------------------------------------
													Initialisation		
	   -------------------------------------------------------------------------------------------------- */	
	
	/* For computational time calculation				 */
	time1 = clock();		

	/* odd number of paths */
	iNumPaths = (int) (iNumPaths / 2) * 2 + 1;
					
	matrix = dmatrix (0, iNumPaths - 1, 0, 2 * (iNbTime - 1) - 1);
	res_evt = dvector (0, iNbProduct - 1);
	path_payoff = dvector (0, iNbProduct - 1); 
	sum_payoff = dvector (0, iNbProduct - 1);
	sum_payoff2 = dvector (0, iNbProduct - 1);

	/* for precalculations */
	sigf = calloc(iNbTime, sizeof(double));
	sigpsi = calloc(iNbTime, sizeof(double));
	muv1 = calloc(iNbTime, sizeof(double));
	muv2 = calloc(iNbTime, sizeof(double));
	sigv = calloc(iNbTime, sizeof(double));
	coef1 = calloc(iNbTime, sizeof(double));
	coef2 = calloc(iNbTime, sizeof(double));

	if (!matrix || !path_payoff || !sum_payoff || !sum_payoff2 || !res_evt ||
		!sigf || !sigpsi || !muv1 || !muv2 || !sigv || !coef1 || !coef2)
	{
		err = "Memory allocation failure in lgmSV_mc_balsam";
		goto FREE_RETURN;
	}

	if (do_optimisation)
	{
		err = mceb_allocate_savevalues_for_GRFN(iNumPaths,
												iNbEvent,
												params,
												&save_values);

		if (err) goto FREE_RETURN;
	}

	memset (sum_payoff, 0, iNbProduct * sizeof (double));
	memset (sum_payoff2, 0, iNbProduct * sizeof (double));

	/* All the needed precalculations */
	switch (Params.iSchemeOrder)
	{
		case 0:
		{
			for (j=1; j<iNbTime; j++)
			{
				dt = dTime[j] - dTime[j-1];
				sqdt = sqrt(dt);
				
				sigf[j] = dSigma[j] * sqdt;		
				sigpsi[j] = sigf[j] * sigf[j];
				
				sigv[j] = dAlpha[j] * sqdt;
				muv1[j] = 1.0 - dLambdaEps[j] * dt;
				muv2[j] = dLvlEps[j] * dt;				

				coef1[j] = dRho[j];
				coef2[j] = sqrt(1.0 - coef1[j] * coef1[j]);				
			}

			break;
		}

		case 1:
		{
			for (j=1; j<iNbTime; j++)
			{
				dt = dTime[j] - dTime[j-1];
				sqdt = sqrt(dt);
				
				sigf[j] = dSigma[j] * sqdt;		
				sigpsi[j] = sigf[j] * sigf[j];				

				if (fabs(dLambdaEps[j]) > 1.0E-08)
				{
					muv1[j] = exp(-dLambdaEps[j] * dt);
					muv2[j] = dLvlEps[j] / dLambdaEps[j] * (1.0 - muv1[j]);
					sigv[j] = dAlpha[j] * sqrt(0.5 / dLambdaEps[j] * (1.0 - muv1[j] * muv1[j]));
					sigv[j] *= sigv[j];
				}
				else
				{
					muv1[j] = 1.0;
					muv2[j] = dLvlEps[j] * dt;
					sigv[j] = dAlpha[j] * sqdt;
					sigv[j] *= sigv[j];
				}

				coef1[j] = dRho[j];
				coef2[j] = sqrt(1.0 - coef1[j] * coef1[j]);
			}

			break;
		}


	case 2:

		{
			muf = calloc(iNbTime, sizeof(double));
			sigf2 = calloc(iNbTime, sizeof(double));
			sigf3 = calloc(iNbTime, sizeof(double));
			sigf4 = calloc(iNbTime, sizeof(double));
			sigf5 = calloc(iNbTime, sizeof(double));
			sigf6 = calloc(iNbTime, sizeof(double));

			mupsi = calloc(iNbTime, sizeof(double));
			mupsi2 = calloc(iNbTime, sizeof(double));

			sigv2 = calloc(iNbTime, sizeof(double));
			sigv3 = calloc(iNbTime, sizeof(double));

			if (!muf || !sigf2 || !sigf3 || !sigf4 || !sigf5 || !sigf6 || !sigv2 || !sigv3 || !mupsi || !mupsi2)
			{
				err = "Memory allocation failure in lgmSV_mc_balsam";
				goto FREE_RETURN;
			}		

			for (j=1; j<iNbTime; j++)
			{
				dt = dTime[j] - dTime[j-1];
				sqdt = sqrt(dt);
			
				coef1[j] = dRho[j];
				coef2[j] = sqrt(1.0 - coef1[j] * coef1[j]);

				muf[j] = -0.25 * dRho[j] * dAlpha[j] * dSigma[j] * dt;
				sigf[j] = coef2[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;
				sigf2[j] = coef2[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf3[j] = coef1[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;
				sigf4[j] = coef1[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf5[j] = 0.25 * coef1[j] * dAlpha[j] * dSigma[j] * dt;
				sigf6[j] = 0.25 * coef2[j] * dAlpha[j] * dSigma[j] * dt;
							
				mupsi[j] = 0.5 * dLvlEps[j] * dSigma[j] * dSigma[j] * dt * dt;
				mupsi2[j] = dSigma[j] * dSigma[j] * dt * (1.0 - 0.25 * dLambdaEps[j] * dt);
				sigpsi[j] = 0.5 * dAlpha[j] * dSigma[j] * dSigma[j] * dt * sqdt;
										
				muv1[j] = (dLvlEps[j] - 0.25 * dAlpha[j] * dAlpha[j] 
							-0.5 * dLvlEps[j] * dLambdaEps[j] * dt) * dt;			
				muv2[j] = (1.0 - dLambdaEps[j] * dt * (1.0 - 0.5 * dLambdaEps[j] * dt));			
				sigv[j] = dAlpha[j] * (1.0 - 0.75 * dLambdaEps[j] * dt) * sqdt;			
				sigv2[j] = dAlpha[j] / 16.0 * (4.0 * dLvlEps[j] - dAlpha[j] * dAlpha[j]) * dt * sqdt;			
				sigv3[j] = 0.25 * dAlpha[j] * dAlpha[j] * dt;
			}

			break;
		}




// case 3: Volatility Projection 

	
		case 3:
		{
			/* Extra allocations */
			muf = calloc(iNbTime, sizeof(double));
			sigf2 = calloc(iNbTime, sizeof(double));
			sigf3 = calloc(iNbTime, sizeof(double));
			sigf4 = calloc(iNbTime, sizeof(double));
			sigf5 = calloc(iNbTime, sizeof(double));
			sigf6 = calloc(iNbTime, sizeof(double));

			mupsi = calloc(iNbTime, sizeof(double));
			mupsi2 = calloc(iNbTime, sizeof(double));

			sigv2 = calloc(iNbTime, sizeof(double));
			sigv3 = calloc(iNbTime, sizeof(double));
			sigv4 = calloc(iNbTime, sizeof(double));
			sigv5 = calloc(iNbTime, sizeof(double));
			sigv6 = calloc(iNbTime, sizeof(double));
			sigv7 = calloc(iNbTime, sizeof(double));
			sigv8 = calloc(iNbTime, sizeof(double));
			sigv9 = calloc(iNbTime, sizeof(double));
			sigv10 = calloc(iNbTime, sizeof(double));
			sigv11 = calloc(iNbTime, sizeof(double));
			sigv12 = calloc(iNbTime, sizeof(double));
			sigv13 = calloc(iNbTime, sizeof(double));


			if (!muf || !sigf2 || !sigf3 || !sigf4 || !sigf5 || !sigf6 || !sigv2 || !sigv3 || !sigv4 || !sigv5 || !sigv6 || !sigv7 || !sigv8 || !sigv9 || !sigv10 || !sigv11 || !sigv12 || !sigv13  || !mupsi || !mupsi2)
			{
				err = "Memory allocation failure in lgmSV_mc_balsam";
				goto FREE_RETURN;
			}		

			for (j=1; j<iNbTime; j++)
			{
				dt = dTime[j] - dTime[j-1];
				sqdt = sqrt(dt);
			
				coef1[j] = dRho[j];
				coef2[j] = sqrt(1.0 - coef1[j] * coef1[j]);

				muf[j] = -0.25 * dRho[j] * dAlpha[j] * dSigma[j] * dt;
				sigf[j] = coef2[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;
				sigf2[j] = coef2[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf3[j] = coef1[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;
				sigf4[j] = coef1[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf5[j] = 0.25 * coef1[j] * dAlpha[j] * dSigma[j] * dt;
				sigf6[j] = 0.25 * coef2[j] * dAlpha[j] * dSigma[j] * dt;
							
				mupsi[j] = 0.5 * dLvlEps[j] * dSigma[j] * dSigma[j] * dt * dt;
				mupsi2[j] = dSigma[j] * dSigma[j] * dt * (1.0 - 0.25 * dLambdaEps[j] * dt);
				sigpsi[j] = 0.5 * dAlpha[j] * dSigma[j] * dSigma[j] * dt * sqdt;
										
				muv1[j] = (dLvlEps[j] - 0.25 * dAlpha[j] * dAlpha[j] 
							-0.5 * dLvlEps[j] * dLambdaEps[j] * dt) * dt;			
				muv2[j] = (1.0 - dLambdaEps[j] * dt * (1.0 - 0.5 * dLambdaEps[j] * dt));			
				sigv[j] = dAlpha[j] * (1.0 - 0.75 * dLambdaEps[j] * dt) * sqdt;			
				sigv2[j] = dAlpha[j] / 16.0 * (4.0 * dLvlEps[j] - dAlpha[j] * dAlpha[j]) * dt * sqdt;			
				sigv3[j] = 0.25 * dAlpha[j] * dAlpha[j] * dt;


// v2: -1.5	 -0.5	0.5


				sigv6[j] = (-11 * pow(dAlpha[j] , 6) + 56 * pow(dAlpha[j] , 4)  * dLvlEps[j] - 48 * pow (dAlpha[j] * dLvlEps[j] , 2)) * dt * dt / ( dAlpha[j] *  1536 );
				sigv7[j] = (24 * pow(dAlpha[j] , 4)*  dLambdaEps[j] - 96 * pow (dAlpha[j] , 2) *  dLambdaEps[j] * dLvlEps[j]) * dt * dt / ( dAlpha[j] *  1536 ) ;
				sigv8[j] = 464 * pow (dAlpha[j]  * dLambdaEps[j] , 2) * dt * dt / ( dAlpha[j] *  1536 );				


// w1: 0 1 


				sigv4[j] = (( pow (dAlpha[j] , 4) - 4.0 * pow ( dAlpha[j] , 2 )  * dLvlEps[j] ) / 96.0 ) * dt;
				sigv5[j] = (-12.0 * pow ( dAlpha[j] , 2 ) * dLambdaEps[j]  / 96.0 ) * dt; 



// w2:	-2	-1	0	1	2


				sigv9[j] =( 1 / (36864 * pow(dAlpha[j], 2))) *  (576 * dAlpha[j] * dAlpha[j] * dLvlEps[j] * dLvlEps[j] * dLvlEps[j] + 364 * dLvlEps[j] * pow( dLvlEps[j], 6) -944 * dLvlEps[j] * dLvlEps[j] * pow( dAlpha[j], 4) -41 * pow(dAlpha[j],8) ) * dt * dt;
				sigv10[j] =( 1 / (36864 * pow(dAlpha[j], 2))) * (-576 * pow( dAlpha[j] * dLvlEps[j] ,2) * dLambdaEps[j] -132 *  dLambdaEps[j] * pow( dAlpha[j], 6) + 672 * dLvlEps[j] * dLambdaEps[j] *pow( dAlpha[j], 4)) * dt * dt;		
				sigv11[j] =( 1 / (36864 * pow(dAlpha[j], 2))) * (1728 * dLvlEps[j] * pow( dAlpha[j] * dLambdaEps[j] , 2)+ 1104 * pow( dLambdaEps[j] * dAlpha[j] * dAlpha[j],2)) * dt * dt;
				sigv12[j] = ( 1 / (36864 * pow(dAlpha[j], 2))) * (-6336 * dAlpha[j] * dAlpha[j] * pow( dLambdaEps[j], 3) -3072 * dLvlEps[j] *pow(dLambdaEps[j],3) )* dt * dt;
				sigv13[j] = ( 1 / (36864 * pow(dAlpha[j], 2))) * (3072 * pow(dLambdaEps[j], 4) )* dt * dt;
		
			}

			break;
			}

//Vol Projection + Expectation of Integral Volatility as a Bessel Function	 

				case 4:
		{
			/* Extra allocations */
			muf		= calloc(iNbTime, sizeof(double));
			sigf2	= calloc(iNbTime, sizeof(double));
			sigf3	= calloc(iNbTime, sizeof(double));
			sigf4	= calloc(iNbTime, sizeof(double));
			sigf5	= calloc(iNbTime, sizeof(double));
			sigf6	= calloc(iNbTime, sizeof(double));
			mupsi	= calloc(iNbTime, sizeof(double));
			mupsi2	= calloc(iNbTime, sizeof(double));
			sigv2	= calloc(iNbTime, sizeof(double));
			sigv3	= calloc(iNbTime, sizeof(double));


			b1= calloc(iNbTime, sizeof(double));
			b2= calloc(iNbTime, sizeof(double));
			b3= calloc(iNbTime, sizeof(double));
			c2= calloc(iNbTime, sizeof(double));
			c3= calloc(iNbTime, sizeof(double));
			c1= calloc(iNbTime, sizeof(double));
			c4= calloc(iNbTime, sizeof(double));
			c5= calloc(iNbTime, sizeof(double));
			c6= calloc(iNbTime, sizeof(double));
			c7= calloc(iNbTime, sizeof(double));

			nu= calloc(iNbTime, sizeof(double));
			xAux = calloc(iNbTime, sizeof(double));

			if (!muf || !sigf2 || !sigf3 || !sigf4 || !sigf5 || !sigf6 || !sigv2 || !sigv3 || !mupsi || !b1 || !b2 || !b3 || !c1 || !c2 || !c3 || !c4 || !c5 || !c6 || !c7 || !nu || !xAux)
			{
				err = "Memory allocation failure in lgmSV_mc_balsam";
				goto FREE_RETURN;
			}		

			for (j=1; j<iNbTime; j++)
			{
				dt = dTime[j] - dTime[j-1];
				sqdt = sqrt(dt);
			
				coef1[j] = dRho[j];
				coef2[j] = sqrt(1.0 - coef1[j] * coef1[j]);

				muf[j] = -0.25 * dRho[j] * dAlpha[j] * dSigma[j] * dt;
				sigf[j] = coef2[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;
				sigf2[j] = coef2[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf3[j] = coef1[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;
				sigf4[j] = coef1[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf5[j] = 0.25 * coef1[j] * dAlpha[j] * dSigma[j] * dt;
				sigf6[j] = 0.25 * coef2[j] * dAlpha[j] * dSigma[j] * dt;
							
				mupsi[j] = 0.5 * dLvlEps[j] * dSigma[j] * dSigma[j] * dt * dt;
				mupsi2[j] = dSigma[j] * dSigma[j] * dt * (1.0 - 0.25 * dLambdaEps[j] * dt);
				sigpsi[j] = 0.5 * dAlpha[j] * dSigma[j] * dSigma[j] * dt * sqdt;
										
				muv1[j] = (dLvlEps[j] - 0.25 * dAlpha[j] * dAlpha[j] 
							-0.5 * dLvlEps[j] * dLambdaEps[j] * dt) * dt;			
				muv2[j] = (1.0 - dLambdaEps[j] * dt * (1.0 - 0.5 * dLambdaEps[j] * dt));			
				sigv[j] = dAlpha[j] * (1.0 - 0.75 * dLambdaEps[j] * dt) * sqdt;			
				sigv2[j] = dAlpha[j] / 16.0 * (4.0 * dLvlEps[j] - dAlpha[j] * dAlpha[j]) * dt * sqdt;			
				sigv3[j] = 0.25 * dAlpha[j] * dAlpha[j] * dt;

				b1[j] =( - 1 + exp( - dLambdaEps[j] * dt) + 0.5 * dLambdaEps[j] * dt * ( 1+ exp( - dLambdaEps[j] * dt)) );
				b2[j] =dAlpha[j] * dAlpha[j] / (( 1- exp( - dLambdaEps[j] * dt)) *  dLambdaEps[j] * dLambdaEps[j]);
				b3[j] = 4 * exp(-0.5 * dLambdaEps[j] * dt) / (dLambdaEps[j] * pow( - 1+ exp(-dLambdaEps[j] * dt),2));



				c2[j]= dSigma[j]  * dRho[j] / dAlpha[j];
				c3[j]= c2[j] * dLambdaEps[j];
				c1[j]= - c2[j] * dLvlEps[j] * dt;
				c4[j]= dSigma[j] * sqrt(1-dRho[j] * dRho[j]);
				c5[j]= b1[j]* b2[j];
				c6[j]= b1[j]* b3[j];
				c7[j]= - ( 2 * dLambdaEps[j] *dt * exp( - dLambdaEps[j] * dt ) - 1 + exp( - 2 * dLambdaEps[j] * dt)) / (dLambdaEps[j] * pow( 1- exp( - dLambdaEps[j] * dt),2));
				nu[j]= 2 * dLvlEps[j] / pow( dAlpha[j], 2) - 1;
				xAux[j]= 4 * dLambdaEps[j] * exp( - 0.5 * dLambdaEps[j] * dt) / ( pow( dAlpha[j], 2) * ( 1 - exp( - dLambdaEps[j] * dt )) );

			}

			break;
		}


		case 5:
		{
			
			/* Extra allocations */
			
			muf		= calloc(iNbTime, sizeof(double));
			sigf2	= calloc(iNbTime, sizeof(double));
			sigf3	= calloc(iNbTime, sizeof(double));
			sigf4	= calloc(iNbTime, sizeof(double));
			sigf5	= calloc(iNbTime, sizeof(double));
			sigf6	= calloc(iNbTime, sizeof(double));
			mupsi	= calloc(iNbTime, sizeof(double));
			mupsi2	= calloc(iNbTime, sizeof(double));
			sigv2	= calloc(iNbTime, sizeof(double));
			sigv3	= calloc(iNbTime, sizeof(double));


			b1= calloc(iNbTime, sizeof(double));
			b2= calloc(iNbTime, sizeof(double));
			b3= calloc(iNbTime, sizeof(double));
			c2= calloc(iNbTime, sizeof(double));
			c3= calloc(iNbTime, sizeof(double));
			c1= calloc(iNbTime, sizeof(double));
			c4= calloc(iNbTime, sizeof(double));
			c5= calloc(iNbTime, sizeof(double));
			c6= calloc(iNbTime, sizeof(double));
			c7= calloc(iNbTime, sizeof(double));
			coefBesselApprox1 = dmatrix(0, iNbTime-1, 0,besselApproxOrder-1);
			coefBesselApprox2 = dmatrix(0, iNbTime-1, 0,besselApproxOrder-1);

			nu= calloc(iNbTime, sizeof(double));
			xAux = calloc(iNbTime, sizeof(double));

			if (!muf || !sigf2 || !sigf3 || !sigf4 || !sigf5 || !sigf6 || !sigv2 || !sigv3 || !mupsi || !b1 || !b2 || !b3 || !c1 || !c2 || !c3 || !c4 || !c5 || !c6 || !c7 || !nu || !xAux || !coefBesselApprox1 || !coefBesselApprox2 )
			{
				err = "Memory allocation failure in lgmSV_mc_balsam";
				goto FREE_RETURN;
			}		

			for (j=1; j<iNbTime; j++)
			{
				dt = dTime[j] - dTime[j-1];
				sqdt = sqrt(dt);
			
				coef1[j] = dRho[j];
				coef2[j] = sqrt(1.0 - coef1[j] * coef1[j]);

				muf[j] = -0.25 * dRho[j] * dAlpha[j] * dSigma[j] * dt;
				sigf[j] = coef2[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;
				sigf2[j] = coef2[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf3[j] = coef1[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;
				sigf4[j] = coef1[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf5[j] = 0.25 * coef1[j] * dAlpha[j] * dSigma[j] * dt;
				sigf6[j] = 0.25 * coef2[j] * dAlpha[j] * dSigma[j] * dt;
							
				mupsi[j] = 0.5 * dLvlEps[j] * dSigma[j] * dSigma[j] * dt * dt;
				mupsi2[j] = dSigma[j] * dSigma[j] * dt * (1.0 - 0.25 * dLambdaEps[j] * dt);
				sigpsi[j] = 0.5 * dAlpha[j] * dSigma[j] * dSigma[j] * dt * sqdt;
										
				muv1[j] = (dLvlEps[j] - 0.25 * dAlpha[j] * dAlpha[j] 
							-0.5 * dLvlEps[j] * dLambdaEps[j] * dt) * dt;			
				muv2[j] = (1.0 - dLambdaEps[j] * dt * (1.0 - 0.5 * dLambdaEps[j] * dt));			
				sigv[j] = dAlpha[j] * (1.0 - 0.75 * dLambdaEps[j] * dt) * sqdt;			
				sigv2[j] = dAlpha[j] / 16.0 * (4.0 * dLvlEps[j] - dAlpha[j] * dAlpha[j]) * dt * sqdt;			
				sigv3[j] = 0.25 * dAlpha[j] * dAlpha[j] * dt;

				b1[j] = ( - 1 + exp( - dLambdaEps[j] * dt) + 0.5 * dLambdaEps[j] * dt * ( 1+ exp( - dLambdaEps[j] * dt)) );
				b2[j] = dAlpha[j] * dAlpha[j] / (( 1- exp( - dLambdaEps[j] * dt)) *  dLambdaEps[j] * dLambdaEps[j]);
				b3[j] = 4 * exp(-0.5 * dLambdaEps[j] * dt) / (dLambdaEps[j] * pow( - 1+ exp(-dLambdaEps[j] * dt),2));



				c2[j]= dSigma[j]  * dRho[j] / dAlpha[j];
				c3[j]= c2[j] * dLambdaEps[j];
				c1[j]= - c2[j] * dLvlEps[j] * dt;
				c4[j]= dSigma[j] * sqrt(1-dRho[j] * dRho[j]);
				c5[j]= b1[j]* b2[j];
				c6[j]= b1[j]* b3[j];
				c7[j]= - ( 2 * dLambdaEps[j] *dt * exp( - dLambdaEps[j] * dt ) - 1 + exp( - 2 * dLambdaEps[j] * dt)) / (dLambdaEps[j] * pow( 1- exp( - dLambdaEps[j] * dt),2));
				nu[j]= 2 * dLvlEps[j] / pow( dAlpha[j], 2) - 1;
				xAux[j]= 4 * dLambdaEps[j] * exp( - 0.5 * dLambdaEps[j] * dt) / ( pow( dAlpha[j], 2) * ( 1 - exp( - dLambdaEps[j] * dt )) );

				for( iBesselApprox =0; iBesselApprox < besselApproxOrder;  iBesselApprox++)
				
				{
					argGamma1 = nu[j]+ 1 + iBesselApprox + 1;
					coefBesselApprox1[j][iBesselApprox] = pow(0.5, nu[j] + 1) * pow( 0.25, iBesselApprox) / ( fact(iBesselApprox) * gamma_( &argGamma1));
					argGamma2 = nu[j]+ iBesselApprox + 1;
					coefBesselApprox2[j][iBesselApprox] = pow(0.5, nu[j]) * pow( 0.25, iBesselApprox) / ( fact(iBesselApprox) * gamma_( &argGamma2 ));
				}


			}

			break;
		}


		case 6:
		{

				/* Extra allocations */
			muf		= calloc(iNbTime, sizeof(double));
			sigf2	= calloc(iNbTime, sizeof(double));
			sigf3	= calloc(iNbTime, sizeof(double));
			sigf4	= calloc(iNbTime, sizeof(double));
			sigf5	= calloc(iNbTime, sizeof(double));
			sigf6	= calloc(iNbTime, sizeof(double));
			mupsi	= calloc(iNbTime, sizeof(double));
			mupsi2	= calloc(iNbTime, sizeof(double));
			sigv2	= calloc(iNbTime, sizeof(double));
			sigv3	= calloc(iNbTime, sizeof(double));
			sigv4 = calloc(iNbTime, sizeof(double));
			sigv5 = calloc(iNbTime, sizeof(double));
			sigv6 = calloc(iNbTime, sizeof(double));
			sigv7 = calloc(iNbTime, sizeof(double));
			sigv8 = calloc(iNbTime, sizeof(double));
			sigv9 = calloc(iNbTime, sizeof(double));
			sigv10 = calloc(iNbTime, sizeof(double));
			sigv11 = calloc(iNbTime, sizeof(double));
			sigv12 = calloc(iNbTime, sizeof(double));
			sigv13 = calloc(iNbTime, sizeof(double));
			b1= calloc(iNbTime, sizeof(double));
			b2= calloc(iNbTime, sizeof(double));
			b3= calloc(iNbTime, sizeof(double));
			c2= calloc(iNbTime, sizeof(double));
			c3= calloc(iNbTime, sizeof(double));
			c1= calloc(iNbTime, sizeof(double));
			c4= calloc(iNbTime, sizeof(double));
			c5= calloc(iNbTime, sizeof(double));
			c6= calloc(iNbTime, sizeof(double));
			c7= calloc(iNbTime, sizeof(double));
			nu= calloc(iNbTime, sizeof(double));
			xAux = calloc(iNbTime, sizeof(double));

			if	(xValues)									free_dmatrix (xValues, 0, iNbTime-1, 1,NbTabCount+1);xValues=0;												      
			if	(TabulatedBesselFunction)					free_dmatrix (TabulatedBesselFunction, 0, iNbTime-1, 1,NbTabCount+1);TabulatedBesselFunction=0;						
			if	(TabulatedBesselFunctionSecondDerivate)		free_dmatrix (TabulatedBesselFunctionSecondDerivate, 0, iNbTime-1, 1,NbTabCount+1);TabulatedBesselFunctionSecondDerivate=0;	

			

			xValues = dmatrix(0, iNbTime-1, 1,NbTabCount+1);	
			TabulatedBesselFunction = dmatrix(0, iNbTime-1, 1,NbTabCount+1);			
			TabulatedBesselFunctionSecondDerivate = dmatrix(0, iNbTime-1, 1,NbTabCount+1);


			if (!muf || !sigf2 || !sigf3 || !sigf4 || !sigf5 || !sigf6 || !sigv2 || !sigv3  || !sigv4 || !sigv5 || !sigv6 || !sigv7 || !sigv8 || !sigv9 || !sigv10 || !sigv11 || !sigv12 || !sigv13   || !mupsi  || !mupsi2 || !b1 || !b2 || !b3 || !c1 || !c2 || !c3 || !c4 || !c5 || !c6 || !c7 || !nu || !xAux)

			{
				err = "Memory allocation failure in lgmSV_mc_balsam";
				goto FREE_RETURN;
			}		

			for (j=1; j<iNbTime; j++)
			{
				dt = dTime[j] - dTime[j-1];
				sqdt = sqrt(dt);
			
				coef1[j] = dRho[j];
				coef2[j] = sqrt(1.0 - coef1[j] * coef1[j]);

				muf[j] = -0.25 * dRho[j] * dAlpha[j] * dSigma[j] * dt;
				sigf[j] = coef2[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;

				sigf2[j] = coef2[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf3[j] = coef1[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;
				sigf4[j] = coef1[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf5[j] = 0.25 * coef1[j] * dAlpha[j] * dSigma[j] * dt;
				sigf6[j] = 0.25 * coef2[j] * dAlpha[j] * dSigma[j] * dt;
							
				mupsi[j] = 0.5 * dLvlEps[j] * dSigma[j] * dSigma[j] * dt * dt;
				mupsi2[j] = dSigma[j] * dSigma[j] * dt * (1.0 - 0.25 * dLambdaEps[j] * dt);
				sigpsi[j] = 0.5 * dAlpha[j] * dSigma[j] * dSigma[j] * dt * sqdt;
										
				muv1[j] = (dLvlEps[j] - 0.25 * dAlpha[j] * dAlpha[j] 
							-0.5 * dLvlEps[j] * dLambdaEps[j] * dt) * dt;			
				muv2[j] = (1.0 - dLambdaEps[j] * dt * (1.0 - 0.5 * dLambdaEps[j] * dt));			
				sigv[j] = dAlpha[j] * (1.0 - 0.75 * dLambdaEps[j] * dt) * sqdt;			
				sigv2[j] = dAlpha[j] / 16.0 * (4.0 * dLvlEps[j] - dAlpha[j] * dAlpha[j]) * dt * sqdt;			
				sigv3[j] = 0.25 * dAlpha[j] * dAlpha[j] * dt;

// v2: -1.5	 -0.5	0.5


				sigv6[j] = (-11 * pow(dAlpha[j] , 6) + 56 * pow(dAlpha[j] , 4)  * dLvlEps[j] - 48 * pow (dAlpha[j] * dLvlEps[j] , 2)) * dt * dt / ( dAlpha[j] *  1536 );
				sigv7[j] = (24 * pow(dAlpha[j] , 4)*  dLambdaEps[j] - 96 * pow (dAlpha[j] , 2) *  dLambdaEps[j] * dLvlEps[j]) * dt * dt / ( dAlpha[j] *  1536 ) ;
				sigv8[j] = 464 * pow (dAlpha[j]  * dLambdaEps[j] , 2) * dt * dt / ( dAlpha[j] *  1536 );				


// w1: 0 1 


				sigv4[j] = (( pow (dAlpha[j] , 4) - 4.0 * pow ( dAlpha[j] , 2 )  * dLvlEps[j] ) / 96.0 ) * dt;
				sigv5[j] = (-12.0 * pow ( dAlpha[j] , 2 ) * dLambdaEps[j]  / 96.0 ) * dt; 



// w2:	-2	-1	0	1	2


				sigv9[j] =( 1 / (36864 * pow(dAlpha[j], 2))) *  (576 * dAlpha[j] * dAlpha[j] * dLvlEps[j] * dLvlEps[j] * dLvlEps[j] + 364 * dLvlEps[j] * pow( dLvlEps[j], 6) -944 * dLvlEps[j] * dLvlEps[j] * pow( dAlpha[j], 4) -41 * pow(dAlpha[j],8) ) * dt * dt;
				sigv10[j] =( 1 / (36864 * pow(dAlpha[j], 2))) * (-576 * pow( dAlpha[j] * dLvlEps[j] ,2) * dLambdaEps[j] -132 *  dLambdaEps[j] * pow( dAlpha[j], 6) + 672 * dLvlEps[j] * dLambdaEps[j] *pow( dAlpha[j], 4)) * dt * dt;		
				sigv11[j] =( 1 / (36864 * pow(dAlpha[j], 2))) * (1728 * dLvlEps[j] * pow( dAlpha[j] * dLambdaEps[j] , 2)+ 1104 * pow( dLambdaEps[j] * dAlpha[j] * dAlpha[j],2)) * dt * dt;
				sigv12[j] = ( 1 / (36864 * pow(dAlpha[j], 2))) * (-6336 * dAlpha[j] * dAlpha[j] * pow( dLambdaEps[j], 3) -3072 * dLvlEps[j] *pow(dLambdaEps[j],3) )* dt * dt;
				sigv13[j] = ( 1 / (36864 * pow(dAlpha[j], 2))) * (3072 * pow(dLambdaEps[j], 4) )* dt * dt;
				muv2[j] = (1.0 - dLambdaEps[j] * dt * (1.0 - 0.5 * dLambdaEps[j] * dt));			
				sigv[j] = dAlpha[j] * (1.0 - 0.75 * dLambdaEps[j] * dt) * sqdt;			
				sigv2[j] = dAlpha[j] / 16.0 * (4.0 * dLvlEps[j] - dAlpha[j] * dAlpha[j]) * dt * sqdt;			
				sigv3[j] = 0.25 * dAlpha[j] * dAlpha[j] * dt;

				b1[j] =( - 1 + exp( - dLambdaEps[j] * dt) + 0.5 * dLambdaEps[j] * dt * ( 1+ exp( - dLambdaEps[j] * dt)) );
				b2[j] =dAlpha[j] * dAlpha[j] / (( 1- exp( - dLambdaEps[j] * dt)) *  dLambdaEps[j] * dLambdaEps[j]);
				b3[j] = 4 * exp(-0.5 * dLambdaEps[j] * dt) / (dLambdaEps[j] * pow( - 1+ exp(-dLambdaEps[j] * dt),2));

				c2[j]= dSigma[j]  * dRho[j] / dAlpha[j];
				c3[j]= c2[j] * dLambdaEps[j];
				c1[j]= - c2[j] * dLvlEps[j] * dt;
				c4[j]= dSigma[j] * sqrt(1-dRho[j] * dRho[j]);
				c5[j]= b1[j]* b2[j];
				c6[j]= b1[j]* b3[j];
				c7[j]= - ( 2 * dLambdaEps[j] *dt * exp( - dLambdaEps[j] * dt ) - 1 + exp( - 2 * dLambdaEps[j] * dt)) / (dLambdaEps[j] * pow( 1- exp( - dLambdaEps[j] * dt),2));
				nu[j]= 2 * dLvlEps[j] / pow( dAlpha[j], 2) - 1;
				xAux[j]= 4 * dLambdaEps[j] * exp( - 0.5 * dLambdaEps[j] * dt) / ( pow( dAlpha[j], 2) * ( 1 - exp( - dLambdaEps[j] * dt )) );
				
				sizeNuMain = (long)floor(nu[j]);				
				sizeNuMain= (int)abs(sizeNuMain);

				sizeNuMain = (long)floor(nu[j]);				
				sizeNuMain= (int)abs(sizeNuMain);
					
				if(nu[j]>=0)
				ModifiedBesselFunctionValues = dmatrix(0, iNbTime-1, 0,sizeNuMain+1);
				if(nu[j]<0)
				ModifiedBesselFunctionValues = dmatrix(0, iNbTime-1, 0,sizeNuMain);				
						
						for(iTabCount = 1;	iTabCount <= NbTabCount+1; iTabCount++)
				{
			
					xValues[j][iTabCount] = xInfBound + ( (xSupBound - xInfBound) / NbTabCount) * (iTabCount-1);
					BesselArrays( nu[j] , xValues[j][iTabCount], ModifiedBesselFunctionValues[j]);					
					
					if(nu[j]>=0)
					TabulatedBesselFunction[j][iTabCount] = ModifiedBesselFunctionValues[j][sizeNuMain+1] / ModifiedBesselFunctionValues[j][sizeNuMain];

					if(nu[j]<0)
					TabulatedBesselFunction[j][iTabCount] = ModifiedBesselFunctionValues[j][sizeNuMain-1] / ModifiedBesselFunctionValues[j][sizeNuMain];


				}
						
				spline(xValues[j], TabulatedBesselFunction[j], NbTabCount + 1, 1.0e30, 1.0e30, TabulatedBesselFunctionSecondDerivate[j]);

			
				if(nu[j]>=0)
				free_dmatrix( ModifiedBesselFunctionValues, 0, iNbTime-1, 0,sizeNuMain+1);
				if(nu[j]<0)
				free_dmatrix( ModifiedBesselFunctionValues, 0, iNbTime-1, 0,sizeNuMain);


			
			}

			break;
		}

// order 3  Volatility Projection + Bessel Approximation by cubic splines

		case 7:
		{

				/* Extra allocations */
			muf		= calloc(iNbTime, sizeof(double));
			sigf2	= calloc(iNbTime, sizeof(double));
			sigf3	= calloc(iNbTime, sizeof(double));
			sigf4	= calloc(iNbTime, sizeof(double));
			sigf5	= calloc(iNbTime, sizeof(double));
			sigf6	= calloc(iNbTime, sizeof(double));
			mupsi	= calloc(iNbTime, sizeof(double));
			mupsi2	= calloc(iNbTime, sizeof(double));
			sigv2	= calloc(iNbTime, sizeof(double));
			sigv3	= calloc(iNbTime, sizeof(double));
			sigv4 = calloc(iNbTime, sizeof(double));
			sigv5 = calloc(iNbTime, sizeof(double));
			sigv6 = calloc(iNbTime, sizeof(double));
			sigv7 = calloc(iNbTime, sizeof(double));
			sigv8 = calloc(iNbTime, sizeof(double));
			sigv9 = calloc(iNbTime, sizeof(double));
			sigv10 = calloc(iNbTime, sizeof(double));
			sigv11 = calloc(iNbTime, sizeof(double));
			sigv12 = calloc(iNbTime, sizeof(double));
			sigv13 = calloc(iNbTime, sizeof(double));
			b1= calloc(iNbTime, sizeof(double));
			b2= calloc(iNbTime, sizeof(double));
			b3= calloc(iNbTime, sizeof(double));
			c2= calloc(iNbTime, sizeof(double));
			c3= calloc(iNbTime, sizeof(double));
			c1= calloc(iNbTime, sizeof(double));
			c4= calloc(iNbTime, sizeof(double));
			c5= calloc(iNbTime, sizeof(double));
			c6= calloc(iNbTime, sizeof(double));
			c7= calloc(iNbTime, sizeof(double));
			nu= calloc(iNbTime, sizeof(double));
			xAux = calloc(iNbTime, sizeof(double));

			if	(xValues)									free_dmatrix (xValues, 0, iNbTime-1, 1,NbTabCount+1);xValues=0;												      
			if	(TabulatedBesselFunction)					free_dmatrix (TabulatedBesselFunction, 0, iNbTime-1, 1,NbTabCount+1);TabulatedBesselFunction=0;						
			if	(TabulatedBesselFunctionSecondDerivate)		free_dmatrix (TabulatedBesselFunctionSecondDerivate, 0, iNbTime-1, 1,NbTabCount+1);TabulatedBesselFunctionSecondDerivate=0;	

			xValues = dmatrix(0, iNbTime-1, 1,NbTabCount+1);	
			TabulatedBesselFunction = dmatrix(0, iNbTime-1, 1,NbTabCount+1);			
			TabulatedBesselFunctionSecondDerivate = dmatrix(0, iNbTime-1, 1,NbTabCount+1);

			if (!muf || !sigf2 || !sigf3 || !sigf4 || !sigf5 || !sigf6 || !sigv2 || !sigv3  || !sigv4 || !sigv5 || !sigv6 || !sigv7 || !sigv8 || !sigv9 || !sigv10 || !sigv11 || !sigv12 || !sigv13   || !mupsi  || !mupsi2 || !b1 || !b2 || !b3 || !c1 || !c2 || !c3 || !c4 || !c5 || !c6 || !c7 || !nu || !xAux)

			{
				err = "Memory allocation failure in lgmSV_mc_balsam";
				goto FREE_RETURN;
			}		

			for (j=1; j<iNbTime; j++)
			{
				dt = dTime[j] - dTime[j-1];
				sqdt = sqrt(dt);
			
				coef1[j] = dRho[j];
				coef2[j] = sqrt(1.0 - coef1[j] * coef1[j]);

				muf[j] = -0.25 * dRho[j] * dAlpha[j] * dSigma[j] * dt;
				sigf[j] = coef2[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;

				sigf2[j] = coef2[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf3[j] = coef1[j] * dSigma[j] * (1.0 - 0.25 * dLambdaEps[j] * dt) * sqdt;
				sigf4[j] = coef1[j] * dSigma[j] * 0.25 * dt * (dLvlEps[j] - dAlpha[j] * dAlpha[j] / 4.0) * sqdt;
				sigf5[j] = 0.25 * coef1[j] * dAlpha[j] * dSigma[j] * dt;
				sigf6[j] = 0.25 * coef2[j] * dAlpha[j] * dSigma[j] * dt;
							
				mupsi[j] = 0.5 * dLvlEps[j] * dSigma[j] * dSigma[j] * dt * dt;
				mupsi2[j] = dSigma[j] * dSigma[j] * dt * (1.0 - 0.25 * dLambdaEps[j] * dt);
				sigpsi[j] = 0.5 * dAlpha[j] * dSigma[j] * dSigma[j] * dt * sqdt;
										
				muv1[j] = (dLvlEps[j] - 0.25 * dAlpha[j] * dAlpha[j] 
							-0.5 * dLvlEps[j] * dLambdaEps[j] * dt) * dt;			
				muv2[j] = (1.0 - dLambdaEps[j] * dt * (1.0 - 0.5 * dLambdaEps[j] * dt));			
				sigv[j] = dAlpha[j] * (1.0 - 0.75 * dLambdaEps[j] * dt) * sqdt;			
				sigv2[j] = dAlpha[j] / 16.0 * (4.0 * dLvlEps[j] - dAlpha[j] * dAlpha[j]) * dt * sqdt;			
				sigv3[j] = 0.25 * dAlpha[j] * dAlpha[j] * dt;

// v2: -1.5	 -0.5	0.5


				sigv6[j] = (-11 * pow(dAlpha[j] , 6) + 56 * pow(dAlpha[j] , 4)  * dLvlEps[j] - 48 * pow (dAlpha[j] * dLvlEps[j] , 2)) * dt * dt / ( dAlpha[j] *  1536 );
				sigv7[j] = (24 * pow(dAlpha[j] , 4)*  dLambdaEps[j] - 96 * pow (dAlpha[j] , 2) *  dLambdaEps[j] * dLvlEps[j]) * dt * dt / ( dAlpha[j] *  1536 ) ;
				sigv8[j] = 464 * pow (dAlpha[j]  * dLambdaEps[j] , 2) * dt * dt / ( dAlpha[j] *  1536 );				


// w1: 0 1 


				sigv4[j] = (( pow (dAlpha[j] , 4) - 4.0 * pow ( dAlpha[j] , 2 )  * dLvlEps[j] ) / 96.0 ) * dt;
				sigv5[j] = (-12.0 * pow ( dAlpha[j] , 2 ) * dLambdaEps[j]  / 96.0 ) * dt; 



// w2:	-2	-1	0	1	2


				sigv9[j] =( 1 / (36864 * pow(dAlpha[j], 2))) *  (576 * dAlpha[j] * dAlpha[j] * dLvlEps[j] * dLvlEps[j] * dLvlEps[j] + 364 * dLvlEps[j] * pow( dLvlEps[j], 6) -944 * dLvlEps[j] * dLvlEps[j] * pow( dAlpha[j], 4) -41 * pow(dAlpha[j],8) ) * dt * dt;
				sigv10[j] =( 1 / (36864 * pow(dAlpha[j], 2))) * (-576 * pow( dAlpha[j] * dLvlEps[j] ,2) * dLambdaEps[j] -132 *  dLambdaEps[j] * pow( dAlpha[j], 6) + 672 * dLvlEps[j] * dLambdaEps[j] *pow( dAlpha[j], 4)) * dt * dt;		
				sigv11[j] =( 1 / (36864 * pow(dAlpha[j], 2))) * (1728 * dLvlEps[j] * pow( dAlpha[j] * dLambdaEps[j] , 2)+ 1104 * pow( dLambdaEps[j] * dAlpha[j] * dAlpha[j],2)) * dt * dt;
				sigv12[j] = ( 1 / (36864 * pow(dAlpha[j], 2))) * (-6336 * dAlpha[j] * dAlpha[j] * pow( dLambdaEps[j], 3) -3072 * dLvlEps[j] *pow(dLambdaEps[j],3) )* dt * dt;
				sigv13[j] = ( 1 / (36864 * pow(dAlpha[j], 2))) * (3072 * pow(dLambdaEps[j], 4) )* dt * dt;
				muv2[j] = (1.0 - dLambdaEps[j] * dt * (1.0 - 0.5 * dLambdaEps[j] * dt));			
				sigv[j] = dAlpha[j] * (1.0 - 0.75 * dLambdaEps[j] * dt) * sqdt;			
				sigv2[j] = dAlpha[j] / 16.0 * (4.0 * dLvlEps[j] - dAlpha[j] * dAlpha[j]) * dt * sqdt;			
				sigv3[j] = 0.25 * dAlpha[j] * dAlpha[j] * dt;

				b1[j] =( - 1 + exp( - dLambdaEps[j] * dt) + 0.5 * dLambdaEps[j] * dt * ( 1+ exp( - dLambdaEps[j] * dt)) );
				b2[j] =dAlpha[j] * dAlpha[j] / (( 1- exp( - dLambdaEps[j] * dt)) *  dLambdaEps[j] * dLambdaEps[j]);
				b3[j] = 4 * exp(-0.5 * dLambdaEps[j] * dt) / (dLambdaEps[j] * pow( - 1+ exp(-dLambdaEps[j] * dt),2));

				c2[j]= dSigma[j]  * dRho[j] / dAlpha[j];
				c3[j]= c2[j] * dLambdaEps[j];
				c1[j]= - c2[j] * dLvlEps[j] * dt;
				c4[j]= dSigma[j] * sqrt(1-dRho[j] * dRho[j]);
				c5[j]= b1[j]* b2[j];
				c6[j]= b1[j]* b3[j];
				c7[j]= - ( 2 * dLambdaEps[j] *dt * exp( - dLambdaEps[j] * dt ) - 1 + exp( - 2 * dLambdaEps[j] * dt)) / (dLambdaEps[j] * pow( 1- exp( - dLambdaEps[j] * dt),2));
				
				nu[j]= 2 * dLvlEps[j] / pow( dAlpha[j], 2) - 1;
				xAux[j]= 4 * dLambdaEps[j] * exp( - 0.5 * dLambdaEps[j] * dt) / ( pow( dAlpha[j], 2) * ( 1 - exp( - dLambdaEps[j] * dt )) );
				
				sizeNuMain = (long)floor(nu[j]);				
				sizeNuMain= (int)abs(sizeNuMain);

				sizeNuMain = (long)floor(nu[j]);				
				sizeNuMain= (int)abs(sizeNuMain);
					
				if(nu[j]>=0)
				ModifiedBesselFunctionValues = dmatrix(0, iNbTime-1, 0,sizeNuMain+1);
				
				if(nu[j]<0)
				ModifiedBesselFunctionValues = dmatrix(0, iNbTime-1, 0,sizeNuMain);				
						
				for(iTabCount = 1;	iTabCount <= NbTabCount+1; iTabCount++)
				{
			
					xValues[j][iTabCount] = xInfBound + ( (xSupBound - xInfBound) / NbTabCount) * (iTabCount-1);
					BesselArrays( nu[j] , xValues[j][iTabCount], ModifiedBesselFunctionValues[j]);					
					
					if(nu[j]>=0)
					TabulatedBesselFunction[j][iTabCount] = ModifiedBesselFunctionValues[j][sizeNuMain+1] / ModifiedBesselFunctionValues[j][sizeNuMain];

					if(nu[j]<0)
					TabulatedBesselFunction[j][iTabCount] = ModifiedBesselFunctionValues[j][sizeNuMain-1] / ModifiedBesselFunctionValues[j][sizeNuMain];

				}
						
				spline(xValues[j], TabulatedBesselFunction[j], NbTabCount + 1, 1.0e30, 1.0e30, TabulatedBesselFunctionSecondDerivate[j]);
			
				if(nu[j]>=0)
				free_dmatrix( ModifiedBesselFunctionValues, 0, iNbTime-1, 0,sizeNuMain+1);

				if(nu[j]<0)
				free_dmatrix( ModifiedBesselFunctionValues, 0, iNbTime-1, 0,sizeNuMain);
			
			}

			break;
		}

	}

	/* fill the Brownian matrix */
	err = balsam_generation (iNumPaths, 2 * (iNbTime - 1), matrix);
	if (err)
	{ 
		goto FREE_RETURN;
	}
	time2 = clock();
	smessage ("Phase 1 -BalSam generation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);

	/* End of initialisation treatment */
	
	/* Initialisation time display */
	time2 = clock();
	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);	
	

	/* --------------------------------------------------------------------------------------------------
													Convolution		
	   -------------------------------------------------------------------------------------------------- */

	for (i=0; i<iNumPaths; i++)
	{
		/* initialisation */
		evtindex = 0;
		stop_path = 0;
		if ( init_func )
			(*init_func)();

		ft = 0;
		v = 1.0;
		sqv = 1.0;
		psi = 0;
		
		matrixi = matrix[i];		
		memset (path_payoff, 0, iNbProduct * sizeof (double));

		if (func_parm_tab[0] && EvalEvent[0])
		{			
			/* Event at time 0.0 */
			err = payoff_func(	0,
								dDate[0],
								dTime[0],
								func_parm_tab[0],
								ft,
								psi,
								v,
								iNbProduct,
								res_evt,
								&stop_path);

			if (err) goto FREE_RETURN;

			df = exp (dff_star[0] + gam_star[0] * ft + gam2_star[0] * psi);
			
			for (k=0; k<iNbProduct; k++)
			{
				path_payoff[k] += res_evt[k] / df;
			}

			if (do_optimisation)
			{				
				mceb_fill_savevalues_from_GRFN(	save_values[evtindex],
												res_evt,
												i,
												df,
												params);
			}

			evtindex++;
		}

		m = 0;

		for (j=1; stop_path == 0 && j<iNbTime; j++)
		{	
			switch (Params.iSchemeOrder)
			{
				case 0:
				{
					brow1 = matrixi[m];
					brow2 = coef1[j] * brow1 + coef2[j] * matrixi[m+1];
					
					sqv = sqrt(v);
					ft += sigf[j] * sqv * brow1;
					psi += sigpsi[j] * v;										
					v = muv2[j] + muv1[j] * v + sigv[j] * sqv * brow2;				

					if (v < 0.0)
					{
						v = 1.0E-08;
					}

					break;
				}

				case 1:
				{
					brow1 = matrixi[m];
					brow2 = coef1[j] * brow1 + coef2[j] * matrixi[m+1];
					
					sqv = sqrt(v);
					ft += sigf[j] * sqv * brow1;
					psi += sigpsi[j] * v;
					log_drift = muv2[j] + muv1[j] * v;

					if (log_drift > 0.0)
					{
						log_var = log(1.0 + sigv[j] * v / log_drift / log_drift);
						v = log_drift * exp(-0.5 * log_var + sqrt(log_var) * brow2);
					}
					else
					{
						v = 1.0E-08;
					}

					break;
				}

				case 2:
				{
					brow1 = matrixi[m];
					brow2 = matrixi[m+1];

					sqv = sqrt(v);
					
					ft += muf[j] + (sigf[j] * sqv + sigf2[j] / sqv) * brow1
						+ (sigf3[j] * sqv + sigf4[j] / sqv) * brow2
						+ (sigf5[j] * brow2 + sigf6[j] * brow1) * brow2;

					psi += mupsi[j] + mupsi2[j] * v + sigpsi[j] * sqv * brow2;

					v = muv1[j] + muv2[j] * v;
					v += (sigv[j] * sqv + sigv2[j] / sqv + sigv3[j] * brow2) * brow2;

					if (v < 0.0)
					{
						v = 0.001;
					}

					break;
				}

// case 3: Volatility Projection 

				case 3:

				{
				
					brow1 = matrixi[m];
					brow2 = matrixi[m+1];

					sqv = sqrt(v);
					
					ft += muf[j] + (sigf[j] * sqv + sigf2[j] / sqv) * brow1
						+ (sigf3[j] * sqv + sigf4[j] / sqv) * brow2
						+ (sigf5[j] * brow2 + sigf6[j] * brow1) * brow2;

					psi += mupsi[j] + mupsi2[j] * v + sigpsi[j] * sqv * brow2;
		
					v0 = v;
					v = muv1[j] + muv2[j] * v0 - (sigv4[j] / v0 + sigv5[j]  + sigv9[j] / (v0 * v0)+ sigv10[j] /v0 + sigv11[j] + sigv12[j] * v0 + sigv13[j] * v0 * v0  ) * dt;
					v += (sigv[j] * sqv + sigv2[j] / sqv + sigv3[j] * brow2 + sqdt * ( sigv6[j] / ( v0 * sqv )  + sigv7[j] / sqv + sigv8[j] *  sqv ) + ( sqdt * (sigv4[j] / v0 + sigv5[j]) + dt * ( sigv9[j]/(v0 * v0) + sigv10[j]/v0 + sigv11[j] + sigv12[j] * v0 + sigv13[j] * v0 * v0 ) ) * brow2) * brow2;


					if (v < 0.0)
					{
						v = 0.001;
					}


					break;				
				
				}				
	
//Vol Projection + Expectation of Integral Volatility as a Bessel Function	 

				case 4:
				{

					brow1 = matrixi[m];
					brow2 = matrixi[m+1];

					sqv = sqrt(v);
					
					v0 = v;
					v = muv1[j] + muv2[j] * v;
					v += (sigv[j] * sqv + sigv2[j] / sqv + sigv3[j] * brow2) * brow2;

					if (v < 0.0)
					{
						v = 0.001;
					}

					sqrtV = sqrt(v*v0);
					sumV = v + v0;
					difV = v - v0;


					x = xAux[j] * sqrtV;
					SimpleBesselFunction = (ModifiedBesselFunction(nu[j]+1,x) / ModifiedBesselFunction(nu[j],x)	+ nu[j] / x);
					BigBesselFunction = c5[j] + c6[j] * sqrtV * SimpleBesselFunction + c7[j] * sumV;
					
					ft +=c1[j] + c2[j] * difV + BigBesselFunction * c3[j] + c4[j] * sqrt(BigBesselFunction) * brow1;  					

					psi += dSigma[j] * dSigma[j] * BigBesselFunction;

					break;



				}
	


				case 5:
				{

					brow1 = matrixi[m];
					brow2 = matrixi[m+1];

					sqv = sqrt(v);
					
					v0 = v;
					v = muv1[j] + muv2[j] * v;
					v += (sigv[j] * sqv + sigv2[j] / sqv + sigv3[j] * brow2) * brow2;

					if (v < 0.0)
					{
						v = 0.001;
					}

					sqrtV = sqrt(v*v0);
					sumV = v + v0;
					difV = v - v0;
					
					besselApproxNum = 0;
					besselApproxDen = 0;

					x = xAux[j] * sqrtV;
					
					for( iBesselApprox = 0; iBesselApprox < besselApproxOrder; iBesselApprox++)
					{
					
					besselApproxNum += coefBesselApprox1[j][iBesselApprox] *  pow(x, nu[j] + 1 + 2 * iBesselApprox);
					besselApproxDen += coefBesselApprox2[j][iBesselApprox] *  pow(x, nu[j]+ 2 * iBesselApprox);
					}
					SimpleBesselFunction = (ModifiedBesselFunction(nu[j]+1,x) / ModifiedBesselFunction(nu[j],x)	+  nu[j]/x );
					SimpleBesselFunction = besselApproxNum / besselApproxDen + nu[j]/x;

					BigBesselFunction = c5[j] + c6[j] * sqrtV * SimpleBesselFunction + c7[j] * sumV;
					
					ft +=c1[j] + c2[j] * difV + BigBesselFunction * c3[j] + c4[j] * sqrt(BigBesselFunction) * brow1;  					

					psi += dSigma[j] * dSigma[j] * BigBesselFunction;

					break;



				}





// order 3  Volatility Projection + Bessel Approximation by cubic splines
				case 6:
				{

					brow1 = matrixi[m];
					brow2 = matrixi[m+1];

					sqv = sqrt(v);
					
					v0 = v;
					v = muv1[j] + muv2[j] * v;
					v += (sigv[j] * sqv + sigv2[j] / sqv + sigv3[j] * brow2) * brow2;

					if (v < 2.0e-308)
					{
						v = 0.001;
					}

					sqrtV = sqrt(v*v0);
					sumV = v + v0;
					difV = v - v0;


					x = xAux[j] * sqrtV;
		  			splint(xValues[j], TabulatedBesselFunction[j], TabulatedBesselFunctionSecondDerivate[j], NbTabCount + 1, x, &yInterpolated);
					SimpleBesselFunction = (yInterpolated	+ nu[j] / x);
					BigBesselFunction = c5[j] + c6[j] * sqrtV * SimpleBesselFunction + c7[j] * sumV;
					
					ft +=c1[j] + c2[j] * difV + BigBesselFunction * c3[j] + c4[j] * sqrt(BigBesselFunction) * brow1;  					

					psi += dSigma[j] * dSigma[j] * BigBesselFunction;

					break;



				}
			

// order 3  Volatility Projection + Bessel Approximation by cubic splines				
				case 7:
				{

					brow1 = matrixi[m];
					brow2 = matrixi[m+1];

					sqv = sqrt(v);
					
					v0 = v;
					v = muv1[j] + muv2[j] * v0 - (sigv4[j] / v0 + sigv5[j]  + sigv9[j] / (v0 * v0)+ sigv10[j] /v0 + sigv11[j] + sigv12[j] * v0 + sigv13[j] * v0 * v0  ) * dt;
					v += (sigv[j] * sqv + sigv2[j] / sqv + sigv3[j] * brow2 + sqdt * ( sigv6[j] / ( v0 * sqv )  + sigv7[j] / sqv + sigv8[j] *  sqv ) + ( sqdt * (sigv4[j] / v0 + sigv5[j]) + dt * ( sigv9[j]/(v0 * v0) + sigv10[j]/v0 + sigv11[j] + sigv12[j] * v0 + sigv13[j] * v0 * v0 ) ) * brow2) * brow2;

					if (v < 0.001)
					{
						v = 0.001;
					}

					sqrtV = sqrt(v*v0);
					sumV = v + v0;
					difV = v - v0;


					x = xAux[j] * sqrtV;
		  			splint(xValues[j], TabulatedBesselFunction[j], TabulatedBesselFunctionSecondDerivate[j], NbTabCount + 1, x, &yInterpolated);
					SimpleBesselFunction = (yInterpolated	+ nu[j] / x);
					BigBesselFunction = c5[j] + c6[j] * sqrtV * SimpleBesselFunction + c7[j] * sumV;
					
					ft +=c1[j] + c2[j] * difV + BigBesselFunction * c3[j] + c4[j] * sqrt(BigBesselFunction) * brow1;  					

					psi += dSigma[j] * dSigma[j] * BigBesselFunction;

					break;



				}


		}		
			
			/*	Case of evaluation events */
			if (EvalEvent[j])
			{
				/* Modification of the Payoff at t */			
				err = payoff_func(	i,
									dDate[j],
									dTime[j],
									func_parm_tab[j],
									ft,
									psi,
									v,
									iNbProduct,
									res_evt,
									&stop_path);
				
				if (err) goto FREE_RETURN;

				df = exp (dff_star[evtindex] + gam_star[evtindex] * ft + gam2_star[evtindex] * psi);

				for (k=0; k<iNbProduct; k++)
				{
					path_payoff[k] += res_evt[k] / df;
				}

				if (do_optimisation)
				{
					mceb_fill_savevalues_from_GRFN(	save_values[evtindex],
													res_evt,
													i,
													df,
													params);
				}

				evtindex++;
			}

			m += 2;
		}

		for (k=0; k<iNbProduct; k++)
		{
			sum_payoff[k] += path_payoff[k] / iNumPaths;
			sum_payoff2[k] += path_payoff[k] * path_payoff[k] / iNumPaths;
		}

		if (do_optimisation && params->iKnockInCol)
		{
			/* we recopy in the col pay the pv of the column */
			for (j=0; j<iNbEvent; j++)
			{
				if (optimise[j])
				{
					save_values[j][params->iNbIndex][i] = path_payoff[(int) (save_values[j][params->iNbIndex][i] + 0.5)];
				}
			}
		}
	}	

	for (k=0; k<iNbProduct; k++)
	{	
		res[k][0] = sum_payoff[k];
		res[k][1] = (sum_payoff2[k] - sum_payoff[k] * sum_payoff[k]) / iNumPaths;

		if (res[k][1] > 0.0)
		{
			res[k][1] = sqrt(res[k][1]);
		}
		else
		{
			res[k][1] = 0.0;
		}
	}	

	/* Convolution time display */
	time1 = clock();
	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (time1 - time2) / CLOCKS_PER_SEC);

	if (do_optimisation)
	{
		/* Free the big matrix of memory first */
		if (matrix) free_dmatrix (matrix, 0, iNumPaths - 1, 0, 2 * (iNbTime - 1) - 1);
		matrix = NULL;

		time1 = clock();

		if (dTime[0] < 1.0E-08 && EvalEvent[0] && optimise[0])
		{
			optimise[0] = 0;
			optim_today = 1;
		}
		else
		{
			optim_today = 0;
		}

		err = find_and_optimise_boundary(	save_values,
											iNbEvent,
											iNumPaths,
											optimise,
											params,
											&(res[iNbProduct][0]),
											&(res[iNbProduct][1]));

		if (err) goto FREE_RETURN;

		if (optim_today)
		{
			if (params->iIsKO)
			{
				if (res[iNbProduct][0] < 0.0)
				{
					res[iNbProduct][0] = 0.0;
					res[iNbProduct][1] = 0.0;
				}
			}
			else
			{
				if (save_values[0][params->iNbIndex][0] > res[iNbProduct][0])
				{
					res[iNbProduct][0] = save_values[0][params->iNbIndex][0];
					res[iNbProduct][1] = 0.0;
				}
			}

			optimise[0] = 1;
		}		

		time2 = clock();
		smessage ("Phase 3 -optimisation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	}
	
FREE_RETURN:

	if (matrix) free_dmatrix (matrix, 0, iNumPaths - 1, 0, 2 * (iNbTime - 1) - 1);
	if (res_evt) free_dvector (res_evt, 0, iNbProduct - 1);
	if (path_payoff) free_dvector (path_payoff, 0, iNbProduct - 1);
	if (sum_payoff) free_dvector (sum_payoff, 0, iNbProduct - 1);
	if (sum_payoff2) free_dvector (sum_payoff2, 0, iNbProduct - 1);

	if (sigf) free (sigf);
	if (sigpsi) free (sigpsi);
	if (muv1) free (muv1);
	if (muv2) free (muv2);
	if (sigv) free (sigv);
	if (coef1) free (coef1);
	if (coef2) free (coef2);

	if (muf) free (muf);
	if (sigf2) free (sigf2);
	if (sigf3) free (sigf3);
	if (sigf4) free (sigf4);
	if (sigf5) free (sigf5);
	
	if (sigv2) free (sigv2);
	if (sigv3) free (sigv3);
	if (sigv4) free (sigv4);
	if (sigv5) free (sigv5);	
	if (sigf6) free (sigf6);
	if (sigv7) free (sigv7);
	if (sigv8) free (sigv8);
	if (sigv9) free (sigv9);
	if (sigv10) free (sigv10);
	if (sigv11) free (sigv11);
	if (sigv12) free (sigv12);
	if (sigv13) free (sigv13);

	if (b1) free (b1);
	if (b2) free (b2);
	if (b3) free (b3);
	if (c1) free (c1);
	if (c2) free (c2);
	if (c3) free (c3);
	if (c4) free (c4);
	if (c5) free (c5);
	if (c6) free (c6);
	if (c7) free (c7);

	if (d1) free (d1);
	if (n0) free (n0);
	if (n2) free (n2);
	if (n4) free (n4);

    if (coefBesselApprox1) free_dmatrix( coefBesselApprox1, 0, iNbTime-1, 0,besselApproxOrder-1);coefBesselApprox1 = 0;    
	if (coefBesselApprox2) free_dmatrix( coefBesselApprox2, 0, iNbTime-1, 0,besselApproxOrder-1);coefBesselApprox2 = 0;

	if (cumulant1) free_dmatrix (cumulant1, 0, iNbTime-1, 0,NbCumulants - 1);cumulant1=0;
	if (cumulant2) free_dmatrix (cumulant2, 0, iNbTime-1, 0,NbCumulants - 1);cumulant2=0;
	if (cumulant3) free_dmatrix (cumulant3, 0, iNbTime-1, 0,NbCumulants - 1);cumulant3=0;	

	if	(xValues)									free_dmatrix (xValues, 0, iNbTime-1, 1,NbTabCount+1);xValues=0;												      
	if	(TabulatedBesselFunction)					free_dmatrix (TabulatedBesselFunction, 0, iNbTime-1, 1,NbTabCount+1);TabulatedBesselFunction=0;						
	if	(TabulatedBesselFunctionSecondDerivate)		free_dmatrix (TabulatedBesselFunctionSecondDerivate, 0, iNbTime-1, 1,NbTabCount+1);TabulatedBesselFunctionSecondDerivate=0;	



	if (gamma) free_dmatrix(gamma, 0, iNbTime-1, 0,NbCumulants - 1);gamma=0;

	if (nu) free (nu);
	if (xAux) free (xAux);

	if (mupsi2) free (mupsi2);
				
	if (save_values) free_f3tensor(save_values, 0, iNbEvent - 1, 0, params->iNbIndex + params->iHasNumeraire + params->iHasFees, 0, iNumPaths - 1);		

	/* Return the error message */
	return err;
}


Err	 lgmSV2F_mc_balsam(	
					/*	Time Information  */
					int			iNbTime,
					int			iNbEvent,
					double		*dTime,
					double		*dDate,
					
					int			iNumPaths,					
						
					/*	Model data Information	*/
					double		dLambdaX1,
					double		dLambdaX2,

					double		*dSigma,
					double		*dAlphaLGM,
					double		*dRhoLGM,
										
					double		*dAlpha,
					double		*dLambdaEps,
					double		*dLvlEps,
					double		*dRho,
					double		*dRho2,

					/* Parameters for DF(t,T*) reconstruction */
					double		*dff_star,
					double		*gam1_star,
					double		*gam2_star,
					double		*gam1_2_star,
					double		*gam2_2_star,
					double		*gam12_star,

					/* Parameters */
					LGMSVParam	Params,

					/*	Product data */
					void		**func_parm_tab, 
					int			*EvalEvent,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					
					/*	Initialisation function to be called at the beggining of each path or NULL if none */
					void		(*init_func)(),
					
					/*	Payoff function */
					Err (*payoff_func)(	long	path_index,
										double	evt_date,
										double	evt_time,
										void	*func_parm,

										/* Model data	*/
										double	ft1,
										double	ft2,
										double	phi1,
										double	phi2,
										double	phi12,
										double	v,
																
										/* Vector of results to be updated */
										int		nprod,
										/* Result	*/
										double	*prod_val,
										int		*stop_path),
					/*	Result */
					int			iNbProduct, 
					double		**res)
{
Err		err = NULL;
clock_t	time1, time2;
long	i, j, k, m;

double	*sum_payoff			= NULL,
		*sum_payoff2		= NULL,
		*path_payoff		= NULL,
		*res_evt			= NULL,				
		**matrix			= NULL,
		*matrixi			= NULL,

		*sigf1				= NULL,
		*sigf2				= NULL,
		*sigpsi1			= NULL,
		*sigpsi2			= NULL,
		*sigpsi12			= NULL,
		*muv1				= NULL,
		*muv2				= NULL,
		*sigv				= NULL,
		*coef1				= NULL,
		*coef2				= NULL,
		*coef3				= NULL,
		*coef4				= NULL,
		*coef5				= NULL,

		***save_values		= NULL;

int		stop_path;
int		evtindex;

double	dt, sqdt;
double	ft1, ft2, v, sqv, psi1, psi2, psi12;
double	brow1, brow2, brow3;
double	log_drift, log_var;
double	df;
int		optim_today;
long	seed = -123456789;

	/* --------------------------------------------------------------------------------------------------
													Initialisation		
	   -------------------------------------------------------------------------------------------------- */	
	
	/* For computational time calculation				 */
	time1 = clock();		

	/* odd number of paths */
	iNumPaths = (int) (iNumPaths / 2) * 2 + 1;
					
	matrix = dmatrix (0, iNumPaths - 1, 0, 3 * (iNbTime - 1) - 1);
	res_evt = dvector (0, iNbProduct - 1);
	path_payoff = dvector (0, iNbProduct - 1); 
	sum_payoff = dvector (0, iNbProduct - 1);
	sum_payoff2 = dvector (0, iNbProduct - 1);

	/* for precalculations */
	sigf1 = dvector(0, iNbTime - 1);
	sigf2 = dvector(0, iNbTime - 1);
	sigpsi1 = dvector(0, iNbTime - 1);
	sigpsi2 = dvector(0, iNbTime - 1);
	sigpsi12 = dvector(0, iNbTime - 1);
	muv1 = dvector(0, iNbTime - 1);
	muv2 = dvector(0, iNbTime - 1);	
	sigv = dvector(0, iNbTime - 1);	
	coef1 = dvector(0, iNbTime - 1);	
	coef2 = dvector(0, iNbTime - 1);
	coef3 = dvector(0, iNbTime - 1);
	coef4 = dvector(0, iNbTime - 1);
	coef5 = dvector(0, iNbTime - 1);

	if (!matrix || !path_payoff || !sum_payoff || !sum_payoff2 || !res_evt ||
		!sigf1 || !sigf2 || !sigpsi1 || !sigpsi2 || !sigpsi12 || !muv1 || !muv2 || !sigv ||
		!coef1 || !coef2 || !coef3 || !coef4 || !coef5)
	{
		err = "Memory allocation failure in lgmSV_mc_balsam";
		goto FREE_RETURN;
	}

	if (do_optimisation)
	{
		save_values = f3tensor(0, iNbEvent - 1, 0, params->iNbIndex + params->iHasNumeraire + params->iHasFees, 0, iNumPaths - 1);
		
		if (!save_values)
		{
			err = "Memory allocation (2) failure in lgmSV_mc_balsam";
			goto FREE_RETURN;
		}
	}

	memset (sum_payoff, 0, iNbProduct * sizeof (double));
	memset (sum_payoff2, 0, iNbProduct * sizeof (double));

	/* All the needed precalculations */
	for (j=1; j<iNbTime; j++)
	{
		dt = dTime[j] - dTime[j-1];
		sqdt = sqrt(dt);
		
		sigf1[j] = dSigma[j] * sqdt;		
		sigf2[j] = sigf1[j] * dAlphaLGM[j];
		sigpsi1[j] = sigf1[j] * sigf1[j];
		sigpsi2[j] = sigf2[j] * sigf2[j];
		sigpsi12[j] = dRhoLGM[j] * sigf1[j] * sigf2[j];

		switch (Params.iSchemeOrder)
		{
			case 0:
			{
				sigv[j] = dAlpha[j] * sqdt;
				muv1[j] = 1.0 - dLambdaEps[j] * dt;
				muv2[j] = dLvlEps[j] * dt;

				break;
			}

			case 1:
			{
				if (fabs(dLambdaEps[j]) > 1.0E-08)
				{
					muv1[j] = exp(-dLambdaEps[j] * dt);
					muv2[j] = dLvlEps[j] / dLambdaEps[j] * (1.0 - muv1[j]);
					sigv[j] = dAlpha[j] * sqrt(0.5 / dLambdaEps[j] * (1.0 - muv1[j] * muv1[j]));
					sigv[j] *= sigv[j];
				}
				else
				{
					muv1[j] = 1.0;
					muv2[j] = dLvlEps[j] * dt;
					sigv[j] = dAlpha[j] * sqdt;
					sigv[j] *= sigv[j];
				}

				break;
			}
		}

		coef1[j] = dRhoLGM[j];
		coef2[j] = sqrt(1.0 - coef1[j] * coef1[j]);
		coef3[j] = dRho[j];
		coef4[j] = (dRho2[j] - coef1[j] * coef3[j]) / coef2[j];
		coef5[j] = 1.0 - coef3[j] * coef3[j] - coef4[j] * coef4[j];

		if (coef5[j] < 0.0)
		{
			err = "Correlation Matrix in LGMSV 2F is not positive definite";
			goto FREE_RETURN;
		}

		coef5[j] = sqrt(coef5[j]);
	}

	/* fill the Brownian matrix */
	err = balsam_generation (iNumPaths, 3 * (iNbTime - 1), matrix);
	if (err)
	{ 
		goto FREE_RETURN;
	}
	time2 = clock();
	smessage ("Phase 1 -BalSam generation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);

	/* End of initialisation treatment */
	
	/* Initialisation time display */
	time2 = clock();
	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);	
	

	/* --------------------------------------------------------------------------------------------------
													Convolution		
	   -------------------------------------------------------------------------------------------------- */

	for (i=0; i<iNumPaths; i++)
	{
		/* initialisation */
		evtindex = 0;
		stop_path = 0;
		if ( init_func )
			(*init_func)();

		ft1 = 0;
		ft2 = 0;
		v = 1;
		psi1 = 0;
		psi2 = 0;
		psi12 = 0;
		
		matrixi = matrix[i];		
		memset (path_payoff, 0, iNbProduct * sizeof (double));

		if (func_parm_tab[0] && EvalEvent[0])
		{			
			/* Event at time 0.0 */
			err = payoff_func(	0,
								dDate[0],
								dTime[0],
								func_parm_tab[0],
								ft1,
								ft2,
								psi1,
								psi2,
								psi12,
								v,
								iNbProduct,
								res_evt,
								&stop_path);

			if (err) goto FREE_RETURN;

			df = exp (dff_star[0] + gam1_star[0] * ft1+ gam2_star[0] * ft2
				+ gam1_2_star[0] * psi1 + gam2_2_star[0] * psi2 + gam12_star[0] * psi12);
			
			for (k=0; k<iNbProduct; k++)
			{
				path_payoff[k] += res_evt[k] / df;
			}

			if (do_optimisation)
			{				
				mceb_fill_savevalues_from_GRFN(	save_values[evtindex],
													res_evt,
													i,
													df,
													params);
			}

			evtindex++;
		}

		m = 0;

		for (j=1; stop_path == 0 && j<iNbTime; j++)
		{																
			brow1 = matrixi[m];
			brow2 = coef1[j] * brow1 + coef2[j] * matrixi[m+1];
			brow3 = coef3[j] * brow1 + coef4[j] * matrixi[m+1] + coef5[j] * matrixi[m+2];
			
			sqv = sqrt(v);
			ft1 += sigf1[j] * sqv * brow1;
			ft2 += sigf2[j] * sqv * brow2;
			psi1 += sigpsi1[j] * v;
			psi2 += sigpsi2[j] * v;
			psi12 += sigpsi12[j] * v;

			switch (Params.iSchemeOrder)
			{
				case 0:
				{
					v = max(muv2[j] + muv1[j] * v + sigv[j] * sqv * brow3, 0.0);
					break;
				}

				case 1:
				{
					log_drift = muv2[j] + muv1[j] * v;

					if (log_drift > 0.0)
					{
						log_var = log(1.0 + sigv[j] * v / log_drift / log_drift);
						v = log_drift * exp(-0.5 * log_var + sqrt(log_var) * brow3);
					}
					else
					{
						v = 0.0;
					}

					break;
				}
			}

			/*	Case of evaluation events */
			if (EvalEvent[j])
			{
				/* Modification of the Payoff at t */			
				err = payoff_func(	i,
									dDate[j],
									dTime[j],
									func_parm_tab[j],
									ft1,
									ft2,
									psi1,
									psi2,
									psi12,
									v,
									iNbProduct,
									res_evt,
									&stop_path);
				
				if (err) goto FREE_RETURN;

				df = exp (dff_star[evtindex] + gam1_star[evtindex] * ft1 + gam2_star[evtindex] * ft2
						+ gam1_2_star[evtindex] * psi1 + gam2_2_star[evtindex] * psi2 + gam12_star[evtindex] * psi12);

				for (k=0; k<iNbProduct; k++)
				{
					path_payoff[k] += res_evt[k] / df;
				}

				if (do_optimisation)
				{
					mceb_fill_savevalues_from_GRFN(	save_values[evtindex],
													res_evt,
													i,
													df,
													params);
				}

				evtindex++;
			}

			m += 3;
		}

		for (k=0; k<iNbProduct; k++)
		{
			sum_payoff[k] += path_payoff[k] / iNumPaths;
			sum_payoff2[k] += path_payoff[k] * path_payoff[k] / iNumPaths;
		}

		if (do_optimisation && params->iKnockInCol)
		{
			/* we recopy in the col pay the pv of the column */
			for (j=0; j<iNbEvent; j++)
			{
				if (optimise[j])
				{
					save_values[j][params->iNbIndex][i] = path_payoff[(int) (save_values[j][params->iNbIndex][i] + 0.5)];
				}
			}
		}
	}	

	for (k=0; k<iNbProduct; k++)
	{	
		res[k][0] = sum_payoff[k];
		res[k][1] = (sum_payoff2[k] - sum_payoff[k] * sum_payoff[k]) / iNumPaths;

		if (res[k][1] > 0.0)
		{
			res[k][1] = sqrt(res[k][1]);
		}
		else
		{
			res[k][1] = 0.0;
		}
	}	

	/* Convolution time display */
	time1 = clock();
	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (time1 - time2) / CLOCKS_PER_SEC);

	if (do_optimisation)
	{
		/* Free the big matrix of memory first */
		if (matrix) free_dmatrix (matrix, 0, iNumPaths - 1, 0, 3 * (iNbTime - 1) - 1);
		matrix = NULL;

		time1 = clock();

		if (dTime[0] < 1.0E-08 && EvalEvent[0] && optimise[0])
		{
			optimise[0] = 0;
			optim_today = 1;
		}
		else
		{
			optim_today = 0;
		}

		err = find_and_optimise_boundary(	save_values,
											iNbEvent,
											iNumPaths,
											optimise,
											params,
											&(res[iNbProduct][0]),
											&(res[iNbProduct][1]));

		if (err) goto FREE_RETURN;

		if (optim_today)
		{
			if (params->iIsKO)
			{
				if (res[iNbProduct][0] < 0.0)
				{
					res[iNbProduct][0] = 0.0;
					res[iNbProduct][1] = 0.0;
				}
			}
			else
			{
				if (save_values[0][params->iNbIndex][0] > res[iNbProduct][0])
				{
					res[iNbProduct][0] = save_values[0][params->iNbIndex][0];
					res[iNbProduct][1] = 0.0;
				}
			}

			optimise[0] = 1;
		}		

		time2 = clock();
		smessage ("Phase 3 -optimisation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	}
	
FREE_RETURN:

	if (matrix) free_dmatrix (matrix, 0, iNumPaths - 1, 0, 3 * (iNbTime - 1) - 1);
	if (res_evt) free_dvector (res_evt, 0, iNbProduct - 1);
	if (path_payoff) free_dvector (path_payoff, 0, iNbProduct - 1);
	if (sum_payoff) free_dvector (sum_payoff, 0, iNbProduct - 1);
	if (sum_payoff2) free_dvector (sum_payoff2, 0, iNbProduct - 1);

	if (sigf1) free_dvector (sigf1, 0, iNbTime - 1);
	if (sigf2) free_dvector (sigf2, 0, iNbTime - 1);
	if (sigpsi1) free_dvector (sigpsi1, 0, iNbTime - 1);
	if (sigpsi2) free_dvector (sigpsi2, 0, iNbTime - 1);
	if (sigpsi12) free_dvector (sigpsi12, 0, iNbTime - 1);
	if (muv1) free_dvector (muv1, 0, iNbTime - 1);
	if (muv2) free_dvector (muv2, 0, iNbTime - 1);
	if (sigv) free_dvector (sigv, 0, iNbTime - 1);
	if (coef1) free_dvector (coef1, 0, iNbTime - 1);
	if (coef2) free_dvector (coef2, 0, iNbTime - 1);
	if (coef3) free_dvector (coef3, 0, iNbTime - 1);
	if (coef4) free_dvector (coef4, 0, iNbTime - 1);
	if (coef5) free_dvector (coef5, 0, iNbTime - 1);

	mceb_free_savevalues_for_GRFN(save_values, iNumPaths, iNbEvent, params);

	/* Return the error message */
	return err;
}

Err	 lgmSV_mc_balsam_rev(	
					/*	Time Information  */
					int			iNbTime,
					int			iNbEvent,
					double		*dTime,
					double		*dDate,
					
					int			iNumPaths,					
						
					/*	Model data Information	*/
					double		dLambdaX,
					double		*dSigma,

					double		*dAlpha,
					double		*dLambdaEps,
					double		*dLvlEps,
					double		*dRho,

					/* Parameters for DF(t,T*) reconstruction */
					double		*dff_star,
					double		*gam_star,
					double		*gam2_star,
					
					/* Parameters */
					LGMSVPARAM	Params,

					/*	Product data */
					void		**func_parm_tab, 
					int			*EvalEvent,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,

					/*	Initialisation function to be called at the beggining of each path or NULL if none */
					void		(*init_func)(),
					
					/*	Payoff function */
					Err (*payoff_func)(	long	path_index,
										double	evt_date,
										double	evt_time,
										void	*func_parm,
										
										double	ft,
										double	psi,
										double	v,
																
										/* Vector of results to be updated */
										int		nprod,
										/* Result	*/
										double	*prod_val,
										int		*stop_path),

					Err	(*payoff_adjust_function) (	long		event_index,
													long		npaths,
													long		nprod,
													double		***saved_values,
													void		*func_parm,
													MCEBPARAMS	mcebparams),

					/*	Result */
					int			iNbProduct, 
					double		**res)
{
Err		err = NULL;
clock_t	time1, time2;
long	i, j, k;
double	**payoff			= NULL,
		*res_evt			= NULL,
		*init_gauss			= NULL,
		*ft					= NULL,
		*v					= NULL,
		*psi				= NULL,
		**matrix			= NULL,
		*matrixk			= NULL,
		***save_values		= NULL,
		**save_valuesj		= NULL;

int		*stop_path			= NULL;

long	seed = -123456789;
double	sqv, brow1, brow2;
double	sigf, sigpsi, sigv, muv1, muv2;
double	log_drift, log_var;
double	coef1, coef2;
double	dt, sqdt;
int		evtindex;
double	df;
int		optim_today;
double	step, prob;
int		rand;

	/* --------------------------------------------------------------------------------------------------
													Initialisation		
	   -------------------------------------------------------------------------------------------------- */	
	
	/* For computational time calculation				 */
	time1 = clock();		

	/* odd number of paths */
	iNumPaths = (int) (iNumPaths / 2) * 2 + 1;
				
	init_gauss = dvector(0, iNumPaths - 1);
	matrix = dmatrix(0, 1, 0, iNumPaths - 1);

	ft = dvector(0, iNumPaths - 1);
	v = dvector(0, iNumPaths - 1);
	psi = dvector(0, iNumPaths - 1);

	payoff = dmatrix (0, iNumPaths - 1, 0, iNbProduct - 1);
	res_evt = dvector (0, iNbProduct - 1);

	stop_path = ivector(0, iNumPaths - 1);

	if (!payoff || !res_evt || !init_gauss || !ft || !v || !psi || !matrix || !stop_path)
	{
		err = "Memory allocation failure in lgmSV_mc_balsam";
		goto FREE_RETURN;
	}

	if (do_optimisation)
	{
		err = mceb_allocate_savevalues_for_GRFN(iNumPaths,
												iNbEvent,
												params,
												&save_values);

		if (err) goto FREE_RETURN;
	}

	/* Gauss initialisation */
	iNumPaths -= 1;
	iNumPaths /= 2;	
	step = 0.5 / (iNumPaths + 1);
	prob = step;

	/* Generation of the fractiles of the gaussian */
	for (i=0; i<iNumPaths; i++)
	{
		init_gauss[i] = inv_cumnorm_fast(prob);
		init_gauss[iNumPaths + i + 1] = -init_gauss[i];
		prob += step;
	}

	init_gauss[iNumPaths] = 0.0;	
	iNumPaths *= 2;
	iNumPaths += 1;

	/* End of initialisation treatment */
	
	/* Initialisation time display */
	time2 = clock();
	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);	

	/* --------------------------------------------------------------------------------------------------
													Convolution		
	   -------------------------------------------------------------------------------------------------- */

	/* Initialisation */
	for (i=0; i<iNumPaths-1; i++)
	{
		ft[i] = 0.0;
		psi[i] = 0.0;
		v[i] = 1.0;
		stop_path[i] = 0;
	}

	if (init_func)
	{
		init_func();
	}

	/* Event at time 0.0 */

	evtindex = 0;

	if (func_parm_tab[0] && EvalEvent[0])
	{			
		err = payoff_func(	0,
							dDate[0],
							dTime[0],
							func_parm_tab[0],
							ft[0],
							psi[0],
							v[0],
							iNbProduct,
							res_evt,
							&(stop_path[0]));

		if (err) goto FREE_RETURN;

		if (dff_star)
		{
			df = exp (dff_star[0] + gam_star[0] * ft[0] + gam2_star[0] * psi[0]);
		}
		else
		{
			df = 1.0;
		}

		for (i=0; i<iNumPaths; i++)
		{
			for (k=0; k<iNbProduct; k++)
			{
				payoff[i][k] += res_evt[k] / df;
			}

			stop_path[i] = stop_path[0];
		}

		if (do_optimisation)
		{
			for (i=0; i<iNumPaths; i++)
			{
				mceb_fill_savevalues_from_GRFN(	save_values[evtindex],
													res_evt,
													i,
													df,
													params);
			}
		}

		evtindex++;
	}

	for (j=1; j<iNbTime; j++)
	{
		/* Generate the random numbers */
		for (k=0; k<2; k++)
		{
			matrixk = matrix[k];

			for (i=0; i<iNumPaths-1; i++)
			{
				/* rand = random_int(nbPaths-1-i, &seed) + i; */
				rand = i + (int) ((iNumPaths-i) * uniform (&seed));
				matrixk[i] = init_gauss[rand];
				init_gauss[rand] = init_gauss[i];
				init_gauss[i] = matrixk[i];
			}

			matrixk[iNumPaths-1] = init_gauss[iNumPaths-1];
		}

		/* Precalculations */
		dt = dTime[j] - dTime[j-1];
		sqdt = sqrt(dt);
		
		sigf = dSigma[j] * sqdt;		
		sigpsi = sigf * sigf;

		switch (Params->iSchemeOrder)
		{
			case 0:
			{
				sigv = dAlpha[j] * sqdt;
				muv1 = 1.0 - dLambdaEps[j] * dt;
				muv2 = dLvlEps[j] * dt;
				break;
			}

			case 1:
			{
				if (fabs(dLambdaEps[j]) > 1.0E-08)
				{
					muv1 = exp(-dLambdaEps[j] * dt);
					muv2 = dLvlEps[j] / dLambdaEps[j] * (1.0 - muv1);
					sigv = dAlpha[j] * sqrt(0.5 / dLambdaEps[j] * (1.0 - muv1 * muv1));
					sigv *= sigv;
				}
				else
				{
					muv1 = 1.0;
					muv2 = dLvlEps[j] * dt;
					sigv = dAlpha[j] * sqdt;
					sigv *= sigv;
				}

				break;
			}
		}

		coef1 = dRho[j];
		coef2 = sqrt(1.0 - coef1 * coef1);

		if (do_optimisation)
		{
			save_valuesj = save_values[evtindex];
		}

		for (i=0; i<iNumPaths; i++)
		{
			if (!stop_path[i])
			{
				brow1 = matrix[0][i];
				brow2 = coef1 * brow1 + coef2 * matrix[1][i];
				
				sqv = sqrt(v[i]);
				ft[i] += sigf * sqv * brow1;
				psi[i] += sigpsi * v[i];

				switch (Params->iSchemeOrder)
				{
					case 0:
					{
						v[i] = max(muv2 + muv1 * v[i] + sigv * sqv * brow2, 0.0);
						break;
					}

					case 1:
					{
						log_drift = muv2 + muv1 * v[i];

						if (log_drift > 0.0)
						{
							log_var = log(1.0 + sigv * v[i] / log_drift / log_drift);
							v[i] = log_drift * exp(-0.5 * log_var + sqrt(log_var) * brow2);
						}
						else
						{
							v[i] = 0.0;
						}

						break;
					}
				}					

				/*	Case of evaluation events */
				if (EvalEvent[j])
				{
					/* Modification of the Payoff at t */			
					err = payoff_func(	i,
										dDate[j],
										dTime[j],
										func_parm_tab[j],
										ft[i],
										psi[i],
										v[i],
										iNbProduct,
										res_evt,
										&(stop_path[i]));
					
					if (err) goto FREE_RETURN;

					if (dff_star)
					{
						df = exp (dff_star[evtindex] + gam_star[evtindex] * ft[i] + gam2_star[evtindex] * psi[i]);

						for (k=0; k<iNbProduct; k++)
						{
							payoff[i][k] += res_evt[k] / df;
						}
					}
					else
					{
						df = 1.0;

						for (k=0; k<iNbProduct; k++)
						{
							payoff[i][k] += res_evt[k];
						}
					}

					if (do_optimisation)
					{
						mceb_fill_savevalues_from_GRFN(	save_valuesj,
													res_evt,
													i,
													df,
													params);
					}
				}
			}
		}
		
		if (EvalEvent[j] && payoff_adjust_function)
		{
			err = payoff_adjust_function(	evtindex,
											iNumPaths,
											iNbProduct,
											save_values,
											func_parm_tab[j],
											params);

			if (err) goto FREE_RETURN;
		}											

		if (EvalEvent[j])
		{
			evtindex = min(evtindex + 1, iNbEvent - 1);
		}			
	}

	for (k=0; k<iNbProduct; k++)
	{
		res[k][0] = 0.0;
		res[k][1] = 0.0;
	}

	for (i=0; i<iNumPaths; i++)
	{
		for (k=0; k<iNbProduct; k++)
		{
			res[k][0] += payoff[i][k] / iNumPaths;
			res[k][1] += payoff[i][k] * payoff[i][k] / iNumPaths;
		}
	}

	for (k=0; k<iNbProduct; k++)
	{	
		res[k][1] = (res[k][1] - res[k][0] * res[k][0]) / iNumPaths;

		if (res[k][1] > 0.0)
		{
			res[k][1] = sqrt(res[k][1]);
		}
		else
		{
			res[k][1] = 0.0;
		}
	}	

	/* Convolution time display */
	time1 = clock();
	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (time1 - time2) / CLOCKS_PER_SEC);

	if (do_optimisation)
	{
		time1 = clock();

		if (dTime[0] < 1.0E-08 && EvalEvent[0] && optimise[0])
		{
			optimise[0] = 0;
			optim_today = 1;
		}
		else
		{
			optim_today = 0;
		}

		err = find_and_optimise_boundary(	save_values,
											iNbEvent,
											iNumPaths,
											optimise,
											params,
											&(res[iNbProduct][0]),
											&(res[iNbProduct][1]));

		if (err) goto FREE_RETURN;

		if (optim_today)
		{
			if (params->iIsKO)
			{
				if (res[iNbProduct][0] < 0.0)
				{
					res[iNbProduct][0] = 0.0;
					res[iNbProduct][1] = 0.0;
				}
			}
			else
			{
				if (save_values[0][params->iNbIndex][0] > res[iNbProduct][0])
				{
					res[iNbProduct][0] = save_values[0][params->iNbIndex][0];
					res[iNbProduct][1] = 0.0;
				}
			}

			optimise[0] = 1;
		}
		
		time2 = clock();
		smessage ("Phase 3 -optimisation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	}
	
FREE_RETURN:

	if (payoff) free_dmatrix (payoff, 0, iNumPaths - 1, 0, iNbProduct - 1);

	if (res_evt) free_dvector (res_evt, 0, iNbProduct - 1);

	if (init_gauss) free_dvector(init_gauss, 0, iNumPaths - 1);

	if (ft) free_dvector(ft, 0, iNumPaths - 1);

	if (v) free_dvector(v, 0, iNumPaths - 1);

	if (psi) free_dvector(psi, 0, iNumPaths - 1);

	if (matrix) free_dmatrix (matrix, 0, 1, 0, iNumPaths - 1);

	if (stop_path) free_ivector(stop_path, 0, iNumPaths - 1);

	mceb_free_savevalues_for_GRFN(save_values, iNumPaths, iNbEvent, params);

	/* Return the error message */
	return err;
}

Err	 lgmSV2F_mc_balsam_rev(	
					/*	Time Information  */
					int			iNbTime,
					int			iNbEvent,
					double		*dTime,
					double		*dDate,
					
					int			iNumPaths,					
						
					/*	Model data Information	*/
					double		dLambdaX1,
					double		dLambdaX2,

					double		*dSigma,
					double		*dAlphaLGM,
					double		*dRhoLGM,
										
					double		*dAlpha,
					double		*dLambdaEps,
					double		*dLvlEps,
					double		*dRho,
					double		*dRho2,

					/* Parameters for DF(t,T*) reconstruction */
					double		*dff_star,
					double		*gam1_star,
					double		*gam2_star,
					double		*gam1_2_star,
					double		*gam2_2_star,
					double		*gam12_star,

					/* Parameters */
					LGMSVPARAM	Params,

					/*	Product data */
					void		**func_parm_tab, 
					int			*EvalEvent,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					
					/*	Initialisation function to be called at the beggining of each path or NULL if none */
					void		(*init_func)(),
					
					/*	Payoff function */
					Err (*payoff_func)(	long	path_index,
										double	evt_date,
										double	evt_time,
										void	*func_parm,

										/* Model data	*/
										double	ft1,
										double	ft2,
										double	phi1,
										double	phi2,
										double	phi12,
										double	v,
																
										/* Vector of results to be updated */
										int		nprod,
										/* Result	*/
										double	*prod_val,
										int		*stop_path),

					Err	(*payoff_adjust_function) (	long		event_index,
													long		npaths,
													long		nprod,
													double		***saved_values,
													void		*func_parm,
													MCEBPARAMS	mcebparams),
													
					/*	Result */
					int			iNbProduct, 
					double		**res)
{
Err		err = NULL;
clock_t	time1, time2;
long	i, j, k;
double	**payoff		= NULL,
		*res_evt		= NULL,
		*init_gauss		= NULL,
		*ft1			= NULL,
		*ft2			= NULL,
		*v				= NULL,
		*psi1			= NULL,
		*psi2			= NULL,
		*psi12			= NULL,
		**matrix		= NULL,
		*matrixk		= NULL,
		***save_values	= NULL,
		**save_valuesj	= NULL;

int		*stop_path	= NULL;

long	seed = -123456789;
double	dt, sqdt;
double	sqv, brow1, brow2, brow3;
double	sigf1, sigf2, sigpsi1, sigpsi2, sigpsi12, sigv, muv1, muv2;
double	coef1, coef2, coef3, coef4, coef5;
double	log_drift, log_var;
int		optim_today;
double	step, prob;
int		rand;
int		evtindex;
double	df;

	/* --------------------------------------------------------------------------------------------------
													Initialisation		
	   -------------------------------------------------------------------------------------------------- */	
	
	/* For computational time calculation				 */
	time1 = clock();		

	/* odd number of paths */
	iNumPaths = (int) (iNumPaths / 2) * 2 + 1;
				
	init_gauss = dvector(0, iNumPaths - 1);
	matrix = dmatrix(0, 2, 0, iNumPaths - 1);

	ft1 = dvector(0, iNumPaths - 1);
	ft2 = dvector(0, iNumPaths - 1);
	v = dvector(0, iNumPaths - 1);
	psi1 = dvector(0, iNumPaths - 1);
	psi2 = dvector(0, iNumPaths - 1);
	psi12 = dvector(0, iNumPaths - 1);

	payoff = dmatrix (0, iNumPaths - 1, 0, iNbProduct - 1);
	res_evt = dvector (0, iNbProduct - 1);

	stop_path = ivector(0, iNumPaths - 1);	

	if (!payoff || !res_evt || !init_gauss || !ft1 || !ft2 || !v || !psi1 || !psi2 || !psi12 || !matrix || !stop_path)
	{
		err = "Memory allocation failure in lgmSV_mc_balsam";
		goto FREE_RETURN;
	}

	if (do_optimisation)
	{
		err = mceb_allocate_savevalues_for_GRFN(iNumPaths,
												iNbEvent,
												params,
												&save_values);

		if (err) goto FREE_RETURN;
	}	

	/* Gauss initialisation */
	iNumPaths -= 1;
	iNumPaths /= 2;	
	step = 0.5 / (iNumPaths + 1);
	prob = step;

	/* Generation of the fractiles of the gaussian */
	for (i=0; i<iNumPaths; i++)
	{
		init_gauss[i] = inv_cumnorm_fast(prob);
		init_gauss[iNumPaths + i + 1] = -init_gauss[i];
		prob += step;
	}

	init_gauss[iNumPaths] = 0.0;	
	iNumPaths *= 2;
	iNumPaths += 1;

	/* End of initialisation treatment */
	
	/* Initialisation time display */
	time2 = clock();
	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);	

	/* --------------------------------------------------------------------------------------------------
													Convolution		
	   -------------------------------------------------------------------------------------------------- */

	/* Initialisation */
	for (i=0; i<iNumPaths; i++)
	{
		ft1[i] = 0.0;
		ft2[i] = 0.0;
		psi1[i] = 0.0;
		psi2[i] = 0.0;
		psi12[i] = 0.0;
		v[i] = 1.0;
		stop_path[i] = 0;
	}

	if (init_func)
	{
		init_func();
	}

	evtindex = 0;

	/* Event at time 0.0 */
	if (func_parm_tab[0] && EvalEvent[0])
	{			
		err = payoff_func(	0,
							dDate[0],
							dTime[0],
							func_parm_tab[0],
							ft1[0],
							ft2[0],
							psi1[0],
							psi2[0],
							psi12[0],
							v[0],
							iNbProduct,
							res_evt,
							&(stop_path[0]));

		if (err) goto FREE_RETURN;

		if (dff_star)
		{
			df = exp (dff_star[0] + gam1_star[0] * ft1[0] + gam2_star[0] * ft2[0]
					+ gam1_2_star[0] * psi1[0] + gam2_2_star[0] * psi2[0] + gam12_star[0] * psi12[0]);
		}
		else
		{
			df = 1.0;
		}
	
		for (i=0; i<iNumPaths; i++)
		{
			for (k=0; k<iNbProduct; k++)
			{
				payoff[i][k] = res_evt[k] / df;
			}

			stop_path[i] = stop_path[0];
		}

		if (do_optimisation)
		{
			for (i=0; i<iNumPaths; i++)
			{
				mceb_fill_savevalues_from_GRFN(	save_values[evtindex],
													res_evt,
													i,
													df,
													params);
			}
		}

		if (payoff_adjust_function)
		{
			err = payoff_adjust_function(	0,
											iNumPaths,
											iNbProduct,
											save_values,
											func_parm_tab[0],
											params);

			if (err) goto FREE_RETURN;
		}

		evtindex++;
	}

	for (j=1; j<iNbTime; j++)
	{
		/* Generate the random numbers */
		for (k=0; k<3; k++)
		{
			matrixk = matrix[k];

			for (i=0; i<iNumPaths-1; i++)
			{
				/* rand = random_int(nbPaths-1-i, &seed) + i; */
				rand = i + (int) ((iNumPaths-i) * uniform (&seed));
				matrixk[i] = init_gauss[rand];
				init_gauss[rand] = init_gauss[i];
				init_gauss[i] = matrixk[i];
			}

			matrixk[iNumPaths-1] = init_gauss[iNumPaths-1];
		}

		/* Precalculations */
		dt = dTime[j] - dTime[j-1];
		sqdt = sqrt(dt);
		
		sigf1 = dSigma[j] * sqdt;
		sigf2 = sigf1 * dAlphaLGM[j];
		sigpsi1 = sigf1 * sigf1;
		sigpsi2 = sigf2 * sigf2;
		sigpsi12 = dRhoLGM[j] * sigf1 * sigf2;

		switch (Params->iSchemeOrder)
		{
			case 0:
			{
				sigv = dAlpha[j] * sqdt;
				muv1 = 1.0 - dLambdaEps[j] * dt;
				muv2 = dLvlEps[j] * dt;
				break;
			}

			case 1:
			{
				if (fabs(dLambdaEps[j]) > 1.0E-08)
				{
					muv1 = exp(-dLambdaEps[j] * dt);
					muv2 = dLvlEps[j] / dLambdaEps[j] * (1.0 - muv1);
					sigv = dAlpha[j] * sqrt(0.5 / dLambdaEps[j] * (1.0 - muv1 * muv1));
					sigv *= sigv;
				}
				else
				{
					muv1 = 1.0;
					muv2 = dLvlEps[j] * dt;
					sigv = dAlpha[j] * sqdt;
					sigv *= sigv;
				}

				break;
			}
		}

		coef1 = dRhoLGM[j];
		coef2 = sqrt(1.0 - coef1 * coef1);
		coef3 = dRho[j];
		coef4 = (dRho2[j] - coef1 * coef3) / coef2;
		coef5 = 1.0 - coef3 * coef3 - coef4 * coef4;
		
		if (coef5 < 0.0)
		{
			err = "Correlation Matrix in LGMSV 2F is not positive definite";
			goto FREE_RETURN;
		}

		coef5 = sqrt(coef5);

		if (do_optimisation)
		{
			save_valuesj = save_values[evtindex];
		}

		for (i=0; i<iNumPaths; i++)
		{
			if (!stop_path[i])
			{
				brow1 = matrix[0][i];
				brow2 = coef1 * brow1 + coef2 * matrix[1][i];
				brow3 = coef3 * brow1 + coef4 * matrix[1][i] + coef5 * matrix[2][i];
				
				sqv = sqrt(v[i]);
				ft1[i] += sigf1 * sqv * brow1;
				ft2[i] += sigf2 * sqv * brow2;
				
				psi1[i] += sigpsi1 * v[i];
				psi2[i] += sigpsi2 * v[i];
				psi12[i] += sigpsi12 * v[i];

				switch (Params->iSchemeOrder)
				{
					case 0:
					{
						v[i] = max(muv2 + muv1 * v[i] + sigv * sqv * brow3, 0.0);
						break;
					}

					case 1:
					{
						log_drift = muv2 + muv1 * v[i];

						if (log_drift > 0.0)
						{
							log_var = log(1.0 + sigv * v[i] / log_drift / log_drift);
							v[i] = log_drift * exp(-0.5 * log_var + sqrt(log_var) * brow3);
						}
						else
						{
							v[i] = 0.0;
						}

						break;
					}
				}				

				/*	Case of evaluation events */
				if (EvalEvent[j])
				{
					/* Modification of the Payoff at t */			
					err = payoff_func(	i,
										dDate[j],
										dTime[j],
										func_parm_tab[j],
										ft1[i],
										ft2[i],
										psi1[i],
										psi2[i],
										psi12[i],
										v[i],
										iNbProduct,
										res_evt,
										&(stop_path[i]));
					
					if (err) goto FREE_RETURN;

					if (dff_star)
					{
						df = exp (dff_star[evtindex] + gam1_star[evtindex] * ft1[i] + gam2_star[evtindex] * ft2[i]
							+ gam1_2_star[evtindex] * psi1[i] + gam2_2_star[evtindex] * psi2[i] + gam12_star[evtindex] * psi12[i]);
					}
					else
					{
						df = 1.0;
					}
					
					for (k=0; k<iNbProduct; k++)
					{
						payoff[i][k] += res_evt[k] / df;
					}
				}

				if (do_optimisation)
				{
					mceb_fill_savevalues_from_GRFN(	save_valuesj,
													res_evt,
													i,
													df,
													params);
				}
			}
		}

		if (EvalEvent[j] && payoff_adjust_function)
		{
			err = payoff_adjust_function(	evtindex,
											iNumPaths,
											iNbProduct,
											save_values,
											func_parm_tab[j],
											params);

			if (err) goto FREE_RETURN;
		}											
		
		if (EvalEvent[j])
		{
			evtindex = min(evtindex + 1, iNbEvent - 1);
		}
	}

	for (k=0; k<iNbProduct; k++)
	{
		res[k][0] = 0.0;
		res[k][1] = 0.0;
	}

	for (i=0; i<iNumPaths; i++)
	{
		for (k=0; k<iNbProduct; k++)
		{
			res[k][0] += payoff[i][k] / iNumPaths;
			res[k][1] += payoff[i][k] * payoff[i][k] / iNumPaths;
		}
	}

	for (k=0; k<iNbProduct; k++)
	{		
		res[k][1] = (res[k][1] - res[k][0] * res[k][0]) / iNumPaths;

		if (res[k][1] > 0.0)
		{
			res[k][1] = sqrt(res[k][1]);
		}
		else
		{
			res[k][1] = 0.0;
		}
	}	

	/* Convolution time display */
	time1 = clock();
	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (time1 - time2) / CLOCKS_PER_SEC);

	if (do_optimisation)
	{
		time1 = clock();

		if (dTime[0] < 1.0E-08 && EvalEvent[0] && optimise[0])
		{
			optimise[0] = 0;
			optim_today = 1;
		}
		else
		{
			optim_today = 0;
		}
		
		err = find_and_optimise_boundary(	save_values,
											iNbEvent,
											iNumPaths,
											optimise,
											params,
											&(res[iNbProduct][0]),
											&(res[iNbProduct][1]));

		if (err) goto FREE_RETURN;

		if (optim_today)
		{
			if (params->iIsKO)
			{
				if (res[iNbProduct][0] < 0.0)
				{
					res[iNbProduct][0] = 0.0;
					res[iNbProduct][1] = 0.0;
				}
			}
			else
			{
				if (save_values[0][params->iNbIndex][0] > res[iNbProduct][0])
				{
					res[iNbProduct][0] = save_values[0][params->iNbIndex][0];
					res[iNbProduct][1] = 0.0;
				}
			}

			optimise[0] = 1;
		}
		
		time2 = clock();
		smessage ("Phase 3 -optimisation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	}
	
FREE_RETURN:

	if (payoff) free_dmatrix (payoff, 0, iNumPaths - 1, 0, iNbProduct - 1);

	if (res_evt) free_dvector (res_evt, 0, iNbProduct - 1);

	if (init_gauss) free_dvector(init_gauss, 0, iNumPaths - 1);

	if (ft1) free_dvector(ft1, 0, iNumPaths - 1);

	if (ft2) free_dvector(ft2, 0, iNumPaths - 1);

	if (v) free_dvector(v, 0, iNumPaths - 1);

	if (psi1) free_dvector(psi1, 0, iNumPaths - 1);
	
	if (psi2) free_dvector(psi2, 0, iNumPaths - 1);

	if (psi12) free_dvector(psi12, 0, iNumPaths - 1);

	if (matrix) free_dmatrix (matrix, 0, 2, 0, iNumPaths - 1);

	if (stop_path) free_ivector(stop_path, 0, iNumPaths - 1);

	mceb_free_savevalues_for_GRFN(save_values, iNumPaths, iNbEvent, params);
	
	/* Return the error message */
	return err;
}

Err	 lgmSV2F_mc_balsam_optim_mem(	
					/*	Time Information  */
					int			iNbTime,
					int			iNbEvent,
					double		*dTime,
					double		*dDate,
					
					int			iNumPaths,					
						
					/*	Model data Information	*/
					double		dLambdaX1,
					double		dLambdaX2,

					double		*dSigma,
					double		*dAlphaLGM,
					double		*dRhoLGM,
										
					double		*dAlpha,
					double		*dLambdaEps,
					double		*dLvlEps,
					double		*dRho,
					double		*dRho2,

					/* Parameters for DF(t,T*) reconstruction */
					double		*dff_star,
					double		*gam1_star,
					double		*gam2_star,
					double		*gam1_2_star,
					double		*gam2_2_star,
					double		*gam12_star,

					/* Parameters */
					LGMSVPARAM	Params,

					/*	Product data */
					void		**func_parm_tab, 
					int			*EvalEvent,

					/* for Optimisation of exercise boundary */
					int			do_optimisation,
					int			*optimise,
					MCEBPARAMS	params,
					
					/*	Initialisation function to be called at the beggining of each path or NULL if none */
					void		(*init_func)(),
					
					/*	Payoff function */
					Err (*payoff_func)(	long	path_index,
										double	evt_date,
										double	evt_time,
										void	*func_parm,

										/* Model data	*/
										double	ft1,
										double	ft2,
										double	phi1,
										double	phi2,
										double	phi12,
										double	v,
																
										/* Vector of results to be updated */
										int		nprod,
										/* Result	*/
										double	*prod_val,
										int		*stop_path),
					/*	Result */
					int			iNbProduct, 
					double		**res)
{
Err		err = NULL;
clock_t	time1, time2;
long	i, j, k;
double	*sum_payoff		= NULL,
		*sum_payoff2	= NULL,
		*res_evt		= NULL,
		*path_payoff	= NULL,
		*init_gauss		= NULL,
		**ft1			= NULL,
		**ft2			= NULL,
		**v				= NULL,
		**psi1			= NULL,
		**psi2			= NULL,
		**psi12			= NULL,
		**matrix		= NULL,
		*matrixk		= NULL,
		***save_values	= NULL,
		**save_valuesj	= NULL;

int		*index_event	= NULL;

int		stop_path;
int		optim_today;
long	seed = -123456789;
double	dt, sqdt;
double	var, sqv, brow1, brow2, brow3;
double	sigf1, sigf2, sigpsi1, sigpsi2, sigpsi12, sigv, muv1, muv2;
double	coef1, coef2, coef3, coef4, coef5;
double	log_drift, log_var;

double	step, prob;
int		rand;
int		evtindex;
double	df;

	/* --------------------------------------------------------------------------------------------------
													Initialisation		
	   -------------------------------------------------------------------------------------------------- */	
	
	/* For computational time calculation				 */
	time1 = clock();		

	/* odd number of paths */
	iNumPaths = (int) (iNumPaths / 2) * 2 + 1;
				
	init_gauss = calloc(iNumPaths, sizeof(double));
	matrix = dmatrix(0, 2, 0, iNumPaths - 1);

	ft1 = dmatrix(0, iNbEvent - 1, 0, iNumPaths - 1);
	ft2 = dmatrix(0, iNbEvent - 1, 0, iNumPaths - 1);
	v = dmatrix(0, iNbEvent - 1, 0, iNumPaths - 1);
	psi1 = dmatrix(0, iNbEvent - 1, 0, iNumPaths - 1);
	psi2 = dmatrix(0, iNbEvent - 1, 0, iNumPaths - 1);
	psi12 = dmatrix(0, iNbEvent - 1, 0, iNumPaths - 1);
	
	res_evt = calloc(iNbProduct, sizeof(double));
	path_payoff = calloc(iNbProduct, sizeof(double));	
	sum_payoff = calloc(iNbProduct, sizeof(double));	
	sum_payoff2 = calloc(iNbProduct, sizeof(double));	
		
	index_event = calloc(iNbEvent, sizeof(int));

	if (!init_gauss || !matrix || !sum_payoff | !sum_payoff2 || !res_evt || !path_payoff || 
		!ft1 || !ft2 || !v || !psi1 || !psi2 || !psi12 || !index_event)
	{
		err = "Memory allocation failure in lgmSV2F_mc_balsam_optim_mem";
		goto FREE_RETURN;
	}

	if (do_optimisation)
	{
		save_values = f3tensor(0, iNbEvent - 1, 0, params->iNbIndex + params->iHasNumeraire + params->iHasFees, 0, iNumPaths - 1);
		
		if (!save_values)
		{
			err = "Memory allocation (2) failure in lgmSV2F_mc_balsam_optim_mem";
			goto FREE_RETURN;
		}
	}

	/* --------------------------------------------------------------------------------------------------
													Path Generation	
	   -------------------------------------------------------------------------------------------------- */

	/* Gauss initialisation */
	iNumPaths -= 1;
	iNumPaths /= 2;	
	step = 0.5 / (iNumPaths + 1);
	prob = step;

	/* Generation of the fractiles of the gaussian */
	for (i=0; i<iNumPaths; i++)
	{
		init_gauss[i] = inv_cumnorm_fast(prob);
		init_gauss[iNumPaths + i + 1] = -init_gauss[i];
		prob += step;
	}

	init_gauss[iNumPaths] = 0.0;	
	iNumPaths *= 2;
	iNumPaths += 1;	

	/* Initialisation */
	evtindex = 0;

	for (i=0; i<iNumPaths; i++)
	{
		ft1[evtindex][i] = 0.0;
		ft2[evtindex][i] = 0.0;
		psi1[evtindex][i] = 0.0;
		psi2[evtindex][i] = 0.0;
		psi12[evtindex][i] = 0.0;
		v[evtindex][i] = 1.0;
	}	

	if (func_parm_tab[0] && EvalEvent[0])
	{
		index_event[evtindex] = 0;

		if (evtindex < iNbEvent-1)
		{
			/* Update next step */
			memcpy(ft1[evtindex+1], ft1[evtindex], iNumPaths * sizeof(double));
			memcpy(ft2[evtindex+1], ft2[evtindex], iNumPaths * sizeof(double));
			memcpy(psi1[evtindex+1], psi1[evtindex], iNumPaths * sizeof(double));
			memcpy(psi2[evtindex+1], psi2[evtindex], iNumPaths * sizeof(double));
			memcpy(psi12[evtindex+1], psi12[evtindex], iNumPaths * sizeof(double));
			memcpy(v[evtindex+1], v[evtindex], iNumPaths * sizeof(double));					
		}

		evtindex++;
	}
	
	/* Generate all the paths */
	for (j=1; j<iNbTime; j++)
	{
		/* Generate the random numbers */
		for (k=0; k<3; k++)
		{
			matrixk = matrix[k];

			for (i=0; i<iNumPaths-1; i++)
			{
				/* rand = random_int(nbPaths-1-i, &seed) + i; */
				rand = i + (int) ((iNumPaths-i) * uniform (&seed));
				matrixk[i] = init_gauss[rand];
				init_gauss[rand] = init_gauss[i];
				init_gauss[i] = matrixk[i];
			}

			matrixk[iNumPaths-1] = init_gauss[iNumPaths-1];
		}

		/* Precalculations */
		dt = dTime[j] - dTime[j-1];
		sqdt = sqrt(dt);
		
		sigf1 = dSigma[j] * sqdt;
		sigf2 = sigf1 * dAlphaLGM[j];
		sigpsi1 = sigf1 * sigf1;
		sigpsi2 = sigf2 * sigf2;
		sigpsi12 = dRhoLGM[j] * sigf1 * sigf2;

		switch (Params->iSchemeOrder)
		{
			case 0:
			{
				sigv = dAlpha[j] * sqdt;
				muv1 = 1.0 - dLambdaEps[j] * dt;
				muv2 = dLvlEps[j] * dt;
				break;
			}

			case 1:
			{
				if (fabs(dLambdaEps[j]) > 1.0E-08)
				{
					muv1 = exp(-dLambdaEps[j] * dt);
					muv2 = dLvlEps[j] / dLambdaEps[j] * (1.0 - muv1);
					sigv = dAlpha[j] * sqrt(0.5 / dLambdaEps[j] * (1.0 - muv1 * muv1));
					sigv *= sigv;
				}
				else
				{
					muv1 = 1.0;
					muv2 = dLvlEps[j] * dt;
					sigv = dAlpha[j] * sqdt;
					sigv *= sigv;
				}

				break;
			}
		}

		coef1 = dRhoLGM[j];
		coef2 = sqrt(1.0 - coef1 * coef1);
		coef3 = dRho[j];
		coef4 = (dRho2[j] - coef1 * coef3) / coef2;
		coef5 = 1.0 - coef3 * coef3 - coef4 * coef4;
		
		if (coef5 < 0.0)
		{
			err = "Correlation Matrix in LGMSV 2F is not positive definite";
			goto FREE_RETURN;
		}

		coef5 = sqrt(coef5);	

		for (i=0; i<iNumPaths; i++)
		{			
			brow1 = matrix[0][i];
			brow2 = coef1 * brow1 + coef2 * matrix[1][i];
			brow3 = coef3 * brow1 + coef4 * matrix[1][i] + coef5 * matrix[2][i];

			var = v[evtindex][i];
			sqv = sqrt(var);
			
			ft1[evtindex][i] += sigf1 * sqv * brow1;
			ft2[evtindex][i] += sigf2 * sqv * brow2;
			
			psi1[evtindex][i] += sigpsi1 * var;
			psi2[evtindex][i] += sigpsi2 * var;
			psi12[evtindex][i] += sigpsi12 * var;

			switch (Params->iSchemeOrder)
			{
				case 0:
				{
					v[evtindex][i] = max(muv2 + muv1 * var + sigv * sqv * brow3, 0.0);
					break;
				}

				case 1:
				{
					log_drift = muv2 + muv1 * var;

					if (log_drift > 0.0)
					{
						log_var = log(1.0 + sigv * var / log_drift / log_drift);
						v[evtindex][i] = log_drift * exp(-0.5 * log_var + sqrt(log_var) * brow3);
					}
					else
					{
						v[evtindex][i] = 0.0;
					}

					break;
				}
			}			
		}

		/*	Case of evaluation events */
		if (EvalEvent[j])
		{
			index_event[evtindex] = j;

			if (evtindex < iNbEvent-1)
			{
				/* Update next step */
				memcpy(ft1[evtindex+1], ft1[evtindex], iNumPaths * sizeof(double));
				memcpy(ft2[evtindex+1], ft2[evtindex], iNumPaths * sizeof(double));
				memcpy(psi1[evtindex+1], psi1[evtindex], iNumPaths * sizeof(double));
				memcpy(psi2[evtindex+1], psi2[evtindex], iNumPaths * sizeof(double));
				memcpy(psi12[evtindex+1], psi12[evtindex], iNumPaths * sizeof(double));
				memcpy(v[evtindex+1], v[evtindex], iNumPaths * sizeof(double));				
			}

			evtindex = min(evtindex + 1, iNbEvent - 1);
		}
	}

	/* End of initialisation treatment */
	
	/* Initialisation time display */
	time2 = clock();
	smessage ("Phase 1 -Path Generation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);		

	/* --------------------------------------------------------------------------------------------------
													Convolution		
	   -------------------------------------------------------------------------------------------------- */

	/* Launch the MC */
	for (i=0; i<iNumPaths; i++)
	{
		/* initialisation */
		evtindex = 0;
		stop_path = 0;
		memset (path_payoff, 0, iNbProduct * sizeof (double));

		if ( init_func )
		{
			(*init_func)();
		}					

		if (func_parm_tab[0] && EvalEvent[0])
		{			
			/* Event at time 0.0 */
			err = payoff_func(	i,
								dDate[index_event[evtindex]],
								dTime[index_event[evtindex]],
								func_parm_tab[index_event[evtindex]],
								ft1[evtindex][i],
								ft2[evtindex][i],
								psi1[evtindex][i],
								psi2[evtindex][i],
								psi12[evtindex][i],
								v[evtindex][i],
								iNbProduct,
								res_evt,
								&stop_path);

			if (err) goto FREE_RETURN;

			df = exp (dff_star[0] 
					+ gam1_star[0] * ft1[evtindex][i]
					+ gam2_star[0] * ft2[evtindex][i]
					+ gam1_2_star[0] * psi1[evtindex][i]
					+ gam2_2_star[0] * psi2[evtindex][i] 
					+ gam12_star[0] * psi12[evtindex][i]);
			
			for (k=0; k<iNbProduct; k++)
			{
				path_payoff[k] += res_evt[k] / df;
			}

			if (do_optimisation)
			{				
				mceb_fill_savevalues_from_GRFN(	save_values[evtindex],
												res_evt,
												i,
												df,
												params);
			}

			evtindex++;
		}

		for (evtindex; stop_path == 0 && evtindex<iNbEvent; evtindex++)
		{																									
			/* Modification of the Payoff at t */			
			err = payoff_func(	i,
								dDate[index_event[evtindex]],
								dTime[index_event[evtindex]],
								func_parm_tab[index_event[evtindex]],
								ft1[evtindex][i],
								ft2[evtindex][i],
								psi1[evtindex][i],
								psi2[evtindex][i],
								psi12[evtindex][i],
								v[evtindex][i],
								iNbProduct,
								res_evt,
								&stop_path);
			
			if (err) goto FREE_RETURN;

			df = exp (dff_star[evtindex]
					+ gam1_star[evtindex] * ft1[evtindex][i]
					+ gam2_star[evtindex] * ft2[evtindex][i]
					+ gam1_2_star[evtindex] * psi1[evtindex][i]
					+ gam2_2_star[evtindex] * psi2[evtindex][i] 
					+ gam12_star[evtindex] * psi12[evtindex][i]);

			for (k=0; k<iNbProduct; k++)
			{
				path_payoff[k] += res_evt[k] / df;
			}

			if (do_optimisation)
			{
				mceb_fill_savevalues_from_GRFN(	save_values[evtindex],
													res_evt,
													i,
													df,
													params);
			}
		}

		for (k=0; k<iNbProduct; k++)
		{
			sum_payoff[k] += path_payoff[k] / iNumPaths;
			sum_payoff2[k] += path_payoff[k] * path_payoff[k] / iNumPaths;
		}

		if (do_optimisation && params->iKnockInCol)
		{
			/* we recopy in the col pay the pv of the column */
			for (j=0; j<iNbEvent; j++)
			{
				if (optimise[j])
				{
					save_values[j][params->iNbIndex][i] = path_payoff[(int) (save_values[j][params->iNbIndex][i] + 0.5)];
				}
			}
		}
	}

	for (k=0; k<iNbProduct; k++)
	{	
		res[k][0] = sum_payoff[k];
		res[k][1] = (sum_payoff2[k] - sum_payoff[k] * sum_payoff[k]) / iNumPaths;

		if (res[k][1] > 0.0)
		{
			res[k][1] = sqrt(res[k][1]);
		}
		else
		{
			res[k][1] = 0.0;
		}
	}

	/* Convolution time display */
	time1 = clock();
	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (time1 - time2) / CLOCKS_PER_SEC);

	if (do_optimisation)
	{
		time1 = clock();

		if (dTime[0] < 1.0E-08 && EvalEvent[0] && optimise[0])
		{
			optimise[0] = 0;
			optim_today = 1;
		}
		else
		{
			optim_today = 0;
		}

		err = find_and_optimise_boundary(	save_values,
											iNbEvent,
											iNumPaths,
											optimise,
											params,
											&(res[iNbProduct][0]),
											&(res[iNbProduct][1]));

		if (err) goto FREE_RETURN;

		if (optim_today)
		{
			if (params->iIsKO)
			{
				if (res[iNbProduct][0] < 0.0)
				{
					res[iNbProduct][0] = 0.0;
					res[iNbProduct][1] = 0.0;
				}
			}
			else
			{
				if (save_values[0][params->iNbIndex][0] > res[iNbProduct][0])
				{
					res[iNbProduct][0] = save_values[0][params->iNbIndex][0];
					res[iNbProduct][1] = 0.0;
				}
			}

			optimise[0] = 1;
		}
		
		time2 = clock();
		smessage ("Phase 3 -optimisation, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);
	}
	
FREE_RETURN:

	if (init_gauss) free(init_gauss);	

	if (matrix) free_dmatrix (matrix, 0, 2, 0, iNumPaths - 1);

	
	if (ft1) free_dmatrix(ft1, 0, iNbEvent - 1, 0, iNumPaths - 1);

	if (ft2) free_dmatrix(ft2, 0, iNbEvent - 1, 0, iNumPaths - 1);

	if (v) free_dmatrix(v, 0, iNbEvent - 1, 0, iNumPaths - 1);

	if (psi1) free_dmatrix(psi1, 0, iNbEvent - 1, 0, iNumPaths - 1);
	
	if (psi2) free_dmatrix(psi2, 0, iNbEvent - 1, 0, iNumPaths - 1);

	if (psi12) free_dmatrix(psi12, 0, iNbEvent - 1, 0, iNumPaths - 1);
		

	if (res_evt) free(res_evt);

	if (path_payoff) free (path_payoff);

	if (sum_payoff) free(sum_payoff);

	if (sum_payoff2) free(sum_payoff2);

	if (index_event) free(index_event);

	mceb_free_savevalues_for_GRFN(save_values, iNumPaths, iNbEvent, params);
	
	/* Return the error message */
	return err;
}

Err	 lgmSV_mc(	/*	Time Information  */
				int			iNbTime,					
				double		*dTime,
				double		*dDate,
				
				int			iNumPaths,					
					
				/*	Model data Information	*/
				double		dLambdaX,
				double		*Sig,
				double		*SigPsi,					
				double		*dAlpha,
				double		*dLambdaEps,
				double		*dLvlEps,
				double		*dRho,
				double		*dRho2,
				
				/* Parameters */
				LGMSVParam	Params,

				/*	Product data */
				void		**func_parm_tab, 
				int			*EvalEvent,
				
				/*	Payoff function */
				Err (*payoff_func)(
									double	evt_date,
									double	evt_time,
									void	*func_parm,
									
									double	ft,
									double	psi,
															
									/* Vector of results to be updated */
									int		nprod,
									/* Result	*/
									double	*prod_val),
				/*	Result */
				int			iNbProduct, 
				double		**res)
{
Err		err = NULL;
clock_t	time1, time2;
double	ft, psi, v;
int		i, j, k;
double	*temp		= NULL;
double	*res_evt	= NULL;
double	*sum_price	= NULL;
double	*sum_2price	= NULL;

long	seed = -123456789;
double	sqv, brow1, brow2;
FILE	*stream = NULL;
double	ortho_fact;

	/* --------------------------------------------------------------------------------------------------
													Initialisation		
	   -------------------------------------------------------------------------------------------------- */

	/* Transformation of the alpha and Lamdba on Eps	*/
	/* to the equivalent alpha on V = Eps^2				*/	

	if (Params.VerifExpectation = LGMSV_TRUE)
	{
		stream = fopen("C:\\MCDistribution.txt","w+");

		if (!stream)
		{
			err = "Cannot open file in LGMSVMC";
			goto FREE_RETURN;
		}

		ortho_fact = dRho[iNbTime - 2] / dAlpha[iNbTime - 2] * Sig[iNbTime - 2];
		
		fprintf(stream, "\n");
		fprintf(stream, "Maturity	%f\n", dTime[iNbTime - 1]);
		fprintf(stream, "Vol	%f\n", Sig[iNbTime - 2] / sqrt(dTime[iNbTime - 1] - dTime[iNbTime - 2]));
		fprintf(stream, "Alpha	%f\n", dAlpha[iNbTime - 2] / sqrt(dTime[iNbTime - 1] - dTime[iNbTime - 2]));
		fprintf(stream, "Rho	%f\n", dRho[iNbTime - 2]);
		fprintf(stream, "Lambda	%f\n", dLambdaEps[iNbTime - 2]);
		fprintf(stream, "\n");
		fprintf(stream, "Fwd	Vol	LogPhi	Ortho	LnVol\n");
	}

	/* For computational time calculation				 */
	time1 = clock();
		
	temp = dvector (0, iNbProduct - 1);
	res_evt = dvector (0, iNbProduct - 1);
	sum_price = dvector (0, iNbProduct - 1);
	sum_2price = dvector (0, iNbProduct - 1);

	if (!temp || !res_evt || !sum_price || !sum_2price)
	{
		err = "Memory allocation failure in lgmSV_mc";
		goto FREE_RETURN;
	}	

	/* End of initialisation treatment */
	
	/* Initialisation time display */
	time2 = clock();
	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);	

	/* --------------------------------------------------------------------------------------------------
													Convolution		
	   -------------------------------------------------------------------------------------------------- */

	for (i=0; i<iNumPaths; i++)
	{
		/* Initialisation */
		ft = 0.0;
		psi = 0.0;
		v = 1.0;
		
		memset (temp, 0, iNbProduct * sizeof (double));

		if (func_parm_tab[0] && EvalEvent[0])
		{			
			err = payoff_func(	dDate[0],
								dTime[0],
								func_parm_tab[0],
								ft,
								psi,
								iNbProduct,
								temp);

			if (err) goto FREE_RETURN;
		}

		/* For each time t backward */
		for (j=1; j<iNbTime; j++)
		{						
			brow1 = gauss_sample (&seed);
			brow2 = dRho[j-1] * brow1 + dRho2[j-1] * gauss_sample (&seed);
			
			sqv = sqrt(v);
			ft += Sig[j-1] * sqv * brow1;
			psi += SigPsi[j-1] * v;
			v = max(v + dLvlEps[j-1] - dLambdaEps[j-1] * v + dAlpha[j-1] * sqv * brow2, 0.0);
			
			/*	Case of evaluation events */
			if (EvalEvent[j])
			{
				/* Modification of the Payoff at t */			
				err = payoff_func(	dDate[j],
									dTime[j],
									func_parm_tab[j],
									ft,
									psi,
									iNbProduct,
									res_evt);
				
				if (err) goto FREE_RETURN;

				for (k=0; k<iNbProduct; k++)
				{
					temp[k] += res_evt[k];
				}
			}
		}

		if (Params.VerifExpectation = LGMSV_TRUE)
		{
			fprintf(stream, "%f	%f	%f	%f	%f\n", ft, v, log(psi), ft - ortho_fact * (v - 1.0), log(max(v, 1.0E-05)));
		}

		for (k=0; k<iNbProduct; k++)
		{
			sum_price[k] += temp[k] / iNumPaths;
			sum_2price[k] += temp[k] * temp[k] / iNumPaths;
		}
	}
	
	for (k=0; k<iNbProduct; k++)
	{
		res[k][0] = sum_price[k];
		res[k][1] = sqrt ((sum_2price[k] - sum_price[k] * sum_price[k]) / iNumPaths);
	}	

	/* Convolution time display */
	time1 = clock();
	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (time1 - time2) / CLOCKS_PER_SEC);

	if (Params.VerifExpectation = LGMSV_TRUE)
	{
		fclose(stream);
	}

FREE_RETURN:

	if (temp)
	{
		free_dvector (temp, 0, iNbProduct - 1);
	}

	if (res_evt)
	{
		free_dvector (res_evt, 0, iNbProduct - 1);
	}

	if (sum_price)
	{
		free_dvector (sum_price, 0, iNbProduct - 1);
	}

	if (sum_2price)
	{
		free_dvector (sum_2price, 0, iNbProduct - 1);
	}

	/* Return the error message */
	return err;
}

Err	 lgmSV2F_mc(	
					/*	Time Information  */
					int			iNbTime,					
					double		*dTime,
					double		*dDate,
					
					int			iNumPaths,					
						
					/*	Model data Information	*/
					double		dLambdaX1,
					double		dLambdaX2,
					double		*Sig1,
					double		*Sig2,
					double		*SigPsi1,
					double		*SigPsi2,
					double		*SigPsi12,
										
					double		*dAlpha,
					double		*dLambdaEps,
					double		*dLvlEps,
					double		*dRho,
					double		*dRho2,

					double		*Coef_Ortho21,
					double		*Coef_Ortho22,
					double		*Coef_Ortho31,
					double		*Coef_Ortho32,
					double		*Coef_Ortho33,

					/* Parameters */
					LGMSVParam	Params,

					/*	Product data */
					void		**func_parm_tab, 
					int			*EvalEvent,											
					
					/*	Payoff function */
					Err (*payoff_func)(	double	evt_date,
										double	evt_time,
										void	*func_parm,

										/* Model data	*/
										double	ft1,
										double	ft2,
										double	phi1,
										double	phi2,
										double	phi12,
																
										/* Vector of results to be updated */
										int		nprod,
										/* Result	*/
										double	*prod_val),
					/*	Result */
					int			iNbProduct, 
					double		**res)
{
Err		err = NULL;
clock_t	time1, time2;
double	ft1, ft2, psi1, psi2, psi12, v, sqv;
int		i, j, k;
double	*temp				= NULL;
double	*res_evt			= NULL;
double	*sum_price			= NULL;
double	*sum_2price			= NULL;

long	seed = -123456789;
double	brow1, brow2, brow3;

FILE	*stream = NULL;

	/* --------------------------------------------------------------------------------------------------
													Initialisation		
	   -------------------------------------------------------------------------------------------------- */	
	
	/* For computational time calculation				 */
	time1 = clock();		

	temp = dvector (0, iNbProduct - 1);
	res_evt = dvector (0, iNbProduct - 1);
	sum_price = dvector (0, iNbProduct - 1);
	sum_2price = dvector (0, iNbProduct - 1);

	if (!temp || !res_evt || !sum_price || !sum_2price)
	{
		err = "Memory allocation failure in lgmSV_mc";
		goto FREE_RETURN;
	}	

	/* End of initialisation treatment */
	
	/* Initialisation time display */
	time2 = clock();
	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);	

	/* --------------------------------------------------------------------------------------------------
													Convolution		
	   -------------------------------------------------------------------------------------------------- */

	for (i=0; i<iNumPaths; i++)
	{
		/* Initialisation */
		ft1 = 0.0;
		ft2 = 0.0;
		psi1 = 0.0;
		psi2 = 0.0;
		psi12 = 0.0;
		v = 1.0;
		
		memset (temp, 0, iNbProduct * sizeof (double));

		if (func_parm_tab[0] && EvalEvent[0])
		{			
			err = payoff_func(	dDate[0],
								dTime[0],
								func_parm_tab[0],								
								ft1,
								ft2,
								psi1,
								psi2,
								psi12,
								iNbProduct,
								temp);

			if (err)
			{
				goto FREE_RETURN;
			}								
		}

		/* For each time t backward */
		for (j=1; j<iNbTime; j++)
		{			
			brow1 = gauss_sample (&seed);
			brow2 = Coef_Ortho21[j-1] * brow1 + Coef_Ortho22[j-1] * gauss_sample (&seed);
			brow3 = Coef_Ortho31[j-1] * brow1 + Coef_Ortho32[j-1] * brow2 + Coef_Ortho33[j-1] * gauss_sample (&seed);
		
			sqv = sqrt(v);
			ft1 += Sig1[j-1] * sqv * brow1;
			ft2 += Sig2[j-1] * sqv * brow2;

			psi1 += SigPsi1[j-1] * v;
			psi2 += SigPsi2[j-1] * v;
			psi12 += SigPsi12[j-1] * v;

			v = max(v + dLvlEps[j-1] - dLambdaEps[j-1] * v + dAlpha[j-1] * sqv * brow3, 0.0);
					
			/*	Case of evaluation events */
			if (EvalEvent[j])
			{
				/* Modification of the Payoff at t */			
				err = payoff_func(	dDate[j],
									dTime[j],
									func_parm_tab[j],
									ft1,
									ft2,
									psi1,
									psi2,
									psi12,
									iNbProduct,
									res_evt);
				
				if (err) goto FREE_RETURN;

				for (k=0; k<iNbProduct; k++)
				{
					temp[k] += res_evt[k];
				}
			}
		}
		
		for (k=0; k<iNbProduct; k++)
		{
			sum_price[k] += temp[k] / iNumPaths;
			sum_2price[k] += temp[k] * temp[k] / iNumPaths;
		}
	}
	
	for (k=0; k<iNbProduct; k++)
	{
		res[k][0] = sum_price[k];
		res[k][1] = sqrt ((sum_2price[k] - sum_price[k] * sum_price[k]) / iNumPaths);
	}	

	/* Convolution time display */
	time1 = clock();
	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (time1 - time2) / CLOCKS_PER_SEC);	

FREE_RETURN:

	if (temp)
	{
		free_dvector (temp, 0, iNbProduct - 1);
	}

	if (res_evt)
	{
		free_dvector (res_evt, 0, iNbProduct - 1);
	}

	if (sum_price)
	{
		free_dvector (sum_price, 0, iNbProduct - 1);
	}

	if (sum_2price)
	{
		free_dvector (sum_2price, 0, iNbProduct - 1);
	}

	/* Return the error message */
	return err;
}

Err	 lgmSV_mc_cv(	
					/*	Time Information  */
					int			iNbTime,					
					double		*dTime,
					double		*dDate,
					
					int			iNumPaths,					
						
					/*	Model data Information	*/
					double		dLambdaX,
					double		*Sig,
					double		*SigPsi,
										
					double		*dAlpha,
					double		*dLambdaEps,
					double		*dLvlEps,
					double		*dRho,
					double		*dRho2,					

					/* Parameters */
					LGMSVParam	Params,

					/*	Product data */
					void		**func_parm_tab, 
					int			*EvalEvent,										
					
					/*	Payoff function */
					Err (*payoff_func)(
										double	evt_date,
										double	evt_time,
										void	*func_parm, 										
										
										double	ft,
										double	psi,
																
										/* Vector of results to be updated */
										int		nprod,
										/* Result	*/
										double	*prod_val
										),
					/*	Result */
					int			iNbProduct, 
					double		**res)
{
Err				err = NULL;
clock_t			time1, time2;
double			ft, psi, v;
double			ft_bs, psi_bs;
int				i, j, k;
double	*temp				= NULL;
double	*res_evt			= NULL;
double	*sum_price			= NULL;
double	*sum_2price			= NULL;

double	*temp_bs			= NULL;
double	*res_evt_bs			= NULL;
double	*sum_price_bs		= NULL;
double	*sum_2price_bs		= NULL;
double	*cov_price_bs		= NULL;

long	seed = -123456789;
double	sqv, brow1, brow2;

	/* --------------------------------------------------------------------------------------------------
													Initialisation		
	   -------------------------------------------------------------------------------------------------- */

	/* For computational time calculation				 */
	time1 = clock();

	/* Transformation of the alpha and Lamdba on Eps	*/
	/* to the equivalent alpha on V = Eps^2				*/	
	
	temp = dvector (0, iNbProduct - 1);
	res_evt = dvector (0, iNbProduct - 1);
	sum_price = dvector (0, iNbProduct - 1);
	sum_2price = dvector (0, iNbProduct - 1);

	temp_bs = dvector (0, iNbProduct - 1);
	res_evt_bs = dvector (0, iNbProduct - 1);
	sum_price_bs = dvector (0, iNbProduct - 1);
	sum_2price_bs = dvector (0, iNbProduct - 1);
	cov_price_bs = dvector (0, iNbProduct - 1);

	if (!temp || !res_evt || !sum_price || !sum_2price || !temp_bs || !res_evt_bs || !sum_price_bs || !sum_2price_bs || !cov_price_bs)
	{
		err = "Memory allocation failure in lgmSV_mc_bs";
		goto FREE_RETURN;
	}	

	/* End of initialisation treatment */
	
	/* Initialisation time display */
	time2 = clock();
	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (time2 - time1) / CLOCKS_PER_SEC);	

	/* --------------------------------------------------------------------------------------------------
													Convolution		
	   -------------------------------------------------------------------------------------------------- */

	for (i=0; i<iNumPaths; i++)
	{
		/* Initialisation */
		ft = 0.0;
		psi = 0.0;
		v = 1.0;

		ft_bs = 0.0;
		psi_bs = 0.0;

		memset (temp, 0, iNbProduct * sizeof (double));
		memset (temp_bs, 0, iNbProduct * sizeof (double));

		if (func_parm_tab[0] && EvalEvent[0])
		{			
			err = payoff_func(	dDate[0],
								dTime[0],
								func_parm_tab[0],								
								ft,
								psi,
								iNbProduct,
								temp);

			if (err)
			{
				goto FREE_RETURN;
			}								
		}

		/* For each time t backward */
		for (j=1; j<iNbTime; j++)
		{						
			brow1 = gauss_sample (&seed);
			brow2 = dRho[j] * brow1 + dRho2[j] * gauss_sample (&seed);
			
			sqv = sqrt(v);
			ft += Sig[j-1] * sqv * brow1;
			psi += SigPsi[j-1] * v;
			v = max(v - dLambdaEps[j] * (v - 1.0) + dAlpha[j] * sqv * brow2, 0.0);

			ft_bs += Sig[j-1] * brow1;
			psi_bs += SigPsi[j-1];
			
			/*	Case of evaluation events */
			if (EvalEvent[j])
			{
				/* Modification of the Payoff at t */			
				err = payoff_func(	dDate[j],
									dTime[j],
									func_parm_tab[j],
									ft,
									psi,
									iNbProduct,
									res_evt);

				/* Modification of the Payoff at t */			
				err = payoff_func(	dDate[j],
									dTime[j],
									func_parm_tab[j],									
									ft_bs,
									psi_bs,
									iNbProduct,
									res_evt_bs);

				if (err) goto FREE_RETURN;

				for (k=0; k<iNbProduct; k++)
				{
					temp[k] += res_evt[k];
					temp_bs[k] += res_evt_bs[k];
				}
			}
		}

		for (k=0; k<iNbProduct; k++)
		{
			sum_price[k] += temp[k] / iNumPaths;
			sum_2price[k] += temp[k] * temp[k] / iNumPaths;

			sum_price_bs[k] += temp_bs[k] / iNumPaths;
			sum_2price_bs[k] += temp_bs[k] * temp_bs[k] / iNumPaths;
			cov_price_bs[k] += temp[k] * temp_bs[k] / iNumPaths;
		}
	}
	
	for (k=0; k<iNbProduct; k++)
	{
		res[k][0] = sum_price[k];
		res[k][1] = sqrt ((sum_2price[k] - sum_price[k] * sum_price[k]) / iNumPaths);

		res[k+iNbProduct][0] = sum_price_bs[k];
		res[k+iNbProduct][1] = sqrt ((sum_2price_bs[k] - sum_price_bs[k] * sum_price_bs[k]) / iNumPaths);

		res[k+2*iNbProduct][0] = (cov_price_bs[k] - res[k][0] * res[k+iNbProduct][0]) / iNumPaths;
		res[k+2*iNbProduct][1] = res[k+iNbProduct][1] * res[k+iNbProduct][1];
	}	

	/* Convolution time display */
	time1 = clock();
	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (time1 - time2) / CLOCKS_PER_SEC);


FREE_RETURN:

	if (temp)
	{
		free_dvector (temp, 0, iNbProduct - 1);
	}

	if (res_evt)
	{
		free_dvector (res_evt, 0, iNbProduct - 1);
	}

	if (sum_price)
	{
		free_dvector (sum_price, 0, iNbProduct - 1);
	}

	if (sum_2price)
	{
		free_dvector (sum_2price, 0, iNbProduct - 1);
	}

	/* Return the error message */
	return err;
}

Err LGMSVOptionMC(
				char			*und_name,						/*	Name of the underlying */
				char			*yc_name,						/*	Name of the yield curve */
				char			*ref_rate_name,					/*	Name of the reference rate */
				char			*swaption_freq,					/*	Frequency and basis of underlying swaptions */
				char			*swaption_basis,

				long			lExDate,
				long			lEndDate,
				double			*dStrike,
				int				nb_strike,
				int				pay_rec,						/*	pay:1 rec:-1 */				
				
				long			iNumPaths,
				long			nb_vol,
				int				use_balsam,

				/* Output */
				double			*pSwaptionPrice,
				double			*pStd)
{
	int				i, j, ncpn;
	SrtCompounding	ifreq;
	SrtBasisCode	ibasis;
	long			cpn_date[MAX_CPN];
	long			theo_dates[MAX_CPN];
	double			cpn_time[MAX_CPN],
					cpn_cvg[MAX_CPN];
				
	long			theo_date, act_date, temp_date;
	long			today; 

	double			lExTime;
	double			swp_rte, swp_cash;

	SrtCurvePtr		yc_ptr;
	Err				err				= NULL;

	int				dNbSig;
	double			*dSigTime		= NULL,
					*dSig			= NULL,
					*dAlphaTS		= NULL,
					*dRhoTS			= NULL,
					*dLambdaEpsTS	= NULL,
					**Coupon		= NULL;

	double			dTau, dAlpha, dRho, dLambdaEps, dTStar;
	int				one2F;

	yc_ptr = lookup_curve (yc_name);
	if (!yc_ptr)
	{
		err = "Yield Curve not found";
		goto FREE_RETURN;
	}

	today = get_today_from_curve (yc_ptr);	

	/*	1.)	Setup the bond schedule and its coupons */
	
	/*	Coupons */

	err = interp_compounding (swaption_freq, &ifreq);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = interp_basis (swaption_basis, &ibasis);
	if (err)
	{
		goto FREE_RETURN;
	}

	theo_date = lEndDate;
	act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
	ncpn = 1;

	while (act_date > lExDate)
	{
		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);
		act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		ncpn++;
	}
	ncpn--;

	if (ncpn < 2)
	{
		err = "Not enough coupons";
		goto FREE_RETURN;		
	}

	theo_date = lEndDate;
	act_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
	i = ncpn - 1;

	while (i >= 0)
	{
		cpn_time[i] = (act_date - today) * YEARS_IN_DAY; 
		cpn_date[i] = act_date;
		theo_dates[i] = theo_date;

		theo_date = add_unit (theo_date, - 12 / ifreq, SRT_MONTH, NO_BUSDAY_CONVENTION);

		temp_date = bus_date_method (theo_date, MODIFIED_SUCCEEDING);
		cpn_cvg[i] = coverage (temp_date, act_date, ibasis);
		act_date = temp_date;

		i--;
	}

	Coupon = dmatrix(0, nb_strike-1, 0, ncpn-1);

	if (!Coupon)
	{
		err = "Memory allocation faillure";
		goto FREE_RETURN;
	}


	
	cpn_cvg[0] = pay_rec * swp_f_df(today, cpn_date[0], yc_name);	

	for (j=0; j<nb_strike; j++)
	{
		Coupon[j][0] = cpn_cvg[0];
	}

	for (i=1; i<ncpn; i++)
	{
		/* calculates the basis swap */
		err = swp_f_ForwardRate(
								cpn_date[i-1],
								theo_dates[i],
								swaption_freq,
								swaption_basis,
								yc_name,
								ref_rate_name,
								&swp_rte);

		swp_cash = (swp_f_df(today, cpn_date[i-1], yc_name) - swp_f_df(today, cpn_date[i], yc_name)) / (cpn_cvg[i] * swp_f_df(today, cpn_date[i], yc_name));
	
		for (j=0; j<nb_strike; j++)
		{
			 Coupon[j][i] = -cpn_cvg[i] * pay_rec * (dStrike[j] - (swp_rte - swp_cash)) * swp_f_df(today, cpn_date[i], yc_name);
		}
	}


	for (j=0; j<nb_strike; j++)
	{
		Coupon[j][ncpn-1] -= pay_rec * swp_f_df(today, cpn_date[ncpn-1], yc_name);
	}

	lExTime = (lExDate - today) * YEARS_IN_DAY;

	/* Get the TS */

	err = Get_LGMSV_TermStructure(	und_name,
									&dSigTime,
									&dSig,
									&dAlphaTS,
									&dRhoTS,
									&dLambdaEpsTS,
									&dTStar,
									&dNbSig,
									&dTau,
									&one2F,
									NULL,
									NULL,
									NULL);

	if (err)
	{
		goto FREE_RETURN;
	}

	/* for now no TS */

	dAlpha = dAlphaTS[0];
	dRho = dRhoTS[0],
	dLambdaEps = dLambdaEpsTS[0];

	for (i=1; i<dNbSig; i++)
	{
		if ((fabs(dAlphaTS[i] - dAlpha) > 1.0E-08) || (fabs(dRhoTS[i] - dRho) > 1.0E-08) || (fabs(dLambdaEpsTS[i] - dLambdaEps) > 1.0E-08))
		{
			err = "no alpha, rho, lameps TS allowed for now";
			goto FREE_RETURN;
		}
	}
	
	LGMSVOptionMCts(
					1.0 / dTau,
					dNbSig,		/* Term Structure of g(t) */
					dSigTime,		
					dSig,
					dTStar,			/* Tstar in years from today */
					dAlpha,		/* Alpha of V = Eps^2 */
					dLambdaEps,	/* LambdaEps of V = Eps^2 */
					dRho,
					lExDate,							/* Exercice date of the swaption  */
					lExTime,		/* Exercice of the swaption in years from today */					 
					ncpn,		/* Description of the cashflows */					
					cpn_time,
					cpn_date,
					Coupon,
					nb_strike,
					yc_name,	/* Yield Curve */					 
					iNumPaths,
					nb_vol,
					pSwaptionPrice,
					pStd);



FREE_RETURN:

	if (dSigTime) free (dSigTime);
	if (dSig) free (dSig);
	if (Coupon) free_dmatrix(Coupon, 0, nb_strike-1, 0, ncpn-1);

	return err;
}

Err LGMSVOptionMCts(	
						/* Parameter of diffusion */
						/* BE CAREFULL : Alpha and LambdaEps for V=Eps^2 */
						double	dLambdaX,
						int		iNbSigTime,		/* Term Structure of g(t) */
						double	*SigTime,		
						double	*Sig,
						double	dTStar,			/* Tstar in years from today */						 
						double	dAlpha,
						double	dLambdaEps,
						double	dRho,

						/* Product description */
						long	lExDate,		/* Exercice date of the swaption  */
						double	dExTime,	/* Exercice of the swaption in years from today */
						int		iNbCoupon,		/* Description of the cashflows */
						double	*CouponTime,
						long	*CouponDate,
						double	**Coupon,
						int		nb_strike,
						char	*cYieldCurve,	/* Yield Curve */

						/* Parameter of grids */
						long	iNumPaths,
						long	nb_vol,
												
						/* Outputs */
						double	*Price,
						double	*Std)
{
Err		err = NULL;
int		i, j, k, nbstep;
double	*time		= NULL,
		*vol		= NULL,
		*cpn		= NULL,
		*lamdt		= NULL,
		*dti		= NULL,
		*beta		= NULL,
		*volX		= NULL,
		*volX2		= NULL,
		*sig		= NULL,
		*var		= NULL,
		*drift		= NULL,
		*init_gauss	= NULL,
		*sum		= NULL,
		*sum2		= NULL,
		*covar		= NULL;

long	*nbpT	= NULL;

double	rho2, sqrho2;
double	lgm_price;
double	brown;
long	seed = -123456789;
double	prob, step;
long	endi, rand;
double	MaxTime, dt, time1, time2;
double	check_var, check_vol, check_vol2, varvar, expect_var;
int		save_file = 0;
FILE	*stream;

	/* odd number of paths */
	iNumPaths = (int) (iNumPaths / 2) * 2 + 1;

	/* construct the time discretisation */
	MaxTime = dExTime / nb_vol;
	endi = Get_Index(dExTime, SigTime, iNbSigTime);
	
	vol = dvector(0, endi);
	volX = dvector(0, endi);
	volX2 = dvector(0, endi);
	lamdt = dvector(0, endi);
	dti = dvector(0, endi);
	nbpT = lvector(0, endi);
	cpn = dvector(0, iNbCoupon-1);
	beta = dvector(0, iNbCoupon-1);	

	sum = dvector(0, nb_strike-1);
	sum2 = dvector(0, nb_strike-1);
	covar = dvector(0, nb_strike-1);

	sig = dvector(0, iNumPaths-1);
	var = dvector(0, iNumPaths-1);
	drift = dvector(0, iNumPaths-1);

	init_gauss = dvector(0, iNumPaths-1);

	if (!vol || !lamdt || !nbpT || !cpn || !beta || !volX || !volX2
		|| !dti || !sum || !sum2 || !covar || !init_gauss || !sig || !var || !drift)
	{
		err = "Memory allocation faillure (1) in LGMSVOptionMCts";
		goto FREE_RETURN;
	}

	/* beta coef for reconstruction formula */
	for (i=0; i<iNbCoupon; i++)
	{
		beta[i] = (1.0 - exp (-dLambdaX * (CouponTime[i] - dTStar))) / dLambdaX;
	}

	/* Time discretisation */
	nbstep = 0;
	expect_var = 0.0;

	for (i=endi; i>=0; i--)
	{
		if (i>0)
		{
			time1 = SigTime[i-1];
		}
		else
		{
			/* First part */
			time1 = 0.0;
		}
		
		if (i==endi || endi==0)
		{
			/* Last part */
			time2 = dExTime;
		}
		else
		{
			time2 = SigTime[i];
		}

		nbpT[i] = max((int) ((time2 - time1) / MaxTime + 0.5), 1);
		nbstep += nbpT[i];
		dt = (time2 - time1) / nbpT[i];
		vol[i] = dAlpha * sqrt(dt);
		volX[i] = Sig[i] * sqrt(dt);
		volX2[i] = volX[i] * volX[i];
		lamdt[i] = dLambdaEps * dt;
		dti[i] = dt;
		expect_var += Sig[i] * Sig[i] * (time2 - time1);
	}
	
	/* constant calculation */
	rho2 = dRho * dRho;
	sqrho2 = sqrt(1.0 - rho2);	

	/* Gauss initialisation */
	iNumPaths -= 1;
	iNumPaths /= 2;	
	step = 0.5 / (iNumPaths + 1);
	prob = step;

	/* Generation of the fractiles of the gaussian */

	for (i=0; i<iNumPaths; i++)
	{
		init_gauss[i] = inv_cumnorm_fast(prob);
		init_gauss[iNumPaths+i+1] = -init_gauss[i];
		prob += step;
	}
	init_gauss[iNumPaths] = 0.0;	
	iNumPaths *= 2;
	iNumPaths += 1;

	/* initialisation */

	memset(sum, 0, nb_strike * sizeof(double));
	memset(sum2, 0, nb_strike * sizeof(double));
	memset(covar, 0, nb_strike * sizeof(double));	
	memset(drift, 0, iNumPaths * sizeof(double));
	memset(var, 0, iNumPaths * sizeof(double));

	for (i=0; i<iNumPaths; i++)
	{
		sig[i] = 1.0;
	}

	/* Launch the MC */
	for (j=0; j<=endi; j++)
	{
		for (i=0; i<iNumPaths; i++)
		{
			var[i] += 0.5 * volX2[j];
		}

		for (k=0; k<nbpT[j]-1; k++)
		{
			for (i=0; i<iNumPaths; i++)
			{
				rand = i + (int) ((iNumPaths-i) * uniform (&seed));
				brown = init_gauss[rand];
				init_gauss[rand] = init_gauss[i];
				init_gauss[i] = brown;

				brown *= sqrt(sig[i]);

				drift[i] += brown * volX[j];
				sig[i] += -lamdt[j] * (sig[i] - 1.0) + vol[j] * brown;

				if (sig[i] < 1.0E-08)
				{
					sig[i] = 1.0E-08;
				}				
				
				var[i] += sig[i] * volX2[j];				
			}
		}

		/* last step */
		for (i=0; i<iNumPaths; i++)
		{
			rand = i + (int) ((iNumPaths-i) * uniform (&seed));
			brown = init_gauss[rand];
			init_gauss[rand] = init_gauss[i];
			init_gauss[i] = brown;

			brown *= sqrt(sig[i]);

			drift[i] += brown * volX[j];
			sig[i] += -lamdt[j] * (sig[i] - 1.0) + vol[j] * brown;

			if (sig[i] < 1.0E-08)
			{
				sig[i] = 1.0E-08;
			}				
			
			var[i] += 0.5 * sig[i] * volX2[j];				
		}
	}

	check_vol = 0.0;
	check_vol2 = 0.0;
	check_var = 0.0;
	varvar = 0.0;	

	for (i=0; i<iNumPaths; i++)
	{
		for (j=0; j<nb_strike; j++)
		{
			lgm_price = lgmopval_stochvol(iNbCoupon, Coupon[j], beta, var[i], dRho * drift[i], rho2, sqrho2);
			sum[j] += lgm_price / iNumPaths;
			sum2[j] += lgm_price * lgm_price / iNumPaths;
			covar[j] += var[i] * lgm_price / iNumPaths;
		}

		check_vol += sig[i] / iNumPaths;
		check_vol2 += sig[i] * sig[i] / iNumPaths;
		check_var += var[i] / iNumPaths;
		varvar += var[i] * var[i] / iNumPaths;
	}

	varvar = varvar - check_var * check_var;

	if (save_file)
	{
		stream = fopen("C:\\LGMSVMCsave.txt","w+");

		nbstep = iNumPaths / 65000 + 1;

		i = 0;
		while (i<iNumPaths)
		{
			for (j=0; j<nbstep && i<iNumPaths; j++)
			{
				fprintf(stream, "%f	%f	%f	", sig[i], drift[i], var[i]);
				i++;
			}

			fprintf(stream, "\n");
		}

		fclose(stream);
	}

	for (j=0; j<nb_strike; j++)
	{		
		Price[j] = sum[j];
		Std[j] = sqrt((sum2[j] - sum[j] * sum[j]) / iNumPaths);
		covar[j] = covar[j] - check_var * sum[j];
		Std[j] = sum[j] - covar[j] / varvar * (check_var - expect_var);
	}

FREE_RETURN:

	if (vol) free_dvector(vol, 0, endi);	
	if (volX) free_dvector(volX, 0, endi);
	if (volX2) free_dvector(volX2, 0, endi);
	if (lamdt) free_dvector(lamdt, 0, endi);
	if (dti) free_dvector(dti, 0, endi);
	if (nbpT) free_lvector(nbpT, 0, endi);
	if (cpn) free_dvector(cpn, 0, iNbCoupon-1);
	if (beta) free_dvector(beta, 0, iNbCoupon-1);	

	if (sum) free_dvector(sum, 0, nb_strike-1);
	if (sum2) free_dvector(sum2, 0, nb_strike-1);
	if (covar) free_dvector(covar, 0, nb_strike-1);

	if (sig) free_dvector(sig, 0, nb_strike-1);
	if (var) free_dvector(var, 0, nb_strike-1);
	if (drift) free_dvector(drift, 0, nb_strike-1);

	if (init_gauss) free_dvector(init_gauss, 0, nb_strike-1);

	return err;
}

/* Put points between original points to ensure time between two points does not exceed Max Time */
Err	fill_time_vector_max_time(	int		iInitNbTimes,
								double	*dInitTimes,							  
								double	dMaxTime,
								int		*iNewNbTimes,
								double	**dNewTimes)
{
Err	err = NULL;
int	i, j, k;
int	nstp;
double	local_min_time;
	
	/* First calculate the new number of steps */
	*iNewNbTimes = 1;

	for (i=0; i<iInitNbTimes-1; i++)
	{
		/* Fill between times[i] and times[i+1] */
		*iNewNbTimes += max((int) ((dInitTimes[i+1] - dInitTimes[i]) / dMaxTime + 0.5), 1);
	}

	/* Then allocate memory */
	*dNewTimes = calloc(*iNewNbTimes, sizeof(double));

	if (!*dNewTimes)
	{
		return "Memory allocation faillure in fill_time_vector_max_time";
	}

	/* Eventually fill the new vector */
	*dNewTimes[0] = dInitTimes[0];
	j = 1;

	for (i=0; i<iInitNbTimes-1; i++)
	{
		/* Fill between times[i] and times[i+1] */
		nstp = max((int) ((dInitTimes[i+1] - dInitTimes[i]) / dMaxTime + 0.5), 1);
		local_min_time  = (dInitTimes[i+1] - dInitTimes[i]) / (nstp * 1.0);
		
		for (k=0; k<nstp-1; k++)
		{
			(*dNewTimes)[j] = (*dNewTimes)[j-1] + local_min_time;
			j++;
		}
		
		(*dNewTimes)[j] = dInitTimes[i+1];
		j++;
	}

	return NULL;
}



