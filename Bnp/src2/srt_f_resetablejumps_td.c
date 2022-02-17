/******************************************************************************
Computes the price of a resetable range floater using a Market model with jumps 

Authors: Ezra Nahum
		 Yann Samuelides



********************************************************************************/



#include "math.h"
#include "srt_h_all.h"
#include "num_h_allhdr.h"
#include "num_h_levenberg.h"
#include "num_h_hermite.h"
#include  "swp_h_cms.h"
#include  <utallhdr.H>
#include  <opfnctns.H>
#include		<swp_h_all.h"
#include		<swp_h_external_fct.h"
#include		<swp_h_vol.h"
#include  <srt_h_resetable.h"

#define NMAXSTRIKES 10

/*Free Memory functions*/

Err free_resetableperiod_td(double **paramperiod,double *periodmat, double *DRStsStart,
			   double *DRSts, double *isfriday, double nfixdays, int  n_resetdates, double **mat, double **paramaftereset)
{
	if (paramperiod) free_dmatrix (paramperiod, 1, (long) nfixdays, 1, 5);
	paramperiod = NULL;
	
	if (periodmat) free_dvector (periodmat, 1, (long) nfixdays);
	periodmat = NULL;
	
	if (DRSts) free_dvector (DRSts, 1, (long) nfixdays);
	DRSts = NULL;

	if (DRStsStart) free_dvector (DRSts, 1, (long) nfixdays);
	DRStsStart = NULL;

	if (isfriday) free_dvector (isfriday, 1, (long) nfixdays);
	isfriday = NULL;

	if (mat) free_dmatrix (mat, 1, (long) nfixdays, 1, n_resetdates);
	mat = NULL;
	
	if (paramaftereset) free_dmatrix (paramaftereset, 1, (long) nfixdays, 1, 5);
	paramaftereset = NULL;
	
	return NULL;
}	

Err free_resetablegeneral_td(double *weights, double *points, double *wn1n2, 
						  double **paramatrix, int nherm,int n_resetdates, double *caliberror)

{
	if (weights) free_dvector (weights, 0, 2*nherm);
	weights = NULL;
	
	if (points) free_dvector (points, 0, 2*nherm);
	points = NULL;
	
	if (wn1n2) free_dvector (wn1n2, 1, 3);
	wn1n2 = NULL;

	if (paramatrix) free_dmatrix (paramatrix, 1, n_resetdates,1,5);
	paramatrix = NULL;

	if (caliberror) free_dvector(caliberror,1,n_resetdates);
	caliberror = NULL;
	
	return NULL;
}




/* Function which calculates the value of a Fra */

Err get_fra(long Fradate, int spotlag, char *cRefRateCode, char *szYieldCurveName, double *dFra);




/* Function which calibrates the parameters of the jump model to a caplet smile
corresponding to a specific matirity date */

Err calibtosmile_td(  long today,
				   long dDate,
				   char *cMarketId,
				   char *cRefRateCode,
				   char *szVolCurveName,	
				   char *szYieldCurveName,
				   int  nperiods,
				   double *param,
				   double *calerr,
				   double *mat,
				   long   *ia)
{
Err err=NULL;
int i;
int spotlag;
double atmvol;
double lognorm=1.0;
long  dDatespot,dDatespot3;
double Strikeinf, Strikesup;
double totalmat;
double dFra;
double Strikespread;
double *strikes;
double *vols;
double *calibvols;
double s;
int  num_of_months;
SrtCrvPtr yldcrv;
SrtBasisCode   float_basis;
SrtCompounding  float_compounding;

yldcrv = lookup_curve(szYieldCurveName);
spotlag = get_spotlag_from_curve(yldcrv);

/*Memory Allocation*/

strikes = dvector(1,NMAXSTRIKES+1);
vols = dvector(1,NMAXSTRIKES+1);
calibvols = dvector(1,NMAXSTRIKES+1);


/* get info from the RefRate Code*/
if (err=swp_f_get_ref_rate_details(cRefRateCode, &float_basis, &float_compounding))
{
	smessage("Error if swp_f_get_ref_rate_details");
}

num_of_months = 12/float_compounding;


/* We now get the strikes and the vols of the smile we will calibrate on*/

err=get_fra(dDate, spotlag, cRefRateCode, szYieldCurveName, &dFra);


dDatespot = add_unit(dDate,spotlag,SRT_BDAY,MODIFIED_SUCCEEDING);


dDatespot3 = add_unit(dDatespot, num_of_months, SRT_MONTH, MODIFIED_SUCCEEDING);


if (err = swp_f_vol(szVolCurveName,(double) dDatespot, (double) dDatespot3,dFra,&atmvol,&lognorm))
{
	smessage("Error in swp_f_vol");
}

totalmat = (dDate-today)/365.0;

Strikeinf = dFra-2*atmvol*dFra*sqrt(totalmat);
Strikesup = dFra+2*atmvol*dFra*sqrt(totalmat);
Strikespread = (Strikesup-Strikeinf)/NMAXSTRIKES+0.00015;


s = Strikeinf-0.001;
i=1;
strikes[i] = s;
if (err = swp_f_vol(szVolCurveName,(double) dDatespot, (double) dDatespot3,strikes[i],&(vols[i]),&lognorm))
{
	smessage("Error in swp_f_vol");
}
i=2;


while (s <= Strikesup)
{
	s = s+Strikespread;
	strikes[i] = s;


if (err = swp_f_vol(szVolCurveName,(double) dDatespot, (double) dDatespot3,strikes[i],&(vols[i]),&lognorm))
{
	smessage("Error in swp_f_vol");
}

	i++;

}


/* call the function which calibrates the Merton model to a given smile */
err = optcalibmertontimedep(dFra,NMAXSTRIKES,nperiods,strikes,vols,param,mat,"Lognormal",calerr,ia,calibvols);



/*free memory*/

free_dvector(strikes,1,NMAXSTRIKES+1);
free_dvector(vols,1,NMAXSTRIKES+1);
free_dvector(calibvols,1,NMAXSTRIKES+1);



return err;
}


/* General calibration function which calibrates the parameters to the smiles of the resetdates
of the resetable as well as the smile of the end date of the product 
(going from the longer maturity to the shorter one)*/



Err resetablecalib_td(long today, int n_resetdates, 
				   double *resetdates,
				   char *cMarketId,
				   char *szVolCurveName,
				   char *cRefRateCode,
				   char *szYieldCurveName,
				   double *paraminitguess,
				   double **paramatrix,
				   double *caliberr)

{
Err err=NULL;
int i,j;
long   *freeze;
double *mat;
double *startparam; /*first: calibration to the longest maturity*/

/*memory allocation*/
freeze = lngvector(1,5*n_resetdates);
mat = dvector(1,n_resetdates);
startparam = dvector(1, 5);

mat[1] = (resetdates[1]-today)/365.0;
for (i = 2; i<=n_resetdates;i++)
{
mat[i] =(resetdates[i]-resetdates[i-1])/365.0; 
}


for (i=1; i<=5;i++)
{
	freeze[i] = 1;
	startparam[i] = paraminitguess[i];
}


/*freeze the intensities of the jumps and the parameters that are determined by the stationarity assumption*/
freeze[3] = 0;
freeze[5] = 0;
for (i=6;i<=5*n_resetdates;i++)
{
	freeze[i] = 0;
}


/* now calibrates the  first smile */

for (i=1;i<=5;i++)
{
	paramatrix[1][i] = startparam[i];
}

err = calibtosmile_td(today, 
					(long) resetdates[1],
				   cMarketId,
				   cRefRateCode,
				   szVolCurveName,
				   szYieldCurveName,
				   1,
				   paramatrix[1],
				   &(caliberr[1]),
				   mat,
				   freeze);

for (j = 2; j<=n_resetdates;j++)
{

	for (i=1;i<=5;i++)
	{
		paramatrix[j][i] = startparam[i];
	}

	for (i = 6; i<=5*j;i++)
	{
		paramatrix[j][i] = paramatrix[j-1][i-5];	

	}


err = calibtosmile_td(today, 
					(long) resetdates[j],
				   cMarketId,
				   cRefRateCode,
				   szVolCurveName,
				   szYieldCurveName,
				   j,
				   paramatrix[j],
				   &(caliberr[j]),
				   mat,
				   freeze);


}

/* free memory */

free_lngvector(freeze,1,5*n_resetdates);
free_dvector(mat,1,n_resetdates);
free_dvector(startparam,1,5);
return err;

}


/* function which linearly interpolates the parameters obtained from the smiles of the beginning
and end of the period , for a maturity in between*/

Err get_param_td(double Fradate, int  indexStart, double *resetdates, 
			  double **paramatrix, double *paramfra)
{
Err err=NULL;
int j;
double datespread;

/*linear interpolation of the parameters*/

datespread = (resetdates[indexStart+1]-resetdates[indexStart]);

paramfra[1] = (paramatrix[indexStart+1][1]-paramatrix[indexStart][1])
*(Fradate-resetdates[indexStart])/datespread+paramatrix[indexStart][1];

paramfra[2] = (paramatrix[indexStart+1][2]-paramatrix[indexStart][2])
*(Fradate-resetdates[indexStart])/datespread+paramatrix[indexStart][2];


paramfra[3] = paramatrix[indexStart][3];



paramfra[4] = (paramatrix[indexStart+1][4]-paramatrix[indexStart][4])
*(Fradate-resetdates[indexStart])/datespread+paramatrix[indexStart][4];


paramfra[5] = paramatrix[indexStart][5];

for (j=6;j<=5*(indexStart+1);j++)
{
	paramfra[j] = paramatrix[indexStart+1][j];
}

return err;
}


/* gets the number of fixing days within a period */

Err get_NumFixingDays (int indexStart, double *resetdates, double *nfixdays);

/* Gets the maturity date of each of the fixing days within a period and assigns the value 1 
in the vector isfriday if the day is a friday or a 0 otherwise */


Err get_periodmatandisfriday(int indexStart, double *resetdates, double nfixdays, double *periodmat, double *isfriday);



/* Gets the DRS value for each fixing day within a period */

Err get_DRSts_td(long today, long Paydate, int nfixdays, double *periodmat, 
				 char *cRefRateCode, char *szYieldCurveName,
					double *DRSts)
{
Err err=NULL;
int index;
long FraSpot;
double  num_of_months;
double maturity;
int spotlag;
double dFra;
SrtBasisCode	float_basis;
SrtCompounding  float_compounding;
SrtCrvPtr yldcrv;

yldcrv = lookup_curve(szYieldCurveName);
spotlag = get_spotlag_from_curve(yldcrv);

for (index=1;index<=nfixdays;index++)
{


/* start date of the Fra*/
FraSpot = add_unit((long) periodmat[index],spotlag,SRT_BDAY,MODIFIED_SUCCEEDING);

/*Get details for the index Fra*/
if (err=swp_f_get_ref_rate_details(cRefRateCode, &float_basis, &float_compounding))
{
	smessage("Error in swp_f_get_ref_rate_details");
}




num_of_months = 12/float_compounding;



maturity = coverage(today,(long) periodmat[index],BASIS_ACT_365);



/* compute the fra */
		
		err = get_fra((long) periodmat[index], spotlag, cRefRateCode, szYieldCurveName, &dFra);
	
			DRSts[index] = dFra;
	
	
}

return err;

}



/* Using the Merton diffusion with fixed parameters, computes the value of a DRs corresponding to 
the solution of the SDE it satisfies*/

Err DRSCalc_td(	  double *wn1n2,
				  int nfixdays,
				  double resetmat,
				  double **periodparam,
				  double **mat,
				  int  nperiods,
				  double *DRSts
				  )

{

Err err=NULL;

int n,j;
double U1ave;
double U2ave;
double *sig;
double *weightsmat;

/*memory allocation*/
weightsmat = dvector(1,nperiods);
sig  = dvector(1,nfixdays);


for (j=1;j<=nfixdays;j++)
{
	sig[j] = 0;
}


for(j=1;j<=nfixdays;j++)
{
	U1ave = 0;
	U2ave = 0;
	weightsmat[1] = mat[j][1];
	
	for (n=2;n<=nperiods-1;n++)
	{
		weightsmat[n] = mat[j][n];
	}
	
	weightsmat[nperiods] = mat[j][nperiods]-mat[j][1];
	

	if (j==1)
	{
		for (n=1;n<=nperiods-1;n++)
		{	
			U1ave +=periodparam[j][2]*weightsmat[n]/resetmat;
			U2ave +=periodparam[j][4]*weightsmat[n]/resetmat;
			sig[j]+=periodparam[j][5*(n-1)+1]*periodparam[j][5*(n-1)+1]*weightsmat[n];
		}
	

	DRSts[j]*=exp(-sig[j]/2+sqrt(sig[j])*wn1n2[1]);

	DRSts[j]*=exp(wn1n2[2]*log(1+U1ave)-periodparam[j][3]*U1ave*resetmat)*
						exp(wn1n2[3]*log(1+U2ave)-periodparam[j][5]*U2ave*resetmat);

	}

	else
	{

		for (n=1;n<=nperiods;n++)
	{
		U1ave +=periodparam[j][2]*weightsmat[n]/resetmat;
		U2ave +=periodparam[j][4]*weightsmat[n]/resetmat;
		sig[j]+=periodparam[j][5*(n-1)+1]*periodparam[j][5*(n-1)+1]*weightsmat[n];

	}


	DRSts[j]*=exp(-sig[j]/2+sqrt(sig[j])*wn1n2[1]);

	DRSts[j]*=exp(wn1n2[2]*log(1+U1ave)-periodparam[j][3]*U1ave*resetmat)*
						exp(wn1n2[3]*log(1+U2ave)-periodparam[j][5]*U2ave*resetmat);

	}
}

/*free memory*/
free_dvector(weightsmat,1,nperiods);
free_dvector(sig,1,nfixdays);

return err;
}







/* Main function*/


Err srt_f_resetable_td(
					int    n_resetdates,
					double *resetdates,
					double *initparam,
					char   *FloatCoupon,
					double margin,
					double width,
					double callspread,
					char *ResetType,
					char   *cRefRateCode,
					char   *cMarketId,
					char   *szVolCurveName,
					char   *szYieldCurveName,
					double  **pv
							)
{

Err err=NULL;
int perindex,j,i;
double  nfixdays;
double **paramatrix;
double *periodmat;
double *periodmatinyears;
double *isfriday;
double *DRSts,*DRStsStart;
double **paramperiod;
double **paramaftereset;
double prix;
double eps;
double probajumps, probainf;
int n,n1,n2;
int nherm;
int N1max,N2max;
double *weights;
double *points;
double *wn1n2;
double resetmat;
double ntotaldays;
double *caliberror;
double **mat;
int spotlag;
/* double *discfact;*/
long today;
SrtBasisCode	float_basis;
SrtCompounding  float_compounding;
SrtCrvPtr yldcrv;

yldcrv = lookup_curve(szYieldCurveName);
today = get_today_from_curve(yldcrv);
spotlag = get_spotlag_from_curve(yldcrv);

/*test if resetdates happen after today*/
if (resetdates[1] <= today)
{
	smessage("The reset dates must be in the future!!");
	return err;
}


/*Get details for the index Fra*/
if (err=swp_f_get_ref_rate_details(cRefRateCode, &float_basis, &float_compounding))
{
	smessage("Error in swp_f_get_ref_rate_details");
	return err;
}


probainf = 0.0001;
nherm = 4;

eps = 0;

/*memory allocation*/
paramatrix = dmatrix(1,n_resetdates,1,5*n_resetdates);
weights = dvector(0,2*nherm);
points = dvector(0,2*nherm);
wn1n2 = dvector(1,3);
caliberror = dvector(1,n_resetdates);

/* gets the weights and points form a Hermite Polynomial for the numerical integration made further down*/
err = HermiteStandard(points,weights,2*nherm);
if (err)
{
free_resetablegeneral_td(weights, points, wn1n2, paramatrix, nherm,n_resetdates, caliberror);
return err;
}

/* calibrate the parameters to the caplets smiles of the resetdates and end date */

err = resetablecalib_td(today,n_resetdates, 
						resetdates,
						cMarketId,
						szVolCurveName,
						cRefRateCode,
						szYieldCurveName,
						initparam,
						paramatrix,
						caliberror);
if (err)
{
free_resetablegeneral_td(weights, points, wn1n2, paramatrix, nherm,n_resetdates, caliberror);
return err;
}

for (i=1;i<=n_resetdates;i++)
{
	pv[2][i] = caliberror[i];
}

/* Big loop begins over the number of periods in the resetable */

for (perindex=1;perindex<n_resetdates;perindex++)
{

	err = get_NumFixingDays(perindex, resetdates, &nfixdays);
	if (err)
	{
		free_resetablegeneral_td(weights, points, wn1n2, paramatrix, nherm,n_resetdates, caliberror);
		return err;
	}

	ntotaldays = resetdates[perindex+1]-resetdates[perindex];
	pv[1][perindex] = 0;

/*memory allocation*/
periodmat = dvector(1,(long) nfixdays);
periodmatinyears = dvector(1,(long) nfixdays);
mat = dmatrix(1,(long) nfixdays,1,n_resetdates);
DRStsStart = dvector(1,(long) nfixdays);	
isfriday = dvector(1,(long) nfixdays);	
paramperiod = dmatrix(1,(long) nfixdays,1,5*n_resetdates);
DRSts = dvector(1,(long) nfixdays);	
paramaftereset = dmatrix(1,(long) nfixdays, 1, 5);


	err =  get_periodmatandisfriday(perindex, resetdates, nfixdays, periodmat, isfriday);
	if (err)
	{
		free_resetableperiod_td(paramperiod,periodmat,DRStsStart,DRSts, isfriday,nfixdays, n_resetdates, mat, paramaftereset);
		free_resetablegeneral_td(weights, points, wn1n2, paramatrix, nherm,n_resetdates, caliberror);
		return err;
	}

	
	err =  get_DRSts_td(today, (long) resetdates[perindex+1], (int) nfixdays, periodmat, 
							cRefRateCode, szYieldCurveName,
								DRStsStart)	;																									
	if (err)
	{
		free_resetableperiod_td(paramperiod,periodmat,DRStsStart,DRSts, isfriday,nfixdays, n_resetdates, mat, paramaftereset);
		free_resetablegeneral_td(weights, points, wn1n2, paramatrix, nherm,n_resetdates, caliberror);
		return err;
	}


/* here, we get the parameters corresponding to each fixing day within the period and put them in
 paramperiod which is a matrix with nfixdays rows and 5 columns (for the 5 parameters)*/

		for (j=1;j<=nfixdays;j++)
		{
			err=get_param_td(periodmat[j], perindex, resetdates, 
			  paramatrix, paramperiod[j]);
		}
		
	resetmat = (periodmat[1]-today)/365.0;

/* get the maturities of the different periods for each fixing day*/

		for (j=1;j<=nfixdays;j++)
		{
			mat[j][1] = (periodmat[j]-resetdates[perindex])/365.0;
			mat[j][2] = (resetdates[1]-today)/365.0;
			for (i = 3;i<=n_resetdates;i++)
			{
				mat[j][i] = (resetdates[i-1]-resetdates[i-2])/365.0;
			}
		}


/* maximum values for the number of positive and negative jumps are computed for the resetdate (start date)
of the period in question */

N1max = (int) (paramperiod[1][3]*resetmat+10*sqrt(paramperiod[1][3]*resetmat));
N2max = (int) (paramperiod[1][5]*resetmat+10*sqrt(paramperiod[1][5]*resetmat));


for (n1=0;n1<=N1max;n1++)
{
	for (n2 = 0;n2<=N2max;n2++)
	{

		probajumps = exp(-paramperiod[1][3]*resetmat-paramperiod[1][5]*resetmat+
							n1*log(paramperiod[1][3]*resetmat)+
							n2*log(paramperiod[1][5]*resetmat))/(fact(n1)*fact(n2));


		if (probajumps>probainf)
		{
			
				for (j=1;j<=nfixdays;j++)
				{
					DRSts[j] = DRStsStart[j];
				}


				for (n = 1;n<=2*nherm;n++)
				{

					
				wn1n2[1] = points[n];
				wn1n2[2] = n1;
				wn1n2[3] = n2;
				
				for (j=1;j<=nfixdays;j++)
				{
					DRSts[j] = DRStsStart[j];
				}


				err = DRSCalc_td(wn1n2,(int) nfixdays,resetmat, paramperiod,mat,perindex+1,DRSts);
                
				if (err)
				{
					free_resetableperiod_td(paramperiod,periodmat,DRStsStart, DRSts,isfriday,nfixdays, n_resetdates, mat, paramaftereset);
					free_resetablegeneral_td(weights, points, wn1n2, paramatrix, nherm,n_resetdates, caliberror);
					return err;
				}
				

				for (j=1;j<=nfixdays;j++)
				{
					for (i=1;i<=5;i++)
					{

						paramaftereset[j][i] = paramperiod[j][5*perindex+i];
					}
				}


				
				if (strcmp(ResetType,"CENTERED") == 0)
				{
					for (i=1;i<(int) nfixdays;i++)
					{
						periodmatinyears[i] = (periodmat[i]-periodmat[1])/365.0;
					}
					
					prix = PrixTimeSwap(width/2,DRSts, (int) nfixdays, periodmatinyears, isfriday,paramaftereset,width, callspread);
						
				}
				
				else
				{
				err = NewtonJumps((int) nfixdays,DRSts,periodmat,isfriday,paramaftereset,width, callspread,&eps, &prix);
				}

				if (err)
				{
					free_resetableperiod_td(paramperiod,periodmat,DRStsStart, DRSts,isfriday,nfixdays, n_resetdates, mat, paramaftereset);
					free_resetablegeneral_td(weights, points, wn1n2, paramatrix, nherm,n_resetdates, caliberror);
					return err;
				}
				
				if (strcmp(FloatCoupon,"YES")==0)
				{
					pv[1][perindex]+=weights[n]*probajumps*(DRSts[1]+margin)*prix/ntotaldays;
				}

				else
				{
					pv[1][perindex]+=weights[n]*probajumps*(margin)*prix/ntotaldays;
				}

			}
		}
	}
}

free_resetableperiod_td(paramperiod,periodmat,DRStsStart, DRSts,isfriday,nfixdays, n_resetdates, mat, paramaftereset);
free_dvector(periodmatinyears,1,(long) nfixdays);
}

free_resetablegeneral_td(weights, points, wn1n2, paramatrix, nherm,n_resetdates, caliberror);

return err;
}


