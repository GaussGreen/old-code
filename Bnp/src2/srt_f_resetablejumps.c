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

Err free_resetableperiod(double **paramperiod,double *periodmat, double *DRStsStart,
			   double *DRSts, double *isfriday, double nfixdays, double *driftint)
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

	if (driftint) free_dvector (driftint, 1, (long) nfixdays);
	isfriday = NULL;

	return NULL;
}	

Err free_resetablegeneral(double *weights, double *points, double *wn1n2, 
						  double **paramatrix, int nherm,int n_resetdates)

{
	if (weights) free_dvector (weights, 0, 2*nherm);
	weights = NULL;
	
	if (points) free_dvector (points, 0, 2*nherm);
	points = NULL;
	
	if (wn1n2) free_dvector (wn1n2, 1, 3);
	wn1n2 = NULL;

	if (paramatrix) free_dmatrix (paramatrix, 1, n_resetdates,1,5);
	paramatrix = NULL;

	return NULL;
}




/* Function which calculates the value of a Fra */
/*FraDate is te fixing date*/
Err get_fra(long Fradate, int spotlag, char *cRefRateCode, char *szYieldCurveName, double *dFra)
{
Err err=NULL;
SwapDP  sdpFra;
double num_of_months;
long FraSpot;
long enddate;
long theoenddate;
SrtBasisCode   float_basis;
SrtCompounding  float_compounding;

FraSpot = add_unit(Fradate,spotlag,SRT_BDAY,MODIFIED_SUCCEEDING);

/* get info from the RefRate Code*/
if (err=swp_f_get_ref_rate_details(cRefRateCode, &float_basis, &float_compounding))
{
	smessage("Error if swp_f_get_ref_rate_details");
}

num_of_months = 12/float_compounding;


/* enddate of the fra and theoenddate*/
theoenddate = add_unit(  FraSpot, (int) num_of_months, SRT_MONTH, NO_BUSDAY_CONVENTION);
enddate = add_unit(theoenddate,0, SRT_BDAY, MODIFIED_SUCCEEDING);



/* Compute the Fra */
/* first initialise the sdpFra */
err = swp_f_setSwapDP(FraSpot, theoenddate, float_compounding, float_basis, &sdpFra);
if (err)
	return err;
/* input the spot lag  */
sdpFra.spot_lag = spotlag;

/* computation of the Fra */
err = swp_f_ForwardRate_SwapDP( 
						&sdpFra, 
						szYieldCurveName, 
						cRefRateCode, 
						dFra);
if (err)
	return err;


return err;
}





/* Function which calibrates the parameters of the jump model to a caplet smile
corresponding to a specific matirity date */

Err calibtosmile(  long today,
				   long dDate,
				   char *cMarketId,
				   char *szVolCurveName,
				   char *cRefRateCode,
				   char *szYieldCurveName,
				   double *param,
				   double *dFra,
				   long *ia)
{
Err err=NULL;
int i;
int spotlag;
double atmvol;
double lognorm=1.0;
long  dDatespot,dDatespot3;
double Strikeinf, Strikesup;
double mat;
double Strikespread;
double *strikes;
double *vols;
double *calibvols;
double s;
int  num_of_months;
double chisq;
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

err=get_fra(dDate, spotlag, cRefRateCode, szYieldCurveName, dFra);


dDatespot = add_unit(dDate,spotlag,SRT_BDAY,MODIFIED_SUCCEEDING);


dDatespot3 = add_unit(dDatespot, num_of_months, SRT_MONTH, MODIFIED_SUCCEEDING);


if (err = swp_f_vol(szVolCurveName,(double) dDatespot, (double) dDatespot3,dFra[0],&atmvol,&lognorm))
{
	smessage("Error in swp_f_vol");
}

mat = (dDate-today)/365.0;

Strikeinf = dFra[0]-2*atmvol*dFra[0]*sqrt(mat);
Strikesup = dFra[0]+2*atmvol*dFra[0]*sqrt(mat);
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
err = optcalibmerton(dFra[0],NMAXSTRIKES,strikes,vols,param,mat,"Lognormal",&chisq,ia,calibvols);



/*free memory*/

free_dvector(strikes,1,NMAXSTRIKES+1);
free_dvector(vols,1,NMAXSTRIKES+1);
free_dvector(calibvols,1,NMAXSTRIKES+1);



return err;
}


/* General calibration function which calibrates the parameters to the smiles of the resetdates
of the resetable as well as the smile of the end date of the product 
(going from the longer maturity to the shorter one)*/



Err resetablecalib(long today, int n_resetdates, 
				   double *resetdates,
				   char *cMarketId,
				   char *szVolCurveName,
				   char *cRefRateCode,
				   char *szYieldCurveName,
				   double *paraminitguess,
				   double **paramatrix)

{
Err err=NULL;
int i,j;
long *freeze;
double dFra;
double cvg;
/*memory allocation*/
freeze = lngvector(1,5);


for (i=1; i<=5;i++)
{
	freeze[i] = 1;
	paramatrix[n_resetdates][i] = paraminitguess[i];
}


/* First we calibrate the parameters to the caplet smile whose maturity corresponds to the end date
of the resetable*/

err = calibtosmile( today,(long) resetdates[n_resetdates],
				   cMarketId,
				   szVolCurveName,
				   cRefRateCode,
				   szYieldCurveName,
				   paramatrix[n_resetdates],
				   &dFra,
				   freeze);


/*freeze the intensities of the jumps*/
freeze[3] = 0;
freeze[5] = 0;


/* now calibrates to the other smiles */


for (j = n_resetdates-1; j>=1;j--)
{
		cvg = (resetdates[j+1]-resetdates[j])/365;

		paramatrix[j][1] = paraminitguess[1];	
		paramatrix[j][2] = paraminitguess[2];
		paramatrix[j][3] = paramatrix[j+1][3]/(1+paramatrix[j+1][2]*cvg*dFra/(1+cvg*dFra));
		paramatrix[j][4] = paraminitguess[4];
		paramatrix[j][5] = paramatrix[j+1][5]/(1+paramatrix[j+1][4]*cvg*dFra/(1+cvg*dFra));



err = calibtosmile(today, (long) resetdates[j],cMarketId,
				   szVolCurveName,
				   cRefRateCode,
				   szYieldCurveName,
				   paramatrix[j],
				   &dFra,
				   freeze);


}

/* free memory */

free_lngvector(freeze,1,5);

return err;

}


/* function which linearly interpolates the parameters obtained from the smiles of the beginning
and end of the period , for a maturity in between*/

Err get_param(double Fradate, int  indexStart, double *resetdates, 
			  double **paramatrix, double *paramfra)
{
Err err=NULL;
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
return err;
}


/* gets the number of fixing days within a period */

Err get_NumFixingDays (int indexStart, double *resetdates, double *nfixdays)
{
Err err=NULL;
long d;

nfixdays[0] = 0;

d = (long) resetdates[indexStart];

while (((double) d) < resetdates[indexStart+1])
{
   if (week_day(d) == FRIDAY)
	{
		d = d+3;
	}

	else 
	{
		d++;
	}

	nfixdays[0] = nfixdays[0]+1;
}

return err;
}


/* Gets the maturity date of each of the fixing days within a period and assigns the value 1 
in the vector isfriday if the day is a friday or a 0 otherwise */


Err get_periodmatandisfriday(int indexStart, double *resetdates, double nfixdays, double *periodmat, double *isfriday)
{
Err err=NULL;
double d;
int j,index;
double firstfriday=0;
double NumFixingDays = nfixdays;
d = resetdates[indexStart];



/*Getting the first friday*/

for (j=0;j<5;j++)
{
	if (week_day((long)(d+j)) == FRIDAY)
	{
		firstfriday = d+j;
	}

}


/* fills in periodmat */

index = 1;

while (d<resetdates[indexStart+1])
{


	periodmat[index] = d;

	if (week_day((long)d) == FRIDAY)
	{
		d = d+3;
	}

	else 
	{
		d++;
	}
	
	index++;
}


/*Initialization of the vectors*/

for (j = 1;j<= NumFixingDays;j++)
{
	isfriday[j]=0;
}


/*fill in isfriday*/

j=  (int) (firstfriday-resetdates[indexStart]+1);

while (j <= NumFixingDays)
{
	
	isfriday[j] = 1;

	j=j+5;
}

return err;
}



/* Gets the DRS value for each fixing day within a period */

Err get_DRSts(long today, long Paydate, int nfixdays, 
			  double *periodmat, char *cRefRateCode, 
			  char *szYieldCurveName,double *DRSts)
{
Err err=NULL;
int index;
long FraSpot;
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


/* compute the fra */
		
		err = get_fra((long) periodmat[index], spotlag, cRefRateCode, szYieldCurveName, &dFra);
	
		DRSts[index] = dFra;

}

return err;

}



/* Using the Merton diffusion with fixed parameters, computes the value of a DRs corresponding to 
the solution of the SDE it satisfies*/

Err DRSCalc(	  double *wn1n2,
				  int nfixdays,
				  double resetmat,
				  double **periodparam,
				  double *drift,
				  double *DRSts
				  )

{

Err err=NULL;

int j;
double T;
double *paramfra;
double sig,U1,U2,l1,l2;

/*memory allocation*/
paramfra = dvector(1,5);


for(j=1;j<=nfixdays;j++)
{

	T=resetmat;

	sig = periodparam[j][1];
	U1  = periodparam[j][2];
	l1  = periodparam[j][3];
	U2  = periodparam[j][4];
	l2  = periodparam[j][5];
	
	

	DRSts[j] = DRSts[j]*exp(-sig*sig*T/2+sig*wn1n2[1]*sqrt(T))*
						exp(wn1n2[2]*log(1+U1)-l1*U1*T)*
						exp(wn1n2[3]*log(1+U2)-l2*U2*T)*exp(drift[j]*T);

}

/*free memory*/
free_dvector(paramfra,1,5);


return err;
}







/* Main function*/



Err srt_f_resetable(
					int    n_resetdates,
					double *resetdates,
					double *initparam,
					char   *FloatCoupon,
					double margin,
					double width,
					char   *cRefRateCode,
					char   *cMarketId,
					char   *szYieldCurveName,
					char   *szVolCurveName,
					double  *pv
							)
{

Err err=NULL;
int perindex,j;
double  nfixdays;
double **paramatrix;
double *periodmat;
double *isfriday;
double *DRSts,*DRStsStart;
double **paramperiod;
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
double drift;
double *driftinterp;
double cvg;
double dFradrift;
int spotlag;
double callspread;
/* double *discfact;*/
long today;
SrtBasisCode	float_basis;
SrtCompounding  float_compounding;
SrtCrvPtr yldcrv;

yldcrv = lookup_curve(szYieldCurveName);
today = get_today_from_curve(yldcrv);
spotlag = get_spotlag_from_curve(yldcrv);


callspread = 0.001;
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
}


probainf = 0.0001;
nherm = 4;

eps = 0;

/*memory allocation*/
paramatrix = dmatrix(1,n_resetdates,1,5);
weights = dvector(0,2*nherm);
points = dvector(0,2*nherm);
wn1n2 = dvector(1,3);
/*disfact = dvector(1,n_resetdates-1);*/



/* gets the weights and points form a Hermite Polynomial for the numerical integration made further down*/
err = HermiteStandard(points,weights,2*nherm);
if (err)
{
free_resetablegeneral(weights, points, wn1n2, paramatrix, nherm,n_resetdates);
return err;
}

/* calibrate the parameters to the caplets smiles of the resetdates and end date */

err = resetablecalib(today,n_resetdates, 
						resetdates,
						cMarketId,
						szVolCurveName,
						cRefRateCode,
						szYieldCurveName,
						initparam,
						paramatrix);
if (err)
{
free_resetablegeneral(weights, points, wn1n2, paramatrix, nherm,n_resetdates);
return err;
}


/* Big loop begins over the number of periods in the resetable */

for (perindex=1;perindex<n_resetdates;perindex++)
{

	err = get_NumFixingDays(perindex, resetdates, &nfixdays);
	if (err)
	{
		free_resetablegeneral(weights, points, wn1n2, paramatrix, nherm,n_resetdates);
		return err;
	}

	ntotaldays = resetdates[perindex+1]-resetdates[perindex];
	pv[perindex] = 0;

/*memory allocation*/
periodmat = dvector(1,(long) nfixdays);
DRStsStart = dvector(1,(long) nfixdays);	
isfriday = dvector(1,(long) nfixdays);	
paramperiod = dmatrix(1,(long) nfixdays,1,5);
DRSts = dvector(1,(long) nfixdays);	
driftinterp = dvector(1, (long) nfixdays);

	err =  get_periodmatandisfriday(perindex, resetdates, nfixdays, periodmat, isfriday);
	if (err)
	{
		free_resetableperiod(paramperiod,periodmat,DRStsStart,DRSts, isfriday,nfixdays, driftinterp);
		free_resetablegeneral(weights, points, wn1n2, paramatrix, nherm,n_resetdates);
		return err;
	}

	
	err =  get_DRSts(today, (long) resetdates[perindex+1], (int) nfixdays, periodmat, cRefRateCode, szYieldCurveName,DRStsStart)	;																									
	if (err)
	{
		free_resetableperiod(paramperiod,periodmat,DRStsStart,DRSts, isfriday,nfixdays, driftinterp);
		free_resetablegeneral(weights, points, wn1n2, paramatrix, nherm,n_resetdates);
		return err;
	}


cvg = coverage((long) resetdates[perindex+1], add_unit((long) resetdates[perindex+1], (int) (12/float_compounding), SRT_MONTH, MODIFIED_SUCCEEDING),float_basis);
/* here, we get the parameters corresponding to each fixing day within the period and put them in
 paramperiod which is a matrix with nfixdays rows and 5 columns (for the 5 parameters)*/


		err=get_fra((long) resetdates[perindex+1], spotlag, cRefRateCode, szYieldCurveName, &dFradrift);

		drift = paramatrix[perindex][1]*cvg*dFradrift/(1+cvg*dFradrift);
		
		for (j=1;j<=nfixdays;j++)
		{
		err=get_param(periodmat[j], perindex, resetdates, 
			  paramatrix, paramperiod[j]);

		driftinterp[j] = drift*(periodmat[j]-resetdates[perindex])/(resetdates[perindex+1]-resetdates[perindex]);

			if (err)
			{
			free_resetableperiod(paramperiod,periodmat,DRStsStart,DRSts, isfriday,nfixdays, driftinterp);
			free_resetablegeneral(weights, points, wn1n2, paramatrix, nherm,n_resetdates);
			return err;
			}
		}
		
	

resetmat = (periodmat[1]-today)/365.0;

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

			for (n = 1;n<=2*nherm;n++)
			{

					
				wn1n2[1] = points[n];
				wn1n2[2] = n1;
				wn1n2[3] = n2;
				
				for (j=1;j<=nfixdays;j++)
				{
					DRSts[j] = DRStsStart[j];
				}


				err = DRSCalc(wn1n2,(int) nfixdays,resetmat,paramperiod,driftinterp,DRSts);
                
				if (err)
				{
					free_resetableperiod(paramperiod,periodmat,DRStsStart, DRSts,isfriday,nfixdays, driftinterp);
					free_resetablegeneral(weights, points, wn1n2, paramatrix, nherm,n_resetdates);
					return err;
				}
				
				
				
				
				
				err = NewtonJumps((int) nfixdays,DRSts,periodmat,isfriday,paramperiod,width,callspread, &eps, &prix);
			
				if (err)
				{
					free_resetableperiod(paramperiod,periodmat,DRStsStart, DRSts,isfriday,nfixdays, driftinterp);
					free_resetablegeneral(weights, points, wn1n2, paramatrix, nherm,n_resetdates);
					return err;
				}
				
				if (strcmp(FloatCoupon,"YES")==0)
				{
					pv[perindex]+=weights[n]*probajumps*(DRSts[1]+margin)*(-prix)/ntotaldays;
				}

				else
				{
					pv[perindex]+=weights[n]*probajumps*(margin)*(-prix)/ntotaldays;
				}

			}
		}
	}
}

free_resetableperiod(paramperiod,periodmat,DRStsStart, DRSts,isfriday,nfixdays, driftinterp);
				
/*
discfact[perindex] = swp_f_df(today,resetdates[perindex+1],szYieldCurveName);

resetprice +=disfact[perindex]*price[perindex];
*/

}

free_resetablegeneral(weights, points, wn1n2, paramatrix, nherm,n_resetdates);

return err;
}


Err get_atm_vols(int nfixings,
				 double *DRSts, 
				 double *Spotdates,
				 double resetmat,
				 char *cRefRateCode, 
				 char *szVolCurveName, 
				 char *szYieldCurveName,
				 double *atmvols)
{
Err err = NULL;
int i;
long EndDate;
double lognorm = 1.0;
SrtBasisCode float_basis;
SrtCompounding float_compounding;
SrtCrvPtr yldcrv;
long today, spotlag, spotdate;
yldcrv = lookup_curve(szYieldCurveName);
today = get_today_from_curve(yldcrv);
spotlag = get_spotlag_from_curve(yldcrv);
spotdate = add_unit(today, spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);



/*Get details for the index Fra*/
err=swp_f_get_ref_rate_details(cRefRateCode, &float_basis, &float_compounding);
if (err)
{
	return err;
}

	for (i=1;i<=nfixings;i++)
	{
		EndDate = add_unit((long) Spotdates[i],12/float_compounding, SRT_MONTH,MODIFIED_SUCCEEDING);
		err = swp_f_vol(szVolCurveName, (double) Spotdates[i], (double) EndDate, DRSts[i], &(atmvols[i]),&lognorm);		
	}

return err;

}

Err get_atm_normal_vols(int nfixings,
				 double *DRSts, 
				 double *Spotdates,
				 char *cRefRateCode, 
				 char *szVolCurveName, 
				 double *atmvols)
{
Err err = NULL;
int i;
long EndDate;
double lognorm = 1.0;
SrtBasisCode float_basis;
SrtCompounding float_compounding;
/*Get details for the index Fra*/
err=swp_f_get_ref_rate_details(cRefRateCode, &float_basis, &float_compounding);
if (err)
{
	return err;
}

	for (i=1;i<=nfixings;i++)
	{
		EndDate = add_unit((long) Spotdates[i],12/float_compounding, SRT_MONTH,MODIFIED_SUCCEEDING);
		err = swp_f_vol(szVolCurveName, (double) Spotdates[i], (double) EndDate, DRSts[i], &(atmvols[i]),&lognorm);	
		atmvols[i]*=DRSts[i];
	}

return err;

}



Err SimulateDRS(double w,
				int nfixings, 
				double resetmat, 
				double *vols, 
				double *DRSts)
{
Err err = NULL;
int i;


	for (i = 1;i<=nfixings;i++)
	{
	DRSts[i] = DRSts[i]*exp(vols[i]*w*sqrt(resetmat)-vols[i]*vols[i]*resetmat/2);
	}

return err;
}


double Compute_Price(int nfixings, 
					 double *DRSts,
					 double width,
					 double lowerbarrier, 
					 double callspread, 
					 long spotdate, 
					 double *periodmat, 
					 double *isfriday,
					 char *FloatCoupon,
					 double margin,
					 SrtCompounding frequency,
					 SrtCompounding ResetFrequency,
					 char *szVolCurveName,
					 char *szYieldCurveName,
					 char *cRefRateCode,
					 char *VolType, 
					 double *VolsInputs, 
					 double VolShift)
{
int i,j;
double *corridorvols;
double lognorm;
long EndDate;
long periodspot;
double price;
double frequencyratio;
int paymentindex;
double call_lower, call_lower_eps, call_upper, call_upper_eps;
double *cvgs;
long *dates;
SrtBasisCode float_basis;
SrtCompounding float_compounding;
int previouspaymentindex = 1;
int drsindex;
SrtCrvPtr yldcrv;
long startdate;
long enddate;
double ATMVol,NormalVol;
double Fwd;
long today;
int spotlag;
double stdevactual;
double stdevoriginal;
yldcrv = lookup_curve(szYieldCurveName);
today = get_today_from_curve(yldcrv);
spotlag = get_spotlag_from_curve(yldcrv);


/*Get details for the index Fra*/
swp_f_get_ref_rate_details(cRefRateCode, &float_basis, &float_compounding);



frequencyratio = (double) (frequency/ResetFrequency);

/*Memory allocation*/
corridorvols = dvector(1,4);
cvgs = dvector(1,(long) frequencyratio);
dates = lngvector(1,(long) frequencyratio+1);


dates[1] = (long) periodmat[1];

for (i = 1;i<frequencyratio+1;i++)
{
dates[i+1] = add_unit(dates[i], 12/frequency, SRT_MONTH, MODIFIED_SUCCEEDING);
cvgs[i] = coverage(dates[i], dates[i+1], float_basis);
}

price = 0.0;
drsindex = 1;

for (i = 1;i<=nfixings;i++)
		{	
			periodspot =	add_unit((long) periodmat[i],2, SRT_BDAY,MODIFIED_SUCCEEDING); 
		
			EndDate = add_unit(periodspot,12/frequency, SRT_MONTH,MODIFIED_SUCCEEDING);
			
			if (strcmp(VolType,"Flat") == 0)
			{

				swp_f_vol(szVolCurveName,periodspot,EndDate,lowerbarrier-callspread, &(corridorvols[1]),&lognorm);  
			
				swp_f_vol(szVolCurveName,periodspot,EndDate,lowerbarrier , &(corridorvols[2]),&lognorm);
			
				swp_f_vol(szVolCurveName,periodspot,EndDate,lowerbarrier+width, &(corridorvols[3]),&lognorm);
			
				swp_f_vol(szVolCurveName,periodspot,EndDate,lowerbarrier+width+callspread, &(corridorvols[4]),&lognorm);
			}

			else if (strcmp(VolType,"Input") == 0)
			{
				corridorvols[1] = VolsInputs[1];
				corridorvols[2] = VolsInputs[2];
				corridorvols[3] = VolsInputs[3];
				corridorvols[4] = VolsInputs[4];
			}


			else if (strcmp(VolType, "LogSliding") == 0)

			{
				startdate = (long) (spotdate+periodmat[i]-periodmat[1]);
				enddate = add_unit((long) (spotdate+periodmat[i]-periodmat[1]),12/frequency,SRT_MONTH, MODIFIED_SUCCEEDING);

				get_fra(add_unit(startdate,-spotlag,SRT_BDAY,MODIFIED_SUCCEEDING),spotlag,cRefRateCode, szYieldCurveName,&Fwd);

				swp_f_vol(szVolCurveName,startdate,enddate,Fwd, &(ATMVol),&lognorm);  
				
				stdevoriginal = ATMVol*Fwd*sqrt((add_unit(startdate,-spotlag,SRT_BDAY,MODIFIED_SUCCEEDING)-today)/365.25);
				stdevactual = ATMVol*DRSts[i]*sqrt((add_unit(startdate,-spotlag,SRT_BDAY,MODIFIED_SUCCEEDING)-today)/365.25);

				swp_f_vol(szVolCurveName,startdate,enddate,((lowerbarrier-callspread)-DRSts[i])*stdevoriginal/stdevactual+Fwd, &(corridorvols[1]),&lognorm);  
			
				swp_f_vol(szVolCurveName,startdate,enddate,(lowerbarrier-DRSts[i])*stdevoriginal/stdevactual+Fwd, &(corridorvols[2]),&lognorm);
			
				swp_f_vol(szVolCurveName,startdate,enddate,((lowerbarrier+width)-DRSts[i])*stdevoriginal/stdevactual+Fwd, &(corridorvols[3]),&lognorm);
			
				swp_f_vol(szVolCurveName,startdate,enddate,((lowerbarrier+width+callspread)-DRSts[i])*stdevoriginal/stdevactual+Fwd, &(corridorvols[4]),&lognorm);
			}

			else /* NormSliding*/
			{
				startdate = (long) (spotdate+periodmat[i]-periodmat[1]);
				enddate = add_unit((long) (spotdate+periodmat[i]-periodmat[1]),12/frequency,SRT_MONTH, MODIFIED_SUCCEEDING);

				get_fra( add_unit(startdate,-spotlag,SRT_BDAY,MODIFIED_SUCCEEDING),spotlag,cRefRateCode, szYieldCurveName,&Fwd);

				swp_f_vol(szVolCurveName,startdate,enddate,Fwd, &(ATMVol),&lognorm);  
				
				NormalVol = ATMVol*Fwd;

				swp_f_vol(szVolCurveName,startdate,enddate,((lowerbarrier-callspread)-DRSts[i])+Fwd, &(corridorvols[1]),&lognorm);  

				corridorvols[1] = NormalVol*corridorvols[1]/(DRSts[i]*ATMVol);
			
				swp_f_vol(szVolCurveName,startdate,enddate,(lowerbarrier-DRSts[i])+Fwd, &(corridorvols[2]),&lognorm);
			
				corridorvols[2] = NormalVol*corridorvols[2]/(DRSts[i]*ATMVol);

				swp_f_vol(szVolCurveName,startdate,enddate,((lowerbarrier+width)-DRSts[i])+Fwd, &(corridorvols[3]),&lognorm);
			
				corridorvols[3] = NormalVol*corridorvols[3]/(DRSts[i]*ATMVol);

				swp_f_vol(szVolCurveName,startdate,enddate,((lowerbarrier+width+callspread)-DRSts[i])+Fwd, &(corridorvols[4]),&lognorm);

				corridorvols[4] = NormalVol*corridorvols[4]/(DRSts[i]*ATMVol);
			}
				corridorvols[1] += VolShift;
				corridorvols[2] -= VolShift;
				corridorvols[3] -= VolShift;
				corridorvols[4] += VolShift;

			call_lower_eps = srt_f_optblksch(DRSts[i], lowerbarrier-callspread, corridorvols[1], (periodmat[i]-periodmat[1])/365.25, 1.0, SRT_CALL, PREMIUM);

			call_lower = srt_f_optblksch(DRSts[i], lowerbarrier, corridorvols[2], (periodmat[i]-periodmat[1])/365.25, 1.0, SRT_CALL, PREMIUM);

			call_upper = srt_f_optblksch(DRSts[i], lowerbarrier+width, corridorvols[3], (periodmat[i]-periodmat[1])/365.25, 1.0, SRT_CALL, PREMIUM);

			call_upper_eps = srt_f_optblksch(DRSts[i], lowerbarrier+width+callspread, corridorvols[4], (periodmat[i]-periodmat[1])/365.25, 1.0, SRT_CALL, PREMIUM);
			
			paymentindex = 1;
			
			for (j = 2;j<frequencyratio+1;j++)
			{
				paymentindex +=(int) (periodmat[i]/dates[j]);
			}
			
			
			if (strcmp(FloatCoupon,"YES")==0)
			{
				if (paymentindex > previouspaymentindex)
				{
					drsindex = i;
				}
				previouspaymentindex = paymentindex;

				price+=frequencyratio*cvgs[paymentindex]*(DRSts[drsindex]+margin)*(1+2*isfriday[i])*((call_lower_eps-call_lower)/callspread-(call_upper-call_upper_eps)/callspread);						
			}
			else
			{
				price+=frequencyratio*cvgs[paymentindex]*margin*(1+2*isfriday[i])*((call_lower_eps-call_lower)/callspread-(call_upper-call_upper_eps)/callspread);
			}
		}

/*Free Memory*/
free_dvector(corridorvols,1,4);
free_lngvector(dates,1,(long) frequencyratio+1);
free_dvector(cvgs,1,(long) frequencyratio);

return price;

}


double Optimise_Time_Swap_Price(int nfixings, 
					 double *DRSts,
					 double width,
					 double lowermin,
					 double lowermax,
					 double *lowerbarrier,
					 double callspread, 
					 long spodate, 
					 double *periodmat, 
					 double *isfriday,
					 char *FloatCoupon,
					 double margin,
					 SrtCompounding frequency,
					 SrtCompounding ResetFrequency, 
					 char *szVolCurveName, 
					 char *szYieldCurveName,
					 char *cRefRateCode,
					 char *VolType, 
					 double *VolsInputs,
					 double VolShift)
{
int iter = 0;
double pleft, pright, pmid;
double xleft, xright, xmid;
double xmidleft,xmidright;
double pmidleft,pmidright;
double precision;
int maxiter = 15;

xleft = lowermin;
xright = lowermax;
xmid = (lowermin+lowermax)/2;
precision = 1.0;
pleft = Compute_Price(nfixings, DRSts,width,xleft, callspread, spodate, periodmat,isfriday, FloatCoupon, margin, frequency,ResetFrequency,szVolCurveName,szYieldCurveName,cRefRateCode, VolType, VolsInputs, VolShift);
pright = Compute_Price(nfixings, DRSts,width,xright, callspread, spodate, periodmat, isfriday,FloatCoupon, margin, frequency,ResetFrequency,szVolCurveName,szYieldCurveName,cRefRateCode, VolType, VolsInputs, VolShift);




	while ((iter < maxiter)&&(precision>=0.0001))
	{

		pmid = Compute_Price(nfixings, DRSts,width,xmid, callspread, spodate, periodmat, isfriday,FloatCoupon, margin, frequency,ResetFrequency,szVolCurveName,szYieldCurveName,cRefRateCode, VolType, VolsInputs, VolShift);

		xmidleft = (xleft+xmid)/2;
		xmidright = (xright+xmid)/2;

		pmidleft = Compute_Price(nfixings, DRSts,width,xmidleft, callspread, spodate, periodmat, isfriday, FloatCoupon, margin, frequency,ResetFrequency, szVolCurveName,szYieldCurveName,cRefRateCode, VolType, VolsInputs, VolShift);
		pmidright = Compute_Price(nfixings, DRSts,width,xmidright, callspread, spodate, periodmat,isfriday, FloatCoupon, margin, frequency, ResetFrequency, szVolCurveName,szYieldCurveName,cRefRateCode, VolType, VolsInputs, VolShift);
		
		if ((pmid>pmidleft)&&(pmid>pmidright))
		{
			xleft = xmidleft;
			xright = xmidright;
		}
		else if ((pmidleft>pmid)&&(pmidleft>pmidright))
		{
			xright = xmid;
		}

		else
		{
			xleft = xmid;
		}

		precision = fabs(pmidleft-pmidright);
		xmid = (xleft+xright)/2;
		iter++;

		
	}

return pmid;

}

double  Price_Time_Swap(int nfixings,
						double *DRSts, 
						double width,
						double callspread, 
						long spotdate, 
						double *periodmat, 
						double *isfriday,
						char *FloatCoupon,
						double margin,
						SrtCompounding float_frequency,
						SrtCompounding ResetFrequency,
						char *szVolCurveName, 
						char *szYieldCurveName,
						char *cRefRateCode,
						char *ResetType, 
						char *VolType, 
						double *VolsInputs,
						double VolShift,
						double lowerbarrier)/*returns the optimal barrier if OPTIMISED or contains the spread to lower barrier if input*/
{
double price;
double lowermin,lowermax;
double eps;

	lowermin = DRSts[1]-2*width;
	lowermax = DRSts[1]+width;

	if (strcmp(ResetType,"NORESET") == 0)
	{
		price = Compute_Price(nfixings, DRSts,width,lowerbarrier, callspread, spotdate, periodmat, isfriday, FloatCoupon, margin, float_frequency, ResetFrequency,szVolCurveName,szYieldCurveName, cRefRateCode,VolType, VolsInputs, VolShift);
	}

	else if (strcmp(ResetType,"CENTERED") == 0)
	{
		price = Compute_Price(nfixings, DRSts,width,DRSts[1]-width/2, callspread, spotdate, periodmat, isfriday, FloatCoupon,margin,float_frequency, ResetFrequency,szVolCurveName,szYieldCurveName, cRefRateCode,VolType, VolsInputs, VolShift);
	}
	
	else
	{
		price = Optimise_Time_Swap_Price(nfixings, DRSts,width,lowermin, lowermax, &eps, callspread, spotdate, periodmat,isfriday, FloatCoupon, margin, float_frequency, ResetFrequency, szVolCurveName,szYieldCurveName, cRefRateCode,VolType, VolsInputs, VolShift);
	}


return price;

}








Err srt_f_simpleresetable(
					double StartDate,
					double  EndDate,
					char   *VolType,/*Flat or LogSliding or Input or NormSliding (that's the default)*/
					double *VolsInput,/*If above is "Input", this contains 4 vols:lowercallspread, lower, upper, uppercallspread*/
					char   *FloatCoupon,
					SrtCompounding ResetFrequency,
					double margin,
					double width,
					double callspread,
					double VolShift,
					char   *ResetType,
					double lowerbarrier,
					char   *cRefRateCode,
					char   *cMarketId,
					char   *szVolCurveName,
					char   *szYieldCurveName,
					double  *pv)
{

Err err=NULL;
int j;
double  nfixdays;
double *periodmat;
double *isfriday;
double *DRSts,*DRStsStart;
int n;
int nherm;
double *weights;
double *points;
double resetmat;
double ntotaldays;
double *atmvols = NULL;
int spotlag;
double w;
double price;
double *dates;

/* double *discfact;*/
long today, spotdate;
SrtBasisCode	float_basis;
SrtCompounding  float_compounding;
SrtCrvPtr yldcrv;

yldcrv = lookup_curve(szYieldCurveName);
today = get_today_from_curve(yldcrv);
spotlag = get_spotlag_from_curve(yldcrv);
spotdate = add_unit(today, spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);


/*Get details for the index Fra*/
if (err=swp_f_get_ref_rate_details(cRefRateCode, &float_basis, &float_compounding))
{
	return err;
}


nherm = 6;



/*memory allocation*/
weights = dvector(0,2*nherm);
points = dvector(0,2*nherm);
dates = dvector(0,1);


	/* gets the weights and points from a Hermite Polynomial for the numerical integration made further down*/
	err = HermiteStandard(points,weights,2*nherm);
	
	dates[0] = StartDate;
	dates[1] = EndDate;

	err = get_NumFixingDays(0, dates, &nfixdays);
	
	ntotaldays = EndDate-StartDate;
	*pv = 0.0;

/*memory allocation*/
atmvols = dvector(1,(long) nfixdays);
periodmat = dvector(1,(long) nfixdays);
DRStsStart = dvector(1,(long) nfixdays);	
isfriday = dvector(1,(long) nfixdays);	
DRSts = dvector(1,(long) nfixdays);	


	/*Gets the dates of fixing within the period and whether or not the date is a friday or not*/
	err =  get_periodmatandisfriday(0, dates, nfixdays, periodmat, isfriday);
	
	
	err =  get_DRSts(today, (long) EndDate, (int) nfixdays, periodmat, 
							cRefRateCode, szYieldCurveName,
								DRStsStart)	;																									

	resetmat = (periodmat[1]-today)/365.25;

    err = get_atm_vols((int) nfixdays,DRStsStart, periodmat,resetmat, cRefRateCode, szVolCurveName, szYieldCurveName,atmvols);

		

	

		for (n = 1;n<=2*nherm;n++)
		{
			w = points[n];
				
				for (j=1;j<=nfixdays;j++)
				{
					DRSts[j] = DRStsStart[j];
				}

				err = SimulateDRS(w,(int) nfixdays,resetmat, atmvols, DRSts);
				price = Price_Time_Swap((int) nfixdays,DRSts, width, callspread, spotdate, periodmat, isfriday,FloatCoupon, margin,float_compounding,ResetFrequency,szVolCurveName, szYieldCurveName,cRefRateCode,ResetType, VolType, VolsInput, VolShift, lowerbarrier);
						
				*pv+=weights[n]*price/ntotaldays;
				
		}
		
/*free memory*/
free_dvector(atmvols,1,(long) nfixdays);
free_dvector(periodmat,1,(long) nfixdays);
free_dvector(DRStsStart,1,(long) nfixdays);	
free_dvector(isfriday,1,(long) nfixdays);	
free_dvector(DRSts,1,(long) nfixdays);	
free_dvector(weights,0,2*nherm);
free_dvector(points,0,2*nherm);
free_dvector(dates,0,1);


return err;
}



double Price_Caplet(long StrikeDate,
					long StartDate, 
					long EndDate, 
					double dFra, 
					double dStrike,
					char *szVolCurveName, 
					char *szYieldCurveName, 
					char *cRefRateCode,
					char *VolType, 
					double VolInput,
					SrtCallPutType CallPut)
{
Err err = NULL;
double price;
long today, spotdate;
int spotlag;
double flatvol;
double forwardvol;
double atmvol;
double mat;
double lognorm = 1.0;
double dfirstfra;
double stdevoriginal,stdevactual;
long start,end,fixing;
SrtBasisCode	float_basis;
SrtCompounding  float_compounding;
SrtCrvPtr yldcrv;

err=swp_f_get_ref_rate_details(cRefRateCode, &float_basis, &float_compounding);
yldcrv = lookup_curve(szYieldCurveName);
today = get_today_from_curve(yldcrv);
spotlag = get_spotlag_from_curve(yldcrv);
spotdate = add_unit(today, spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);

mat = (StartDate-StrikeDate)/365.25;

if (strcmp(VolType,"Input")==0)
{
	if (CallPut == SRT_STRADDLE)
		price = srt_f_optblksch(dFra, dStrike, VolInput, mat, 1.0,SRT_CALL , PREMIUM)
			+ srt_f_optblksch(dFra, dStrike, VolInput, mat, 1.0,SRT_PUT , PREMIUM);
	else 
	{
		price = srt_f_optblksch(dFra, dStrike, VolInput, mat, 1.0,CallPut , PREMIUM);
	}
}

else if (strcmp(VolType,"Flat")==0)
{
	err = swp_f_vol(szVolCurveName, StartDate, EndDate,dStrike, &flatvol,&lognorm);

	if (CallPut == SRT_STRADDLE)
		price = srt_f_optblksch(dFra, dStrike, flatvol, mat, 1.0,SRT_CALL , PREMIUM)
			+ srt_f_optblksch(dFra, dStrike, flatvol, mat, 1.0,SRT_PUT , PREMIUM);
	else price = srt_f_optblksch(dFra, dStrike, flatvol, mat, 1.0,CallPut , PREMIUM);
}
else if (strcmp(VolType,"LogSliding")==0)
{
	start = add_unit(spotdate, 12/float_compounding, SRT_MONTH, MODIFIED_SUCCEEDING);
	end  = add_unit(start, 12/float_compounding, SRT_MONTH, MODIFIED_SUCCEEDING);
	fixing = add_unit(start, -spotlag, SRT_BDAY,MODIFIED_SUCCEEDING);
	err = get_fra(fixing, spotlag,  cRefRateCode, szYieldCurveName, &dfirstfra);
	err = swp_f_vol(szVolCurveName, start, end,dfirstfra, &atmvol,&lognorm);
	stdevoriginal = dfirstfra*atmvol*sqrt(mat);
	stdevactual = dFra*atmvol*sqrt(mat);
	
	err = swp_f_vol(szVolCurveName, start, end, dfirstfra+(dStrike-dFra)*stdevoriginal/stdevactual, &forwardvol,&lognorm);

	if (CallPut == SRT_STRADDLE)
		price = srt_f_optblksch(dFra, dStrike,forwardvol , mat, 1.0,SRT_CALL , PREMIUM)
			+ srt_f_optblksch(dFra, dStrike, forwardvol, mat, 1.0,SRT_PUT , PREMIUM);
	else price = srt_f_optblksch(dFra, dStrike, forwardvol, mat, 1.0,CallPut , PREMIUM);
	
}
else
{
	start = add_unit(spotdate, 12/float_compounding, SRT_MONTH, MODIFIED_SUCCEEDING);
	end  = add_unit(start, 12/float_compounding, SRT_MONTH, MODIFIED_SUCCEEDING);
	fixing = add_unit(start, -spotlag, SRT_BDAY,MODIFIED_SUCCEEDING);
	err = get_fra(fixing, spotlag,  cRefRateCode, szYieldCurveName, &dfirstfra);

	err = swp_f_vol(szVolCurveName, start, end, dfirstfra+(dStrike-dFra), &forwardvol,&lognorm);

	if (CallPut == SRT_STRADDLE)
		price = srt_f_optblksch(dFra, dStrike,forwardvol , mat, 1.0,SRT_CALL , PREMIUM)
			+ srt_f_optblksch(dFra, dStrike, forwardvol, mat, 1.0,SRT_PUT , PREMIUM);
	else price = srt_f_optblksch(dFra, dStrike, forwardvol, mat, 1.0,CallPut , PREMIUM);
}


return price;
}



Err srt_f_simpleresetcap(
					long StrikeDate,
					long StartDate,
					long  EndDate,
					char   *VolType,/*Flat or NormSliding or LogSliding or Input*/
					double VolInput,/*Contains a vol input buy the user used only if above is input*/
					char   *cRefRateCode,
					char   *cMarketId,
					char   *szVolCurveName,
					char   *szYieldCurveName,
					SrtCallPutType CallPut,
					double  *pv)
{

Err err=NULL;
int n;
int nherm;
double lognorm = 1.0;
double *weights;
double *points;
double resetmat;
double atmvolforward, atmvolstrike;
int spotlag;
double w;
long Fixing;
long StrikeFixing;
long today, spotdate;
double dForward;
double price;
double dFwd,Strike, StrikeForward;
SrtBasisCode	float_basis;
SrtCompounding  float_compounding;
SrtCrvPtr yldcrv;

yldcrv = lookup_curve(szYieldCurveName);
today = get_today_from_curve(yldcrv);
spotlag = get_spotlag_from_curve(yldcrv);
spotdate = add_unit(today, spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);


/*Get details for the index Fra*/
err=swp_f_get_ref_rate_details(cRefRateCode, &float_basis, &float_compounding);
if (err)
{
	return err;
}


nherm = 4;



/*memory allocation*/
weights = dvector(0,2*nherm);
points = dvector(0,2*nherm);


	/* gets the weights and points from a Hermite Polynomial for the numerical integration made further down*/
	err = HermiteStandard(points,weights,2*nherm);
	
	*pv = 0.0;



	StrikeFixing = add_unit(StrikeDate, -spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);

	
	
	if (StrikeFixing<today)
		StrikeFixing = today;

	err =  get_fra(StrikeFixing, spotlag,  cRefRateCode, szYieldCurveName, &StrikeForward);																									

	
	resetmat = (
		StrikeFixing-today)/365.25;

    err = swp_f_vol(szVolCurveName,StrikeDate, StartDate, StrikeForward,&atmvolstrike,&lognorm);
	
	Fixing = add_unit(StartDate, -spotlag, SRT_BDAY, MODIFIED_SUCCEEDING);
	
	err = get_fra(Fixing, spotlag,  cRefRateCode, szYieldCurveName, &dForward);
	
	err = swp_f_vol(szVolCurveName,StartDate, EndDate, StrikeForward,&atmvolforward,&lognorm);



		for (n = 1;n<=2*nherm;n++)
		{
			w = points[n];
			Strike = StrikeForward*exp(atmvolstrike*w*sqrt(resetmat)-atmvolstrike*atmvolstrike*resetmat/2);
			dFwd = dForward+Strike-StrikeForward;
			price = Price_Caplet(StrikeDate,StartDate, EndDate, dFwd, Strike, szVolCurveName, szYieldCurveName, cRefRateCode, VolType, VolInput, CallPut);

			*pv+=weights[n]*price;
		}
		
/*free memory*/
free_dvector(weights,0,2*nherm);
free_dvector(points,0,2*nherm);


return err;
}




 














