/* Calibrates with Levenberg the parameters of a simple Merton model*/ 

#include "srt_h_all.h"
#include "num_h_allhdr.h"
#include "num_h_levenberg.h"

#include  <utallhdr.H>
#include  <opfnctns.H>
#include  <srt_h_resetable.h"
#include  <math.h"

/* Free Memory */
void free_levenberg_memory2 (double **covar,double **alpha,int nparam)
{
if (covar) free_dmatrix (covar, 1, nparam, 1, nparam);
	covar = NULL;
	if (alpha) free_dmatrix (alpha, 1, nparam, 1, nparam);
	alpha = NULL;	
}	

#define n_opt 5     /* number of parameters in the model */
#define shift 0.001   /* shift used to compute gradients */

/* Global Variables declaration */
static double *strikes = NULL;
static double *dvols = NULL;
static double mat,dFwd;
static int n_strikes;
static double *vols = NULL;
static double **gradient = NULL;
static char *Lognorm;

Err SmileLinearRegression(double dFwd, int length, double *strikes, double *volats,  double **abc)
{
Err err=NULL;
int i,j;
double **xtx;
double **xty;
double *x;

/*Memory allocation*/
x = dvector(1, length);
xtx = dmatrix(1,3,1,3);
xty = dmatrix(1,3,1,1);

for(j=1;j<=length;j++)
{
	x[j] = log(dFwd/strikes[j]);
}

for(i=1;i<=3;i++)
{	
	for(j=1;j<=3;j++)
	{
		xtx[i][j]=0;
	}
	xty[i][1]=0;
}	


for(j=1;j<=length;j++)
{
xtx[1][1]+=1;
xtx[1][2]+=x[j];
xtx[1][3]+=x[j]*x[j];
xtx[2][1]+=x[j];
xtx[2][2]+=x[j]*x[j];
xtx[2][3]+=x[j]*x[j]*x[j];
xtx[3][1]+=x[j]*x[j];
xtx[3][2]+=x[j]*x[j]*x[j];
xtx[3][3]+=x[j]*x[j]*x[j]*x[j];
xty[1][1]+=volats[j];
xty[2][1]+=volats[j]*x[j];
xty[3][1]+=volats[j]*x[j]*x[j];
}
 
err = gaussj(xtx,3,xty,1);
if (err)
{
	free_dvector(x,1,length);
	free_dmatrix(xtx,1,3,1,3);
	free_dmatrix(xty,1,3,1,1);
}

for(i=1;i<=3;i++)
{
	abc[i][1] = xty[i][1];
}

	free_dvector(x,1,length);
	free_dmatrix(xtx,1,3,1,3);
	free_dmatrix(xty,1,3,1,1);

	return err;
}



Err ComputeabcfromMerton(double dFwd, double T, double *param, double *abc)
{
Err err=NULL;
int n1,n2;
double sigma;
double U1;
double l1;
double U2;
double l2;
double probainf,probajumps,N1max,N2max;
double Sk,Wk;
double phi;



sigma = param[1];
U1	=	param[2];
l1	=	param[3];
U2	=	param[4];
l2	=	param[5];


probainf = 0.000001;
Sk = 0;
Wk = 0;
srt_f_optimpvol(optmertonpremium(dFwd, dFwd, sigma, U1, l1, U2, l2, T, "Lognormal"),
				dFwd,dFwd,T,1.0,SRT_CALL, SRT_LOGNORMAL,&(abc[1]));


N1max = l1*T+10*sqrt(l1*T);
N2max = l2*T+10*sqrt(l2*T);

phi = exp(-abc[1]*abc[1]*T/8)/(sqrt(2*SRT_PI));


for (n1 = 0;n1<N1max;n1++)
{
	for (n2 = 0;n2<N2max;n2++)
	{

		
		probajumps = exp(-(l1+l2)*T)*pow(l1*T,n1)*pow(l2*T,n2)/(fact(n1)*fact(n2));

		if (probajumps > probainf)
		{
			Sk+=probajumps*
				norm((n1*log(1+U1)+n2*log(1+U2)-(l1*U1+l2*U2)*T-sigma*sigma*T/2)/(sigma*sqrt(T)));			
			
			Wk+=probajumps*exp(-(n1*log(1+U1)+n2*log(1+U2)-(l1*U1+l2*U2)*T-sigma*sigma*T/2)*(n1*log(1+U1)+n2*log(1+U2)-(l1*U1+l2*U2)*T-sigma*sigma*T/2)/(2*sigma*sigma*T))/sqrt(2*SRT_PI);

		
		}
	}
}

abc[2] = -(norm(-abc[1]*sqrt(T)/2)-Sk)/(sqrt(T)*phi);
abc[3] = 0.5*(Wk-1/(abc[1]*T)+0.25*abc[1]*abc[2]*abc[2]*T);

return err;
}


Err ComputeVolsfromabc(double dFwd, int length, double *strikes, double *abc, double *volats)
{
Err err=NULL;
int j;

for (j = 1;j<=length;j++)
{
	volats[j] = abc[1]+abc[2]*log(dFwd/strikes[j])+abc[3]*log(dFwd/strikes[j])*log(dFwd/strikes[j]);
}


return err;
}


/* Computes the gradient of the price */
Err priceforcalibnew(double str,double param[],
					 double *value,double *Deriv,int n_param)
{
	
	int i,j,k;
	double *value_shiftee;
	double *param_shiftes;
	double *abc;
	int out_range=0;
	Err err=NULL;
	
	if (str == strikes[1])	
	{						
		
		if ((param[1] <0)
			|| (param[1] > 0.5) 
			|| (param[2] < 0) 
			|| (param[2] > 1)
			|| (param[3] < 0)
			|| (param[3] > 5) 
			|| (param[4] > 0) 
			|| (param[4] < -1)
			|| (param[5] < 0)
			|| (param[5] > 5)) out_range = 1; 
	
		if (out_range == 1)
		{
			param[1] = 0.1;
			param[2] = 0.15;
			param[3] = 1;
			param[4] = -0.15;
			param[5] = 1;
		}
			
		/*memory allocation for abc*/
		abc = dvector(1,3);

		/* Computes the implied vols for all strikes */
		err=ComputeabcfromMerton(dFwd, mat, param, abc);
		if (err)
		{
			smessage("Error in ComputeabcfromMerton");
			free_dvector(abc,1,3);
			return err;
		}
		
		err=ComputeVolsfromabc(dFwd, n_strikes, strikes, abc, vols);
		if (err)
		{
			smessage("Error in ComputeVolsfromabc");
			free_dvector(abc,1,3);
			return err;
		}
		
		/* Calculates the derivatives: optimization parameters */
		param_shiftes = dvector(1,n_param);
		value_shiftee = dvector(1,n_strikes);
		for (j=1;j<(n_param+1);++j)
		{
			for (k=1;k<(n_param+1);k++) param_shiftes[k] = param[k];
			param_shiftes[j] += param[j]*shift;
		
		/* Computes the implied vols for all strikes */

		err=ComputeabcfromMerton(dFwd, mat, param_shiftes, abc);
		if (err)
		{
			smessage("Error in ComputeabcfromMerton");
			free_dvector(param_shiftes,1,n_param);
			free_dvector(value_shiftee,1,n_strikes);
			free_dvector(abc,1,3);
			return err;
		}
		
		err=ComputeVolsfromabc(dFwd, n_strikes, strikes, abc, value_shiftee);
		if (err)
		{
			smessage("Error in ComputeVolsfromabc");
			free_dvector(param_shiftes,1,n_param);
			free_dvector(value_shiftee,1,n_strikes);
			free_dvector(abc,1,3);
			return err;
		}
				
			

			
			
			for (i=1;i<(n_strikes+1);++i)
				gradient[i][j] = (value_shiftee[i]-vols[i])/(shift*param[j]);
	}
		
		out_range = 0;
		
	
		if ((param[1] <0)
			|| (param[1] > 0.5) 
			|| (param[2] < 0) 
			|| (param[2] > 1)
			|| (param[3] < 0)
			|| (param[3] > 5) 
			|| (param[4] > 0) 
			|| (param[4] < -1)
			|| (param[5] < 0)
			|| (param[5] > 5)) out_range = 1; 
	
		if (out_range==1)
			for (i=1;i<(n_strikes+1);++i) 
			{
				vols[i] = -10000;
				for (j=1;j<(n_param+1);++j) gradient[i][j] = 10000;
			}
		free_dvector(param_shiftes,1,n_param);
		free_dvector(value_shiftee,1,n_strikes);
		free_dvector(abc,1,3);
	}
	
	/* Returns the elements of the vector after each calls of the function*/
	for (i=1;i<(n_strikes+1);i++)		
	{
		if (str == strikes[i])
		{
			*value = vols[i];
			for (j=1;j<(n_param+1);j++) Deriv[j] = gradient[i][j];
			break;
		}
	}   
	return err;

}








/* Computes the gradient of the price */
Err priceforcalib(double str,double param[],
					 double *value,double *Deriv,int n_param)
{
	
	int i,j,k;
	double *value_shiftee;
	double *param_shiftes;
	int out_range=0;
	Err err=NULL;
	
	if (str == strikes[1])	
	{						
		
		if (strcmp(Lognorm,"Normal")==0)
		{
			if ((param[1] <0)
				|| (param[1] > 0.1) 
				|| (param[2] < 0) 
				|| (param[2] > 0.02)
				|| (param[3] < 0)
				|| (param[3] > 10) 
				|| (param[4] > 0) 
				|| (param[4] < -0.02)
				|| (param[5] < 0)
				|| (param[5] > 10)) out_range = 1; 
		
			if (out_range == 1)
			{
				param[1] = 0.01;
				param[2] = 0.0025;
				param[3] = 1;
				param[4] = -0.0025;
				param[5] = 1;
			}
		}
	else
	{
		if ((param[1] <0)
			|| (param[1] > 0.5) 
			|| (param[2] < 0) 
			|| (param[2] > 1)
			|| (param[3] < 0)
			|| (param[3] > 5) 
			|| (param[4] > 0) 
			|| (param[4] < -1)
			|| (param[5] < 0)
			|| (param[5] > 5)) out_range = 1; 
	
		if (out_range == 1)
		{
			param[1] = 0.1;
			param[2] = 0.15;
			param[3] = 1;
			param[4] = -0.15;
			param[5] = 1;
		}
	}		
		/* Computes the implied vols for all strikes */
		err=optmertonsmile(
			dFwd,
			n_strikes,
			strikes,
			param,
			mat,
			Lognorm,
			vols /*implied vols*/
			);
		if (err)
		{ 
			smessage("Error in optmertonsmile");
			return err;
		}
		
		
		/* Calculates the derivatives: optimization parameters */
		param_shiftes = dvector(1,n_param);
		value_shiftee = dvector(1,n_strikes);
		for (j=1;j<(n_param+1);++j)
		{
			for (k=1;k<(n_param+1);k++) param_shiftes[k] = param[k];
			param_shiftes[j] += param[j]*shift;
		
			err=optmertonsmile(
				dFwd,
				n_strikes,
				strikes,
				param_shiftes,
				mat,
				Lognorm,
				value_shiftee/*resulting implied vols*/
				);
		if (err)
			{
			smessage("Error in optmertonsmile");
			free_dvector(param_shiftes,1,n_param);
			free_dvector(value_shiftee,1,n_strikes);
			return err;	
			}
		

			
			
			for (i=1;i<(n_strikes+1);++i)
				gradient[i][j] = (value_shiftee[i]-vols[i])/(shift*param[j]);
		}
		
		out_range = 0;
		
		if (strcmp(Lognorm,"Normal")==0)
		{

		if ((param[1] <0)
			|| (param[1] > 0.1) 
			|| (param[2] < 0) 
			|| (param[2] > 0.02)
			|| (param[3] < 0)
			|| (param[3] > 10) 
			|| (param[4] > 0) 
			|| (param[4] < -0.02)
			|| (param[5] < 0)
			|| (param[5] > 10)) out_range = 1; 
			
		}

		else
		{

		if ((param[1] <0)
			|| (param[1] > 0.5) 
			|| (param[2] < 0) 
			|| (param[2] > 1)
			|| (param[3] < 0)
			|| (param[3] > 5) 
			|| (param[4] > 0) 
			|| (param[4] < -1)
			|| (param[5] < 0)
			|| (param[5] > 5)) out_range = 1; 
		}
		if (out_range==1)
			for (i=1;i<(n_strikes+1);++i) 
			{
				vols[i] = -10000;
				for (j=1;j<(n_param+1);++j) gradient[i][j] = 10000;
			}
		free_dvector(param_shiftes,1,n_param);
		free_dvector(value_shiftee,1,n_strikes);
	}
	
	/* Returns the elements of the vector after each calls of the function*/
	for (i=1;i<(n_strikes+1);i++)		
	{
		if (str == strikes[i])
		{
			*value = vols[i];
			for (j=1;j<(n_param+1);j++) Deriv[j] = gradient[i][j];
			break;
		}
	}   
	return err;

}

/* Functions which finds the parameters which fit the market vols*/
Err optcalibmerton(
					  double dFwd_loc,
					  int n_strikes_loc,
	                  double *strikes_loc,
					  double *market_vols,
					  double *param,
					  double maturity,
					  char   *Logornorm,
					  double *chisq, 
					  long   *ia,
					  double *calibvols
					  )
{
double *sig;
double Strike_inf,Strike_sup;
int i,i1,i2;
int nparam;
Err err=NULL;

	
	Lognorm = strdup(Logornorm);
	nparam = 5;
	dFwd = dFwd_loc;
	n_strikes = n_strikes_loc;
	mat = maturity;
	strikes = dvector(1,n_strikes);
	dvols = dvector(1,n_strikes);
	vols = dvector(1,n_strikes);
	gradient = dmatrix(1,n_strikes,1,nparam);
	for (i=1;i<(n_strikes+1);i++) strikes[i] = strikes_loc[i];
	for (i=1;i<(n_strikes+1);i++) dvols[i] = market_vols[i];
	sig = dvector(1,n_strikes);
	


	/* Assignment of weights to the different strikes. More weight is given to strikes 
	closer to the ATM. Also, basically no weight is given to strikes further away than
	2 sd's from the money. 
	*/
	i=1;

	while (strikes[i] < dFwd) i++;
	i--;

	if (strcmp(Lognorm,"Normal")==0)
	{
		Strike_inf = dFwd - 2*market_vols[i]*sqrt(mat);
		Strike_sup = dFwd + 2*market_vols[i]*sqrt(mat);
	}
	else
	{
		Strike_inf = dFwd *(1- 2*market_vols[i]*sqrt(mat));
		Strike_sup = dFwd *(1+ 2*market_vols[i]*sqrt(mat));
	}
	i1 = 1;
	i2 = n_strikes;

	for (i=1;i<n_strikes-1;i++)
	{
		if ((strikes[i]<Strike_inf)&&(strikes[i+1]>Strike_inf)) i1 = i;
		if ((strikes[i]<Strike_sup)&&(strikes[i+1]>Strike_sup)) i2 = i;
	}
	
	for (i=1;i<i1;i++) sig[i] = 100;
	for (i=i1;i<i2+1;i++) sig[i] = 10*fabs(strikes[i]-dFwd)+0.0001;
	for (i=i2+1;i<n_strikes+1;i++) sig[i] = 100;


	nparam = 5;							
	
	err = levenberg_marquardt_select(
					strikes,
					dvols,
	 				sig,
					n_strikes_loc,
					param,
					ia,
					nparam,
					25,
					priceforcalib,
					chisq);
	if (err)
	{
		smessage("Error in levenberg_marquardt_select");
	} else {
		err= optmertonsmile(dFwd,n_strikes, strikes, 
									 param, mat,Logornorm,calibvols);
		if (err)
			smessage("Error in optmertonsmile");
	}


	free_dvector(strikes,1,n_strikes);
	free_dvector(dvols,1,n_strikes);
	free_dvector(sig,1,n_strikes);
	free_dvector(vols,1,n_strikes);
	free_dmatrix(gradient,1,n_strikes,1,nparam);
	free(Lognorm);
	return err;
}




/* Functions which finds the parameters which fit the market vols*/
Err optcalibmertonnew(
					  double dFwd_loc,
					  int n_strikes_loc,
	                  double *strikes_loc,
					  double *market_vols,
					  double *param,
					  double maturity,
					  char   *Logornorm,
					  double *chisq, 
					  long   *ia,
					  double *calibvols
					  )
{
double *sig;
double Strike_inf,Strike_sup;
int i,i1,i2;
int nparam;
Err err=NULL;

	
	Lognorm = strdup(Logornorm);
	nparam = 5;
	dFwd = dFwd_loc;
	n_strikes = n_strikes_loc;
	mat = maturity;
	strikes = dvector(1,n_strikes);
	dvols = dvector(1,n_strikes);
	vols = dvector(1,n_strikes);
	gradient = dmatrix(1,n_strikes,1,nparam);
	for (i=1;i<(n_strikes+1);i++) strikes[i] = strikes_loc[i];
	for (i=1;i<(n_strikes+1);i++) dvols[i] = market_vols[i];
	sig = dvector(1,n_strikes);
	


	/* Assignment of weights to the different strikes. More weight is given to strikes 
	closer to the ATM. Also, basically no weight is given to strikes further away than
	2 sd's from the money. 
	*/
	i=1;

	while (strikes[i] < dFwd) i++;


	if (strcmp(Lognorm,"Normal")==0)
	{
		Strike_inf = dFwd - 2*market_vols[i]*sqrt(mat);
		Strike_sup = dFwd + 2*market_vols[i]*sqrt(mat);
	}
	else
	{
		Strike_inf = dFwd *(1- 2*market_vols[i]*sqrt(mat));
		Strike_sup = dFwd *(1+ 2*market_vols[i]*sqrt(mat));
	}
	i1 = 1;
	i2 = n_strikes;

	for (i=1;i<n_strikes-1;i++)
	{
		if ((strikes[i]<Strike_inf)&&(strikes[i+1]>Strike_inf)) i1 =i;
		if ((strikes[i]<Strike_sup)&&(strikes[i+1]>Strike_sup)) i2 = i;
	}
	
	for (i=1;i<i1;i++) sig[i] = 100;
	for (i=i1;i<i2+1;i++) sig[i] = 10*fabs(strikes[i]-dFwd)+0.0001;
	for (i=i2+1;i<n_strikes+1;i++) sig[i] = 100;


	nparam = 5;							
	
	err = levenberg_marquardt_select(
					strikes,
					dvols,
					sig,
					n_strikes_loc,
					param,
					ia,
					nparam,
					25,
					priceforcalib,
					chisq);
	if (err)
	{
		smessage("Error in levenberg_marquardt_select");
	} else {
		err= optmertonsmile(dFwd,n_strikes, strikes, 
									 param, mat,Logornorm,calibvols);
		if (err)
			smessage("Error in optmertonsmile");
	}


	free_dvector(strikes,1,n_strikes);
	free_dvector(dvols,1,n_strikes);
	free_dvector(sig,1,n_strikes);
	free_dvector(vols,1,n_strikes);
	free_dmatrix(gradient,1,n_strikes,1,nparam);
	free(Lognorm);
	return err;
}
