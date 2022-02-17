/* ==========================================================================
   FILE_NAME:	DoubleLGM1FQuanto_pde.c

   PURPOSE:		PDE implementation of the double 1 Factor LGM quanto model.
				Discretisation is LOD.

   DATE:		02/09/01
   
   AUTHOR:		B.B. (from L.C.'s LGM2Fpde.c)
   ========================================================================== */

#define	NSTD_LGM		7.0
#define	NB_CHECK		10

#include "DoubleLGM1FQuanto_pde.h"
#include "math.h"

static double Ji0_func(double lam, double t1, double t2)
{
	return (exp(lam*t2) - exp(lam*t1)) / lam;
}

static double Ji1_func(double lam, double t1, double t2)
{
	return ( (t2*exp(lam*t2) - t1*exp(lam*t1)) - Ji0_func(lam, t1, t2) ) / lam ;
}

static double Ji2_func(double lam, double t1, double t2)
{
	return ( (t2*t2*exp(lam*t2) - t1*t1*exp(lam*t1)) - 2 * Ji1_func(lam, t1, t2) ) / lam ;
}

static double Ji3_func(double lam, double t1, double t2)
{
	return ( (t2*t2*t2*exp(lam*t2) - t1*t1*t1*exp(lam*t1)) - 3 * Ji2_func(lam, t1, t2) ) / lam ;
}

static double Ji4_func(double lam, double t1, double t2)
{
	return ( (t2*t2*t2*t2*exp(lam*t2) - t1*t1*t1*t1*exp(lam*t1)) - 4 * Ji3_func(lam, t1, t2) ) / lam ;
}


Err convertVol(
					double lam,
					int nstept, 
					double *time, 
					int nb_sig, 
					double *sig_time, 
					double *sig, 
					double *sigtilda,
					double *contsig,
					double **ablinsig)
{
	int i,j;
	double t1, t2;
	char *err = NULL;
	double **ablinsigloc = NULL;
	double e, A, B, C, D, E, F, Delta;

	ablinsigloc = dmatrix(0, 1, 0, nb_sig);
	if (!ablinsigloc)
	{
		err = "Memory allocation error (1) in convertVol";
		goto FREE_RETURN;
	}

	sigtilda[0] = sig[0];
	ablinsigloc[0][0] = sig[0];
	ablinsigloc[1][0] = 0.0;

	for(i=0;i<nb_sig-1;++i)
	{
		e = sig_time[i+1]/(sig_time[i+1]-sig_time[i]);
		A = Ji0_func(2*lam, sig_time[i], sig_time[i+1]);
		B = 2*Ji1_func(2*lam, sig_time[i], sig_time[i+1])/(sig_time[i+1]-sig_time[i]);
		C = Ji2_func(2*lam, sig_time[i], sig_time[i+1])/((sig_time[i+1]-sig_time[i])*(sig_time[i+1]-sig_time[i]));
		D = A*(1-e)*(1-e) + B*(1-e) + C;
		E = 2*A*e*(1-e)*sigtilda[i] - B*(1-2*e)*sigtilda[i] - 2*C*sigtilda[i];
		F = (A*e*e - B*e + C)*sigtilda[i]*sigtilda[i] - sig[i]*sig[i]*A;
		Delta = E*E - 4*D*F;
		if(Delta<0)
		{
			err = "Cannot convert PW constant vol in continuous PW affine vol";
			goto FREE_RETURN;
		}

		sigtilda[i+1] = (-E+sqrt(Delta))/(2*D);
		if(sigtilda[i+1]<0)
		{
			err = "continuous PW affine vol is negative";
			goto FREE_RETURN;
		}
	}

	for(i=1;i<nb_sig;++i)
	{
		ablinsigloc[0][i] = sigtilda[i-1] - (sigtilda[i]-sigtilda[i-1])*sig_time[i-1]/(sig_time[i]-sig_time[i-1]);
		ablinsigloc[1][i] = (sigtilda[i]-sigtilda[i-1])/(sig_time[i]-sig_time[i-1]);
	}

	j=0;
	for(i=0;i<nstept;++i)
	{
		t2 = time[i];
		t1 = sig_time[j];

		while (t1 < t2 && j < nb_sig-1)
		{
			j++;
			t1 = sig_time[j];
		}

		ablinsig[0][i] = ablinsigloc[0][j];
		ablinsig[1][i] = ablinsigloc[1][j];

		if(t2>sig_time[nb_sig-1])
		{
			contsig[i] = sigtilda[nb_sig-1];
		}
		else
		{
			contsig[i] = ablinsigloc[0][j] + ablinsigloc[1][j] * t2;
		}
	}

FREE_RETURN:

	if (ablinsigloc) free_dmatrix(ablinsigloc, 0, 1, 0, nb_sig);

	return err;

}


Err convertVol2(
					double lam,
					int nstept, 
					double *time, 
					int nb_sig, 
					double *sig_time, 
					double *sig, 
					double *sigtilda,
					double *contsig,
					double **ablinsig)
{
	int i,j;
	double t1, t2;
	char *err = NULL;
	double **ablinsigloc = NULL;
	double A, B, C, Delta;

	ablinsigloc = dmatrix(0, 1, 0, nb_sig);
	if (!ablinsigloc)
	{
		err = "Memory allocation error (1) in convertVol2";
		goto FREE_RETURN;
	}

	sigtilda[0] = sig[0];
	ablinsigloc[0][0] = sig[0];
	ablinsigloc[1][0] = 0.0;

	for(i=0;i<nb_sig-1;++i)
	{
		if( (sig_time[i+1] - sig_time[i]) > 2.0/365.0 )
		{

			A = Ji0_func(2*lam, sig_time[i], sig_time[i+1])*sig_time[i]*sig_time[i]
				-2*Ji1_func(2*lam, sig_time[i], sig_time[i+1])*sig_time[i]
				+Ji2_func(2*lam, sig_time[i], sig_time[i+1]);

			B = -2*sigtilda[i]*(
				Ji0_func(2*lam, sig_time[i], sig_time[i+1])*sig_time[i+1]*sig_time[i]
				-Ji1_func(2*lam, sig_time[i], sig_time[i+1])*(sig_time[i+1]+sig_time[i])
				+Ji2_func(2*lam, sig_time[i], sig_time[i+1])
							);

			C = sigtilda[i]*sigtilda[i]*(
				Ji0_func(2*lam, sig_time[i], sig_time[i+1])*sig_time[i+1]*sig_time[i+1]
				-2*Ji1_func(2*lam, sig_time[i], sig_time[i+1])*sig_time[i+1]
				+Ji2_func(2*lam, sig_time[i], sig_time[i+1])
									)
				- sig[i+1]*sig[i+1]*Ji0_func(2*lam, sig_time[i], sig_time[i+1])*(sig_time[i+1]-sig_time[i])*(sig_time[i+1]-sig_time[i]);

			Delta = B*B - 4*A*C;
			if(Delta<0)
			{
				err = "Cannot convert PW constant vol in continuous PW affine vol";
				goto FREE_RETURN;
			}
		
			sigtilda[i+1] = (-B+sqrt(Delta))/(2*A);
			if(sigtilda[i+1]<0)
			{
				err = "continuous PW affine vol is negative";
				goto FREE_RETURN;
			}
		}
		else
		{
			sigtilda[i+1] = sigtilda[i];
		}

	}

	for(i=1;i<nb_sig;++i)
	{
		ablinsigloc[0][i] = sigtilda[i-1] - (sigtilda[i]-sigtilda[i-1])*sig_time[i-1]/(sig_time[i]-sig_time[i-1]);
		ablinsigloc[1][i] = (sigtilda[i]-sigtilda[i-1])/(sig_time[i]-sig_time[i-1]);
	}

	j=0;
	for(i=0;i<nstept;++i)
	{
		t2 = time[i];
		t1 = sig_time[j];

		while (t1 < t2 && j < nb_sig-1)
		{
			j++;
			t1 = sig_time[j];
		}

		ablinsig[0][i] = ablinsigloc[0][j];
		ablinsig[1][i] = ablinsigloc[1][j];

		if(t2>sig_time[nb_sig-1])
		{
			contsig[i] = sigtilda[nb_sig-1];
		}
		else
		{
			contsig[i] = ablinsigloc[0][j] + ablinsigloc[1][j] * t2;
		}
	}

FREE_RETURN:

	if (ablinsigloc) free_dmatrix(ablinsigloc, 0, 1, 0, nb_sig);

	return err;

}


Err convertVol3(
					double lam,
					int nstept, 
					double *time, 
					int nb_sig, 
					double *sig_time, 
					double *sig, 
					double *sigtilda,
					double *contsig,
					double **ablinsig)
{
	int i,j;
	double t1, t2;
	char *err = NULL;
	double **ablinsigloc = NULL;

	ablinsigloc = dmatrix(0, 1, 0, nb_sig);
	if (!ablinsigloc)
	{
		err = "Memory allocation error (1) in convertVol2";
		goto FREE_RETURN;
	}

	sigtilda[0] = sig[0];
	ablinsigloc[0][0] = sig[0];
	ablinsigloc[1][0] = 0.0;

	for(i=0;i<nb_sig;++i)
	{
		sigtilda[i] = sig[i];
	}

	for(i=1;i<nb_sig;++i)
	{
		ablinsigloc[0][i] = sigtilda[i-1] - (sigtilda[i]-sigtilda[i-1])*sig_time[i-1]/(sig_time[i]-sig_time[i-1]);
		ablinsigloc[1][i] = (sigtilda[i]-sigtilda[i-1])/(sig_time[i]-sig_time[i-1]);
	}

	j=0;
	for(i=0;i<nstept;++i)
	{
		t2 = time[i];
		t1 = sig_time[j];

		while (t1 < t2 && j < nb_sig-1)
		{
			j++;
			t1 = sig_time[j];
		}

		ablinsig[0][i] = ablinsigloc[0][j];
		ablinsig[1][i] = ablinsigloc[1][j];

		if(t2>sig_time[nb_sig-1])
		{
			contsig[i] = sigtilda[nb_sig-1];
		}
		else
		{
			contsig[i] = ablinsigloc[0][j] + ablinsigloc[1][j] * t2;
		}
	}

FREE_RETURN:

	if (ablinsigloc) free_dmatrix(ablinsigloc, 0, 1, 0, nb_sig);

	return err;

}



Err convertCorrelationWithTwoAffineVols(
					double domlam,
					double forlam,
					int nstept, 
					double *time, 
					int nb_sig, 
					double *sig_time, 
					double *correl, 
					double *domsig, 
					double *forsig, 
					double *domsigtilda,
					double *forsigtilda,
					double *correltilda,
					double **ablincorrel,
					double *contcorrel)
{
	int i,j;
	double t1, t2;
	char *err = NULL;
	double **ablincorrelloc = NULL;
	double A, B;
	double DeltaT;
	double domsiga, domsigb;
	double forsiga, forsigb;
	double rho;
	double J0,J1,J2,J3;


	ablincorrelloc = dmatrix(0, 1, 0, nb_sig);
	if (!ablincorrelloc)
	{
		err = "Memory allocation error (1) in convertVol2";
		goto FREE_RETURN;
	}

	correltilda[0] = correl[0];
	ablincorrelloc[0][0] = correl[0];
	ablincorrelloc[1][0] = 0.0;

	for(i=0;i<nb_sig-1;++i)
	{
		DeltaT = sig_time[i+1]-sig_time[i];

		if(DeltaT > 30.0/365.0)
		{
			J0 = Ji0_func(domlam + forlam, sig_time[i], sig_time[i+1]); 
			J1 = Ji1_func(domlam + forlam, sig_time[i], sig_time[i+1]); 
			J2 = Ji2_func(domlam + forlam, sig_time[i], sig_time[i+1]); 
			J3 = Ji3_func(domlam + forlam, sig_time[i], sig_time[i+1]); 

			domsiga = (domsigtilda[i+1]-domsigtilda[i])/(sig_time[i+1]-sig_time[i]);
			domsigb = domsigtilda[i] - domsiga * sig_time[i];

			forsiga = (forsigtilda[i+1]-forsigtilda[i])/(sig_time[i+1]-sig_time[i]);
			forsigb = forsigtilda[i] - forsiga * sig_time[i];

			A =	(
				J3 * domsiga * forsiga
				+ J2 * (domsiga*forsigb+forsiga*domsigb-sig_time[i]*domsiga*forsiga)
				+ J1 * (domsigb*forsigb-sig_time[i]*(domsiga*forsigb+forsiga*domsigb))
				- J0 * sig_time[i]*domsigb*forsigb
				)/DeltaT;
		
			B = (
				- J3 * domsiga * forsiga
				+ J2 * (sig_time[i+1]*domsiga*forsiga-(domsiga*forsigb+forsiga*domsigb))
				+ J1 * ((domsiga*forsigb+forsiga*domsigb)*sig_time[i+1]-domsigb*forsigb)
				+ J0 * domsigb * forsigb * sig_time[i+1]
				)/DeltaT;

			rho = (correl[i]*domsig[i]*forsig[i]*J0
								- B*correltilda[i]) / A;

			if((rho<-1)||(rho>1))
			{
				smessage ("Cannot convert PW Const correl in continuous PW affine correl, Abs(correl)>1");
//				err = "Cannot convert PW Const correl in continuous PW affine correl, Abs(correl)>1";
//				goto FREE_RETURN;
			}

			correltilda[i+1] = max(-0.99,min(0.99,rho));
//			correltilda[i+1] = correl[i+1];
		}
		else
		{
			correltilda[i+1] = correltilda[i];
//			correltilda[i+1] = correl[i+1];
		}
	}

	for(i=1;i<nb_sig;++i)
	{
		ablincorrelloc[0][i] = correltilda[i-1] - (correltilda[i]-correltilda[i-1])*sig_time[i-1]/(sig_time[i]-sig_time[i-1]);
		ablincorrelloc[1][i] = (correltilda[i]-correltilda[i-1])/(sig_time[i]-sig_time[i-1]);
	}

	j=0;
	for(i=0;i<nstept;++i)
	{
		t2 = time[i];
		t1 = sig_time[j];

		while (t1 < t2 && j < nb_sig-1)
		{
			j++;
			t1 = sig_time[j];
		}

		ablincorrel[0][i] = ablincorrelloc[0][j];
		ablincorrel[1][i] = ablincorrelloc[1][j];

		if(t2>sig_time[nb_sig-1])
		{
			contcorrel[i] = correltilda[nb_sig-1];
		}
		else
		{
			contcorrel[i] = ablincorrelloc[0][j] + ablincorrelloc[1][j] * t2;
		}
	}

FREE_RETURN:

	if (ablincorrelloc) free_dmatrix(ablincorrelloc, 0, 1, 0, nb_sig);

	return err;

}


Err convertCorrelationWithOneAffineVol(
					double forlam,
					int nstept, 
					double *time, 
					int nb_sig, 
					double *sig_time, 
					double *correl, 
					double *forsig, 
					double *forsigtilda,
					double *correltilda,
					double **ablincorrel,
					double *contcorrel)
{
	int i,j;
	double t1, t2;
	char *err = NULL;
	double **ablincorrelloc = NULL;
	double A, B;
	double DeltaT;
	double forsiga, forsigb;
	double rho;
	double J0,J1,J2;


	ablincorrelloc = dmatrix(0, 1, 0, nb_sig);
	if (!ablincorrelloc)
	{
		err = "Memory allocation error (1) in convertVol2";
		goto FREE_RETURN;
	}

	correltilda[0] = correl[0];
	ablincorrelloc[0][0] = correl[0];
	ablincorrelloc[1][0] = 0.0;

	for(i=0;i<nb_sig-1;++i)
	{
		DeltaT = sig_time[i+1]-sig_time[i];

		if(DeltaT>30.0/365.0)
		{
			J0 = Ji0_func(forlam, sig_time[i], sig_time[i+1]); 
			J1 = Ji1_func(forlam, sig_time[i], sig_time[i+1]); 
			J2 = Ji2_func(forlam, sig_time[i], sig_time[i+1]); 

			forsiga = (forsigtilda[i+1]-forsigtilda[i])/(sig_time[i+1]-sig_time[i]);
			forsigb = forsigtilda[i] - forsiga * sig_time[i];

			A =	forsiga * J2 + forsigb * J1 
				- forsiga * sig_time[i] * J1 - forsigb * sig_time[i] * J0;
		
			B = forsiga * J1 * sig_time[i+1] - forsiga * J2 
				+ forsigb * J0 * sig_time[i+1] - forsigb * J1;

			rho = (correl[i]*forsig[i]*J0*DeltaT
							- B*correltilda[i]) / A;

			if(fabs(rho)>1)
			{
//				err = "Cannot convert PW Const correl in continuous PW affine correl, Abs(correl)>1";
//				goto FREE_RETURN;
			}
			correltilda[i+1] = max(-0.99,min(0.99,rho));
//			correltilda[i+1] = correl[i+1];
		}
		else
		{
			correltilda[i+1] = correltilda[i];
//			correltilda[i+1] = correl[i+1];
		}

	}

	for(i=1;i<nb_sig;++i)
	{
		ablincorrelloc[0][i] = correltilda[i-1] - (correltilda[i]-correltilda[i-1])*sig_time[i-1]/(sig_time[i]-sig_time[i-1]);
		ablincorrelloc[1][i] = (correltilda[i]-correltilda[i-1])/(sig_time[i]-sig_time[i-1]);
	}

	j=0;
	for(i=0;i<nstept;++i)
	{
		t2 = time[i];
		t1 = sig_time[j];

		while (t1 < t2 && j < nb_sig-1)
		{
			j++;
			t1 = sig_time[j];
		}

		ablincorrel[0][i] = ablincorrelloc[0][j];
		ablincorrel[1][i] = ablincorrelloc[1][j];

		if(t2>sig_time[nb_sig-1])
		{
			contcorrel[i] = correltilda[nb_sig-1];
		}
		else
		{
			contcorrel[i] = ablincorrelloc[0][j] + ablincorrelloc[1][j] * t2;
		}
	}

FREE_RETURN:

	if (ablincorrelloc) free_dmatrix(ablincorrelloc, 0, 1, 0, nb_sig);

	return err;

}


/*----------------------------------------------------------------------*/
/*--Function to convert a PW constant vol to a continuous PW quadratic--*/
/*------------keeping the factor's variance unchanged ------------------*/
/*----------------------------------------------------------------------*/
Err convertVolQuadratic(
					double lam,
					int nstept, 
					double *time, 
					int nb_sig, 
					double *sig_time, 
					double *sig, 
					double *alpha, 
					double *contsig,
					double **abclinsig)
{
	int i,j;
	double t1, t2;
	char *err = NULL;
	double **abclinsigloc = NULL;
	double A, B, C, D, E, F, G, Delta;

	abclinsigloc = dmatrix(0, 2, 0, nb_sig-1);
	if (!abclinsigloc)
	{
		err = "Memory allocation error (1) in convertVolQuadratic";
		goto FREE_RETURN;
	}

	abclinsigloc[0][0] = sig[0];
	abclinsigloc[1][0] = 0.0;
	abclinsigloc[2][0] = 0.0;
	alpha[0] = 0.0;

	for(i=0;i<nb_sig-1;++i)
	{
		A = (sig[i+1]-sig[i])/(sig_time[i+1]-sig_time[i]);
		B = sig_time[i+1] + sig_time[i];
		C = sig_time[i+1] * sig_time[i];
		D = (sig[i]*sig_time[i+1]-sig[i+1]*sig_time[i])/(sig_time[i+1]-sig_time[i]);

		E = Ji4_func(2*lam, sig_time[i], sig_time[i+1])
			- 2*B*Ji3_func(2*lam, sig_time[i], sig_time[i+1])
			+ (2*C+B*B)*Ji2_func(2*lam, sig_time[i], sig_time[i+1])
			- 2*B*C*Ji1_func(2*lam, sig_time[i], sig_time[i+1])
			+ C*C*Ji0_func(2*lam, sig_time[i], sig_time[i+1]);
		F = A*Ji3_func(2*lam, sig_time[i], sig_time[i+1])
			+ D*Ji2_func(2*lam, sig_time[i], sig_time[i+1])
			- A*B*Ji2_func(2*lam, sig_time[i], sig_time[i+1])
			+ (A*C-B*D)*Ji1_func(2*lam, sig_time[i], sig_time[i+1])
			+ C*D*Ji0_func(2*lam, sig_time[i], sig_time[i+1]);
		G = A*A*Ji2_func(2*lam, sig_time[i], sig_time[i+1])
			+ 2*A*D*Ji1_func(2*lam, sig_time[i], sig_time[i+1])
			+ (D*D-sig[i+1]*sig[i+1])*Ji0_func(2*lam, sig_time[i], sig_time[i+1]);

		Delta = F*F-E*G;
		if(Delta<0)
		{
			err = "Cannot convert PW constant vol in continuous PW quadratic vol";
			goto FREE_RETURN;
		}

		alpha[i+1] = (-F-sqrt(Delta))/E;
	}

	for(i=1;i<nb_sig;++i)
	{
		abclinsigloc[0][i] = sig[i-1] + alpha[i] * sig_time[i] * sig_time[i-1]
								- sig_time[i-1] * (sig[i]-sig[i-1])/(sig_time[i]-sig_time[i-1]);
		abclinsigloc[1][i] = (sig[i]-sig[i-1])/(sig_time[i]-sig_time[i-1])
								- alpha[i] * (sig_time[i]+sig_time[i-1]);
		abclinsigloc[2][i] = alpha[i];
	}

	j=0;
	for(i=0;i<nstept;++i)
	{
		t2 = time[i];
		t1 = sig_time[j];

		while (t1 < t2 && j < nb_sig-1)
		{
			j++;
			t1 = sig_time[j];
		}

		abclinsig[0][i] = abclinsigloc[0][j];
		abclinsig[1][i] = abclinsigloc[1][j];
		abclinsig[2][i] = abclinsigloc[2][j];

		if(t2>sig_time[nb_sig-1])
		{
			contsig[i] = sig[nb_sig-1];
		}
		else
		{
			contsig[i] = abclinsigloc[0][j] + abclinsigloc[1][j] * t2 + abclinsigloc[2][j] * t2 * t2;
		}
	}

FREE_RETURN:

	if (abclinsigloc) free_dmatrix(abclinsigloc, 0, 2, 0, nb_sig-1);

	return err;

}
/*-------------------------------------------------------------------*/
/*------------------End Of convertVolquadratic function--------------*/
/*-------------------------------------------------------------------*/


void fillVol(
					int nstept, 
					double *time, 
					int nb_sig, 
					double *sig_time, 
					double *sig, 
					double *sigall)
{
	int i,j;
	double t1, t2;

	j=0;
	for(i=0;i<nstept;++i)
	{
		t2 = time[i];
		t1 = sig_time[j];

		while (t1 < t2 && j < nb_sig-1)
		{
			j++;
			t1 = sig_time[j];
		}

		sigall[i] = sig[j];
	}
}


void fillVolAndCorrel(
					int nstept, 
					double *time, 
					int nb_sig, 
					double *sig_time, 
					double *domsig, 
					double *forsig, 
					double *fxsig, 
					double *domforrho, 
					double *quantorho, 
					double *domsigall,
					double *forsigall,
					double *fxsigall,
					double *domforrhoall,
					double *quantorhoall)
{
	int i,j;
	double t1, t2;

	j=0;
	for(i=0;i<nstept;++i)
	{
		t2 = time[i];
		t1 = sig_time[j];

		while (t1 < t2 && j < nb_sig-1)
		{
			j++;
			t1 = sig_time[j];
		}

		domsigall[i] = domsig[j];
		forsigall[i] = forsig[j];
		fxsigall[i] = fxsig[j];
		domforrhoall[i] = domforrho[j];
		quantorhoall[i] = quantorho[j];
	}
}


/*	Function to evaluate the expectations	
	of the variables r1_dim1 and r3				
	This function suppose that all the vol dates
	are included in the time discretisation			*/


static void DoubleLGM1FQuantoExpectations(
						int			nstept,
						double		*time,
						double		domlam,
						double		forlam,
						double		*sig_time,	//	domsig, forsig and fxsig must have 
						double		*domsig,	//	the same sig_time
						double		*forsig,
						double		*fxsig,
						int			nb_sig,
						double		quantorho,
						double		domforrho,
						double		*fwd1,
						double		*fwd3,
						double		*var1,
						double		*var2,
						double		*domphi,
						double		*forphi,
						double		*domforphi)
{
	double	t1, t2, ta, tb;
	int		i, j, nb_sig_minus1;
	double domI1, domI2;
	double forI1, forI2;
	double domforI2;
	double domH, forH;
	double doma, fora;
	double domb, forb;
	double QAdj;
	double std1;

	nb_sig_minus1 = nb_sig - 1;

	// initialisation 
	t1 = 0.0;
	j = 0;

	doma = domsig[0];
	domb = 0;
	fora = forsig[0];
	forb = 0;

	domphi[0] = 0.0;
	forphi[0] = 0.0;
	domforphi[0] = 0.0;
	
	fwd1[0] = 0;
	fwd3[0] = 0;
	domI1 = 0;
	domI2 = 0;
	forI1 = 0;
	forI2 = 0;
	domforI2 = 0;
	domH = 0;
	forH = 0;
	QAdj = 0;

	for (i=1; i<nstept; i++)
	{		
		t2 = time[i];
		ta = t1;
		tb = sig_time[j];

		std1 = 0;

		while (tb < t2 && j < nb_sig_minus1)
		{
			if(j>0)
			{
				doma = domsig[j-1] - (domsig[j] - domsig[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
				domb = (domsig[j] - domsig[j-1]) / (sig_time[j] - sig_time[j-1]);
			
				fora = forsig[j-1] - (forsig[j] - forsig[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
				forb = (forsig[j] - forsig[j-1]) / (sig_time[j] - sig_time[j-1]);
			}

			domI1 += doma * doma * Ji0_func(domlam, ta, tb)
						+ 2 * doma * domb * Ji1_func(domlam, ta, tb)
						+ domb * domb * Ji2_func(domlam, ta, tb);

			domI2 += doma * doma * Ji0_func(2*domlam, ta, tb)
						+ 2 * doma * domb * Ji1_func(2*domlam, ta, tb)
						+ domb * domb * Ji2_func(2*domlam, ta, tb);

			forI1 += fora * fora * Ji0_func(forlam, ta, tb)
						+ 2 * fora * forb * Ji1_func(forlam, ta, tb)
						+ forb * forb * Ji2_func(forlam, ta, tb);

			forI2 += fora * fora * Ji0_func(2*forlam, ta, tb)
						+ 2 * fora * forb * Ji1_func(2*forlam, ta, tb)
						+ forb * forb * Ji2_func(2*forlam, ta, tb);

			domforI2 += doma * fora * Ji0_func(domlam+forlam, ta, tb)
						+ (doma * forb + fora * domb) * Ji1_func(domlam+forlam, ta, tb)
						+ domb * forb * Ji2_func(domlam+forlam, ta, tb);

			QAdj += quantorho * fxsig[j]
				* ( fora * Ji0_func(forlam, ta, tb) + forb * Ji1_func(forlam, ta, tb) );

			std1 += doma*doma*(tb-ta) + domb*domb*(tb*tb*tb-ta*ta*ta)/3.0 + doma*domb*(tb*tb-ta*ta);

			j++;
			ta = tb;
			tb = sig_time[j];
		}

		if(j>0)
		{
			doma = domsig[j-1] - (domsig[j] - domsig[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
			domb = (domsig[j] - domsig[j-1]) / (sig_time[j] - sig_time[j-1]);
			
			fora = forsig[j-1] - (forsig[j] - forsig[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
			forb = (forsig[j] - forsig[j-1]) / (sig_time[j] - sig_time[j-1]);
		}

		domI1 += doma * doma * Ji0_func(domlam, ta, t2)
					+ 2 * doma * domb * Ji1_func(domlam, ta, t2)
					+ domb * domb * Ji2_func(domlam, ta, t2);

		domI2 += doma * doma * Ji0_func(2*domlam, ta, t2)
					+ 2 * doma * domb * Ji1_func(2*domlam, ta, t2)
					+ domb * domb * Ji2_func(2*domlam, ta, t2);

		forI1 += fora * fora * Ji0_func(forlam, ta, t2)
					+ 2 * fora * forb * Ji1_func(forlam, ta, t2)
					+ forb * forb * Ji2_func(forlam, ta, t2);

		forI2 += fora * fora * Ji0_func(2*forlam, ta, t2)
					+ 2 * fora * forb * Ji1_func(2*forlam, ta, t2)
					+ forb * forb * Ji2_func(2*forlam, ta, t2);

		domforI2 += doma * fora * Ji0_func(domlam+forlam, ta, t2)
					+ (doma * forb + fora * domb) * Ji1_func(domlam+forlam, ta, t2)
					+ domb * forb * Ji2_func(domlam+forlam, ta, t2);

		QAdj += quantorho * fxsig[j]
				* ( fora * Ji0_func(forlam, ta, t2) + forb * Ji1_func(forlam, ta, t2) );

		std1 += doma*doma*(t2-ta) + domb*domb*(t2*t2*t2-ta*ta*ta)/3.0 + doma*domb*(t2*t2-ta*ta);

		fwd1[i] = (exp(-domlam*t2) * domI1 - exp(-2*domlam*t2) * domI2)/domlam;

		fwd3[i] = (exp(-forlam*t2) * forI1 - exp(-2*forlam*t2) * forI2)/forlam
					- exp(-forlam*t2) * QAdj;

		domphi[i] = exp(-2*domlam*t2) * domI2;
		forphi[i] = exp(-2*forlam*t2) * forI2;
		domforphi[i] = domforrho * exp(-(domlam+forlam)*t2) * domforI2;

		var1[i-1] = std1;
		var2[i-1] = std1;

		t1 = t2;
	}
}


static void DoubleLGM1FQuantoExpectationsWithCorrelTS(
						int			nstept,
						double		*time,
						double		domlam,
						double		forlam,
						double		*sig_time,	//	domsig, forsig and fxsig must have 
						double		*domsig,	//	the same sig_time
						double		*forsig,
						double		*fxsig,
						int			nb_sig,
						double		*quantorho,
						double		*domforrho,
						double		*fwd1,
						double		*fwd3,
						double		*var1,
						double		*var2,
						double		*domphi,
						double		*forphi,
						double		*domforphi)
{
	double	t1, t2, ta, tb;
	int		i, j, nb_sig_minus1;
	double domI1, domI2;
	double forI1, forI2;
	double domforI2;
	double domH, forH;
	double doma, fora, quantorhoa, domforrhoa;
	double domb, forb, quantorhob, domforrhob;
	double QAdj;
	double std1;

	nb_sig_minus1 = nb_sig - 1;

	// initialisation 
	t1 = 0.0;
	j = 0;

	doma = domsig[0];
	domb = 0;
	fora = forsig[0];
	forb = 0;
	quantorhoa = quantorho[0];
	quantorhob = 0;
	domforrhoa = domforrho[0];
	domforrhob = 0;

	domphi[0] = 0.0;
	forphi[0] = 0.0;
	domforphi[0] = 0.0;
	
	fwd1[0] = 0;
	fwd3[0] = 0;
	domI1 = 0;
	domI2 = 0;
	forI1 = 0;
	forI2 = 0;
	domforI2 = 0;
	domH = 0;
	forH = 0;
	QAdj = 0;

	for (i=1; i<nstept; i++)
	{		
		t2 = time[i];
		ta = t1;
		tb = sig_time[j];

		std1 = 0;

		while (tb < t2 && j < nb_sig_minus1)
		{
			if(j>0)
			{
				doma = domsig[j-1] - (domsig[j] - domsig[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
				domb = (domsig[j] - domsig[j-1]) / (sig_time[j] - sig_time[j-1]);
			
				fora = forsig[j-1] - (forsig[j] - forsig[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
				forb = (forsig[j] - forsig[j-1]) / (sig_time[j] - sig_time[j-1]);

				quantorhoa = quantorho[j-1] - (quantorho[j] - quantorho[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
				quantorhob = (quantorho[j] - quantorho[j-1]) / (sig_time[j] - sig_time[j-1]);

				domforrhoa = domforrho[j-1] - (domforrho[j] - domforrho[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
				domforrhob = (domforrho[j] - domforrho[j-1]) / (sig_time[j] - sig_time[j-1]);
			}

			domI1 += doma * doma * Ji0_func(domlam, ta, tb)
						+ 2 * doma * domb * Ji1_func(domlam, ta, tb)
						+ domb * domb * Ji2_func(domlam, ta, tb);

			domI2 += doma * doma * Ji0_func(2*domlam, ta, tb)
						+ 2 * doma * domb * Ji1_func(2*domlam, ta, tb)
						+ domb * domb * Ji2_func(2*domlam, ta, tb);

			forI1 += fora * fora * Ji0_func(forlam, ta, tb)
						+ 2 * fora * forb * Ji1_func(forlam, ta, tb)
						+ forb * forb * Ji2_func(forlam, ta, tb);

			forI2 += fora * fora * Ji0_func(2*forlam, ta, tb)
						+ 2 * fora * forb * Ji1_func(2*forlam, ta, tb)
						+ forb * forb * Ji2_func(2*forlam, ta, tb);

			domforI2 += domforrhoa * doma * fora * Ji0_func(domlam+forlam, ta, tb)
						+ (domforrhoa * doma * forb + domforrhoa * fora * domb + domforrhob * doma * fora) * Ji1_func(domlam+forlam, ta, tb)
						+ (domforrhoa * domb * forb + domforrhob * doma * forb + domforrhob * domb * fora) * Ji2_func(domlam+forlam, ta, tb)
						+ domforrhob * domb * forb * Ji3_func(domlam+forlam, ta, tb);

			QAdj += fxsig[j] * ( 
								quantorhoa * fora * Ji0_func(forlam, ta, tb) 
								+ (quantorhoa * forb + quantorhob * fora) * Ji1_func(forlam, ta, tb)
								+ quantorhob * forb * Ji2_func(forlam, ta, tb)
								);

			std1 += doma*doma*(tb-ta) + domb*domb*(tb*tb*tb-ta*ta*ta)/3.0 + doma*domb*(tb*tb-ta*ta);

			j++;
			ta = tb;
			tb = sig_time[j];
		}

		if(j>0)
		{
			doma = domsig[j-1] - (domsig[j] - domsig[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
			domb = (domsig[j] - domsig[j-1]) / (sig_time[j] - sig_time[j-1]);
			
			fora = forsig[j-1] - (forsig[j] - forsig[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
			forb = (forsig[j] - forsig[j-1]) / (sig_time[j] - sig_time[j-1]);

			quantorhoa = quantorho[j-1] - (quantorho[j] - quantorho[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
			quantorhob = (quantorho[j] - quantorho[j-1]) / (sig_time[j] - sig_time[j-1]);

			domforrhoa = domforrho[j-1] - (domforrho[j] - domforrho[j-1]) * sig_time[j-1] / (sig_time[j] - sig_time[j-1]);
			domforrhob = (domforrho[j] - domforrho[j-1]) / (sig_time[j] - sig_time[j-1]);
		}

		domI1 += doma * doma * Ji0_func(domlam, ta, t2)
					+ 2 * doma * domb * Ji1_func(domlam, ta, t2)
					+ domb * domb * Ji2_func(domlam, ta, t2);

		domI2 += doma * doma * Ji0_func(2*domlam, ta, t2)
					+ 2 * doma * domb * Ji1_func(2*domlam, ta, t2)
					+ domb * domb * Ji2_func(2*domlam, ta, t2);

		forI1 += fora * fora * Ji0_func(forlam, ta, t2)
					+ 2 * fora * forb * Ji1_func(forlam, ta, t2)
					+ forb * forb * Ji2_func(forlam, ta, t2);

		forI2 += fora * fora * Ji0_func(2*forlam, ta, t2)
					+ 2 * fora * forb * Ji1_func(2*forlam, ta, t2)
					+ forb * forb * Ji2_func(2*forlam, ta, t2);

		domforI2 += domforrhoa * doma * fora * Ji0_func(domlam+forlam, ta, t2)
						+ (domforrhoa * doma * forb + domforrhoa * fora * domb + domforrhob * doma * fora) * Ji1_func(domlam+forlam, ta, t2)
						+ (domforrhoa * domb * forb + domforrhob * doma * forb + domforrhob * domb * fora) * Ji2_func(domlam+forlam, ta, t2)
						+ domforrhob * domb * forb * Ji3_func(domlam+forlam, ta, t2);

		QAdj += fxsig[j] * ( 
							quantorhoa * fora * Ji0_func(forlam, ta, t2) 
							+ (quantorhoa * forb + quantorhob * fora) * Ji1_func(forlam, ta, t2)
							+ quantorhob * forb * Ji2_func(forlam, ta, t2)
							);

		std1 += doma*doma*(t2-ta) + domb*domb*(t2*t2*t2-ta*ta*ta)/3.0 + doma*domb*(t2*t2-ta*ta);

		fwd1[i] = (exp(-domlam*t2) * domI1 - exp(-2*domlam*t2) * domI2)/domlam;

		fwd3[i] = (exp(-forlam*t2) * forI1 - exp(-2*forlam*t2) * forI2)/forlam
					- exp(-forlam*t2) * QAdj;

		domphi[i] = exp(-2*domlam*t2) * domI2;
		forphi[i] = exp(-2*forlam*t2) * forI2;
		domforphi[i] = exp(-(domlam+forlam)*t2) * domforI2;

		var1[i-1] = std1;
		var2[i-1] = std1;

		t1 = t2;
	}
}



static void DoubleLGM1FQuantoExpectationsWithCorrelationTS(
						int			nstept,
						double		*time,
						double		domlam,
						double		forlam,
						double		*sig_time,	//	domsig, forsig and fxsig must have 
						double		*domsig,	//	the same sig_time
						double		*forsig,
						double		*fxsig,
						int			nb_sig,
						double		*quantorho,
						double		*domforrho,
						double		*fwd1,
						double		*fwd3,
						double		*var1,
						double		*var2,
						double		*domphi,
						double		*forphi,
						double		*domforphi)
{
	double	t1, t2, ta, tb;
	int		i, j, nb_sig_minus1;
	double domI1, domI2;
	double forI1, forI2;
	double domforI2;
	double domH, forH;
	double doma, fora, quantorhoa, domforrhoa;
	double domb, forb, quantorhob, domforrhob;
	double QAdj;
	double std1;

	nb_sig_minus1 = nb_sig - 1;

	// initialisation 
	t1 = 0.0;
	j = 0;

	doma = domsig[0];
	domb = 0;
	fora = forsig[0];
	forb = 0;
	quantorhoa = quantorho[0];
	quantorhob = 0;
	domforrhoa = domforrho[0];
	domforrhob = 0;

	domphi[0] = 0.0;
	forphi[0] = 0.0;
	domforphi[0] = 0.0;
	
	fwd1[0] = 0;
	fwd3[0] = 0;
	domI1 = 0;
	domI2 = 0;
	forI1 = 0;
	forI2 = 0;
	domforI2 = 0;
	domH = 0;
	forH = 0;
	QAdj = 0;

	for (i=1; i<nstept; i++)
	{		
		t2 = time[i];
		ta = t1;
		tb = sig_time[j];

		std1 = 0;

		while (tb < t2 && j < nb_sig_minus1)
		{
 
			doma = domsig[j];
			domb = 0;
			
			fora = forsig[j];
			forb = 0;

			quantorhoa = quantorho[j];
			quantorhob = 0;

			domforrhoa = domforrho[j];
			domforrhob = 0;

			domI1 += doma * doma * Ji0_func(domlam, ta, tb)
						+ 2 * doma * domb * Ji1_func(domlam, ta, tb)
						+ domb * domb * Ji2_func(domlam, ta, tb);

			domI2 += doma * doma * Ji0_func(2*domlam, ta, tb)
						+ 2 * doma * domb * Ji1_func(2*domlam, ta, tb)
						+ domb * domb * Ji2_func(2*domlam, ta, tb);

			forI1 += fora * fora * Ji0_func(forlam, ta, tb)
						+ 2 * fora * forb * Ji1_func(forlam, ta, tb)
						+ forb * forb * Ji2_func(forlam, ta, tb);

			forI2 += fora * fora * Ji0_func(2*forlam, ta, tb)
						+ 2 * fora * forb * Ji1_func(2*forlam, ta, tb)
						+ forb * forb * Ji2_func(2*forlam, ta, tb);

			domforI2 += domforrhoa * doma * fora * Ji0_func(domlam+forlam, ta, tb)
						+ (domforrhoa * doma * forb + domforrhoa * fora * domb + domforrhob * doma * fora) * Ji1_func(domlam+forlam, ta, tb)
						+ (domforrhoa * domb * forb + domforrhob * doma * forb + domforrhob * domb * fora) * Ji2_func(domlam+forlam, ta, tb)
						+ domforrhob * domb * forb * Ji3_func(domlam+forlam, ta, tb);

			QAdj += fxsig[j] * ( 
								quantorhoa * fora * Ji0_func(forlam, ta, tb) 
								+ (quantorhoa * forb + quantorhob * fora) * Ji1_func(forlam, ta, tb)
								+ quantorhob * forb * Ji2_func(forlam, ta, tb)
								);

			std1 += doma*doma*(tb-ta) + domb*domb*(tb*tb*tb-ta*ta*ta)/3.0 + doma*domb*(tb*tb-ta*ta);

			j++;
			ta = tb;
			tb = sig_time[j];
		}

		doma = domsig[j];
		domb = 0;
			
		fora = forsig[j];
		forb = 0;

		quantorhoa = quantorho[j];
		quantorhob = 0;

		domforrhoa = domforrho[j];
		domforrhob = 0;

		domI1 += doma * doma * Ji0_func(domlam, ta, t2)
					+ 2 * doma * domb * Ji1_func(domlam, ta, t2)
					+ domb * domb * Ji2_func(domlam, ta, t2);

		domI2 += doma * doma * Ji0_func(2*domlam, ta, t2)
					+ 2 * doma * domb * Ji1_func(2*domlam, ta, t2)
					+ domb * domb * Ji2_func(2*domlam, ta, t2);

		forI1 += fora * fora * Ji0_func(forlam, ta, t2)
					+ 2 * fora * forb * Ji1_func(forlam, ta, t2)
					+ forb * forb * Ji2_func(forlam, ta, t2);

		forI2 += fora * fora * Ji0_func(2*forlam, ta, t2)
					+ 2 * fora * forb * Ji1_func(2*forlam, ta, t2)
					+ forb * forb * Ji2_func(2*forlam, ta, t2);

		domforI2 += domforrhoa * doma * fora * Ji0_func(domlam+forlam, ta, t2)
						+ (domforrhoa * doma * forb + domforrhoa * fora * domb + domforrhob * doma * fora) * Ji1_func(domlam+forlam, ta, t2)
						+ (domforrhoa * domb * forb + domforrhob * doma * forb + domforrhob * domb * fora) * Ji2_func(domlam+forlam, ta, t2)
						+ domforrhob * domb * forb * Ji3_func(domlam+forlam, ta, t2);

		QAdj += fxsig[j] * ( 
							quantorhoa * fora * Ji0_func(forlam, ta, t2) 
							+ (quantorhoa * forb + quantorhob * fora) * Ji1_func(forlam, ta, t2)
							+ quantorhob * forb * Ji2_func(forlam, ta, t2)
							);

		std1 += doma*doma*(t2-ta) + domb*domb*(t2*t2*t2-ta*ta*ta)/3.0 + doma*domb*(t2*t2-ta*ta);

		fwd1[i] = (exp(-domlam*t2) * domI1 - exp(-2*domlam*t2) * domI2)/domlam;

		fwd3[i] = (exp(-forlam*t2) * forI1 - exp(-2*forlam*t2) * forI2)/forlam
					- exp(-forlam*t2) * QAdj;

		domphi[i] = exp(-2*domlam*t2) * domI2;
		forphi[i] = exp(-2*forlam*t2) * forI2;
		domforphi[i] = exp(-(domlam+forlam)*t2) * domforI2;

		var1[i-1] = std1;
		var2[i-1] = std1;

		t1 = t2;
	}
}


/*----Same function as DoubleLGM1FQuantoExpectations but with quadratic vols---*/

static void DoubleLGM1FQuantoExpectationsQuadratic(
						int			nstept,
						double		*time,
						double		domlam,
						double		forlam,
						double		*sig_time,	//	domsig, forsig and fxsig must have 
						double		*domsig,	//	the same sig_time
						double		*domalpha,
						double		*forsig,
						double		*foralpha,
						double		*fxsig,
						int			nb_sig,
						double		quantorho,
						double		domforrho,
						double		*fwd1,
						double		*fwd3,
						double		*var1,
						double		*var2,
						double		*domphi,
						double		*forphi,
						double		*domforphi)
{
	double	t1, t2, ta, tb;
	int		i, j, nb_sig_minus1;
	double domI1, domI2;
	double forI1, forI2;
	double domforI2;
	double domH, forH;
	double doma, fora;
	double domb, forb;
	double domc, forc;
	double QAdj;
	double std1;

	nb_sig_minus1 = nb_sig - 1;

	// initialisation 
	t1 = 0.0;
	j = 0;

	doma = 0;
	domb = 0;
	domc = domsig[0];
	fora = 0;
	forb = 0;
	forc = forsig[0];

	domphi[0] = 0.0;
	forphi[0] = 0.0;
	domforphi[0] = 0.0;
	
	fwd1[0] = 0;
	fwd3[0] = 0;
	domI1 = 0;
	domI2 = 0;
	forI1 = 0;
	forI2 = 0;
	domforI2 = 0;
	domH = 0;
	forH = 0;
	QAdj = 0;

	for (i=1; i<nstept; i++)
	{		
		t2 = time[i];
		ta = t1;
		tb = sig_time[j];

		std1 = 0;

		while (tb < t2 && j < nb_sig_minus1)
		{
			if(j>0)
			{
				doma = domalpha[j];
				domb = (domsig[j]-domsig[j-1])/(sig_time[j]-sig_time[j-1]) - domalpha[j]*(sig_time[j]+sig_time[j-1]);
				domc = domsig[j-1] + domalpha[j]*sig_time[j]*sig_time[j-1] - sig_time[j-1]*(domsig[j]-domsig[j-1])/(sig_time[j]-sig_time[j-1]);
		
				fora = foralpha[j];
				forb = (forsig[j]-forsig[j-1])/(sig_time[j]-sig_time[j-1]) - foralpha[j]*(sig_time[j]+sig_time[j-1]);
				forc = forsig[j-1] + foralpha[j]*sig_time[j]*sig_time[j-1] - sig_time[j-1]*(forsig[j]-forsig[j-1])/(sig_time[j]-sig_time[j-1]);
			}

			domI1 += doma * doma * Ji4_func(domlam, ta, tb)
						+ 2 * doma * domb * Ji3_func(domlam, ta, tb)
						+ (2 * doma * domc + domb * domb) * Ji2_func(domlam, ta, tb)
						+ 2 * domb * domc * Ji1_func(domlam, ta, tb)
						+ domc * domc * Ji0_func(domlam, ta, tb);

			domI2 += doma * doma * Ji4_func(2*domlam, ta, tb)
						+ 2 * doma * domb * Ji3_func(2*domlam, ta, tb)
						+ (2 * doma * domc + domb * domb) * Ji2_func(2*domlam, ta, tb)
						+ 2 * domb * domc * Ji1_func(2*domlam, ta, tb)
						+ domc * domc * Ji0_func(2*domlam, ta, tb);

			forI1 += fora * fora * Ji4_func(forlam, ta, tb)
						+ 2 * fora * forb * Ji3_func(forlam, ta, tb)
						+ (2 * fora * forc + forb * forb) * Ji2_func(forlam, ta, tb)
						+ 2 * forb * forc * Ji1_func(forlam, ta, tb)
						+ forc * forc * Ji0_func(forlam, ta, tb);

			forI2 += fora * fora * Ji4_func(2*forlam, ta, tb)
						+ 2 * fora * forb * Ji3_func(2*forlam, ta, tb)
						+ (2 * fora * forc + forb * forb) * Ji2_func(2*forlam, ta, tb)
						+ 2 * forb * forc * Ji1_func(2*forlam, ta, tb)
						+ forc * forc * Ji0_func(2*forlam, ta, tb);

			domforI2 += doma * fora * Ji4_func(domlam+forlam, ta, tb)
						+ (doma * forb + domb * fora)* Ji3_func(domlam+forlam, ta, tb)
						+ (doma * forc + domb * forb + domc * fora) * Ji2_func(domlam+forlam, ta, tb)
						+ (domb * forc + domc * forb) * Ji1_func(domlam+forlam, ta, tb)
						+ domc * forc * Ji0_func(domlam+forlam, ta, tb);

			QAdj += quantorho * fxsig[j]
				* ( fora * Ji2_func(forlam, ta, tb)
					+ forb * Ji1_func(forlam, ta, tb)
					+ forc * Ji0_func(forlam, ta, tb));

			std1 += doma*doma*(tb*tb*tb*tb*tb-ta*ta*ta*ta*ta)/5.0
					+2*doma*domb*(tb*tb*tb*tb-ta*ta*ta*ta)/4.0
					+(2*doma*domc+domb*domb)*(tb*tb*tb-ta*ta*ta)/3.0
					+2*domb*domc*(tb*tb-ta*ta)/2.0
					+domc*domc*(tb-ta);

			j++;
			ta = tb;
			tb = sig_time[j];
		}

		if(j>0)
		{
			doma = domalpha[j];
			domb = (domsig[j]-domsig[j-1])/(sig_time[j]-sig_time[j-1]) - domalpha[j]*(sig_time[j]+sig_time[j-1]);
			domc = domsig[j-1] + domalpha[j]*sig_time[j]*sig_time[j-1] - sig_time[j-1]*(domsig[j]-domsig[j-1])/(sig_time[j]-sig_time[j-1]);
		
			fora = foralpha[j];
			forb = (forsig[j]-forsig[j-1])/(sig_time[j]-sig_time[j-1]) - foralpha[j]*(sig_time[j]+sig_time[j-1]);
			forc = forsig[j-1] + foralpha[j]*sig_time[j]*sig_time[j-1] - sig_time[j-1]*(forsig[j]-forsig[j-1])/(sig_time[j]-sig_time[j-1]);
		}

			domI1 += doma * doma * Ji4_func(domlam, ta, t2)
						+ 2 * doma * domb * Ji3_func(domlam, ta, t2)
						+ (2 * doma * domc + domb * domb) * Ji2_func(domlam, ta, t2)
						+ 2 * domb * domc * Ji1_func(domlam, ta, t2)
						+ domc * domc * Ji0_func(domlam, ta, t2);

			domI2 += doma * doma * Ji4_func(2*domlam, ta, t2)
						+ 2 * doma * domb * Ji3_func(2*domlam, ta, t2)
						+ (2 * doma * domc + domb * domb) * Ji2_func(2*domlam, ta, t2)
						+ 2 * domb * domc * Ji1_func(2*domlam, ta, t2)
						+ domc * domc * Ji0_func(2*domlam, ta, t2);

			forI1 += fora * fora * Ji4_func(forlam, ta, t2)
						+ 2 * fora * forb * Ji3_func(forlam, ta, t2)
						+ (2 * fora * forc + forb * forb) * Ji2_func(forlam, ta, t2)
						+ 2 * forb * forc * Ji1_func(forlam, ta, t2)
						+ forc * forc * Ji0_func(forlam, ta, t2);

			forI2 += fora * fora * Ji4_func(2*forlam, ta, t2)
						+ 2 * fora * forb * Ji3_func(2*forlam, ta, t2)
						+ (2 * fora * forc + forb * forb) * Ji2_func(2*forlam, ta, t2)
						+ 2 * forb * forc * Ji1_func(2*forlam, ta, t2)
						+ forc * forc * Ji0_func(2*forlam, ta, t2);

			domforI2 += doma * fora * Ji4_func(domlam+forlam, ta, t2)
						+ (doma * forb + domb * fora)* Ji3_func(domlam+forlam, ta, t2)
						+ (doma * forc + domb * forb + domc * fora) * Ji2_func(domlam+forlam, ta, t2)
						+ (domb * forc + domc * forb) * Ji1_func(domlam+forlam, ta, t2)
						+ domc * forc * Ji0_func(domlam+forlam, ta, t2);

			QAdj += quantorho * fxsig[j]
				* ( fora * Ji2_func(forlam, ta, t2)
					+ forb * Ji1_func(forlam, ta, t2)
					+ forc * Ji0_func(forlam, ta, t2));

			std1 += doma*doma*(t2*t2*t2*t2*t2-ta*ta*ta*ta*ta)/5.0
					+2*doma*domb*(t2*t2*t2*t2-ta*ta*ta*ta)/4.0
					+(2*doma*domc+domb*domb)*(t2*t2*t2-ta*ta*ta)/3.0
					+2*domb*domc*(t2*t2-ta*ta)/2.0
					+domc*domc*(t2-ta);

		fwd1[i] = (exp(-domlam*t2) * domI1 - exp(-2*domlam*t2) * domI2)/domlam;

		fwd3[i] = (exp(-forlam*t2) * forI1 - exp(-2*forlam*t2) * forI2)/forlam
					- exp(-forlam*t2) * QAdj;

		domphi[i] = exp(-2*domlam*t2) * domI2;
		forphi[i] = exp(-2*forlam*t2) * forI2;
		domforphi[i] = domforrho * exp(-(domlam+forlam)*t2) * domforI2;

		var1[i-1] = std1;
		var2[i-1] = std1;

		t1 = t2;
	}
}


/*	---------------------------- */
/*	Main function without Tau TS */
/*	---------------------------- */

Err	 doublelgm1fQuanto_adi(	
						//Time data		
					int			nstept,
					double		*time,
					double		*date,

						//Discretisation	
					int			nsteps,
										
						//Model data		
					double		domlam,
					double		forlam,
					double		*sig_time,	//	domsig, forsig and fxsig must have 
					double		*domsig,	//	the same sig_time
					double		*forsig,
					double		*fxsig,
					int			nb_sig,
					double		quantorho,
					double		domforrho,

						//Product data
					void		**func_parm_tab, 
					int			*eval_evt,
					
						//Market data 
					double		*dom_ifr,
					double		*for_ifr,
					char		*dom_yc,
					char		*for_yc,
					
						//Payoff function 
					Err (*payoff_func)( //Event 
									double	evt_date,
									double	evt_time,
									void	*func_parm, 
					
									/* Market data	*/										
									void	*dom_yc,
									void	*for_yc,
									
									/* Model data	*/
									double	domlam,
									double	forlam,
									double	domphi,
									double	forphi,
									
									/* Grid data	*/
									int		l1,
									int		u1,
									int		l2,
									int		u2,
									double	*r1,
									double	**r2,
															
									/* Vector of results to be updated */
									int		nprod,
									double	***prod_val
									),
						//Result 
					int			nprod, 
					double		*res)
{

Err				err = NULL;
				
int				i, j, step, index_x, index_z;
int				nstepx, nstepz;
double			meshx, meshz, dt;
double			mu_r1i, mu_r3i;
double			r_temp;

double			std1, std3, r2_temp;

double			*r1_bar			= NULL,
				*r3_bar			= NULL,
				*r1_dim1		= NULL,
				**r1_dim2		= NULL,
				**r2			= NULL,
				***values		= NULL,
				***values_p1	= NULL,
				***values_temp	= NULL,
				**mu_r1			= NULL,
				**mu_r3			= NULL,
				**var_r1		= NULL,
				**var_r3		= NULL,
				**r				= NULL,
				**r_init		= NULL,
				*fwd1			= NULL,
				*fwd3			= NULL,
				*var1			= NULL,
				*var3			= NULL,
				*domphi			= NULL,
				*forphi			= NULL,
				*domforphi		= NULL;

int				lx, ux, lz, uz, tem;

double *domsigtilda=NULL;
double *forsigtilda=NULL;
double *contdomsig=NULL;
double *contforsig=NULL;
double *allfxsig=NULL;
double **domabsig=NULL;
double **forabsig=NULL;

clock_t			t1, t2;

CNPDE_TEMP_2D_ADI	pdestr, *pde = &pdestr;

	t1 = clock();

		//Allocations	of time vectors						
	domphi = dvector(0, nstept-1);
	forphi = dvector(0, nstept-1);
	domforphi = dvector(0, nstept-1);
	var1 = dvector(0, nstept-1);
	var3 = dvector(0, nstept-1);
	fwd1 = dvector(0, nstept-1);
	fwd3 = dvector(0, nstept-1);
	domsigtilda = dvector(0, nb_sig-1);
	forsigtilda = dvector(0, nb_sig-1);
	contdomsig = dvector(0, nstept-1);
	contforsig = dvector(0, nstept-1);
	allfxsig = dvector(0, nstept-1);
	domabsig = dmatrix(0, 1, 0, nstept-1);
	forabsig = dmatrix(0, 1, 0, nstept-1);

	if (!pde || !var1 || !var3 || !domphi || !forphi || !domforphi || !fwd1 || !fwd3 || !domsigtilda || !forsigtilda || !allfxsig || !contdomsig || !contforsig || !domabsig || !forabsig)
	{
		err = "Memory allocation error (1) in doublelgm1fQuanto_adi";
		goto FREE_RETURN;
	}

	err = convertVol2(domlam, nstept, time, nb_sig, sig_time, domsig, domsigtilda, contdomsig, domabsig);
	if (err)
	{
		goto FREE_RETURN;
	}
	err = convertVol2(forlam, nstept, time, nb_sig, sig_time, forsig, forsigtilda, contforsig, forabsig);
	if (err)
	{
		goto FREE_RETURN;
	}
	fillVol(nstept, time, nb_sig, sig_time, fxsig, allfxsig);

		//Calculate the expecations of r1 and r3 
	DoubleLGM1FQuantoExpectations(nstept,
								time,
								domlam,
								forlam,
								sig_time,
								domsigtilda,
								forsigtilda,
								fxsig,
								nb_sig,
								quantorho,
								domforrho,
								fwd1,
								fwd3,
								var1,
								var3,
								domphi,
								forphi,
								domforphi);

		//Calculation of the number of steps in each direction: since local volatility			
		//is the same, the mesh has to be the same, but the number of steps has to be adjusted	

	std1 = sqrt(domphi[nstept-1]);
	std3 = sqrt(
					( (contdomsig[nstept-1]/contforsig[nstept-1])*(contdomsig[nstept-1]/contforsig[nstept-1])*forphi[nstept-1]
						+ domforrho*domforrho*domphi[nstept-1]
						- 2*domforrho*(contdomsig[nstept-1]/contforsig[nstept-1])*domforphi[nstept-1]
					)/(1-domforrho*domforrho)
				);

	nstepx = (int) (nsteps * sqrt(std1 / std3) + 0.5);
	nstepz = (int) (nsteps * sqrt(std3 / std1) + 0.5);
		
//	nstepx = nsteps;
//	nstepz = nsteps;

	//nstep has to be a odd nuber			
	nstepx = ((int) (nstepx / 2)) * 2 + 1;
	nstepz = ((int) (nstepz / 2)) * 2 + 1;

		//we want at least three points in each directions										
	if (nstepx < 3)
	{
		nstepx = 3;
		nstepz = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepz = ((int) (nstepz / 2)) * 2 + 1;

		if (nstepz < 3)
		{
			nstepz = 3;		
		}
	}

	if (nstepz < 3)
	{
		nstepz = 3;
		nstepx = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepx = ((int) (nstepx / 2)) * 2 + 1;

		if (nstepx < 3)
		{
			nstepx = 3;		
		}
	}

		//corresponding index to the 0 value of x and z	
	index_x = (nstepx - 1) / 2;
	index_z = (nstepz - 1) / 2;
	
	meshx = 2.0 * NSTD_LGM * std1 / (nstepx - 1);
	meshz = 2.0 * NSTD_LGM * std3 / (nstepz - 1);
		
		//Allocations of space vectors					
	r1_bar = dvector(0, nstepx-1);
	r3_bar = dvector(0, nstepz-1);
	r1_dim1 = dvector(0, nstepx-1);
	values = f3tensor(0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	values_p1 = f3tensor(0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	mu_r1 = dmatrix(0, nstepx-1, 0, nstepz-1);
	mu_r3 = dmatrix(0, nstepx-1, 0, nstepz-1);
	var_r1 = dmatrix(0, nstepx-1, 0, nstepz-1);
	var_r3 = dmatrix(0, nstepx-1, 0, nstepz-1);
	r = dmatrix(0, nstepx-1, 0, nstepz-1);
	r_init = dmatrix(0, nstepx-1, 0, nstepz-1);
	r1_dim2 = dmatrix(0, nstepx-1, 0, nstepz-1);
	r2 = dmatrix(0, nstepx-1, 0, nstepz-1);

	if (!r1_dim2 || !r2 || !values || !values_p1 || !mu_r1 || !mu_r3 || !var_r1 || !var_r3 || !r || !r_init ||
		!r1_bar || !r1_dim1 || !r3_bar)
	{
		err = "Memory allocation error (2) in doublelgm1fQuanto_adi";
		goto FREE_RETURN;
}

	t2 = clock();

	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	smessage ("Phase 2 -convolution, stept: %d stepx: %d stepz: %d", nstept, nstepx, nstepz);

	t1 = clock();

		//Then discretise space in the orthogonal system r1 / r3	
						
	r1_bar[0] = -(nstepx - 1) / 2.0 * meshx;
	r3_bar[0] = -(nstepz - 1) / 2.0 * meshz;

	for (i=1; i<nstepx; i++)
	{
		r1_bar[i] = r1_bar[i-1] + meshx;
	}

	for (j=1; j<nstepz; j++)
	{
		r3_bar[j] = r3_bar[j-1] + meshz;
	}

		//Corresponding dom_r and for_r					
	for (i=0; i<nstepx; i++)
	{
		r1_dim1[i] = r1_bar[i] + fwd1[nstept-1];
		r2_temp = domforrho * r1_bar[i] * (contforsig[nstept-1]/contdomsig[nstept-1]);
		r_temp = domforrho * (1 + (contforsig[nstept-1]/contdomsig[nstept-1])) * r1_bar[i];

		for (j=0; j<nstepz; j++)
		{
			r1_dim2[i][j] = r1_dim1[i];
			r2[i][j] = fwd3[nstept-1] + r2_temp + (contforsig[nstept-1]/contdomsig[nstept-1]) * sqrt(1-domforrho*domforrho) * r3_bar[j];
			r_init[i][j] = r1_bar[i];
		}
	}

		//Final payoff valuation					
	if (!eval_evt[nstept-1])
	{
		err = "No event at last step in lgm2f_pde";
		goto FREE_RETURN;
	}
	
	err = payoff_func (	date[nstept-1],
						time[nstept-1],
						func_parm_tab[nstept-1],
						dom_yc,
						for_yc,
						domlam,
						forlam,
						domphi[nstept-1],
						forphi[nstept-1],
						0,
						nstepx-1,
						0,
						nstepz-1,
						r1_dim1,						
						r2,
						nprod,
						values_p1
						);

	if (err)
	{
		goto FREE_RETURN;
	}

		//Initialize the CNPDE_TEMP_2D		

	num_f_pde_init_2d_adi(	pde,
							nstepx,
							nstepz,
							nprod
							);
	
	if (!pde)
	{
		err = "Memory allocation error (3) in lgm2fpde";
		goto FREE_RETURN;
	}
			
	lx = 0;
	ux = nstepx - 1;
	lz = 0;
	uz = nstepz - 1;

		//now do the backward pde					

	for (step=nstept-2; step>=0; step--)
	{
		dt = time[step+1] - time[step];
						
		r_temp = dom_ifr[step] + fwd1[step];

		for (i=lx; i<=ux; i++)
		{							
			mu_r1i = -domlam * r1_bar[i] * dt;
			mu_r3i = -(forlam - domlam) * domforrho * r1_bar[i] / sqrt(1-domforrho*domforrho);

			for (j=lz; j<=uz; j++)
			{
				mu_r1[i][j] = mu_r1i;
				mu_r3[i][j] = (
								mu_r3i - forlam * r3_bar[j]
								+ (r3_bar[j] + domforrho * r1_bar[i] / sqrt(1-domforrho*domforrho) )
									*(domabsig[1][step]/contdomsig[step] - forabsig[1][step]/contforsig[step])
							  )* dt;
				
				var_r1[i][j] = var1[step];
				var_r3[i][j] = var3[step];
												
				r[i][j] = (r_init[i][j] + r_temp) * dt;
			}
		}

			//convolve
		
		num_f_pde_one_step_backward_2f_adi(	pde,
											nstepx,
											r1_bar,
											nstepz,
											r3_bar,
											0,
											nprod-1,
											values_p1,
											mu_r1,
											mu_r3,
											var_r1,
											var_r3,
											r,
											values,
											lx,
											ux,
											lz,
											uz);
		
			//Eval payoff 
		if (eval_evt[step])
		{

				//Corresponding r1 and r2					
			for (i=lx; i<=ux; i++)
			{
				r1_dim1[i] = r1_bar[i] + fwd1[step];
				r2_temp = domforrho * r1_bar[i] * (contforsig[step]/contdomsig[step]);

				for (j=0; j<nstepz; j++)
				{
					r1_dim2[i][j] = r1_dim1[i];
					r2[i][j] = fwd3[step] + r2_temp + (contforsig[step]/contdomsig[step]) * sqrt(1-domforrho*domforrho) * r3_bar[j];
				}
			}

			err = payoff_func (	date[step],
						time[step],
						func_parm_tab[step],
						dom_yc,
						for_yc,
						domlam,
						forlam,
						domphi[step],
						forphi[step],
						lx,
						ux,
						lz,
						uz,
						r1_dim1,						
						r2,
						nprod,
						values
						);

			if (err)
			{
				goto FREE_RETURN;
			}
		}

		values_temp = values_p1;
		values_p1 = values;
		values = values_temp;

		 //new indexes: we cut the PDE after NSTD_LGM number of standard deviations 
		 //but we need at least three points to do the pde						

		tem = (int) (NSTD_LGM * sqrt(domphi[step]) / meshx + 0.5);
		ux = min(nstepx-1, index_x + tem);
		lx = max(0, index_x - tem);
		
		if (ux - lx < 2)
		{
			tem += 1;
			ux = min(nstepx-1, index_x + tem);
			lx = max(0, index_x - tem);
		}

		tem = (int) (NSTD_LGM * std3 / meshx + 0.5);			
		uz = min(nstepz-1, index_z + tem);
		lz = max(0, index_z - tem);

		if (uz - lz < 2)
		{
			tem += 1;
			uz = min(nstepz-1, index_z + tem);
			lz = max(0, index_z - tem);
		}
	}

	 //copy the result					
	for (i=0; i<nprod; i++)
	{
		res[i] = values_p1[index_x][index_z][i];
	}

	t2 = clock();

	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);


FREE_RETURN:

		//allocation (1)		
	if (domphi) free_dvector(domphi, 0, nstept-1);
	if (forphi) free_dvector(forphi, 0, nstept-1);
	if (domforphi) free_dvector(domforphi, 0, nstept-1);
	if (var1) free_dvector(var1, 0, nstept-1);
	if (var3) free_dvector(var3, 0, nstept-1);
	if (fwd1) free_dvector(fwd1, 0, nstept-1);
	if (fwd3) free_dvector(fwd3, 0, nstept-1);
	if (domsigtilda) free_dvector(domsigtilda, 0, nb_sig-1);
	if (forsigtilda) free_dvector(forsigtilda, 0, nb_sig-1);
	if (contdomsig) free_dvector(contdomsig, 0, nstept-1);
	if (contforsig) free_dvector(contforsig, 0, nstept-1);
	if (allfxsig) free_dvector(allfxsig, 0, nstept-1);
	if (domabsig) free_dmatrix(domabsig, 0, 1, 0, nstept-1);
	if (forabsig) free_dmatrix(forabsig, 0, 1, 0, nstept-1);

		//allocation (2)		
	if (pde) num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

	if (r1_bar) free_dvector(r1_bar, 0, nstepx-1);
	if (r3_bar) free_dvector(r3_bar, 0, nstepz-1);
	if (r1_dim1) free_dvector(r1_dim1, 0, nstepx-1);
	if (r1_dim2) free_dmatrix(r1_dim2, 0, nstepx-1, 0, nstepz-1);
	if (r2) free_dmatrix(r2, 0, nstepx-1, 0, nstepz-1);
	if (values) free_f3tensor(values, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	if (values_p1) free_f3tensor(values_p1, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	if (mu_r1) free_dmatrix(mu_r1, 0, nstepx-1, 0, nstepz-1);
	if (mu_r3) free_dmatrix(mu_r3, 0, nstepx-1, 0, nstepz-1);
	if (var_r1) free_dmatrix(var_r1, 0, nstepx-1, 0, nstepz-1);
	if (var_r3) free_dmatrix(var_r3, 0, nstepx-1, 0, nstepz-1);
	if (r) free_dmatrix(r, 0, nstepx-1, 0, nstepz-1);
	if (r_init) free_dmatrix(r_init, 0, nstepx-1, 0, nstepz-1);
			
	return err;
}

/*---------------------------------------------------------------------*/
/*---Same function as doublelgm1fQuanto_adi but with quadratic vols----*/
/*---------------------------------------------------------------------*/

Err	 doublelgm1fQuanto_adi2(	
						//Time data		
					int			nstept,
					double		*time,
					double		*date,

						//Discretisation	
					int			nsteps,
										
						//Model data		
					double		domlam,
					double		forlam,
					double		*sig_time,	//	domsig, forsig and fxsig must have 
					double		*domsig,	//	the same sig_time
					double		*forsig,
					double		*fxsig,
					int			nb_sig,
					double		quantorho,
					double		domforrho,

						//Product data
					void		**func_parm_tab, 
					int			*eval_evt,
					
						//Market data 
					double		*dom_ifr,
					double		*for_ifr,
					char		*dom_yc,
					char		*for_yc,
					
						//Payoff function 
					Err (*payoff_func)( //Event 
									double	evt_date,
									double	evt_time,
									void	*func_parm, 
					
									/* Market data	*/										
									void	*dom_yc,
									void	*for_yc,
									
									/* Model data	*/
									double	domlam,
									double	forlam,
									double	domphi,
									double	forphi,
									
									/* Grid data	*/
									int		l1,
									int		u1,
									int		l2,
									int		u2,
									double	*r1,
									double	**r2,
															
									/* Vector of results to be updated */
									int		nprod,
									double	***prod_val
									),
						//Result 
					int			nprod, 
					double		*res)
{

Err				err = NULL;
				
int				i, j, step, index_x, index_z;
int				nstepx, nstepz;
double			meshx, meshz, dt;
double			mu_r1i, mu_r3i;
double			r_temp;

double			std1, std3, r2_temp;

double			*r1_bar			= NULL,
				*r3_bar			= NULL,
				*r1_dim1		= NULL,
				**r1_dim2		= NULL,
				**r2			= NULL,
				***values		= NULL,
				***values_p1	= NULL,
				***values_temp	= NULL,
				**mu_r1			= NULL,
				**mu_r3			= NULL,
				**var_r1		= NULL,
				**var_r3		= NULL,
				**r				= NULL,
				**r_init		= NULL,
				*fwd1			= NULL,
				*fwd3			= NULL,
				*var1			= NULL,
				*var3			= NULL,
				*domphi			= NULL,
				*forphi			= NULL,
				*domforphi		= NULL;

int				lx, ux, lz, uz, tem;

double *domalpha=NULL;
double *foralpha=NULL;
double *contdomsig=NULL;
double *contforsig=NULL;
double *allfxsig=NULL;
double **domabcsig=NULL;
double **forabcsig=NULL;

clock_t			t1, t2;

CNPDE_TEMP_2D_ADI	pdestr, *pde = &pdestr;

	t1 = clock();

		//Allocations	of time vectors						
	domphi = dvector(0, nstept-1);
	forphi = dvector(0, nstept-1);
	domforphi = dvector(0, nstept-1);
	var1 = dvector(0, nstept-1);
	var3 = dvector(0, nstept-1);
	fwd1 = dvector(0, nstept-1);
	fwd3 = dvector(0, nstept-1);
	domalpha = dvector(0, nb_sig-1);
	foralpha = dvector(0, nb_sig-1);
	contdomsig = dvector(0, nstept-1);
	contforsig = dvector(0, nstept-1);
	allfxsig = dvector(0, nstept-1);
	domabcsig = dmatrix(0, 2, 0, nstept-1);
	forabcsig = dmatrix(0, 2, 0, nstept-1);

	if (!pde || !var1 || !var3 || !domphi || !forphi || !domforphi || !fwd1 || !fwd3 || !domalpha || !foralpha || !allfxsig || !contdomsig || !contforsig || !domabcsig || !forabcsig)
	{
		err = "Memory allocation error (1) in doublelgm1fQuanto_adi2";
		goto FREE_RETURN;
	}

	err = convertVolQuadratic(domlam, nstept, time, nb_sig, sig_time, domsig, domalpha, contdomsig, domabcsig);
	if (err)
	{
		goto FREE_RETURN;
	}
	err = convertVolQuadratic(forlam, nstept, time, nb_sig, sig_time, forsig, foralpha, contforsig, forabcsig);
	if (err)
	{
		goto FREE_RETURN;
	}
	fillVol(nstept, time, nb_sig, sig_time, fxsig, allfxsig);

		//Calculate the expecations of r1 and r3 
	DoubleLGM1FQuantoExpectationsQuadratic(nstept,
								time,
								domlam,
								forlam,
								sig_time,
								domsig,
								domalpha,
								forsig,
								foralpha,
								fxsig,
								nb_sig,
								quantorho,
								domforrho,
								fwd1,
								fwd3,
								var1,
								var3,
								domphi,
								forphi,
								domforphi);

		//Calculation of the number of steps in each direction: since local volatility			
		//is the same, the mesh has to be the same, but the number of steps has to be adjusted	

	std1 = sqrt(domphi[nstept-1]);
	std3 = sqrt(
					( (contdomsig[nstept-1]/contforsig[nstept-1])*(contdomsig[nstept-1]/contforsig[nstept-1])*forphi[nstept-1]
						+ domforrho*domforrho*domphi[nstept-1]
						- 2*domforrho*(contdomsig[nstept-1]/contforsig[nstept-1])*domforphi[nstept-1]
					)/(1-domforrho*domforrho)
				);

	nstepx = (int) (nsteps * sqrt(std1 / std3) + 0.5);
	nstepz = (int) (nsteps * sqrt(std3 / std1) + 0.5);
		
	//nstep has to be a odd nuber			
	nstepx = ((int) (nstepx / 2)) * 2 + 1;
	nstepz = ((int) (nstepz / 2)) * 2 + 1;

		//we want at least three points in each directions										
	if (nstepx < 3)
	{
		nstepx = 3;
		nstepz = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepz = ((int) (nstepz / 2)) * 2 + 1;

		if (nstepz < 3)
		{
			nstepz = 3;		
		}
	}

	if (nstepz < 3)
	{
		nstepz = 3;
		nstepx = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepx = ((int) (nstepx / 2)) * 2 + 1;

		if (nstepx < 3)
		{
			nstepx = 3;		
		}
	}

		//corresponding index to the 0 value of x and z	
	index_x = (nstepx - 1) / 2;
	index_z = (nstepz - 1) / 2;
	
	meshx = 2.0 * NSTD_LGM * std1 / (nstepx - 1);
	meshz = 2.0 * NSTD_LGM * std3 / (nstepz - 1);
		
		//Allocations of space vectors					
	r1_bar = dvector(0, nstepx-1);
	r3_bar = dvector(0, nstepz-1);
	r1_dim1 = dvector(0, nstepx-1);
	values = f3tensor(0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	values_p1 = f3tensor(0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	mu_r1 = dmatrix(0, nstepx-1, 0, nstepz-1);
	mu_r3 = dmatrix(0, nstepx-1, 0, nstepz-1);
	var_r1 = dmatrix(0, nstepx-1, 0, nstepz-1);
	var_r3 = dmatrix(0, nstepx-1, 0, nstepz-1);
	r = dmatrix(0, nstepx-1, 0, nstepz-1);
	r_init = dmatrix(0, nstepx-1, 0, nstepz-1);
	r1_dim2 = dmatrix(0, nstepx-1, 0, nstepz-1);
	r2 = dmatrix(0, nstepx-1, 0, nstepz-1);

	if (!r1_dim2 || !r2 || !values || !values_p1 || !mu_r1 || !mu_r3 || !var_r1 || !var_r3 || !r || !r_init ||
		!r1_bar || !r1_dim1 || !r3_bar)
	{
		err = "Memory allocation error (2) in doublelgm1fQuanto_adi";
		goto FREE_RETURN;
	}

	t2 = clock();

	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	smessage ("Phase 2 -convolution, stept: %d stepx: %d stepz: %d", nstept, nstepx, nstepz);

	t1 = clock();

		//Then discretise space in the orthogonal system r1 / r3	
						
	r1_bar[0] = -(nstepx - 1) / 2.0 * meshx;
	r3_bar[0] = -(nstepz - 1) / 2.0 * meshz;

	for (i=1; i<nstepx; i++)
	{
		r1_bar[i] = r1_bar[i-1] + meshx;
	}

	for (j=1; j<nstepz; j++)
	{
		r3_bar[j] = r3_bar[j-1] + meshz;
	}

		//Corresponding dom_r and for_r					
	for (i=0; i<nstepx; i++)
	{
		r1_dim1[i] = r1_bar[i] + fwd1[nstept-1];
		r2_temp = domforrho * r1_bar[i] * (contforsig[nstept-1]/contdomsig[nstept-1]);
		r_temp = domforrho * (1 + (contforsig[nstept-1]/contdomsig[nstept-1])) * r1_bar[i];

		for (j=0; j<nstepz; j++)
		{
			r1_dim2[i][j] = r1_dim1[i];
			r2[i][j] = fwd3[nstept-1] + r2_temp + (contforsig[nstept-1]/contdomsig[nstept-1]) * sqrt(1-domforrho*domforrho) * r3_bar[j];
			r_init[i][j] = r1_bar[i];
		}
	}

		//Final payoff valuation					
	if (!eval_evt[nstept-1])
	{
		err = "No event at last step in doublelgm1fQuanto_adi";
		goto FREE_RETURN;
	}
	
	err = payoff_func (	date[nstept-1],
						time[nstept-1],
						func_parm_tab[nstept-1],
						dom_yc,
						for_yc,
						domlam,
						forlam,
						domphi[nstept-1],
						forphi[nstept-1],
						0,
						nstepx-1,
						0,
						nstepz-1,
						r1_dim1,						
						r2,
						nprod,
						values_p1
						);

	if (err)
	{
		goto FREE_RETURN;
	}

		//Initialize the CNPDE_TEMP_2D		

	num_f_pde_init_2d_adi(	pde,
							nstepx,
							nstepz,
							nprod
							);
	
	if (!pde)
	{
		err = "Memory allocation error (3) in doublelgm1fQuanto_adi";
		goto FREE_RETURN;
	}
			
	lx = 0;
	ux = nstepx - 1;
	lz = 0;
	uz = nstepz - 1;

		//now do the backward pde					

	for (step=nstept-2; step>=0; step--)
	{
		dt = time[step+1] - time[step];
						
		r_temp = dom_ifr[step] + fwd1[step];

		for (i=lx; i<=ux; i++)
		{							
			mu_r1i = -domlam * r1_bar[i] * dt;
			mu_r3i = -(forlam - domlam) * domforrho * r1_bar[i] / sqrt(1-domforrho*domforrho);

			for (j=lz; j<=uz; j++)
			{
				mu_r1[i][j] = mu_r1i;
				mu_r3[i][j] = (
								mu_r3i - forlam * r3_bar[j]
								+ (r3_bar[j] + domforrho * r1_bar[i] / sqrt(1-domforrho*domforrho) )
									*((domabcsig[1][step]+2*domabcsig[2][step]*time[step])/contdomsig[step] - (forabcsig[1][step]+2*forabcsig[2][step]*time[step])/contforsig[step])
							  )* dt;
				

				var_r1[i][j] = var1[step];
				var_r3[i][j] = var3[step];
												
				r[i][j] = (r_init[i][j] + r_temp) * dt;
			}
		}

			//convolve
	
		num_f_pde_one_step_backward_2f_adi(	pde,
											nstepx,
											r1_bar,
											nstepz,
											r3_bar,
											0,
											nprod-1,
											values_p1,
											mu_r1,
											mu_r3,
											var_r1,
											var_r3,
											r,
											values,
											lx,
											ux,
											lz,
											uz);
		
			//Eval payoff 
		if (eval_evt[step])
		{

				//Corresponding r1 and r2					
			for (i=lx; i<=ux; i++)
			{
				r1_dim1[i] = r1_bar[i] + fwd1[step];
				r2_temp = domforrho * r1_bar[i] * (contforsig[step]/contdomsig[step]);

				for (j=0; j<nstepz; j++)
				{
					r1_dim2[i][j] = r1_dim1[i];
					r2[i][j] = fwd3[step] + r2_temp + (contforsig[step]/contdomsig[step]) * sqrt(1-domforrho*domforrho) * r3_bar[j];
				}
			}

			err = payoff_func (	date[step],
						time[step],
						func_parm_tab[step],
						dom_yc,
						for_yc,
						domlam,
						forlam,
						domphi[step],
						forphi[step],
						lx,
						ux,
						lz,
						uz,
						r1_dim1,						
						r2,
						nprod,
						values
						);

			if (err)
			{
				goto FREE_RETURN;
			}
		}

		values_temp = values_p1;
		values_p1 = values;
		values = values_temp;

		 //new indexes: we cut the PDE after NSTD_LGM number of standard deviations 
		 //but we need at least three points to do the pde						

		tem = (int) (NSTD_LGM * sqrt(domphi[step]) / meshx + 0.5);
		ux = min(nstepx-1, index_x + tem);
		lx = max(0, index_x - tem);
		
		if (ux - lx < 2)
		{
			tem += 1;
			ux = min(nstepx-1, index_x + tem);
			lx = max(0, index_x - tem);
		}

		tem = (int) (NSTD_LGM * std3 / meshx + 0.5);			
		uz = min(nstepz-1, index_z + tem);
		lz = max(0, index_z - tem);

		if (uz - lz < 2)
		{
			tem += 1;
			uz = min(nstepz-1, index_z + tem);
			lz = max(0, index_z - tem);
		}
	}

	 //copy the result					
	for (i=0; i<nprod; i++)
	{
		res[i] = values_p1[index_x][index_z][i];
	}

	t2 = clock();

	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);


FREE_RETURN:

		//allocation (1)		
	if (domphi) free_dvector(domphi, 0, nstept-1);
	if (forphi) free_dvector(forphi, 0, nstept-1);
	if (domforphi) free_dvector(domforphi, 0, nstept-1);
	if (var1) free_dvector(var1, 0, nstept-1);
	if (var3) free_dvector(var3, 0, nstept-1);
	if (fwd1) free_dvector(fwd1, 0, nstept-1);
	if (fwd3) free_dvector(fwd3, 0, nstept-1);
	if (domalpha) free_dvector(domalpha, 0, nb_sig-1);
	if (foralpha) free_dvector(foralpha, 0, nb_sig-1);
	if (contdomsig) free_dvector(contdomsig, 0, nstept-1);
	if (contforsig) free_dvector(contforsig, 0, nstept-1);
	if (allfxsig) free_dvector(allfxsig, 0, nstept-1);
	if (domabcsig) free_dmatrix(domabcsig, 0, 2, 0, nstept-1);
	if (forabcsig) free_dmatrix(forabcsig, 0, 2, 0, nstept-1);

		//allocation (2)		
	if (pde) num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

	if (r1_bar) free_dvector(r1_bar, 0, nstepx-1);
	if (r3_bar) free_dvector(r3_bar, 0, nstepz-1);
	if (r1_dim1) free_dvector(r1_dim1, 0, nstepx-1);
	if (r1_dim2) free_dmatrix(r1_dim2, 0, nstepx-1, 0, nstepz-1);
	if (r2) free_dmatrix(r2, 0, nstepx-1, 0, nstepz-1);
	if (values) free_f3tensor(values, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	if (values_p1) free_f3tensor(values_p1, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	if (mu_r1) free_dmatrix(mu_r1, 0, nstepx-1, 0, nstepz-1);
	if (mu_r3) free_dmatrix(mu_r3, 0, nstepx-1, 0, nstepz-1);
	if (var_r1) free_dmatrix(var_r1, 0, nstepx-1, 0, nstepz-1);
	if (var_r3) free_dmatrix(var_r3, 0, nstepx-1, 0, nstepz-1);
	if (r) free_dmatrix(r, 0, nstepx-1, 0, nstepz-1);
	if (r_init) free_dmatrix(r_init, 0, nstepx-1, 0, nstepz-1);
			
	return err;
}


/*	---------------------------- */
/*	Main function without Tau TS */
/*	---------------------------- */

Err	 doublelgm1fQuanto_adi_correl(	
						//Time data		
					int			nstept,
					double		*time,
					double		*date,

						//Discretisation	
					int			nsteps,
										
						//Model data		
					double		domlam,
					double		forlam,
					double		*sig_time,	//	domsig, forsig and fxsig must have 
					double		*domsig,	//	the same sig_time
					double		*forsig,
					double		*fxsig,
					int			nb_sig,
					double		*quantorho,
					double		*domforrho,

						//Product data
					void		**func_parm_tab, 
					int			*eval_evt,
					
						//Market data 
					double		*dom_ifr,
					double		*for_ifr,
					char		*dom_yc,
					char		*for_yc,
					
						//Payoff function 
					Err (*payoff_func)( //Event 
									double	evt_date,
									double	evt_time,
									void	*func_parm, 
					
									/* Market data	*/										
									void	*dom_yc,
									void	*for_yc,
									
									/* Model data	*/
									double	domlam,
									double	forlam,
									double	domphi,
									double	forphi,
									
									/* Grid data	*/
									int		l1,
									int		u1,
									int		l2,
									int		u2,
									double	*r1,
									double	**r2,
															
									/* Vector of results to be updated */
									int		nprod,
									double	***prod_val
									),
						//Result 
					int			nprod, 
					double		*res)
{

Err				err = NULL;
				
int				i, j, step, index_x, index_z;
int				nstepx, nstepz;
double			meshx, meshz, dt;
double			mu_r1i, mu_r3i;
double			r_temp;

double			std1, std3, r2_temp;

double			*r1_bar			= NULL,
				*r3_bar			= NULL,
				*r1_dim1		= NULL,
				**r1_dim2		= NULL,
				**r2			= NULL,
				***values		= NULL,
				***values_p1	= NULL,
				***values_temp	= NULL,
				**mu_r1			= NULL,
				**mu_r3			= NULL,
				**var_r1		= NULL,
				**var_r3		= NULL,
				**r				= NULL,
				**r_init		= NULL,
				*fwd1			= NULL,
				*fwd3			= NULL,
				*var1			= NULL,
				*var3			= NULL,
				*domphi			= NULL,
				*forphi			= NULL,
				*domforphi		= NULL;

int				lx, ux, lz, uz, tem;

double *domsigtilda=NULL;
double *forsigtilda=NULL;
double *correltilda=NULL;
double *contdomsig=NULL;
double *contforsig=NULL;

double *contdomforcorrel=NULL;
double *contquantocorrel=NULL;

double *allfxsig=NULL;
double **domabsig=NULL;
double **forabsig=NULL;
double **abcorrel=NULL;

double *domforcorreltilda=NULL;
double **ablindomforcorrel=NULL;
double **ablinquantocorrel=NULL;
double *quantocorreltilda=NULL;

clock_t			t1, t2;

double temp;

CNPDE_TEMP_2D_ADI	pdestr, *pde = &pdestr;

	t1 = clock();

		//Allocations	of time vectors						
	domphi = dvector(0, nstept-1);
	forphi = dvector(0, nstept-1);
	domforphi = dvector(0, nstept-1);
	var1 = dvector(0, nstept-1);
	var3 = dvector(0, nstept-1);
	fwd1 = dvector(0, nstept-1);
	fwd3 = dvector(0, nstept-1);
	domsigtilda = dvector(0, nb_sig-1);
	forsigtilda = dvector(0, nb_sig-1);
	domforcorreltilda = dvector(0, nb_sig-1);
	quantocorreltilda = dvector(0, nb_sig-1);
	contdomsig = dvector(0, nstept-1);
	contforsig = dvector(0, nstept-1);

	contdomforcorrel = dvector(0, nstept-1);
	contquantocorrel = dvector(0, nstept-1);

	allfxsig = dvector(0, nstept-1);
	domabsig = dmatrix(0, 1, 0, nstept-1);
	forabsig = dmatrix(0, 1, 0, nstept-1);
	ablindomforcorrel = dmatrix(0, 1, 0, nstept-1);
	ablinquantocorrel = dmatrix(0, 1, 0, nstept-1);

	if (!pde || !var1 || !var3 || !domphi || !forphi || !domforphi || !fwd1 || !fwd3 || !domsigtilda || !forsigtilda || !domforcorreltilda || !quantocorreltilda || !allfxsig || !contdomsig || !contforsig || !contdomforcorrel || !contquantocorrel || !domabsig || !forabsig || !ablindomforcorrel || !ablinquantocorrel)
	{
		err = "Memory allocation error (1) in doublelgm1fQuanto_adi";
		goto FREE_RETURN;
	}

	err = convertVol2(domlam, nstept, time, nb_sig, sig_time, domsig, domsigtilda, contdomsig, domabsig);
	if (err)
	{
		goto FREE_RETURN;
	}
	err = convertVol2(forlam, nstept, time, nb_sig, sig_time, forsig, forsigtilda, contforsig, forabsig);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = convertCorrelationWithTwoAffineVols(domlam, forlam, nstept, time, nb_sig, sig_time, 
						domforrho, domsig, forsig, 
						domsigtilda, forsigtilda, domforcorreltilda,
						ablindomforcorrel, contdomforcorrel);
	if (err)
	{
		goto FREE_RETURN;
	}

	err = convertCorrelationWithOneAffineVol(forlam, nstept, time, nb_sig, sig_time, 
						quantorho, forsig, forsigtilda, quantocorreltilda,
						ablinquantocorrel, contquantocorrel);
	if (err)
	{
		goto FREE_RETURN;
	}

	fillVol(nstept, time, nb_sig, sig_time, fxsig, allfxsig);

		//Calculate the expecations of r1 and r3 
	DoubleLGM1FQuantoExpectationsWithCorrelTS(nstept,
								time,
								domlam,
								forlam,
								sig_time,
								domsigtilda,
								forsigtilda,
								fxsig,
								nb_sig,
								quantocorreltilda,
								domforcorreltilda,
								fwd1,
								fwd3,
								var1,
								var3,
								domphi,
								forphi,
								domforphi);

		//Calculation of the number of steps in each direction: since local volatility			
		//is the same, the mesh has to be the same, but the number of steps has to be adjusted	

	std1 = sqrt(domphi[nstept-1]);
	std3 = sqrt(
					( (contdomsig[nstept-1]/contforsig[nstept-1])*(contdomsig[nstept-1]/contforsig[nstept-1])*forphi[nstept-1]
						+ contdomforcorrel[nstept-1]*contdomforcorrel[nstept-1]*domphi[nstept-1]
						- 2*contdomforcorrel[nstept-1]*(contdomsig[nstept-1]/contforsig[nstept-1])*domforphi[nstept-1]
					)/(1-contdomforcorrel[nstept-1]*contdomforcorrel[nstept-1])
				);

	nstepx = (int) (nsteps * sqrt(std1 / std3) + 0.5);
	nstepz = (int) (nsteps * sqrt(std3 / std1) + 0.5);
		
	//nstep has to be a odd nuber			
	nstepx = ((int) (nstepx / 2)) * 2 + 1;
	nstepz = ((int) (nstepz / 2)) * 2 + 1;

		//we want at least three points in each directions										
	if (nstepx < 3)
	{
		nstepx = 3;
		nstepz = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepz = ((int) (nstepz / 2)) * 2 + 1;

		if (nstepz < 3)
		{
			nstepz = 3;		
		}
	}

	if (nstepz < 3)
	{
		nstepz = 3;
		nstepx = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepx = ((int) (nstepx / 2)) * 2 + 1;

		if (nstepx < 3)
		{
			nstepx = 3;		
		}
	}

		//corresponding index to the 0 value of x and z	
	index_x = (nstepx - 1) / 2;
	index_z = (nstepz - 1) / 2;
	
	meshx = 2.0 * NSTD_LGM * std1 / (nstepx - 1);
	meshz = 2.0 * NSTD_LGM * std3 / (nstepz - 1);
		
		//Allocations of space vectors					
	r1_bar = dvector(0, nstepx-1);
	r3_bar = dvector(0, nstepz-1);
	r1_dim1 = dvector(0, nstepx-1);
	values = f3tensor(0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	values_p1 = f3tensor(0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	mu_r1 = dmatrix(0, nstepx-1, 0, nstepz-1);
	mu_r3 = dmatrix(0, nstepx-1, 0, nstepz-1);
	var_r1 = dmatrix(0, nstepx-1, 0, nstepz-1);
	var_r3 = dmatrix(0, nstepx-1, 0, nstepz-1);
	r = dmatrix(0, nstepx-1, 0, nstepz-1);
	r_init = dmatrix(0, nstepx-1, 0, nstepz-1);
	r1_dim2 = dmatrix(0, nstepx-1, 0, nstepz-1);
	r2 = dmatrix(0, nstepx-1, 0, nstepz-1);

	if (!r1_dim2 || !r2 || !values || !values_p1 || !mu_r1 || !mu_r3 || !var_r1 || !var_r3 || !r || !r_init ||
		!r1_bar || !r1_dim1 || !r3_bar)
	{
		err = "Memory allocation error (2) in doublelgm1fQuanto_adi";
		goto FREE_RETURN;
	}

	t2 = clock();

	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	smessage ("Phase 2 -convolution, stept: %d stepx: %d stepz: %d", nstept, nstepx, nstepz);

	t1 = clock();

		//Then discretise space in the orthogonal system r1 / r3	
						
	r1_bar[0] = -(nstepx - 1) / 2.0 * meshx;
	r3_bar[0] = -(nstepz - 1) / 2.0 * meshz;

	for (i=1; i<nstepx; i++)
	{
		r1_bar[i] = r1_bar[i-1] + meshx;
	}

	for (j=1; j<nstepz; j++)
	{
		r3_bar[j] = r3_bar[j-1] + meshz;
	}

		//Corresponding dom_r and for_r					
	for (i=0; i<nstepx; i++)
	{
		r1_dim1[i] = r1_bar[i] + fwd1[nstept-1];
		r2_temp = contdomforcorrel[nstept-1] * r1_bar[i] * (contforsig[nstept-1]/contdomsig[nstept-1]);
//		r_temp = contdomforcorrel[nstept-1] * (1 + (contforsig[nstept-1]/contdomsig[nstept-1])) * r1_bar[i];

		for (j=0; j<nstepz; j++)
		{
			r1_dim2[i][j] = r1_dim1[i];
			r2[i][j] = fwd3[nstept-1] + r2_temp + (contforsig[nstept-1]/contdomsig[nstept-1]) * sqrt(1-contdomforcorrel[nstept-1]*contdomforcorrel[nstept-1]) * r3_bar[j];
			r_init[i][j] = r1_bar[i];
		}
	}

		//Final payoff valuation					
	if (!eval_evt[nstept-1])
	{
		err = "No event at last step in doublelgm1fQuanto_adi_correl";
		goto FREE_RETURN;
	}
	
	err = payoff_func (	date[nstept-1],
						time[nstept-1],
						func_parm_tab[nstept-1],
						dom_yc,
						for_yc,
						domlam,
						forlam,
						domphi[nstept-1],
						forphi[nstept-1],
						0,
						nstepx-1,
						0,
						nstepz-1,
						r1_dim1,						
						r2,
						nprod,
						values_p1
						);

	if (err)
	{
		goto FREE_RETURN;
	}

		//Initialize the CNPDE_TEMP_2D		

	num_f_pde_init_2d_adi(	pde,
							nstepx,
							nstepz,
							nprod
							);
	
	if (!pde)
	{
		err = "Memory allocation error (3) in lgm2fpde";
		goto FREE_RETURN;
	}
			
	lx = 0;
	ux = nstepx - 1;
	lz = 0;
	uz = nstepz - 1;

	//now do the backward pde

	for (step=nstept-2; step>=0; step--)
	{
		dt = time[step+1] - time[step];
						
		r_temp = dom_ifr[step] + fwd1[step];

		for (i=lx; i<=ux; i++)
		{							
			mu_r1i = -domlam * r1_bar[i] * dt;
			mu_r3i = -(forlam - domlam) * contdomforcorrel[step] * r1_bar[i] / sqrt(1-contdomforcorrel[step]*contdomforcorrel[step]);

			for (j=lz; j<=uz; j++)
			{
				mu_r1[i][j] = mu_r1i;
				temp = (contdomforcorrel[step] * r3_bar[j] / sqrt(1-contdomforcorrel[step]*contdomforcorrel[step]) - r1_bar[i] )
									*ablindomforcorrel[1][step] / sqrt(1-contdomforcorrel[step]*contdomforcorrel[step]);

				mu_r3[i][j] = (
								mu_r3i - forlam * r3_bar[j]
								+ (r3_bar[j] + contdomforcorrel[step] * r1_bar[i] / sqrt(1-contdomforcorrel[step]*contdomforcorrel[step]) )
									*(domabsig[1][step]/contdomsig[step] - forabsig[1][step]/contforsig[step])
								+ (contdomforcorrel[step] * r3_bar[j] / sqrt(1-contdomforcorrel[step]*contdomforcorrel[step]) - r1_bar[i] )
									*ablindomforcorrel[1][step] / sqrt(1-contdomforcorrel[step]*contdomforcorrel[step])
							  )* dt;
				
				var_r1[i][j] = var1[step];
				var_r3[i][j] = var3[step];
												
				r[i][j] = (r_init[i][j] + r_temp) * dt;
			}
		}

			//convolve
	
		num_f_pde_one_step_backward_2f_adi(	pde,
											nstepx,
											r1_bar,
											nstepz,
											r3_bar,
											0,
											nprod-1,
											values_p1,
											mu_r1,
											mu_r3,
											var_r1,
											var_r3,
											r,
											values,
											lx,
											ux,
											lz,
											uz);
		
			//Eval payoff 
		if (eval_evt[step])
		{

				//Corresponding r1 and r2					
			for (i=lx; i<=ux; i++)
			{
				r1_dim1[i] = r1_bar[i] + fwd1[step];
				r2_temp = contdomforcorrel[step] * r1_bar[i] * (contforsig[step]/contdomsig[step]);

				for (j=0; j<nstepz; j++)
				{
					r1_dim2[i][j] = r1_dim1[i];
					r2[i][j] = fwd3[step] + r2_temp + (contforsig[step]/contdomsig[step]) * sqrt(1-contdomforcorrel[step]*contdomforcorrel[step]) * r3_bar[j];
				}
			}

			err = payoff_func (	date[step],
						time[step],
						func_parm_tab[step],
						dom_yc,
						for_yc,
						domlam,
						forlam,
						domphi[step],
						forphi[step],
						lx,
						ux,
						lz,
						uz,
						r1_dim1,						
						r2,
						nprod,
						values
						);

			if (err)
			{
				goto FREE_RETURN;
			}
		}

		values_temp = values_p1;
		values_p1 = values;
		values = values_temp;

		 //new indexes: we cut the PDE after NSTD_LGM number of standard deviations 
		 //but we need at least three points to do the pde						

		tem = (int) (NSTD_LGM * sqrt(domphi[step]) / meshx + 0.5);
		ux = min(nstepx-1, index_x + tem);
		lx = max(0, index_x - tem);
		
		if (ux - lx < 2)
		{
			tem += 1;
			ux = min(nstepx-1, index_x + tem);
			lx = max(0, index_x - tem);
		}

		tem = (int) (NSTD_LGM * std3 / meshx + 0.5);			
		uz = min(nstepz-1, index_z + tem);
		lz = max(0, index_z - tem);

		if (uz - lz < 2)
		{
			tem += 1;
			uz = min(nstepz-1, index_z + tem);
			lz = max(0, index_z - tem);
		}
	}

	 //copy the result					
	for (i=0; i<nprod; i++)
	{
		res[i] = values_p1[index_x][index_z][i];
	}

	t2 = clock();

	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);


FREE_RETURN:

		//allocation (1)		
	if (domphi) free_dvector(domphi, 0, nstept-1);
	if (forphi) free_dvector(forphi, 0, nstept-1);
	if (domforphi) free_dvector(domforphi, 0, nstept-1);
	if (var1) free_dvector(var1, 0, nstept-1);
	if (var3) free_dvector(var3, 0, nstept-1);
	if (fwd1) free_dvector(fwd1, 0, nstept-1);
	if (fwd3) free_dvector(fwd3, 0, nstept-1);
	if (domsigtilda) free_dvector(domsigtilda, 0, nb_sig-1);
	if (forsigtilda) free_dvector(forsigtilda, 0, nb_sig-1);
	if (contdomsig) free_dvector(contdomsig, 0, nstept-1);
	if (contforsig) free_dvector(contforsig, 0, nstept-1);
	if (allfxsig) free_dvector(allfxsig, 0, nstept-1);
	if (domabsig) free_dmatrix(domabsig, 0, 1, 0, nstept-1);
	if (forabsig) free_dmatrix(forabsig, 0, 1, 0, nstept-1);

		//allocation (2)		
	if (pde) num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

	if (r1_bar) free_dvector(r1_bar, 0, nstepx-1);
	if (r3_bar) free_dvector(r3_bar, 0, nstepz-1);
	if (r1_dim1) free_dvector(r1_dim1, 0, nstepx-1);
	if (r1_dim2) free_dmatrix(r1_dim2, 0, nstepx-1, 0, nstepz-1);
	if (r2) free_dmatrix(r2, 0, nstepx-1, 0, nstepz-1);
	if (values) free_f3tensor(values, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	if (values_p1) free_f3tensor(values_p1, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	if (mu_r1) free_dmatrix(mu_r1, 0, nstepx-1, 0, nstepz-1);
	if (mu_r3) free_dmatrix(mu_r3, 0, nstepx-1, 0, nstepz-1);
	if (var_r1) free_dmatrix(var_r1, 0, nstepx-1, 0, nstepz-1);
	if (var_r3) free_dmatrix(var_r3, 0, nstepx-1, 0, nstepz-1);
	if (r) free_dmatrix(r, 0, nstepx-1, 0, nstepz-1);
	if (r_init) free_dmatrix(r_init, 0, nstepx-1, 0, nstepz-1);
			
	return err;
}


/*	Adjust function */
void AdjustJumpingDrift(
/*	Size of the X array	*/
int							nx,
/*	X values	*/
double						*x,
/*	Size of the Y array	*/
int							ny,
/*	X values	*/
double						*y,
/*	Number of functions	*/
int							mstart,
int							mend,
/*	f values @ t+	*/
/*	i: x state, j: y state, k: function number */
double						***f_t_plus,
/*	Jumps Size*/
double						**Jumpx,
double						**Jumpy,
/*	Result: f values @ t	*/
double						***f_t,
/*	Lower and upper x and y */
int							lx,
int							ux,
int							ly,
int							uy)
{
	static int i, j, k;
	double yderivative, JumpYDividedByDeltaY;
//	double ***ftemp = NULL;

//	ftemp = f3tensor(0, nx-1, 0, ny-1, 0, mend);

	
	for (i=lx; i<=ux; i++)
	{
		for (k=mstart; k<=mend; k++) 
		{
			f_t[i][ly][k] = f_t_plus[i][ly][k];
		}

		for (j=ly+1; j<=uy; j++) 
		{
			for (k=mstart; k<=mend; k++)
			{
				yderivative = (f_t_plus[i][j][k] - f_t_plus[i][j-1][k]) / (y[j]-y[j-1]);
				JumpYDividedByDeltaY = 0.5 * Jumpy[i][j] / (y[j]-y[j-1]);
				f_t[i][j][k] = f_t_plus[i][j][k] + 0.5 * Jumpy[i][j] * yderivative;

//				f_t[i][j][k] = (f_t_plus[i][j][k] 
//								+ 0.5 * Jumpy[i][j] * yderivative
//								- f_t[i][j-1][k] * JumpYDividedByDeltaY
//								) / (1 - JumpYDividedByDeltaY);
			}
		}
	}

//	free_f3tensor(ftemp, 0, nx-1, 0, ny-1, 0, mend);
}




Err	 doublelgm1fQuanto_adi_correl2(	
						//Time data		
					int			nstept,
					double		*time,
					double		*date,

						//Discretisation
					int			nsteps,
										
						//Model data
					double		domlam,
					double		forlam,
					double		*sig_time,	//	domsig, forsig and fxsig must have 
					double		*domsig,	//	the same sig_time
					double		*forsig,
					double		*fxsig,
					int			nb_sig,
					double		*quantorho,
					double		*domforrho,

						//Product data
					void		**func_parm_tab, 
					int			*eval_evt,
					
						//Market data 
					double		*dom_ifr,
					double		*for_ifr,
					char		*dom_yc,
					char		*for_yc,
					
						//Payoff function 
					Err (*payoff_func)( //Event 
									double	evt_date,
									double	evt_time,
									void	*func_parm, 
					
									/* Market data	*/										
									void	*dom_yc,
									void	*for_yc,
									
									/* Model data	*/
									double	domlam,
									double	forlam,
									double	domphi,
									double	forphi,
									
									/* Grid data	*/
									int		l1,
									int		u1,
									int		l2,
									int		u2,
									double	*r1,
									double	**r2,
															
									/* Vector of results to be updated */
									int		nprod,
									double	***prod_val
									),
						//Result 
					int			nprod, 
					double		*res)
{

Err				err = NULL;

int				n = nb_sig-1;
int				i, j, step, index_x, index_z;
int				nstepx, nstepz;
double			meshx, meshz, dt;
double			mu_r1i, mu_r3i;
double			r_temp;

double			std1, std3, r2_temp;

double			*r1_bar			= NULL,
				*r3_bar			= NULL,
				*r1_dim1		= NULL,
				**r1_dim2		= NULL,
				**r2			= NULL,
				***values		= NULL,
				***values_p1	= NULL,
				***values_temp	= NULL,
				**mu_r1			= NULL,
				**mu_r3			= NULL,
				**var_r1		= NULL,
				**var_r3		= NULL,
				**zerovar1		= NULL,
				**zerovar3		= NULL,
				**r				= NULL,
				**zeror			= NULL,
				**r_init		= NULL,
				*fwd1			= NULL,
				*fwd3			= NULL,
				*var1			= NULL,
				*var3			= NULL,
				*domphi			= NULL,
				*forphi			= NULL,
				*domforphi		= NULL;

int				lx, ux, lz, uz, tem;

double *contdomsig=NULL;
double *contforsig=NULL;
double *contdomforcorrel=NULL;
double *contquantocorrel=NULL;

double *allfxsig=NULL;

double **Jump_r1 = NULL;
double **Jump_r3 = NULL;


clock_t			t1, t2;

CNPDE_TEMP_2D_ADI	pdestr, *pde = &pdestr;

	t1 = clock();

		//Allocations	of time vectors						
	domphi = dvector(0, nstept-1);
	forphi = dvector(0, nstept-1);
	domforphi = dvector(0, nstept-1);
	var1 = dvector(0, nstept-1);
	var3 = dvector(0, nstept-1);
	fwd1 = dvector(0, nstept-1);
	fwd3 = dvector(0, nstept-1);

	contdomsig = dvector(0, nstept-1);
	contforsig = dvector(0, nstept-1);
	contdomforcorrel = dvector(0, nstept-1);
	contquantocorrel = dvector(0, nstept-1);
	allfxsig = dvector(0, nstept-1);

	if (!pde || !var1 || !var3 || !domphi || !forphi || !domforphi || !fwd1 || !fwd3 || !allfxsig || !contdomsig || !contforsig || !contdomforcorrel || !contquantocorrel)
	{
		err = "Memory allocation error (1) in doublelgm1fQuanto_adi";
		goto FREE_RETURN;
	}

	fillVolAndCorrel(nstept, time, nb_sig, sig_time, domsig, forsig, fxsig, domforrho, quantorho, contdomsig, contforsig, allfxsig, contdomforcorrel, contquantocorrel);

		//Calculate the expecations of r1 and r3 
	DoubleLGM1FQuantoExpectationsWithCorrelationTS(nstept,
								time,
								domlam,
								forlam,
								sig_time,
								domsig,
								forsig,
								fxsig,
								nb_sig,
								quantorho,
								domforrho,
								fwd1,
								fwd3,
								var1,
								var3,
								domphi,
								forphi,
								domforphi);

		//Calculation of the number of steps in each direction: since local volatility			
		//is the same, the mesh has to be the same, but the number of steps has to be adjusted	

	std1 = sqrt(domphi[nstept-1]);
	std3 = sqrt(
					( (contdomsig[nstept-1]/contforsig[nstept-1])*(contdomsig[nstept-1]/contforsig[nstept-1])*forphi[nstept-1]
						+ contdomforcorrel[nstept-1]*contdomforcorrel[nstept-1]*domphi[nstept-1]
						- 2*contdomforcorrel[nstept-1]*(contdomsig[nstept-1]/contforsig[nstept-1])*domforphi[nstept-1]
					)/(1-contdomforcorrel[nstept-1]*contdomforcorrel[nstept-1])
				);

	nstepx = (int) (nsteps * sqrt(std1 / std3) + 0.5);
	nstepz = (int) (nsteps * sqrt(std3 / std1) + 0.5);
		
	//nstep has to be a odd nuber			
	nstepx = ((int) (nstepx / 2)) * 2 + 1;
	nstepz = ((int) (nstepz / 2)) * 2 + 1;

		//we want at least three points in each directions										
	if (nstepx < 3)
	{
		nstepx = 3;
		nstepz = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepz = ((int) (nstepz / 2)) * 2 + 1;

		if (nstepz < 3)
		{
			nstepz = 3;		
		}
	}

	if (nstepz < 3)
	{
		nstepz = 3;
		nstepx = (int) (nsteps * nsteps / 3.0 + 0.5);
		nstepx = ((int) (nstepx / 2)) * 2 + 1;

		if (nstepx < 3)
		{
			nstepx = 3;		
		}
	}

		//corresponding index to the 0 value of x and z	
	index_x = (nstepx - 1) / 2;
	index_z = (nstepz - 1) / 2;
	
	meshx = 2.0 * NSTD_LGM * std1 / (nstepx - 1);
	meshz = 2.0 * NSTD_LGM * std3 / (nstepz - 1);
		
		//Allocations of space vectors					
	r1_bar = dvector(0, nstepx-1);
	r3_bar = dvector(0, nstepz-1);
	r1_dim1 = dvector(0, nstepx-1);
	values = f3tensor(0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	values_p1 = f3tensor(0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	mu_r1 = dmatrix(0, nstepx-1, 0, nstepz-1);
	mu_r3 = dmatrix(0, nstepx-1, 0, nstepz-1);
	var_r1 = dmatrix(0, nstepx-1, 0, nstepz-1);
	var_r3 = dmatrix(0, nstepx-1, 0, nstepz-1);
	zerovar1 = dmatrix(0, nstepx-1, 0, nstepz-1);
	zerovar3 = dmatrix(0, nstepx-1, 0, nstepz-1);
	r = dmatrix(0, nstepx-1, 0, nstepz-1);
	zeror = dmatrix(0, nstepx-1, 0, nstepz-1);
	r_init = dmatrix(0, nstepx-1, 0, nstepz-1);
	r1_dim2 = dmatrix(0, nstepx-1, 0, nstepz-1);
	r2 = dmatrix(0, nstepx-1, 0, nstepz-1);

	Jump_r1 = dmatrix(0, nstepx-1, 0, nstepz-1);
	Jump_r3 = dmatrix(0, nstepx-1, 0, nstepz-1);

	if (!r1_dim2 || !r2 || !values || !values_p1 || !mu_r1 || !mu_r3 || !zerovar1 || !zerovar3 || !var_r1 || !var_r3 || !zeror || !r || !r_init ||
		!r1_bar || !r1_dim1 || !r3_bar)
	{
		err = "Memory allocation error (2) in doublelgm1fQuanto_adi";
		goto FREE_RETURN;
	}

	t2 = clock();

	smessage ("Phase 1 -preprocessing, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);
	smessage ("Phase 2 -convolution, stept: %d stepx: %d stepz: %d", nstept, nstepx, nstepz);

	t1 = clock();

		//Then discretise space in the orthogonal system r1 / r3	
						
	r1_bar[0] = -(nstepx - 1) / 2.0 * meshx;
	r3_bar[0] = -(nstepz - 1) / 2.0 * meshz;

	for (i=1; i<nstepx; i++)
	{
		r1_bar[i] = r1_bar[i-1] + meshx;
	}

	for (j=1; j<nstepz; j++)
	{
		r3_bar[j] = r3_bar[j-1] + meshz;
	}

		//Corresponding dom_r and for_r					
	for (i=0; i<nstepx; i++)
	{
		r1_dim1[i] = r1_bar[i] + fwd1[nstept-1];
		r2_temp = contdomforcorrel[nstept-1] * r1_bar[i] * (contforsig[nstept-1]/contdomsig[nstept-1]);

		for (j=0; j<nstepz; j++)
		{
			r1_dim2[i][j] = r1_dim1[i];
			r2[i][j] = fwd3[nstept-1] + r2_temp + (contforsig[nstept-1]/contdomsig[nstept-1]) * sqrt(1-contdomforcorrel[nstept-1]*contdomforcorrel[nstept-1]) * r3_bar[j];
			r_init[i][j] = r1_bar[i];
		}
	}

		//Final payoff valuation					
	if (!eval_evt[nstept-1])
	{
		err = "No event at last step in doublelgm1fQuanto_adi_correl";
		goto FREE_RETURN;
	}
	
	err = payoff_func (	date[nstept-1],
						time[nstept-1],
						func_parm_tab[nstept-1],
						dom_yc,
						for_yc,
						domlam,
						forlam,
						domphi[nstept-1],
						forphi[nstept-1],
						0,
						nstepx-1,
						0,
						nstepz-1,
						r1_dim1,						
						r2,
						nprod,
						values_p1
						);

	if (err)
	{
		goto FREE_RETURN;
	}

		//Initialize the CNPDE_TEMP_2D		

	num_f_pde_init_2d_adi(	pde,
							nstepx,
							nstepz,
							nprod
							);
	
	if (!pde)
	{
		err = "Memory allocation error (3) in lgm2fpde";
		goto FREE_RETURN;
	}
			
	lx = 0;
	ux = nstepx - 1;
	lz = 0;
	uz = nstepz - 1;

	//now do the backward pde

	for (step=nstept-2; step>=0; step--)
	{
		dt = time[step+1] - time[step];
						
		r_temp = dom_ifr[step] + fwd1[step];

		for (i=lx; i<=ux; i++)
		{							
			mu_r1i = -domlam * r1_bar[i] * dt;
			mu_r3i = -(forlam - domlam) * contdomforcorrel[step] * r1_bar[i] / sqrt(1-contdomforcorrel[step]*contdomforcorrel[step]);

			for (j=lz; j<=uz; j++)
			{
				mu_r1[i][j] = mu_r1i;

				mu_r3[i][j] = (	mu_r3i - forlam * r3_bar[j] )* dt;
				
				var_r1[i][j] = var1[step];
				var_r3[i][j] = var3[step];

				Jump_r1[i][j] = 0.0;
				Jump_r3[i][j] = (
									(contdomsig[step+1]/(contforsig[step+1]*sqrt(1-contdomforcorrel[step+1]*contdomforcorrel[step+1])))
									-(contdomsig[step]/(contforsig[step]*sqrt(1-contdomforcorrel[step]*contdomforcorrel[step])))
								)
								* (
									(
									0.5 * (contforsig[step+1]/contdomsig[step+1]) * r3_bar[j] * sqrt(1-contdomforcorrel[step+1]*contdomforcorrel[step+1])
									+ 0.5 * (contforsig[step+1]/contdomsig[step+1]) * contdomforcorrel[step+1] * r1_bar[i]
									)
									+
									(
									0.5 * (contforsig[step]/contdomsig[step]) * r3_bar[j] * sqrt(1-contdomforcorrel[step]*contdomforcorrel[step])
									+ 0.5 * (contforsig[step]/contdomsig[step]) * contdomforcorrel[step] * r1_bar[i]
									)
								  )
								-
								(
									contdomforcorrel[step+1] / sqrt(1-contdomforcorrel[step+1]*contdomforcorrel[step+1])
									-contdomforcorrel[step] / sqrt(1-contdomforcorrel[step]*contdomforcorrel[step])
								) * r1_bar[i];

				r[i][j] = (r_init[i][j] + r_temp) * dt;
			}
		}


		if((contforsig[step]!=contforsig[step+1])||(contdomsig[step]!=contdomsig[step+1])||(contdomforcorrel[step]!=contdomforcorrel[step+1]))
		{
			num_f_pde_one_step_backward_2f_adi(	pde,
												nstepx,
												r1_bar,
												nstepz,
												r3_bar,
												0,
												nprod-1,
												values_p1,
												Jump_r1,
												Jump_r3,
												zerovar1,
												zerovar3,
												zeror,
												values,
												lx,
												ux,
												lz,
												uz);
			values_temp = values_p1;
			values_p1 = values;
			values = values_temp;
		}

			//convolve
		num_f_pde_one_step_backward_2f_adi(	pde,
											nstepx,
											r1_bar,
											nstepz,
											r3_bar,
											0,
											nprod-1,
											values_p1,
											mu_r1,
											mu_r3,
											var_r1,
											var_r3,
											r,
											values,
											lx,
											ux,
											lz,
											uz);
		
			//Eval payoff 
		if (eval_evt[step])
		{

				//Corresponding r1 and r2					
			for (i=lx; i<=ux; i++)
			{
				r1_dim1[i] = r1_bar[i] + fwd1[step];
				r2_temp = contdomforcorrel[step] * r1_bar[i] * (contforsig[step]/contdomsig[step]);

				for (j=0; j<nstepz; j++)
				{
					r1_dim2[i][j] = r1_dim1[i];
					r2[i][j] = fwd3[step] + r2_temp + (contforsig[step]/contdomsig[step]) * sqrt(1-contdomforcorrel[step]*contdomforcorrel[step]) * r3_bar[j];
				}
			}

			err = payoff_func (	date[step],
						time[step],
						func_parm_tab[step],
						dom_yc,
						for_yc,
						domlam,
						forlam,
						domphi[step],
						forphi[step],
						lx,
						ux,
						lz,
						uz,
						r1_dim1,						
						r2,
						nprod,
						values
						);

			if (err)
			{
				goto FREE_RETURN;
			}
		}

		values_temp = values_p1;
		values_p1 = values;
		values = values_temp;

		 //new indexes: we cut the PDE after NSTD_LGM number of standard deviations 
		 //but we need at least three points to do the pde						

		tem = (int) (NSTD_LGM * sqrt(domphi[step]) / meshx + 0.5);
		ux = min(nstepx-1, index_x + tem);
		lx = max(0, index_x - tem);
		
		if (ux - lx < 2)
		{
			tem += 1;
			ux = min(nstepx-1, index_x + tem);
			lx = max(0, index_x - tem);
		}

		tem = (int) (NSTD_LGM * std3 / meshx + 0.5);			
		uz = min(nstepz-1, index_z + tem);
		lz = max(0, index_z - tem);

		if (uz - lz < 2)
		{
			tem += 1;
			uz = min(nstepz-1, index_z + tem);
			lz = max(0, index_z - tem);
		}
	}

	 //copy the result					
	for (i=0; i<nprod; i++)
	{
		res[i] = values_p1[index_x][index_z][i];
	}

	t2 = clock();

	smessage ("Phase 2 -convolution, time in sec: %.2f", (double) (t2 - t1) / CLOCKS_PER_SEC);


FREE_RETURN:

		//allocation (1)		
	if (domphi) free_dvector(domphi, 0, nstept-1);
	if (forphi) free_dvector(forphi, 0, nstept-1);
	if (domforphi) free_dvector(domforphi, 0, nstept-1);
	if (var1) free_dvector(var1, 0, nstept-1);
	if (var3) free_dvector(var3, 0, nstept-1);
	if (fwd1) free_dvector(fwd1, 0, nstept-1);
	if (fwd3) free_dvector(fwd3, 0, nstept-1);
	if (contdomsig) free_dvector(contdomsig, 0, nstept-1);
	if (contforsig) free_dvector(contforsig, 0, nstept-1);
	if (contdomforcorrel) free_dvector(contdomforcorrel, 0, nstept-1);
	if (contquantocorrel) free_dvector(contquantocorrel, 0, nstept-1);
	if (allfxsig) free_dvector(allfxsig, 0, nstept-1);

		//allocation (2)		
	if (pde) num_f_pde_free_2d_adi(pde, nstepx, nstepz, nprod);

	if (r1_bar) free_dvector(r1_bar, 0, nstepx-1);
	if (r3_bar) free_dvector(r3_bar, 0, nstepz-1);
	if (r1_dim1) free_dvector(r1_dim1, 0, nstepx-1);
	if (r1_dim2) free_dmatrix(r1_dim2, 0, nstepx-1, 0, nstepz-1);
	if (r2) free_dmatrix(r2, 0, nstepx-1, 0, nstepz-1);
	if (values) free_f3tensor(values, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	if (values_p1) free_f3tensor(values_p1, 0, nstepx-1, 0, nstepz-1, 0, nprod-1);
	if (mu_r1) free_dmatrix(mu_r1, 0, nstepx-1, 0, nstepz-1);
	if (mu_r3) free_dmatrix(mu_r3, 0, nstepx-1, 0, nstepz-1);
	if (Jump_r1) free_dmatrix(Jump_r1, 0, nstepx-1, 0, nstepz-1);
	if (Jump_r3) free_dmatrix(Jump_r3, 0, nstepx-1, 0, nstepz-1);
	if (var_r1) free_dmatrix(var_r1, 0, nstepx-1, 0, nstepz-1);
	if (var_r3) free_dmatrix(var_r3, 0, nstepx-1, 0, nstepz-1);
	if (zerovar1) free_dmatrix(zerovar1, 0, nstepx-1, 0, nstepz-1);
	if (zerovar3) free_dmatrix(zerovar3, 0, nstepx-1, 0, nstepz-1);
	if (r) free_dmatrix(r, 0, nstepx-1, 0, nstepz-1);
	if (zeror) free_dmatrix(zeror, 0, nstepx-1, 0, nstepz-1);
	if (r_init) free_dmatrix(r_init, 0, nstepx-1, 0, nstepz-1);
			
	return err;
}
