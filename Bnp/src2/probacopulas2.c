
// S. Galluccio: 21 November 2000
// Created Source File

/**************************************************************************************
*                                                                                     *
*                    Converts a Rank Correlation matrix into                          *
*                     a linear correlation matrix                                     *
*                                                                                     *
**************************************************************************************/

#include "math.h"
#include "UTALLHDR.H>
#include "num_h_allhdr.h"
#include "num_h_proba.h"
#include "utconst.h"
	
#define NRANSI
#define ITMAX 100
#define EPS_CPL 3.0e-8

double RankToLinCorr  (	 
						const double SpRho,
						const int degree,
						const int acc_degree)
{
	int iter;
	double x1=0.0,x2=0.99,tol=1.e-4;
	double a=x1,b=x2,c=x2,d,e,min1,min2;
	double fa,fb,fc,p,q,r,s,tol1,xm;

	switch (degree) {

	case 0: 
			return 2.*sin(SRT_PI*SpRho/6.);  // if degree=0 returns the gaussian case (NB, it is a convention, since 
										 // the gaussian is recovered for degree -> 00 )
		break;
	
	default:

		 fa= CopulaRhoSpearman(	a,
							SpRho,
							degree,
							acc_degree);

		 fb= CopulaRhoSpearman(	b,
							SpRho,
							degree,
							acc_degree);


		if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) goto end;
	//		nrerror("Root must be bracketed in zbrent");
		fc=fb;
		for (iter=1;iter<=ITMAX;iter++) {
			if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
				c=a;
				fc=fa;
				e=d=b-a;
			}
			if (fabs(fc) < fabs(fb)) {
				a=b;
				b=c;
				c=a;
				fa=fb;
				fb=fc;
				fc=fa;
			}
			tol1=2.0*EPS_CPL*fabs(b)+0.5*tol;
			xm=0.5*(c-b);
			if (fabs(xm) <= tol1 || fb == 0.0) return b;
			if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
				s=fb/fa;
				if (a == c) {
					p=2.0*xm*s;
					q=1.0-s;
				} else {
					q=fa/fc;
					r=fb/fc;
					p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
					q=(q-1.0)*(r-1.0)*(s-1.0);
				}
				if (p > 0.0) q = -q;
				p=fabs(p);
				min1=3.0*xm*q-fabs(tol1*q);
				min2=fabs(e*q);
				if (2.0*p < (min1 < min2 ? min1 : min2)) {
					e=d;
					d=p/q;
				} else {
					d=xm;
					e=d;
				}
			} else {
				d=xm;
				e=d;
			}
			a=b;
			fa=fb;
			if (fabs(d) > tol1)
				b += d;
			else
				b += SIGN(tol1,xm);
			fb=CopulaRhoSpearman(
								b,
								SpRho,
								degree,
								acc_degree);
			}

//	nrerror("Maximum number of iterations exceeded in zbrent");
		end:
		return 0.0;
	break;

	} // end of switch
}
#undef ITMAX
#undef EPS_CPL
#undef NRANSI

/**************************************************************************************
*                                                                                     *
*                    The next four functions are used to evaluate                    *
*                     the integral on the unit square [0,1]^2
*                     that is used to compute the rank correlation                    *
*                                                                                     *
**************************************************************************************/
double CopulaRhoSpearman(
						const double r,
						const double SpRho,
						const int degree,
						const int acc_degree
					)
{
	double Val,Int;

	Int	=	IntegrateCopula(	r,
								degree,
								acc_degree
							);
	
	Val=12.*Int-3.-SpRho;
	return Val;
}

//////////////////////////////////////////////////////////////////////////////////////

#define N_SMALL 10
double IntegrateCopula ( 
						const double r,
						const int n,
						const int acc_degree
					  )
{

	int i,j;
	double *x,*y,*wx,*wy,integrand,Value=0.0,*xin,*yin;

	x=dvector(1,N_SMALL);
	y=dvector(1,N_SMALL);
	wx=dvector(1,N_SMALL);
	wy=dvector(1,N_SMALL);
	xin=dvector(1,N_SMALL);
	yin=dvector(1,N_SMALL);

	GaussLeg(0.0,1.0,x,wx,N_SMALL);
	GaussLeg(0.0,1.0,y,wy,N_SMALL);

	for (i=1;i<=N_SMALL;i++) 	xin[i]=yin[i] = Student_Dis_Inv(
																1.e-7,
																x[i],
																n
																);
	

	for(i=1;i<=N_SMALL;i++) {
		for (j=1;j<=N_SMALL;j++) {
 
			integrand=EvalStudCpl(	
											-100.*acc_degree,
											xin[i],
											-100.*acc_degree,
											yin[j],
											r,
											n,
											acc_degree);

			Value+=wx[i]*wy[j]*integrand;
		}
	}

	///// free the memory 

	free_dvector(x,1,N_SMALL);
	free_dvector(y,1,N_SMALL);
	free_dvector(wx,1,N_SMALL);
	free_dvector(wy,1,N_SMALL);
	free_dvector(xin,1,N_SMALL);
	free_dvector(yin,1,N_SMALL);

	return Value;
	
}
#undef N_SMALL

//////////////////////////////////////////////////////////////////////////////////////

double EvalStudCpl(
					const double		xmin,
					const double		xmax,
					const double		ymin,
					const double		ymax,
					const double		r,
					const int	degree,
					const int	acc_degree)

{
	int i,j,n_gl_pts=100*acc_degree;
	double *x,*y,*wx,*wy,Cpl=0.0,integrand;

	x=dvector(1,n_gl_pts);
	y=dvector(1,n_gl_pts);
	wx=dvector(1,n_gl_pts);
	wy=dvector(1,n_gl_pts);

	///// determines the Gauss-Lengendre points and weights in the 2 dimensions
	GaussLeg(xmin,xmax,x,wx,n_gl_pts);
	GaussLeg(ymin,ymax,y,wy,n_gl_pts);

	///// Comoutes the integral

	for(i=1;i<=n_gl_pts;i++) {
		for (j=1;j<=n_gl_pts;j++) {

			integrand=EvalStudCplIntegrand(x[i],y[j],r,degree);
			Cpl+=wx[i]*wy[j]*integrand;

		}
	}
	///// free the memory 

	free_dvector(x,1,n_gl_pts);
	free_dvector(y,1,n_gl_pts);
	free_dvector(wx,1,n_gl_pts);
	free_dvector(wy,1,n_gl_pts);

	return Cpl;
}

/////////////////////////////////////////////////////////////////////////////
////////////////////  Evaluates the integrand of the 2D Student Copula //////

double EvalStudCplIntegrand(	
						const double x,
						const double y,
						const double r,
						const int n
						)
{
	int i; 
	double A,value=1.;

	A=1./(2.*SRT_PI*sqrt(1.-r*r));
	
	for (i=1;i<=n+2;i++) value*=1.+(x*x-2.*r*x*y+y*y)/n/(1.-r*r);
	value=A/sqrt(value);

	return value;
}




//------------------------------------------------------------------------------------------------------------------------------------
//
// Calculates the linear correlation using either the rank correlation or the Kendall's tau
//
//------------------------------------------------------------------------------------------------------------------------------------
// Given data arrays data1[1..n] and data2[1..n], this program returns Kendall's ½ as tau,
// its number of standard deviations from zero as z, and its two-sided signicance level as prob.
// Small values of prob indicate a signicant correlation (tau positive) or anticorrelation (tau
// negative).  
//@@@@@@ Modified to just return the tau value
double getKendallTau(double data1[], double data2[], unsigned long n )
{
	double erfcc(double x);
	unsigned long n2=0,n1=0,k,j;
	long is=0;
	double aa,a2,a1;
	for (j=1;j<n;j++) {											// Loop over first member of pair,
		for (k=(j+1);k<=n;k++) {								// and second member.
			a1=data1[j]-data1[k];
			a2=data2[j]-data2[k];
			aa=a1*a2;
			if (aa) {											// Neither array has a tie.
			++n1;
			++n2;
			aa > 0.0 ? ++is : --is;
		} else {												// One or both arrays have ties.
		if (a1) ++n1;											// An extra x" event.
		if (a2) ++n2;											// An extra y" event.
		}
		}
	}
	return is/(sqrt((double) n1)*sqrt((double) n2));				// Equation (14.6.8).
}



void initLinearCorrelation( double* dv1Data1, double* dv1Data2, unsigned int uiNum, double* dPar, unsigned long uiRankToLinAcc )
{
	if ( uiRankToLinAcc > 0 )
		*dPar = spear_rho( dv1Data1, dv1Data2, uiNum );
	else
		*dPar = sin( 0.5 * SRT_PI * getKendallTau( dv1Data1, dv1Data2, uiNum )  );
}


double getLinearCorrelation( unsigned int uiCopulaDegree, double dPar, unsigned long uiRankToLinAcc )
{
	if ( uiRankToLinAcc > 0 )
		return RankToLinCorr( dPar, uiCopulaDegree, uiRankToLinAcc );
	else
		return dPar;
}


/**************************************************************************************************************************************/
// 
// Stud_Deg_Fit
//
// Calculates the MLE value of the degree of a bivariate Student copula and the corresponding linear correlation coefficient using an
// input data series with fitted generalized Student parameters
//
// Outputs:
//			m						the MLE value of the Student degree for the bivariate copula
//			linear_corr				the MLE value of the linear correlation parameter for the bivariate copula
//
// Inputs: 
//			data1					input historical data of changes in forward rate 1.  First element is data[1]
//			m1						fitted degree of generalized Student distribution for data1
//			sigma1					fitted variance parameter of generalized Student distribution for data1
//			data2					input historical data of changes in forward rate 2.  First element is data[1]
//			m2						fitted degree of generalized Student distribution for data1
//			sigma2					fitted variance parameter of generalized Student distribution for data1
//			N						number of data points
//			Inv_Stud_acc			the accuracy in the inverse Student distribution calculation
//			Max_Stud_Degree			maximum value of the Student degree tested for the MLE
//			RankToLin_acc			accuracy parameter for calculating the linear correlation from the rank correlation.  Uses the 
//										Kendall tau approximation for the linear correlation if the value is set to 0.
//
// Return Value:
//			none
//								
//
/**************************************************************************************************************************************/
void Stud_Deg_Fit(			 
					unsigned long *m,
					double *linear_corr,				
					double *data1,
					const double m1,
					const double sigma1,
					double *data2,						
					const double m2,
					const double sigma2,
					const unsigned short N,
					double Inv_Stud_acc,
					unsigned long Max_Stud_Degree,		
					unsigned long RankToLin_acc )
{						
// Local variables.  Because some routines start from 0 and some from 1, we need a pointer to the first data element.
	double max_llhd, temp_linear_corr, temp_max_llhd, dLinearCorrelationParameter;
	double *PtrToData1 = &data1[1];
	double *PtrToData2 = &data2[1];

// Calculate the cumulative probability densities for the input data.  It currently only works for integer values, so round.
	double *CumDensity1 = dvector(1,N);
	double *CumDensity2 = dvector(1,N);
	int m1int = (int) floor( m1+0.5 );
	int m2int = (int) floor( m2+0.5 );
	unsigned short i;
	for ( i=1; i<=N; i++ )
	{
		CumDensity1[ i ] = Student_Dis( m1int, data1[i]/sigma1 );
		CumDensity2[ i ] = Student_Dis( m2int, data2[i]/sigma2 );
	}

// Compute the rank correlation.  Then loop over all possible values of the student degree, calculating the implied linear correlation
// and then the modified log-likelihood value.  Keep the maximum log-likelihood values.
	initLinearCorrelation( data1, data2, N, &dLinearCorrelationParameter, RankToLin_acc );
	*linear_corr = getLinearCorrelation( 0, dLinearCorrelationParameter, RankToLin_acc ); 
	*m = 0;
	max_llhd = Stud_cop_llhd( 0, *linear_corr, CumDensity1, CumDensity2, N, Inv_Stud_acc );
	for ( i=1; i<=Max_Stud_Degree; i++ )
	{
		temp_linear_corr = getLinearCorrelation( i, dLinearCorrelationParameter, RankToLin_acc );
		temp_max_llhd = Stud_cop_llhd( i, temp_linear_corr, CumDensity1, CumDensity2, N, Inv_Stud_acc );
		if ( temp_max_llhd > max_llhd )
		{
			*m = i;
			max_llhd = temp_max_llhd;
			*linear_corr = temp_linear_corr;
		}
	}


// Free the memory and call it a day
	free_dvector(CumDensity1,1,N);
	free_dvector(CumDensity2,1,N);

}



/**************************************************************************************************************************************/
// 
// Stud_Bivar_Copula
//
// Calculates all the likelihood values for a range of degrees of a bivariate Student copula with the corresponding linear correlation 
// coefficients determined by a fit to the rank correlation of an
// input data series with fitted Student parameters.
//
// Needs to be altered to allow for a fit with a generalized Student distribution
//
// Outputs:
//			m						the degree of the Student bivariate copula.  First element stored in data[1].  Assumes that memory
//										is declared elsewhere.
//			llhd					the likelihood value for the corresponding degree.  Same memory and storage conventions as m.
//			linear_corr				the linear correlation parameter for the corresponding degree.  Same memory and storage conventions 
//										as m.
//
// Inputs: 
//			data1					input historical data of changes in forward rate 1.  First element is data[1]
//			m1						fitted degree of generalized Student distribution for data1
//			sigma1					fitted variance parameter of generalized Student distribution for data1
//			data2					input historical data of changes in forward rate 2.  First element is data[1]
//			m2						fitted degree of generalized Student distribution for data1
//			sigma2					fitted variance parameter of generalized Student distribution for data1
//			N						number of data points
//			Inv_Stud_acc			the accuracy in the inverse Student distribution calculation
//			Max_Stud_Degree			maximum value of the Student degree tested for the MLE
//			RankToLin_acc			accuracy parameter for calculating the linear correlation from the rank correlation
//
// Return Value:
//			none
//								
//
/**************************************************************************************************************************************/
void Stud_Bivar_Copula(			 
					unsigned long *m,
					double *llhd,
					double *linear_corr,				
					double *data1,
					const double m1,
					const double sigma1,
					double *data2,						
					const double m2,
					const double sigma2,
					const unsigned short N,
					double Inv_Stud_acc,
					unsigned long Max_Stud_Degree,		
					unsigned long RankToLin_acc )
{						
// Local variables.  Because some routines start from 0 and some from 1, we need a pointer to the first data element.
	double rank_corr;
	double *PtrToData1 = &data1[1];
	double *PtrToData2 = &data2[1];

// Calculate the cumulative probability densities for the input data.  It currently only works for integer values, so round.
	double *CumDensity1 = dvector(1,N);
	double *CumDensity2 = dvector(1,N);
	int m1int = (int) floor( m1+0.5 );
	int m2int = (int) floor( m2+0.5 );
	unsigned short i;
	for ( i=1; i<=N; i++ )
	{
		CumDensity1[ i ] = Student_Dis( m1int, data1[i]/sigma1 );
		CumDensity2[ i ] = Student_Dis( m2int, data2[i]/sigma2 );
	}

// Compute the rank correlation.  Then loop over all possible values of the student degree, calculating the implied linear correlation
// and then the modified log-likelihood value.  Store the relevant values. 
	rank_corr = spear_rho( data1, data2, N );
	linear_corr[1] = RankToLinCorr( rank_corr, 0, RankToLin_acc );
	m[1] = 0;
	llhd[1] = Stud_cop_llhd( 0, linear_corr[1], CumDensity1, CumDensity2, N, Inv_Stud_acc );
	for ( i=1; i<=Max_Stud_Degree; i++ )
	{
		m[i+1] = i;
		linear_corr[i+1] = RankToLinCorr( rank_corr, i, RankToLin_acc );
		llhd[i+1] = Stud_cop_llhd( i, linear_corr[i+1], CumDensity1, CumDensity2, N, Inv_Stud_acc );
	}


// Free the memory and call it a day
	free_dvector(CumDensity1,1,N);
	free_dvector(CumDensity2,1,N);

}


//
// Calculates the modified log-likelihood function for a bivariate Student copula
//
double Stud_cop_llhd( unsigned long m, double rho, double *CumDensity1, double *CumDensity2, unsigned short N, double acc )
{
	double sum = 0.0;
	unsigned short i;

// Loop over the data points, calculating the required value
	for ( i=1; i<N; i++ )
		sum += log_Student_copula_dist( CumDensity1[i], CumDensity2[i], m, rho, acc );
	return sum;

}


/**************************************************************************************************************************************/
// 
// log_Student_copula_dist
//
// Returns the log of Student's t bivariate copula
//
// Inputs: u1		uniform variate at which to evaluate the copula 
//		   u2		uniform variate at which to evaluate the copula 
//		   m		degree of the distribution (0 is Gaussian, otherwise can be any positive integer) 
//         rho		correlation parameter,
//		   acc		accuracy parameter for the inverse Student distribution	
//
/**************************************************************************************************************************************/
double log_Student_copula_dist( double u1, double u2, unsigned long m, double rho, double acc )
{
	double x1 = Inverse_Student_Distribution( u1, m, acc );
	double x2 = Inverse_Student_Distribution( u2, m, acc );
	double det = 1-rho*rho;
	if ( m==0 )
		return -rho*0.5*( rho*(x1*x1+x2*x2)  - 2.0*x1*x2  )/det - 0.5*log(det);
	else
		return gammln( 0.5*(m+2.0) ) + gammln( 0.5*m ) - 2.0* gammln( 0.5*(m+1.0) ) - 0.5*log(det) - 
			0.5*(m+2.0)*log( 1+(x1*x1-2.0*rho*x1*x2+x2*x2)/det/m ) + 0.5*(m+1.0)*log( (1+x1*x1/m)*(1+x2*x2/m) );
}


/**************************************************************************************************************************************/
// 
// log_Student_dist
//
// Returns the log of Student's t probability density
//
// Inputs: x		value at which the probability density is to be evaluated (can be any real number)
//		   m		degree of the distribution (0 is Gaussian, otherwise ca be any positive integer) 
//
/**************************************************************************************************************************************/
double log_Student_dist( double x, unsigned long m )
{
	if ( m==0 )
		return -0.5*( x*x + log(2.0*SRT_PI) );
	else
		return gammln( 0.5*(m+1.0) ) - gammln( 0.5*m ) - 0.5*log( SRT_PI*m ) - 0.5*(m+1.0)*log( 1+x*x/m );
}


double Student_Density( double x, double m, double sigma )
{
	double sigma2 = sigma*sigma;
	return exp( gammln( 0.5*(m+1.0) ) - gammln( 0.5*m ) - 0.5*log( SRT_PI*m*sigma2 ) - 0.5*(m+1.0)*log( 1+x*x/m/sigma2 )  );
}


double Std_Student_Density( double x, double m )
{
	return Student_Density( x, m, 1.0 );
}


/**************************************************************************************************************************************/
// 
// Hill_estimator
//
// Calculates the Hill estimator for a set of data
//
// Inputs: 
//			data				the ordered or unordered input data.  The first element is data[1].  Is unchanged by the function
//			k					k+1 is the number of upper order statistics used
//          N					set to 0 if the data is ordered, otherwise is set equal to the number of data elements 
// Outputs:
//
// Return Value:
//			the Hill estimator
//								
//
/**************************************************************************************************************************************/
double Hill_estimator( double *data, unsigned long k, unsigned long N )
{
	double sum = 0.0;
	unsigned long i, *index;

// Check to see if the data is unordered and then create an ordered index into it.
	if ( N )
	{
		index = lvector( 1, N );
		index_data( N, data, index );
		for ( i=1; i<=k; i++ )
			sum += log( data[index[i]] );
	}
	else
	{
		for ( i=1; i<=k; i++ )
		sum += log( data[i] );
	}

// Return the value.
	sum /= k;
	return sum - log( data[k+1] );

}


/**************************************************************************************************************************************/
// 
// Stud_Log_Likelihd
//
// Returns the derivative of the log-likelihood function for the Student distribution and and its derivative
//
// Inputs: 
//			sigma				the standard deviation parameter
//			data				An array of input data.  data[0] is the Student degree.  data[1] is the number of observations, 
//									while data[2] onwards gives the actual data
// Outputs:
//			f					the value of the derivative of the log-likelihood function
//			df					the value of the 2nd derivative of the log-likelihood function
//
// Return Value:
//			Null if no error
//								
//
/**************************************************************************************************************************************/
Err Stud_Log_Likelihd( double sigma, double *f, double *df, double *data )
{
	unsigned long i;
	double m = data[0];
	double msigma2 = m*sigma*sigma;
	double N = data[1]+1;

	double x2;
	double denom;

	*f = 0.0;
	*df = 0.0;
	for ( i = 2; i < N; i++ )
	{
		x2 = data[i]*data[i];
		denom = msigma2 + x2;
		*f += x2/denom;
		*df = x2 /denom / denom;
	}

// Final alterations
	*f *= (m+1.0);
	*f -= floor(data[1]);
	*df = -2.0*m*(m+1.0);

	return NULL;
}


/**************************************************************************************************************************************/
// 
// Stud_Deg_Diff
//
// Returns the difference between the model and historical StDev for a given Student degree.  Also gives the derivative
//
// Inputs: 
//			m					the Student degree parameter parameter
//			data				An array of input data.  data[0] is the historical StDev.  data[1] is the number of 
//									observations, while data[2] onwards gives the actual data
// Outputs:
//			f					the value of the function
//			df					the value of the derivative of the function
//
// Return Value:
//			Null if no error
//								
//
/**************************************************************************************************************************************/
Err Stud_Deg_Diff( double m, double *f, double *df, double *data )
{

// Local variables and switch around the array data, such that it can be used to call the log-likelihood function
	Err err = NULL;
	double HistStDev = data[0];
	double HistSigma;
	double sigma, StDev, StDevP;
	if ( m <= 2.0 )
		m = 2.0001;
	data[0] = m;
	HistSigma = HistStDev * sqrt( (m-2.0) / m );

// Calculate the value of sigma with the current Student degree
	if ( err = rtsafe_with_par(&Stud_Log_Likelihd, HistSigma/3.0, HistSigma*3.0, HistSigma/100000.0, 100, &sigma, data) )
		return err;
	StDev = sigma*sqrt(m/(m-2.0));

// Calculate the derivative by a simple forward differencing
	data[0] = m*1.01;
	if ( err = rtsafe_with_par(&Stud_Log_Likelihd, HistSigma/3.0, HistSigma*3.0, HistSigma/100000.0, 100, &sigma, data) )
		return err;
	StDevP = sigma*sqrt(m*1.01/(m*1.01-2.0));
	
// Restore the data array to its original state and return the values
	data[0] = HistStDev;
	*f = StDev - HistStDev;
	*df = (StDevP - StDev) / m * 100;

	return err;
}


//
// generate two vectors of N data with Student marginals, joined by a Student copula
// The output vectors are assumed to already have memory allocated and run from 1 to N in good NR fashion
// 
void Student_Simluation( double *data1, 
					    double *data2, 
						const unsigned int N, 
						const unsigned int m, 
						const double rho, 
						const unsigned int m1, 
						const double sigma1, 
						const unsigned int m2, 
						const double sigma2, 
						const double acc,
						long seed )
{
	unsigned int i;
	double u1, u2;

// Get the uniform variates from the copula and convert them to t-variates
	for ( i=1; i<=N; i++ )
	{
		BiStudCpl_Variates( &u1, &u2, m, rho, &seed );
		data1[i] = sigma1*Inverse_Student_Distribution( u1, m1, acc );
		data2[i] = sigma2*Inverse_Student_Distribution( u2, m2, acc );
	}

}




void BiStudCpl_Variates( double *u1, double *u2, const unsigned int m, const double rho, long *seed )
{
// Local Variables
	unsigned int i;
	double scale, chi2 = 0.0;

// Generate two correlated Gaussian random variables 
	double W1 = gauss_sample (seed);
	double W2 = rho*W1 + sqrt(1-rho*rho)*gauss_sample(seed);

// Generate the chi^2 random variable by summing m Gaussians

	for ( i=0; i<m; i++ )
		chi2 += pow(gauss_sample(seed),2);
	scale = sqrt(1.*m) / sqrt( chi2 );

// Produce the output variables by using the Student cumulative distribution
	*u1 = Student_Dis( m, W1*scale );
	*u2 = Student_Dis( m, W2*scale );

}


/* ========================================================================
   FUNC: GetStudentCplDev
   DESC: Generates  a matrix (0,p-1)x(0,d-1) that contains a sample of correlated 
		 random variables according to a Student Copula of degree n, and compatible
		 with a set of marginal cumulative distributions

		degree:  degree of the Student Copula (must be <=10)
		mean_v: vector of 0
		corr_mtx: linear correlatio  matrix (extracted from the rank corr. matrix)
		p: number of simulations
		d: total n. of degrees of freedom, i.e. degree + n. of assets
		xa: a matrix (0,n_pts-1)x(0,n_assets-1) containing the abscissas of the marginal
			cumulative distributions
		ya: a matrix (0,n_pts-1)x(0,n_assets-1) containing the ordinates of the marginal
			cumulative distributions
		n_pts: number of abscissas
		n_conv: number of convolutions
		
   MODIFIES:
   DECLARATION:
   ======================================================================== */
/*
	void GetStudentCplDev (
						const int degree,
	  					const double *mean_v,
						double **corr_mtx,   
						const long p,
						const int d, 
						double **xa,
						double **ya,
						const long n_pts,
						int n_conv,
						double **res
					   )	

{
	// NB d is the total numbers of dimensions: d=#degrees + dim space

	long i,idum=-893580627;
	int j,k,eff_d=(int)(d+fmod(d,2));
	double chi2=0.0;
	double **GaussSample = dmatrix(0,2*p-1,0,eff_d-1), 
		   **Unif_dev,arg;
	double **sqrt_corr_mtx,IntermRes,q;
	
 		sqrt_corr_mtx=dmatrix(0,d-degree-1,0,d-degree-1);	
		Unif_dev=dmatrix( 0, p-1, 0, eff_d-1);

		switch (n_conv) {

		case 0:

//		ABSCube( Unif_dev, 0, p-1, 0, eff_d-1, 0, 0, &idum );

		Unif_dev=GetUnifDev(p,
							eff_d);  // get the uniform sample matrix 
		
			break;
		default: 

		Unif_dev=dmatrix(0,p-1,0,eff_d-1);
			for (i=0;i<p;i++) 
			{
				for (j=0;j<eff_d;j++) {
					Unif_dev[i][j]=uniform(&idum);
				}
			}
			
			break;
		}

		CholDec	(&CholAlg,
				 corr_mtx,       // performs Choleski decomposition on the corr. matrix
				 sqrt_corr_mtx,
				 d-degree);

		for (i=0;i<p;i++) {
			for (j=0;j < (int)(eff_d)/2 ;j++) {

				arg=sqrt(-2.*log(Unif_dev[i][2*j]));
				GaussSample[i][2*j] = arg*cos(2.*SRT_PI*Unif_dev[i][2*j+1]); 			
				GaussSample[i][2*j+1] = arg*sin(2.*SRT_PI*Unif_dev[i][2*j+1]);	

				GaussSample[i+p][2*j] = -GaussSample[i][2*j]; 			
				GaussSample[i+p][2*j+1] = -GaussSample[i][2*j+1];	
			
			}
					
						// correlates the first two. 

			for (j=0; j<d-degree; j++){
				res[i][j]=0.0;
				res[i+p][j]=0.0;

				for (k=0; k<d-degree; k++) {
					res[i][j]+=sqrt_corr_mtx[j][k]*GaussSample[i][k];
					res[i+p][j]+=sqrt_corr_mtx[j][k]*GaussSample[i][k];
				}

				res[i][j]+=mean_v[j];
				res[i+p][j]+=mean_v[j];
			}
					
			///////////////////////////////////////////////////////////////////////////////
			// computes the chi2 draw
			///////////////////////////////////////////////////////////////////////////////

			chi2=0.0;
			for (j=d-degree; j<d; j++) 	chi2+=GaussSample[i+p][j]*GaussSample[i+p][j];


			for (j=0; j<d-degree; j++) {
				if (degree > 0) {           // checks if we are in the gauassian copula case (degree = 0) 
						IntermRes=sqrt(degree/chi2)*res[i+p][j];
				} else {                    // we are in the gaussian copula case
						IntermRes=res[i+p][j];
				}

// computes a draw of the Student variable

				IntermRes=Student_Dis   (   degree,
											IntermRes
										);                                

				res[i+p][j]=interp_columns(	ya, 
											xa, 
											n_pts, 
											IntermRes, 
											0, 
											&q,
											j); 


			}
						// computes the chi2 draw
			chi2=0.0;
			for (j=d-degree; j<d; j++) 	chi2+=GaussSample[i][j]*GaussSample[i][j];


			for (j=0; j<d-degree; j++) {
				if (degree > 0) {           // checks if we are in the gauassian copula case (degree = 0) 
						IntermRes=sqrt(degree/chi2)*res[i][j];
				} else {                    // we are in the gaussian copula case
						IntermRes=res[i][j];
				}

// computes a draw of the Student variable

				IntermRes=Student_Dis(	degree,
										IntermRes
										);                                

																		   			
				res[i][j]=interp_columns(	ya, 
											xa, 
											n_pts, 
											IntermRes, 
											0, 
											&q,
											j);
			
			}
		
			 // computes a draw of the Copula
		}
				
	free_dmatrix(Unif_dev,0,p-1,0,eff_d-1);
	free_dmatrix(GaussSample,0,2*p-1,0,eff_d-1);
	free_dmatrix(sqrt_corr_mtx,0,d-degree-1,0,d-degree-1);

}

*/

/* ========================================================================
   FUNC: GetStudentCplDev
   DESC: Generates  a matrix (0,p-1)x(0,d-1) that contains a sample of correlated 
		 random variables according to a Student Copula of degree n, and compatible
		 with a set of marginal cumulative distributions

		degree:  degree of the Student Copula (must be <=10)
		mean_v: vector of 0
		corr_mtx: linear correlatio  matrix (extracted from the rank corr. matrix)
		p: number of simulations
		d: total n. of degrees of freedom, i.e. degree + n. of assets
		xa: a matrix (0,n_pts-1)x(0,n_assets-1) containing the abscissas of the marginal
			cumulative distributions
		ya: a matrix (0,n_pts-1)x(0,n_assets-1) containing the ordinates of the marginal
			cumulative distributions
		n_pts: number of abscissas
		n_conv: number of convolutions
		
   MODIFIES:
   DECLARATION:
   ======================================================================== */

	void GetStudentCplDev (
						const int degree,
	  					const double *mean_v,
						double **corr_mtx,   
						const long p,
						const int d, 
						double **xa,
						double **ya,
						const long n_pts,
						int n_conv,
						SrtMCSamType  MCType,
						double **res
					   )	

{
	// NB d is the total numbers of dimensions: d=#degrees + dim space

	long i,idum=-8935807;
	int j,k,eff_d=d-degree;
	double	chi2=0.0,
			***GaussSample=NULL,**UnifSample=NULL,
			**sqrt_corr_mtx=NULL,
			IntermRes,q;
				
		//////////  Memory allocation
		GaussSample=dcube(0,p,0,d-1,0,0);
		sqrt_corr_mtx=dmatrix(0,eff_d-1,0,eff_d-1);

		//////////  Step1: performs Chol decomposition of the correlation matrix corr_mtx
		CholDec	(&CholAlg,
				 corr_mtx,       
				 sqrt_corr_mtx,
				 eff_d);

		switch(MCType) {
		case ABS:    // uses ABS cube

	 		ABSCube( GaussSample,0,p-1,0,d-1,0,0,&idum);
			//////////  Step3: Generates the matrix of r.v. compatible with the Copula and the marginals
			for (i=0;i<p;i++) {	
					
				for (j=0; j<eff_d; j++){
					res[i][j]=0.0;
					for (k=0; k<eff_d ;k++) {
						res[i][j]+=sqrt_corr_mtx[j][k]*GaussSample[i][k][0];
					}
					res[i][j]+=mean_v[j];
				}
						
					// computes a draw of the Chi2 variable
				chi2=0.0;
				for (j=eff_d; j<d; j++) 	chi2+=GaussSample[i][j][0]*GaussSample[i][j][0];
				
				for (j=0; j<eff_d; j++) {
					if (degree > 0) {           // checks if we are in the gauassian copula case (degree = 0) 
							IntermRes=sqrt(degree/chi2)*res[i][j];
					} else {                    // we are in the gaussian copula case
							IntermRes=res[i][j];
					}		

					// computes a draw of the Student variable	

					IntermRes=Student_Dis(	degree,IntermRes);                                
																		   		
					res[i][j]=interp_columns(	ya, 
												xa, 
												n_pts, 
												IntermRes, 
												0, 
												&q,
												j);			
				}
			}

			break;

		case SOBOLBM:   // uses SOBOL with Box-Muller algorithm
			break;
		
		case RANDOM_GAUSS:

			UnifSample=dmatrix(0,p,0,d-1);
			gauss_matrix(UnifSample, 0,p,0,d-1, &idum); // note that in this case the function returns a gaussian
			                                            // matrix of uncorrelated gaussian rv, even though the name
			                                            // is UnifSample

			//////////  Step3: Generates the matrix of r.v. compatible with the Copula and the marginals
			for (i=0;i<p;i++) {	
					
				
				for (j=0; j<eff_d; j++){
					res[i][j]=0.0;
					for (k=0; k<eff_d ;k++) {
						res[i][j]+=sqrt_corr_mtx[j][k]*UnifSample[i][k];
					}
					res[i][j]+=mean_v[j];
				}
						
					// computes a draw of the Chi2 variable
				chi2=0.0;
				for (j=eff_d; j<d; j++) 	chi2+=UnifSample[i][j]*UnifSample[i][j];
				
				for (j=0; j<eff_d; j++) {
					if (degree > 0) {           // checks if we are in the gauassian copula case (degree = 0) 
							IntermRes=sqrt(degree/chi2)*res[i][j];
					} else {                    // we are in the gaussian copula case
							IntermRes=res[i][j];
					}		

					// computes a draw of the Student variable	

					IntermRes=Student_Dis(	degree,IntermRes);                                
																		   		
					res[i][j]=interp_columns(	ya, 
												xa, 
												n_pts, 
												IntermRes, 
												0, 
												&q,
												j);			
				}
			}


			free_dmatrix(UnifSample,0,p,0,d-1);
			break;

		default:   // by default uses SOBOL with Normal Inversion 

			UnifSample=dmatrix(0,p,0,d-1);
			GetSobolMatrix(p,d,UnifSample);
			
			//////////  Step3: Generates the matrix of r.v. compatible with the Copula and the marginals
			for (i=1;i<p;i++) {	
					
				for (j=0; j<d; j++) GaussSample[i][j][0]=inv_cumnorm_fast(UnifSample[i][j]);
				
				for (j=0; j<eff_d; j++){
					res[i][j]=0.0;
					for (k=0; k<eff_d ;k++) {
						res[i][j]+=sqrt_corr_mtx[j][k]*GaussSample[i][k][0];
					}
					res[i][j]+=mean_v[j];
				}
						
					// computes a draw of the Chi2 variable
				chi2=0.0;
				for (j=eff_d; j<d; j++) 	chi2+=GaussSample[i][j][0]*GaussSample[i][j][0];
				
				for (j=0; j<eff_d; j++) {
					if (degree > 0) {           // checks if we are in the gauassian copula case (degree = 0) 
							IntermRes=sqrt(degree/chi2)*res[i][j];
					} else {                    // we are in the gaussian copula case
							IntermRes=res[i][j];
					}		

					// computes a draw of the Student variable	

					IntermRes=Student_Dis(	degree,IntermRes);                                
																		   		
					res[i][j]=interp_columns(	ya, 
												xa, 
												n_pts, 
												IntermRes, 
												0, 
												&q,
												j);			
				}
			}
			free_dmatrix(UnifSample,0,p,0,d-1);

			break;
		}
				
	free_dcube(GaussSample,0,p,0,d-1,0,0);
	free_dmatrix(sqrt_corr_mtx,0,eff_d-1,0,eff_d-1);

}



/* ========================================================================
   FUNC: GetStudentCplDev
   DESC: Same function as the previous one, modified for GRS applications. At 
		 each call, it generates a vector of r.v. associated to a Student Copula
		 of degree "degree".

		degree:  degree of the Student Copula (must be <=10)
		mean_v: vector of 0
		corr_mtx: linear correlatio  matrix (extracted from the rank corr. matrix)
		p: number of simulations
		d: total n. of degrees of freedom, i.e. degree + n. of assets
		xa: a matrix (0,n_pts-1)x(0,n_assets-1) containing the abscissas of the marginal
			cumulative distributions
		ya: a matrix (0,n_pts-1)x(0,n_assets-1) containing the ordinates of the marginal
			cumulative distributions
		n_pts: number of abscissas
		n_conv: number of convolutions
		
   MODIFIES:
   DECLARATION:
   ======================================================================== */

/* ------------------------------------------------------------------------------------------ */

double *GetStudentCplDev_Rand (
								const int degree,
	  							const double *mean_v,
								double **sqrt_corr_mtx,   
								const long p,
								const int d, 
								double **xa,
								double **ya,
								const long n_pts,
								int n_conv,
								long *idum,
								double *GaussSample
					          )	

{
	// NB d is the total numbers of dimensions: d=#degrees + dim space

	int j,k,eff_d=d-degree;
	double *res = dvector(0,d-degree-1),chi2=0.0;
//	double *GaussSample = dvector(0,d-1), *Unif_dev = dvector(0,d-1);
	double IntermRes,q;	

			for (j=0; j<d; j++) GaussSample[j] = gauss_sample (idum);  					
								
						// correlates the gaussian vector. 

			for (j=0; j<eff_d; j++){
				res[j]=0.0;

				for (k=0; k< eff_d; k++) {
					res[j]+=sqrt_corr_mtx[j][k]*GaussSample[k];
				}

				res[j]+=mean_v[j];
			}
					
			///////////////////////////////////////////////////////////////////////////////
			// computes the chi2 draw
			///////////////////////////////////////////////////////////////////////////////

			chi2=0.0;
			for (j=d-degree; j<d; j++) 	chi2+=GaussSample[j]*GaussSample[j];

			for (j=0; j<d-degree; j++) {
				if (degree > 0) {           // checks if we are in the gauassian copula case (degree = 0) 
						IntermRes=sqrt(degree/chi2)*res[j];
				} else {                    // we are in the gaussian copula case
						IntermRes=res[j];
				}

// computes a draw of the Student variable

				IntermRes=Student_Dis   (   degree,
											IntermRes
										);                                

				res[j]=interp_columns(	ya, 
									    xa, 
									    n_pts, 
									    IntermRes, 
										0, 
										&q, 
										j); 

			}
						// computes the chi2 draw
										
	return (res);
}