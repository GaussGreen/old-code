/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file PDE_Solver.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2005
 */

#include <glob/firsttoinc.h>

#include <cmath>

#include <vector>
#include <iostream>
#include <iomanip>
#include <glob/expt.h>
#include "gpnumlib/gaussiananalytics.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/tridiagonalsolve.h"
#include "gpclosedforms/PDE_Solver.h"



CC_BEGIN_NAMESPACE(ARM)

inline double min(double x, double y) {return ((x<=y) ? x : y);}
inline double max(double x, double y) {return ((x<=y) ? y : x);}

///////////////////////////////////////////////////////////////////////////////////////////
///
///
///  Shema with Jumps , time and space dependent coefficients
///
///
///////////////////////////////////////////////////////////////////////////////////////////

const double EPS_LIMIT    = 1.0e-10;

void  PDE_JTX_Explicit_Schema::backwardcolumn(vector<double>& pn)
{
	double  fk, K0,K1,K2,lambda0, lambdap, Q, dt,df,dfup,dfdown,dfsum,sum,complementary_jumpterms;
	int L=fgrid->size(); /* should be odd */
	int TimeStepsNb=tgrid->size();
	int i,j,k;
	vector<double> pn2(L);
	for(i=0;i<L;i++)
	{
		pn[i]=terminalvalue(i);
	}
	for(k=TimeStepsNb-1;k>0;k--)
	{
		dt=(*tgrid)[k]-(*tgrid)[k-1];
		
		K0=dt*localriskneutralrate(k);
		lambda0=localjump_probability(k)*dt;
		for(i=1;i<L-1;i++)
		{
			dfup=(*fgrid)[i+1]-(*fgrid)[i];
			dfdown=(*fgrid)[i]-(*fgrid)[i-1];
			dfsum=dfup+dfdown;
			df=0.5*dfsum;
			fk=(*fgrid)[i];

/// handling of jumps
			double localvoljumpValue=localvoljump(k,i);
			Q=0.5*df*df/(localvoljumpValue*localvoljumpValue);
			
			lambdap=lambda0*df/(localvoljumpValue*ARM_NumericConstants::ARM_SQRT_2_PI);
			
			if(fabs(lambdap)>EPS_LIMIT)
			{
				if(2*i+1>L) 
				{
					sum=df*exp(-(L-i)*(L-i)*Q)/(2.*Q);
				}
				else
				{
					sum=df*exp(-(-L-i)*(-L-i)*Q)/(2.*Q);
				}
				for(j=-i;j<L-i;j++)
				{
					sum+=(pn[i+j]-pn[i])*exp(-j*j*Q);
				}
				complementary_jumpterms= jumpupintegral(k,i,localvoljumpValue)+jumpdownintegral(k,i,localvoljumpValue);
			}
			else
			{
				sum=0;
				complementary_jumpterms=0;
			}

////
			double localvolValue=localvol(k,i);
			if (fabs(localvolValue)>1.0e-12)
			{	
				K1=dt*localdrift(k,i);
				K2=dt*localvolValue*localvolValue;
				pn2[i]=pn[i]*(1-K0)+
					K2*((pn[i+1]-pn[i])/dfup-(pn[i]-pn[i-1])/dfdown)/dfsum+
					K1*(pn[i+1]-pn[i-1])/dfsum+
					lambdap*sum+lambda0*complementary_jumpterms;		
			}
			else
			{
				K1=dt*localdrift(k,i);
					pn2[i]=pn[i]*(1-K0)+
					K1*(pn[i+1]-pn[i-1])/dfsum+
					lambdap*sum+lambda0*complementary_jumpterms;	
			}
		}

		fk=(*fgrid)[0];i=0;
		K1=dt*localdrift(k,i);
		dfsum=(*fgrid)[1]-fk;
		pn2[0]=pn[0]*(1-K0)+K1*(pn[1]-pn[0])/dfsum;
		fk=(*fgrid)[L-1];i=L-1;
		K1=dt*localdrift(k,i);
		dfsum=fk-(*fgrid)[L-2];
		pn2[L-1]=pn[L-1]*(1-K0)+K1*(pn[L-1]-pn[L-2])/dfsum;
	
		for(i=0;i<L;i++)
		{
			pn[i]=NormalizationService(pn2[i],k,i);
		}
	}
	return ;
}

void  PDE_JTX_Theta_Schema::backwardcolumn(vector<double>& pn)
{
	double  fk, K0,K1,K2,lambda0, lambdap, Q, dt,df,dfup,dfdown,dfsum,sum,complementary_jumpterms;
	int L=fgrid->size(); /* should be odd */
	int TimeStepsNb=tgrid->size();
	int i,j,k;
	double un_moins_theta=1.0-theta;
	vector<double> pn2(L);
	vector<double> diagVec(L);
	vector<double> diagMoinsVec(L-1);
	vector<double> diagPlusVec(L-1);
	vector<double> constantVec(L);
	vector<double> tampon(L);
	double IR_Effect;

	for(i=0;i<L;i++)
	{
		pn[i]=terminalvalue(i);
	}
	for(k=TimeStepsNb-1;k>0;k--)
	{
		dt=(*tgrid)[k]-(*tgrid)[k-1];
		K0=dt*localriskneutralrate(k);
		lambda0=localjump_probability(k)*dt;
		for(i=1;i<L-1;i++)
		{
			dfup=(*fgrid)[i+1]-(*fgrid)[i];
			dfdown=(*fgrid)[i]-(*fgrid)[i-1];
			dfsum=dfup+dfdown;
			df=0.5*dfsum;
			fk=(*fgrid)[i];
			if (fabs(localvol(k,i))>1.0e-12)
			{
				double localvoljumpValue=localvoljump(k,i);
				double localvolValue=localvol(k,i);
				
				K1=dt*localdrift(k,i);
				K2=dt*localvolValue*localvolValue;
				Q=0.5*df*df/(localvoljumpValue*localvoljumpValue);
				lambdap=lambda0*df/(localvoljumpValue*ARM_NumericConstants::ARM_SQRT_2_PI);
				if(fabs(lambdap)>EPS_LIMIT)
				{
					if(2*i+1>L) 
					{
						sum=df*exp(-(L-i)*(L-i)*Q)/(2.*Q);
					}
					else
					{
						sum=df*exp(-(-L-i)*(-L-i)*Q)/(2.*Q);
					}
					for(j=-i;j<L-i;j++)
					{
						sum+=(pn[i+j]-pn[i])*exp(-j*j*Q);
					}
					complementary_jumpterms= jumpupintegral(k,i,localvoljumpValue)+jumpdownintegral(k,i,localvoljumpValue);
				}
				else
				{
					sum=0;
					complementary_jumpterms=0;
				}
				constantVec[i]=pn[i]+
					un_moins_theta*K2*((pn[i+1]-pn[i])/dfup-(pn[i]-pn[i-1])/dfdown)/dfsum+
					un_moins_theta*K1*(pn[i+1]-pn[i-1])/dfsum+
					lambdap*sum+lambda0*complementary_jumpterms;

				diagVec[i]=1.0+theta*(K2/(dfup*dfdown))+K0;
				IR_Effect=theta*K1/dfsum;

				diagPlusVec[i]=-theta*K2/(dfup*dfsum)-IR_Effect;
				diagMoinsVec[i-1]=-theta*K2/(dfdown*dfsum)+IR_Effect;
			}
			else
			{
				diagVec[i]=1.0;
				constantVec[i]=pn[i];
				diagPlusVec[i]=0;
				diagMoinsVec[i-1]=0;
			}
		}
	

		fk=(*fgrid)[0];i=0;
		K1=dt*localdrift(k,i);
		dfsum=(*fgrid)[1]-fk;
		constantVec[0]=pn[0]*(1-K0)+K1*(pn[1]-pn[0])/dfsum;

		fk=(*fgrid)[L-1];i=L-1;
		K1=dt*localdrift(k,i);
		dfsum=fk-(*fgrid)[L-2];
		constantVec[L-1]=pn[L-1]*(1-K0)+K1*(pn[L-1]-pn[L-2])/dfsum;

		diagVec[0]=diagVec[1];
		diagVec[L-1]=diagVec[L-2];

		diagPlusVec[0]=diagPlusVec[1];
		diagMoinsVec[L-2]=diagMoinsVec[L-3];
		tridiagonalsolve(diagMoinsVec,diagVec,diagPlusVec,constantVec,pn2,tampon);
		for(i=0;i<L;i++)
		{
				pn[i]=NormalizationService(pn2[i],k,i);
		}
	}
	return ;

}


void  PDE_JTX_MidPoint_Schema::backwardcolumn(vector<double>& pn)
{
	double  fk, K0,K1,K2,lambda0, lambdap, Q, dt,df,dfup,dfdown,sum,complementary_jumpterms;
	int L=fgrid->size(); /* should be odd */
	int TimeStepsNb=tgrid->size();
	int i,j,k;
	vector<double> pn2(L);
	vector<double> pn21(L);
	vector<double> diagVec(L);
	vector<double> diagMoinsVec(L-1);
	vector<double> diagPlusVec(L-1);
	vector<double> constantVec(L);
	vector<double> tampon(L);
	for(i=0;i<L;i++)
	{
		pn[i]=terminalvalue(i);
	}
	for(k=0;k<TimeStepsNb;k++)
	{
		if(k>0)
		{
			dt=((*tgrid)[k]-(*tgrid)[k-1])/2.;
		}
		else
		{
			dt=(*tgrid)[0]/2.;
		}
		K0=dt*localriskneutralrate(k);
		lambda0=localjump_probability(k)*dt;
		///                                  Midpoint calculation
		for(i=1;i<L-1;i++)
		{
			dfup=(*fgrid)[i+1]-(*fgrid)[i];
			dfdown=(*fgrid)[i]-(*fgrid)[i-1];
			df=0.5*(dfup+dfdown);
			fk=(*fgrid)[i];
			if (fabs(localvol(k,i))>1.0e-12)
			{
				double localvoljumpValue=localvoljump(k,i);
				double localvolValue=localvol(k,i);
				
				K1=dt*localdrift(k,i)/(2*df);
				K2=0.5*dt*localvolValue*localvolValue/(dfup*dfdown);
				Q=0.5*df*df/(localvoljumpValue*localvoljumpValue);
				lambdap=lambda0*df/(localvoljumpValue*ARM_NumericConstants::ARM_SQRT_2_PI);
				if(fabs(lambdap)>EPS_LIMIT)
				{
					if(2*i+1>L) 
					{
						sum=df*exp(-(L-i)*(L-i)*Q)/(2.*Q);
					}
					else
					{
						sum=df*exp(-(-L-i)*(-L-i)*Q)/(2.*Q);
					}
					for(j=-i;j<L-i;j++)
					{
						sum+=(pn[i+j]-pn[i])*exp(-j*j*Q);
					}
					complementary_jumpterms= jumpupintegral(k,i,localvoljumpValue)+jumpdownintegral(k,i,localvoljumpValue);
				}
				else
				{
					complementary_jumpterms=0;
					sum=0;
				}
				constantVec[i]=pn[i]*(1.-K0)+(pn[i+1]-pn[i-1])*K1
					+lambdap*sum+lambda0*complementary_jumpterms;
				diagVec[i]=1.0+2*K2;
				diagPlusVec[i]=-K2;
				diagMoinsVec[i-1]=-K2;
			}
			else
			{
				diagVec[i]=1.0;
				constantVec[i]=pn[i];
				diagPlusVec[i]=0;
				diagMoinsVec[i-1]=0;
			}
		}
		constantVec[0]=pn[0];
		constantVec[L-1]=pn[L-1];
		diagVec[0]=diagVec[1];
		diagVec[L-1]=diagVec[L-2];
		diagPlusVec[0]=diagPlusVec[1];
		diagMoinsVec[L-2]=diagMoinsVec[L-3];
		tridiagonalsolve(diagMoinsVec,diagVec,diagPlusVec,constantVec,pn21,tampon);

		///									Final point Calculation
		dt*=2.;
		K0*=2.;
		lambda0*=2.;
		for(i=1;i<L-1;i++)
		{
			dfup=(*fgrid)[i+1]-(*fgrid)[i];
			dfdown=(*fgrid)[i]-(*fgrid)[i-1];
			df=0.5*(dfup+dfdown);
			fk=(*fgrid)[i];
			if (fabs(localvol(k,i))>1.0e-12)
			{
				double localvoljumpValue=localvoljump(k,i);
				double localvolValue=localvol(k,i);
				
				K1=dt*localdrift(k,i)/(2*df);
				K2=0.5*dt*localvolValue*localvolValue/(dfup*dfdown);
				Q=0.5*df*df/(localvoljumpValue*localvoljumpValue);
				lambdap=lambda0*df/(localvoljumpValue*ARM_NumericConstants::ARM_SQRT_2_PI);
				if(fabs(lambdap)>EPS_LIMIT)
				{
					if(2*i+1>L) 
					{
						sum=df*exp(-(L-i)*(L-i)*Q)/(2.*Q);
					}
					else
					{
						sum=df*exp(-(-L-i)*(-L-i)*Q)/(2.*Q);
					}
					for(j=-i;j<L-i;j++)
					{
						sum+=(pn21[i+j]-pn21[i])*exp(-j*j*Q);
					}
					complementary_jumpterms= jumpupintegral(k,i,localvoljumpValue)+jumpdownintegral(k,i,localvoljumpValue);
				}
				else
				{
					sum=0;
					complementary_jumpterms=0;
				}
				pn2[i]=pn[i]+pn21[i]*(-2*K2-K0)+pn21[i+1]*(K2+K1)+pn21[i-1]*(K2-K1)+lambdap*sum+lambda0*complementary_jumpterms;
				
			}
			else
			{
				pn2[i]=pn21[i];
			}
		}
		pn2[0]=pn21[0];
		pn2[L-1]=pn21[L-1];
		for(i=0;i<L;i++)
		{
			pn[i]=NormalizationService(pn2[i],k,i);
		}
	}
	return ;

}
///////////////////////////////////////////////////////////////////////////////////////////
///
///
///  Shema diffusion only time and space dependent coefficients
///
///
///////////////////////////////////////////////////////////////////////////////////////////


void  PDE_TX_Explicit_Schema::backwardcolumn(vector<double>& pn)
{
	double  fk, K0,K1,K2, dt,df,dfup,dfdown,dfsum;
	int L=fgrid->size(); /* should be odd */
	int TimeStepsNb=tgrid->size();
	int i,k;
	vector<double> pn2(L);
	for(i=0;i<L;i++)
	{
		pn[i]=terminalvalue(i);
	}
	for(k=TimeStepsNb-1;k>0;k--)
	{
		dt=(*tgrid)[k]-(*tgrid)[k-1];
		
		K0=dt*localriskneutralrate(k);
		for(i=1;i<L-1;i++)
		{
			dfup=(*fgrid)[i+1]-(*fgrid)[i];
			dfdown=(*fgrid)[i]-(*fgrid)[i-1];
			dfsum=dfup+dfdown;
			df=0.5*dfsum;
			fk=(*fgrid)[i];
			if (fabs(localvol(k,i))>1.0e-12)
			{
				double localvolValue=localvol(k,i);
				
				K1=dt*localdrift(k,i);
				K2=dt*localvolValue*localvolValue;
	
				pn2[i]=pn[i]*(1-K0)+
					K2*((pn[i+1]-pn[i])/dfup-(pn[i]-pn[i-1])/dfdown)/dfsum+
					K1*(pn[i+1]-pn[i-1])/dfsum;
				
			}
			else
			{
				pn2[i]=pn[i];
			}
		}


		pn2[0]=pn[0];
		pn2[L-1]=pn[L-1];
		for(i=0;i<L;i++)
		{
			pn[i]=NormalizationService(pn2[i],k,i);
		}
	}
	return ;
}

void  PDE_TX_Theta_Schema::backwardcolumn(vector<double>& pn)
{
	double  fk, K0,K1,K2, dt,df,dfup,dfdown,dfsum;
	int L=fgrid->size(); /* should be odd */
	int TimeStepsNb=tgrid->size();
	int i,k;
	double un_moins_theta=1.0-theta;
	vector<double> pn2(L);
	vector<double> diagVec(L);
	vector<double> diagMoinsVec(L-1);
	vector<double> diagPlusVec(L-1);
	vector<double> constantVec(L);
	vector<double> tampon(L);
	double IR_Effect;

	for(i=0;i<L;i++)
	{
		pn[i]=terminalvalue(i);
	}
	for(k=TimeStepsNb-1;k>0;k--)
	{
		dt=(*tgrid)[k]-(*tgrid)[k-1];
		K0=dt*localriskneutralrate(k);
		for(i=1;i<L-1;i++)
		{
			dfup=(*fgrid)[i+1]-(*fgrid)[i];
			dfdown=(*fgrid)[i]-(*fgrid)[i-1];
			dfsum=dfup+dfdown;
			df=0.5*dfsum;
			fk=(*fgrid)[i];
			if (fabs(localvol(k,i))>1.0e-12)
			{
				double localvolValue=localvol(k,i);
				
				K1=dt*localdrift(k,i);
				K2=dt*localvolValue*localvolValue;
			
				constantVec[i]=pn[i]+
					un_moins_theta*K2*((pn[i+1]-pn[i])/dfup-(pn[i]-pn[i-1])/dfdown)/dfsum+
					un_moins_theta*K1*(pn[i+1]-pn[i-1])/dfsum;

				diagVec[i]=1.0+theta*(K2/(dfup*dfdown))+K0;
				IR_Effect=theta*K1/dfsum;

				diagPlusVec[i]=-theta*K2/(dfup*dfsum)-IR_Effect;
				diagMoinsVec[i-1]=-theta*K2/(dfdown*dfsum)+IR_Effect;
			}
			else
			{
				diagVec[i]=1.0;
				constantVec[i]=pn[i];
				diagPlusVec[i]=0;
				diagMoinsVec[i-1]=0;
			}
		}
		constantVec[0]=pn[0];
		constantVec[L-1]=pn[L-1];
		diagVec[0]=diagVec[1];
		diagVec[L-1]=diagVec[L-2];
		diagPlusVec[0]=diagPlusVec[1];
		diagMoinsVec[L-2]=diagMoinsVec[L-3];
		tridiagonalsolve(diagMoinsVec,diagVec,diagPlusVec,constantVec,pn2,tampon);
		for(i=0;i<L;i++)
		{
			pn[i]=NormalizationService(pn2[i],k,i);
		}
	}
	return ;

}


///////////////////////////////////////////////////////////////////////////////////////////
///
///
///  Shema with Diffusion Only , time  dependent coefficients
///
///
///////////////////////////////////////////////////////////////////////////////////////////

void  PDE_T_Explicit_Schema::backwardcolumn(vector<double>& pn)
{
	double  fk, K0,K1,K2, dt,df,dfup,dfdown,dfsum;
	int L=fgrid->size(); /* should be odd */
	int TimeStepsNb=tgrid->size();
	int i,k;
	vector<double> pn2(L);
	for(i=0;i<L;i++)
	{
		pn[i]=terminalvalue(i);
	}
	for(k=TimeStepsNb-1;k>0;k--)
	{
		dt=(*tgrid)[k]-(*tgrid)[k-1];
		
		K0=dt*localriskneutralrate(k);
		for(i=1;i<L-1;i++)
		{
			dfup=(*fgrid)[i+1]-(*fgrid)[i];
			dfdown=(*fgrid)[i]-(*fgrid)[i-1];
			dfsum=dfup+dfdown;
			df=0.5*dfsum;
			fk=(*fgrid)[i];
			if (fabs(localvol(k))>1.0e-12)
			{
				double localvolValue=localvol(k);
				
				K1=dt*localdrift(k);
				K2=dt*localvolValue*localvolValue;
	
				pn2[i]=pn[i]*(1-K0)+
					K2*((pn[i+1]-pn[i])/dfup-(pn[i]-pn[i-1])/dfdown)/dfsum+
					K1*(pn[i+1]-pn[i-1])/dfsum;
				
			}
			else
			{
				pn2[i]=pn[i];
			}
		}
		pn2[0]=pn[0];
		pn2[L-1]=pn[L-1];
		for(i=0;i<L;i++)
		{
			pn[i]=NormalizationService(pn2[i],k);
		}
	}
	return ;
}

void  PDE_T_Theta_Schema::backwardcolumn(vector<double>& pn)
{
	double  fk, K0,K1,K2, dt,df,dfup,dfdown,dfsum;
	int L=fgrid->size(); /* should be odd */
	int TimeStepsNb=tgrid->size();
	int i,k;
	double un_moins_theta=1.0-theta;
	vector<double> pn2(L);
	vector<double> diagVec(L);
	vector<double> diagMoinsVec(L-1);
	vector<double> diagPlusVec(L-1);
	vector<double> constantVec(L);
	vector<double> tampon(L);
	double IR_Effect;

	for(i=0;i<L;i++)
	{
		pn[i]=terminalvalue(i);
	}
	for(k=TimeStepsNb-1;k>0;k--)
	{
		dt=(*tgrid)[k]-(*tgrid)[k-1];
		K0=dt*localriskneutralrate(k);
		for(i=1;i<L-1;i++)
		{
			dfup=(*fgrid)[i+1]-(*fgrid)[i];
			dfdown=(*fgrid)[i]-(*fgrid)[i-1];
			dfsum=dfup+dfdown;
			df=0.5*dfsum;
			fk=(*fgrid)[i];
			if (fabs(localvol(k))>1.0e-12)
			{
				double localvolValue=localvol(k);
				
				K1=dt*localdrift(k);
				K2=dt*localvolValue*localvolValue;
			
				constantVec[i]=pn[i]+
					un_moins_theta*K2*((pn[i+1]-pn[i])/dfup-(pn[i]-pn[i-1])/dfdown)/dfsum+
					un_moins_theta*K1*(pn[i+1]-pn[i-1])/dfsum;

				diagVec[i]=1.0+theta*(K2/(dfup*dfdown))+K0;
				IR_Effect=theta*K1/dfsum;

				diagPlusVec[i]=-theta*K2/(dfup*dfsum)-IR_Effect;
				diagMoinsVec[i-1]=-theta*K2/(dfdown*dfsum)+IR_Effect;
			}
			else
			{
				diagVec[i]=1.0;
				constantVec[i]=pn[i];
				diagPlusVec[i]=0;
				diagMoinsVec[i-1]=0;
			}
		}
		constantVec[0]=pn[0];
		constantVec[L-1]=pn[L-1];
		diagVec[0]=diagVec[1];
		diagVec[L-1]=diagVec[L-2];
		diagPlusVec[0]=diagPlusVec[1];
		diagMoinsVec[L-2]=diagMoinsVec[L-3];
		tridiagonalsolve(diagMoinsVec,diagVec,diagPlusVec,constantVec,pn2,tampon);
		for(i=0;i<L;i++)
		{
			pn[i]=NormalizationService(pn2[i],k);
		}
	}
	return ;

}



///////////////////////////////////////////////////////////////////////////////////////////
///
///
///  Shema diffusion only  space dependent coefficients
///
///
///////////////////////////////////////////////////////////////////////////////////////////

void  PDE_X_Explicit_Schema::backwardcolumn(vector<double>& pn)
{
	double  fk, K0,K1,K2, dt,df,dfup,dfdown,dfsum;
	int L=fgrid->size(); /* should be odd */
	int TimeStepsNb=tgrid->size();
	int i,k;
	vector<double> pn2(L);
	for(i=0;i<L;i++)
	{
		pn[i]=terminalvalue(i);
	}
	for(k=TimeStepsNb-1;k>0;k--)
	{
		dt=(*tgrid)[k]-(*tgrid)[k-1];
		
		K0=dt*localriskneutralrate();
		for(i=1;i<L-1;i++)
		{
			dfup=(*fgrid)[i+1]-(*fgrid)[i];
			dfdown=(*fgrid)[i]-(*fgrid)[i-1];
			dfsum=dfup+dfdown;
			df=0.5*dfsum;
			fk=(*fgrid)[i];
			if (fabs(localvol(i))>1.0e-12)
			{
				double localvolValue=localvol(i);
				
				K1=dt*localdrift(i);
				K2=dt*localvolValue*localvolValue;
	
				pn2[i]=pn[i]*(1-K0)+
					K2*((pn[i+1]-pn[i])/dfup-(pn[i]-pn[i-1])/dfdown)/dfsum+
					K1*(pn[i+1]-pn[i-1])/dfsum;
				
			}
			else
			{
				pn2[i]=pn[i];
			}
		}
		pn2[0]=pn[0];
		pn2[L-1]=pn[L-1];
		for(i=0;i<L;i++)
		{
			pn[i]=NormalizationService(pn2[i],i);
		}
	}
	return ;
}

void  PDE_X_Theta_Schema::backwardcolumn(vector<double>& pn)
{
	double  fk, K0,K1,K2, dt,df,dfup,dfdown,dfsum;
	int L=fgrid->size(); /* should be odd */
	int TimeStepsNb=tgrid->size();
	int i,k;
	double un_moins_theta=1.0-theta;
	vector<double> pn2(L);
	vector<double> diagVec(L);
	vector<double> diagMoinsVec(L-1);
	vector<double> diagPlusVec(L-1);
	vector<double> constantVec(L);
	vector<double> tampon(L);
	double IR_Effect;

	for(i=0;i<L;i++)
	{
		pn[i]=terminalvalue(i);
	}
	for(k=TimeStepsNb-1;k>0;k--)
	{
		dt=(*tgrid)[k]-(*tgrid)[k-1];
		K0=dt*localriskneutralrate();
		for(i=1;i<L-1;i++)
		{
			dfup=(*fgrid)[i+1]-(*fgrid)[i];
			dfdown=(*fgrid)[i]-(*fgrid)[i-1];
			dfsum=dfup+dfdown;
			df=0.5*dfsum;
			fk=(*fgrid)[i];
			if (fabs(localvol(i))>1.0e-12)
			{
				double localvolValue=localvol(i);
				
				K1=dt*localdrift(i);
				K2=dt*localvolValue*localvolValue;
			
				constantVec[i]=pn[i]+
					un_moins_theta*K2*((pn[i+1]-pn[i])/dfup-(pn[i]-pn[i-1])/dfdown)/dfsum+
					un_moins_theta*K1*(pn[i+1]-pn[i-1])/dfsum;

				diagVec[i]=1.0+theta*(K2/(dfup*dfdown))+K0;
				IR_Effect=theta*K1/dfsum;

				diagPlusVec[i]=-theta*K2/(dfup*dfsum)-IR_Effect;
				diagMoinsVec[i-1]=-theta*K2/(dfdown*dfsum)+IR_Effect;
			}
			else
			{
				diagVec[i]=1.0;
				constantVec[i]=pn[i];
				diagPlusVec[i]=0;
				diagMoinsVec[i-1]=0;
			}
		}
		constantVec[0]=pn[0];
		constantVec[L-1]=pn[L-1];
		diagVec[0]=diagVec[1];
		diagVec[L-1]=diagVec[L-2];
		diagPlusVec[0]=diagPlusVec[1];
		diagMoinsVec[L-2]=diagMoinsVec[L-3];
		tridiagonalsolve(diagMoinsVec,diagVec,diagPlusVec,constantVec,pn2,tampon);
		for(i=0;i<L;i++)
		{
			pn[i]=NormalizationService(pn2[i],i);
		}
	}
	return ;

}



CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/