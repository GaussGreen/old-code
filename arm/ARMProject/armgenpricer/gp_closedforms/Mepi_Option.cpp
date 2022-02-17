
/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file merton.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"

#include <cmath>
#include <vector>

#include "gpclosedforms/Mepi_Option.h"
#include "gpclosedforms/vanilla_bs.h"
#include "gpclosedforms/normal.h"
#include "gpbase/numericconstant.h"
#include "gpclosedforms/basic_distributions.h"
#include "gpclosedforms/gaussian_integrals.h"
#include "gpclosedforms/PDE_Solver.h"

using std::vector;



CC_BEGIN_NAMESPACE(ARM)

const double EPS_LIMIT    = 1.0e-10;


inline double min(double x, double y) {return ((x<=y) ? x : y);}
inline double max(double x, double y) {return ((x<=y) ? y : x);}



/*		Exposition: amount of money invested in the exposed asset.
		f00	: Starting portfolio value
		gamma0 : protection at time =0: forexemple 0.8
		gamma1 : protection at time = maturity : for exemple : 1.0
		R : risk factor: for exemple 5
		Emin: percentage of the initial investment that will remain exposed exemeple : 0.3
		Lmax: percentage of the initial investment that can be borrowed to buy asset: exemple :0.7

		 in the constructor, fgrid[i] and tgrid[j] are computed 
		 and exposition[i,j] is computed 

*/


template <typename A>
Mepi_Option_PDE<A>::Mepi_Option_PDE(double f00,double T0,double K0,
												   double R0,double Emin0,double Lmax0,double gamma00,double gamma10,double sig0,
												   double lambda0,double sigJ0,double r0,double s0,double mu0,
												   double fees0,int callput0,
												   int nbSpaceSteps0,double nbSpaceStdDev0,double volDrift0,double volVol0,
												   vector<double>* tgrid0,double theta0)
												   : K(K0),R(R0),Emin(Emin0),
												   Lmax(Lmax0),gamma0(gamma00),gamma1(gamma10),sig(sig0),sigJ(sigJ0),r(r0),s(s0),
												   mu(mu0),lambda(lambda0),fees(fees0),nbSpaceSteps(nbSpaceSteps0),nbSpaceStdDev(nbSpaceStdDev0),
												   volDrift(volDrift0),volVol(volVol0),callput(callput0),timeGrid(tgrid0)
{
				timeStepsNb=timeGrid->size();
				spaceGrid=new vector<double>(2*nbSpaceSteps+1);
				tgrid=timeGrid;
				fgrid=spaceGrid;
				A::Set_Theta(theta0);
				maturity=T0;
				f0=f00;
				exposition=new vector<double>((2*nbSpaceSteps+1)*timeStepsNb);
}

template <typename A>
void Mepi_Option_PDE<A>::Mepi_Option_Determine_SpaceGrid_and_Exposition(double sigma,double limit_divisor)
{
	timeStepsNb=timeGrid->size();
	Set_Sigma(sigma);

	/// fill up of the space grid
	double T=(*timeGrid)[timeGrid->size()-1];
	double multplicator=nbSpaceStdDev*sigma*sqrt(T);
	if(multplicator>80)
	{
		multplicator=80;/// limit of the extend of the log grid : cannot go beyond the operating scale !( with some room)
	}
	double S0=f0*exp(r*T);
	double acc=multplicator/nbSpaceSteps;
	for(int i=-nbSpaceSteps;i<=nbSpaceSteps;i++)
	{
		if(limit_divisor*i>-nbSpaceSteps)
		{
			(*spaceGrid)[nbSpaceSteps+i]=S0*exp(acc*i);
		}
		else
		{
			(*spaceGrid)[nbSpaceSteps+i]=S0*exp(-acc*nbSpaceSteps/limit_divisor)*(1.+(i+nbSpaceSteps/limit_divisor)*(1.-exp(-acc)));
		}
	}
	/// fill up of the precomputed exposition
	for(int k=0;k<timeStepsNb;k++)
				{
		for(int im=-nbSpaceSteps;im<=nbSpaceSteps;im++)
		{
			double fk=(*spaceGrid)[nbSpaceSteps+im];
			double t0=(*timeGrid)[k];
			double f1=(gamma0+(*timeGrid)[k]/maturity*(gamma1-gamma0))*f0;
			double f2=max((fk-(gamma0+(*timeGrid)[k]/maturity*(gamma1-gamma0))*f0)*R,Emin*fk);
			double expo=	min(max((fk-(gamma0+(*timeGrid)[k]/maturity*(gamma1-gamma0))*f0)*R,Emin*fk),fk+Lmax);
			(*exposition)[timeStepsNb*(nbSpaceSteps+im)+k]=
///
/////////  here we describe the investment procedure     /////////////////////////////////////////////
///
///
				min(max((fk-(gamma0+(*timeGrid)[k]/maturity*(gamma1-gamma0))*f0)*R,Emin*fk),fk+Lmax);
///
//////////////////////////////////////////////////////////////////////////////////////////////////////
///

		}
				}
}

template <typename A>
Mepi_Option_PDE<A>::~Mepi_Option_PDE()
{
	delete exposition;
	delete spaceGrid;
}

template <typename A>
inline double Mepi_Option_PDE<A>::terminalvalue(int& i)
{
				if(callput==K_CALL)
				{
					return max((*spaceGrid)[i]-K,0);
				}
				else
				{
					return max(K-(*spaceGrid)[i],0);
				}
				
}


template <typename A>
inline double Mepi_Option_PDE<A>::jumpupintegral(int& k,int& i,double& eta)
//  approx. of the Expectation(jumps above the discret limit): eta = localvoljump(k,i)=expo*sigJ
{
				if((i>0)&&(callput==K_CALL))
				{
					int spaceStepsMax=spaceGrid->size();
					double drift=localdriftjump(k,i);
					double w=((*spaceGrid)[spaceStepsMax]-(*spaceGrid)[i] -drift)/eta;
					if(fabs(drift)<EPS_LIMIT)
					{
						return eta*exp(-w*w/2.)/ARM_NumericConstants::ARM_SQRT_2_PI;
					}
					else
					{
						return eta*exp(-w*w/2.)/ARM_NumericConstants::ARM_SQRT_2_PI+drift*(NormalCDF(w)-0.5);
					}
				}
				else
				{
					return 0.0;
				}
}


template <typename A>
inline double Mepi_Option_PDE<A>::jumpdownintegral(int& k,int& i,double& eta)  //  approx. of the Expectation(jumps below the discret limit)
{
				if((i<0)&&(callput==K_PUT))
				{
					double drift=localdriftjump(k,i);
					double w=((*spaceGrid)[0]-(*spaceGrid)[i] -drift)/eta;
					if(fabs(drift)<EPS_LIMIT)
					{
						return eta*exp(-w*w/2.)/ARM_NumericConstants::ARM_SQRT_2_PI;
					}
					else
					{
						return eta*exp(-w*w/2.)/ARM_NumericConstants::ARM_SQRT_2_PI+drift*(NormalCDF(w)-0.5);
					}
				}
				else
					
				{
					return 0.0;
				}
}


template <typename A>
inline double Mepi_Option_PDE<A>::localdrift(int& k,int& i)
			{
				double expo=(*exposition)[timeStepsNb*i+k];
				double fk=(*spaceGrid)[i];
///
//////////////   Here we describe the return of the non invested liquidities ////////
///
				return r*(fk-expo)+s*min(fk-expo,0)+mu*expo-fees*fk;
///
/////////////////////////////////////////////////////////////////////////////////////
///
			}


template <typename A>
inline double Mepi_Option_PDE<A>::NormalizationService(double fn,int& k, int& i)
			{
			
///
//////////////   Here we describe the normalization of the option between its intrinsic value decreased by the fees and underlying  ////////
///
				return min(max(fn,-(*spaceGrid)[i]),(*spaceGrid)[i]);
///
/////////////////////////////////////////////////////////////////////////////////////
///
			}





//*****************************************************************************
//
//
//			Fonction that solves the PDE
//
//
//*****************************************************************************

/// Function that implemente the space grid 

vector<double>* mepi_spacegrid_create(int nb, double sig, double T ,double S0, double c, double nbSD)
{
	vector<double>* result=new vector<double>(2*nb+1);
	double acc=nbSD*sig*sqrt(T)/nb;
	for(int i=-nb;i<=nb;i++)
	{
		if(c*i>-nb)
		{
			(*result)[nb+i]=S0*exp(acc*i);
		}
		else
		{
			(*result)[nb+i]=S0*exp(-acc*nb/c)*(1.+(i+nb/c)*(1.-exp(-acc)));
		}

	}
	return result;
}

/// Function that implement the time grid

vector<double>* mepi_timegrid_create(int nb, double alpha,double maturity)
{
	vector<double>* result=new vector<double>(nb+1);

	(*result)[0]=0;
	for(int i=1;i<=nb;i++)
	{
// FIXMEFRED: mig.vc8 (23/05/2007 14:25:36):cast
		(*result)[i]=maturity*pow((double(i))/nb,alpha);
	}
	return result;
}



vector<double>* Mepi_EDP_VanillaOption(
							  double P0,
							  double f0,
							  double T,
							  double K,
							  double R,
							  double Emin,
							  double Lmax,
							  double gamma0,
							  double gamma1,  
							  double sig,
							  double lambda,
							  double sigJ,
							  double r,
							  double s,
							  double mu,
							  double fees,
							  int nsSpaceSteps,
							  double nbSpaceStdDev,
							  double volDrift,
							  double volVol,
							  double limit_divisor,
							  double relativeShift,
							  int NbHermiteComponent,
							  int CallPut,
							  vector<double>* 	timegrid,
							  int	schema
							  )
{	
	double sum=0,deltasum=0;
	double price,deltaprice;
	double maturity=(*timegrid)[timegrid->size()-1];
	ReducedGaussHermite_Coefficients c(NbHermiteComponent);
	PDE_JTX* pde1;
	vector<double>* result=new vector<double>(2);
	switch (schema)
	{
	case 1:
		{
			pde1=new Mepi_Option_PDE<PDE_JTX_Explicit_Schema>(P0,T,K,R,Emin,Lmax,gamma0,gamma1,sig,
				lambda,sigJ,r,s,mu,fees,CallPut,nsSpaceSteps, nbSpaceStdDev, volDrift, volVol,timegrid,0);
			break;
		}
	case 2:
		{
			double theta_test=0.5;
			pde1=new Mepi_Option_PDE<PDE_JTX_Theta_Schema>(P0,T,K,R,Emin,Lmax,gamma0,gamma1,sig,
				lambda,sigJ,r,s,mu,fees,CallPut,nsSpaceSteps, nbSpaceStdDev, volDrift, volVol,timegrid,theta_test);
			break;	
		}
		
	default:
		{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Mepi_EDP_VanillaOption : schema  bad input :");
			break;
		}
	}

	double weigthsum=0;
	for( int ivol=0;ivol<NbHermiteComponent;ivol++)
	{
		double currentsigma=sig*exp((volDrift-volVol*volVol/2.)*maturity/2.+c.get_point(ivol)*volVol*sqrt(2.*maturity/3.));
		
		if(currentsigma*sqrt(maturity)*nbSpaceStdDev<10.)
		{	
			weigthsum+=c.get_weight(ivol);
			switch (schema)
			{
			case 1:
				{			
					((Mepi_Option_PDE<PDE_JTX_Explicit_Schema>*)pde1)->Mepi_Option_Determine_SpaceGrid_and_Exposition( currentsigma, limit_divisor);	
					DeltaSolve<PDE_JTX_Explicit_Schema>((Mepi_Option_PDE<PDE_JTX_Explicit_Schema>*)pde1,price,deltaprice,relativeShift);
					sum+=price*c.get_weight(ivol);
					deltasum+=deltaprice*c.get_weight(ivol);
					break;
				}
			case 2:
				{
					((Mepi_Option_PDE<PDE_JTX_Theta_Schema>*)pde1)->Mepi_Option_Determine_SpaceGrid_and_Exposition( currentsigma, limit_divisor);
					DeltaSolve<PDE_JTX_Theta_Schema>((Mepi_Option_PDE<PDE_JTX_Theta_Schema>*)pde1,price,deltaprice,relativeShift);
					sum+=price*c.get_weight(ivol);
					deltasum+=deltaprice*c.get_weight(ivol);
					break;	
				}
				
			default:
				{
					throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Mepi_EDP_VanillaOption : schema  bad input :");
					break;
				}
			}
		}
		
	}

	(*result)[0]=sum/weigthsum;
	(*result)[1]=deltasum/weigthsum*P0/f0*(1.-gamma0)*R;
	return result;
	
}



CC_END_NAMESPACE()



/*---------------------------------------------------------------------------*/
/*---- End of file ----*/

