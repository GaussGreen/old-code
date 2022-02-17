/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file stochasticvol_ln.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
 
#ifndef _GP_CF_MEPI_OPTION_H
#define _GP_CF_MEPI_OPTION_H


#include "firsttoinc.h"
#include "gpbase/port.h"
#include "expt.h"
#include "gpclosedforms/PDE_Solver.h"


CC_BEGIN_NAMESPACE(ARM)

/// Function that implemente the space grid 

vector<double>* mepi_spacegrid_create(int nb, double sig, double T ,double S0, double c, double nbSD);

/// Function that implement the time grid and is used in Export_Mepi_EDP_VanillaOption

vector<double>* mepi_timegrid_create(int nb, double alpha,double maturity);



template <typename A>
class Mepi_Option_PDE : public A
		{
		public:
			int spaceStepsNb;
			int timeStepsNb;
			double K;
			double R;
			double Emin;
			double Lmax;
			double gamma0;
			double gamma1;
			double sig;
			double lambda;
			double sigJ;
			double r;
			double s;
			double mu;
			double fees;
			int nbSpaceSteps;
			double nbSpaceStdDev;
			double volDrift;
			double volVol;
			int callput;
			double currentSigma;
			vector<double>* 	timeGrid;
			vector<double>* 	spaceGrid;
			vector<double>*		exposition;
			
			Mepi_Option_PDE(double f00,double T0,double K0,
				double R0,double Emin0,double Lmax0,double gamma00,double gamma10,double sig0,
				double lambda0,double sigJ0,double r0,double s0,double mu0,double fees0,int callput0,
			int nsSpaceSteps0,double nbSpaceStdDev0,double volDrift0,double volVol0,
				vector<double>* timeGrid0,double theta0);
			
			void Mepi_Option_Determine_SpaceGrid_and_Exposition(double sigma,double limit_divisor);
			
			~Mepi_Option_PDE();
			
			inline double localvol(int& k,int& i)
			{
				return (*exposition)[timeStepsNb*i+k]*currentSigma;
			}
			inline double localvoljump(int& k,int& i)
			{
				return (*exposition)[timeStepsNb*i+k]*sigJ;
			}
			inline double localdrift(int& k,int& i);
			inline double localdriftjump(int& k,int& i)
			{
				return mu;
			}
			inline double localjump_probability(int& k)
			{
				return lambda;
			}
			inline double localriskneutralrate(int& k)
			{
				return r;
			}
			
			inline double terminalvalue(int& i);
			inline double jumpupintegral(int& k,int& i,double& eta);
			inline double jumpdownintegral(int& k,int& i,double& eta);
			inline double NormalizationService(double fn,int& k, int& i);
			void Set_Sigma(double sigma)
			{
				currentSigma=sigma;
			}
};



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
							  int nbHermiteComponent,
							  int CallPut,
							  vector<double>* timegrid,
							  int schema
							  );






CC_END_NAMESPACE()
 



#endif
/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
