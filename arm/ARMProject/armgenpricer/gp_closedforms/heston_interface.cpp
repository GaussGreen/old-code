/*!
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file heston.cpp
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/typedef.h"
#include "gpbase/gpvector.h"

#include "gpclosedforms/Simple_Heston.h"
#include "gpclosedforms/heston.h"
#include "gpclosedforms/normal_heston.h"
#include "gpclosedforms/heston_formula.h"
#include "gpclosedforms/vanille_bs_formula.h"
#include "gpclosedforms/heston_pricer.h"
#include "gpclosedforms/sabrbdiff1.h"
#include "gpclosedforms/heston_interface.h"
#include "gpclosedforms/smile_calibration.h"

#include <math.h>



///////////////////////////////////////////////////////////////////////////////
///
///					Process :
///			dS= S(rdt+V^(1/2) dW1+ Jdq
///			dV=(omega-theta*V)dt +ksi*V^(1/2) dW2 
///			dW1.dW2=rho*dt
///
///			Jump J : probability lambda, volatility sigmaJ, log-size muJ
///
/////////////////////////////////////////////////////////////////////////////:


CC_BEGIN_NAMESPACE(ARM)

#define ARM_CF_PI 3.1415926535897932
#define ARM_CF_EPS 1.0e-14

 ///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



/// callput =  1 (K_CALL) for call
/// callput = -1 (K_PUT) for put

double Export_GHeston_VanillaOption(
						double F,
						double K,
						double V0,
						double t,
						double longtermV,
						double speed,
						double volvol,
						double rho,
						double lambda,
						double muJ, 
						double sigmaJ,
						int callorput,int nb)
{
	ArgumentList a(F,K,V0,t,longtermV,speed,volvol,rho,lambda,muJ,sigmaJ,callorput,nb);

	Power_Expression<ARM_CF_Heston_JumpDiffusion_Formula> y;
	return y(a);
}

///////////////////////////////////////////////////////////////////////
///  
///			1st Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_GHeston_VanillaOption(int i,
						double F,
						double K,
						double V0,
						double t,
						double longtermV,
						double speed,
						double volvol,
						double rho,
						double lambda,
						double muJ, 
						double sigmaJ,
						int callorput,int nb)
{
	ArgumentList a(F,K,V0,t,longtermV,speed,volvol,rho,lambda,muJ,sigmaJ,callorput,nb);
	
	Power_Expression<ARM_CF_Heston_JumpDiffusion_Formula> y;
	return y(i,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_GHeston_VanillaOption(
									int i,
									int j,
									double F,
									double K,
									double V0,
									double t,
									double longtermV,
									double speed,
									double volvol,
									double rho,
									double lambda,
									double muJ, 
									double sigmaJ,
									int callorput,
									int nb)
{
	ArgumentList a(F,K,V0,t,longtermV,speed,volvol,rho,lambda,muJ,sigmaJ,callorput,nb);
	
	Power_Expression<ARM_CF_Heston_JumpDiffusion_Formula> y;
	return y(i,j,a);
}

///////////////////////////////////////////////////////////////////////
///  
///			Simple pricing from a vector of models  
///
///////////////////////////////////////////////////////////////////////

double Export_GHeston_VanillaOption_ModelVector(
										 double K,
										 double T,
										 int callput,
										 int Interpolation_Method,
										 ARM_GP_T_Vector<double>* Maturity_Vec,
										 ARM_GP_T_Vector<double>* F_Vec,
										 ARM_GP_T_Vector<double>* InitialVol_Vec,
										 ARM_GP_T_Vector<double>* longtermV_Vec,
										 ARM_GP_T_Vector<double>* theta_Vec,
										 ARM_GP_T_Vector<double>* ksi_Vec,
										 ARM_GP_T_Vector<double>* rho_Vec,
										 ARM_GP_T_Vector<double>* lambda_Vec,
										 ARM_GP_T_Vector<double>* muJ_Vec,
										 ARM_GP_T_Vector<double>* sigmaJ_Vec,
										 int nbsteps)

{
	/// Check of the size
	int n_model=Maturity_Vec->size();
	if(n_model<=0) 
	{
			throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_GHeston_VanillaOption_ModelVector::empty vector of maturity");
	}
	bool check=(InitialVol_Vec->size()==n_model)&&(longtermV_Vec->size()==n_model)&&(theta_Vec->size()==n_model)&&
		(ksi_Vec->size()==n_model)&&(rho_Vec->size()==n_model)&&(lambda_Vec->size()==n_model)&&(muJ_Vec->size()==n_model)&&(sigmaJ_Vec->size()==n_model);
	if(!check)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_GHeston_VanillaOption_ModelVector::mismatch of the vector's length");
	}
	if(n_model==1) 
	{
		ArgumentList a(F_Vec->Elt(0),K,
			InitialVol_Vec->Elt(0),
			T,
			longtermV_Vec->Elt(0),
			theta_Vec->Elt(0),
			ksi_Vec->Elt(0),
			rho_Vec->Elt(0),
			lambda_Vec->Elt(0),
			muJ_Vec->Elt(0),
			sigmaJ_Vec->Elt(0),
			callput,
			nbsteps);
		Power_Expression<ARM_CF_Heston_JumpDiffusion_Formula> y;
		return y(a);
	}
	/// on suppose que les vecteurs sont ordonnés dans le temps
	/// determination du rang dans les maturité associé a T
	int i=0;
	while ((i<n_model)&&(T > (*Maturity_Vec)[i]) )
		i++;
	///  si  T> Maturity_Vec[n_model-1] ,c'est le cas ou i>=n_model
	if(i>=n_model) 
	{
		ArgumentList a(F_Vec->Elt(n_model-1),K,
			InitialVol_Vec->Elt(n_model-1),
			T,
			longtermV_Vec->Elt(n_model-1),
			theta_Vec->Elt(n_model-1),
			ksi_Vec->Elt(n_model-1),
			rho_Vec->Elt(n_model-1),
			lambda_Vec->Elt(n_model-1),
			muJ_Vec->Elt(n_model-1),
			sigmaJ_Vec->Elt(n_model-1),
			callput,
			nbsteps);
		Power_Expression<ARM_CF_Heston_JumpDiffusion_Formula> y;
		return y(a);
	}

	/// si on est avant la premiere maturity
	if(i==0)
		{
		ArgumentList a(F_Vec->Elt(0),K,
			InitialVol_Vec->Elt(0),
			T,
			longtermV_Vec->Elt(0),
			theta_Vec->Elt(0),
			ksi_Vec->Elt(0),
			rho_Vec->Elt(0),
			lambda_Vec->Elt(0),
			muJ_Vec->Elt(0),
			sigmaJ_Vec->Elt(0),
			callput,
			nbsteps);
		Power_Expression<ARM_CF_Heston_JumpDiffusion_Formula> y;
		return y(a);
	}
	/// maintenant Maturity_Vec[i]<=T<Maturity_Vec[i+1]
	ArgumentList a1((*F_Vec)[i-1],K,
		InitialVol_Vec->Elt(i-1),
		T,
		longtermV_Vec->Elt(i-1),
		theta_Vec->Elt(i-1),
		ksi_Vec->Elt(i-1),
		rho_Vec->Elt(i-1),
		lambda_Vec->Elt(i-1),
		muJ_Vec->Elt(i-1),
		sigmaJ_Vec->Elt(i-1),
		callput,
		nbsteps);
	Power_Expression<ARM_CF_Heston_JumpDiffusion_Formula> y1;
	double c1= y1(a1);
	ArgumentList a2(F_Vec->Elt(i),K,
		InitialVol_Vec->Elt(i),
		T,
		longtermV_Vec->Elt(i),
		theta_Vec->Elt(i),
		ksi_Vec->Elt(i),
		rho_Vec->Elt(i),
		lambda_Vec->Elt(i),
		muJ_Vec->Elt(i),
		sigmaJ_Vec->Elt(i),
		callput,
		nbsteps);
	Power_Expression<ARM_CF_Heston_JumpDiffusion_Formula> y2;
	double c2= y2(a2);
	double t1=Maturity_Vec->Elt(i-1);
	double t2=Maturity_Vec->Elt(i);
	return c1+(c2-c1)*(T-t1)/(t2-t1);
}
										 
double Export_GHeston_Implicit_Volatility_ModelVector(
										 double K,
										 double T,
										 int Interpolation_Method,
										 ARM_GP_T_Vector<double>* Maturity_Vec,
										 ARM_GP_T_Vector<double>* F_Vec,
										 ARM_GP_T_Vector<double>* InitialVol_Vec,
										 ARM_GP_T_Vector<double>* longtermV_Vec,
										 ARM_GP_T_Vector<double>* theta_Vec,
										 ARM_GP_T_Vector<double>* ksi_Vec,
										 ARM_GP_T_Vector<double>* rho_Vec,
										 ARM_GP_T_Vector<double>* lambda_Vec,
										 ARM_GP_T_Vector<double>* muJ_Vec,
										 ARM_GP_T_Vector<double>* sigmaJ_Vec,
										 int nbsteps)
{
	/// computation of the option price

	double opt=Export_GHeston_VanillaOption_ModelVector(
		K,
		T,
		K_CALL,
		Interpolation_Method,
		Maturity_Vec,
		F_Vec,
		InitialVol_Vec,
		longtermV_Vec,
		theta_Vec,
		ksi_Vec,
		rho_Vec,
		lambda_Vec,
		muJ_Vec,
		sigmaJ_Vec,
		nbsteps);

	/// Coputation of the forward price

	double F;
	/// Check of the size
	int n_model=Maturity_Vec->size();
	if(n_model<=0) 
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_GHeston_VanillaOption_ModelVector::empty vector of maturity");
	}
	bool check=(InitialVol_Vec->size()==n_model)&&(longtermV_Vec->size()==n_model)&&(theta_Vec->size()==n_model)&&
		(ksi_Vec->size()==n_model)&&(rho_Vec->size()==n_model)&&(lambda_Vec->size()==n_model)&&(muJ_Vec->size()==n_model)&&(sigmaJ_Vec->size()==n_model);
	if(!check)
	{
		throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,"Export_GHeston_VanillaOption_ModelVector::mismatch of the vector's length");
	}
	if(n_model==1) 
	{
		F=F_Vec->Elt(0);
	}
	/// on suppose que les vecteurs sont ordonnés dans le temps
	/// determination du rang dans les maturité associé a T
	int i=0;
	while ((i<n_model)&&(T>Maturity_Vec->Elt(i))) 
		i++;
	///  si  T> Maturity_Vec[n_model-1] ,c'est le cas ou i>=n_model
	if(i>=n_model) 
	{
		F=F_Vec->Elt(n_model-1);
	}
	
	else
		/// si on est avant la premiere maturity
		if(i==0)
		{
			F=F_Vec->Elt(0);
		}
		/// maintenant Maturity_Vec[i]<=T<Maturity_Vec[i+1]
		else
		{
			double c1= F_Vec->Elt(i-1);
			double c2= F_Vec->Elt(i);
			double t1=Maturity_Vec->Elt(i-1);
			double t2=Maturity_Vec->Elt(i);
			F=c1+(c2-c1)*(T-t1)/(t2-t1);
		}

	///  Computation of the implicit volatility

	return ARM_CF_BS_Formula::callimplicit_totalvolatility(F,1.,
		K,K_CALL,opt,1e-12) / sqrt(T);
}

///////////////////////////////////////////////////////////////////////
///  
///			Shifted Heston 
///
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


double Export_Shifted_Heston_VanillaOption(
						double F,
						double K,
						double V0,
						double t,
						double longtermV,
						double speed,
						double volvol,
						double rho,
						double shift,
						int callorput,int nb1,int nb,int NbStage, int NbOscill,double prec)
{
	if(callorput>0)
		return Shifted_Heston(F,K,V0,t,speed ,longtermV,volvol,rho,shift,nb1,nb,NbStage,NbOscill,prec);
	else
		return (K-F) + Shifted_Heston(F,K,V0,t,speed ,longtermV,volvol,rho,shift,nb1,nb,NbStage,NbOscill,prec);
}

///////////////////////////////////////////////////////////////////////
///  
///			Generalized Heston 
///
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


double Export_Generalized_Heston_VanillaOption(
						double F,
						double K,
						double V0,
						double t,
						double longtermV,
						double speed,
						double volvol,
						double rho,
						double muJ, double sigmaJ, double lambda,
						int callorput,int nb1,int nb,int NbStage, int NbOscill,double prec)
{
	if(callorput>0)
		return GeneralizedHeston2(F,K,V0,t ,longtermV,speed,volvol,rho,muJ,sigmaJ,lambda,nb1,nb,NbStage,NbOscill,prec);
	else
		return (K-F) + GeneralizedHeston2(F,K,V0,t ,longtermV,speed,volvol,rho,muJ,sigmaJ,lambda,nb1,nb,NbStage,NbOscill,prec);
}


///////////////////////////////////////////////////////////////////////
///  
///			SABR Heston 
///
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////



double Export_SABR_Heston_VanillaOption(
						double F,
						double K,
						double V0,
						double t,
						double longtermV,
						double speed,
						double volvol,
						double rho,
						double beta,
						int callorput,int nb1,int nb,int NbStage, int NbOscill,double prec)
{
	return 0;
}

#undef ARM_CF_PI
#undef ARM_CF_EPS


// Resize the gpvector with size. then fill it first with vec and then with the last 
// element of vec.

void fillVector(int size, std::vector<double>& gpvec, const vector<double>& vec)
{
	gpvec.resize(size);

	for (int i = 0; i < size; ++i)
		if (i < vec.size())
			gpvec[i] = vec[i];
		else
			gpvec[i] = vec[vec.size()-1];
}

double Export_Heston_OptionPrice(double resetTime, double forward, double strike, int callorput, double v0, double kappa, 
								 double theta, double vvol, double rho, double shift, 
								 const std::vector<double>& times, const std::vector<double>& levels)
{
	ARM_HestonOptionPricer pricer(resetTime, forward, strike, callorput, 
								  v0, kappa, theta, rho, vvol, shift, times, levels);

	return pricer.price();
}

double Export_Heston2b_OptionPrice(double resetTime, double forward, double strike, int callorput, 
								   double v01, double kappa1, double theta1, double vvol1, double rho1, 
								   double v02, double kappa2, double theta2, double vvol2, double rho2, 
								   double shift, const std::vector<double>& times, const std::vector<double>& levels)
{
	ARM_Heston2BOptionPricer pricer(resetTime, forward, strike, callorput, 
								  v01, kappa1, theta1, rho1, vvol1,
								  v02, kappa2, theta2, rho2, vvol2,
								  shift, times, levels);

	return pricer.price();
}


double Export_MixteHeston_OptionPrice(double resetTime, double forward, double strike, int callorput, 
									  double sigma, double v0, double kappa, 
								 double theta, double vvol, double rho, double shift, 
								 const std::vector<double>& times, const std::vector<double>& levels)
{
	ARM_MixteHestonOptionPricer pricer(resetTime, forward, strike, callorput, sigma,
								  v0, kappa, theta, rho, vvol, shift, times, levels);

	return pricer.price();
}

///////////////////////////////////////////////////////////////////////
///  
///			Normal Heston 
///
///////////////////////////////////////////////////////////////////////
double Export_Normal_Heston_VanillaOption(double rho,double lambdaV,double thetaV,
						double kappaV,double V0,  double S0,double k,double T,
						double lambdaB,int callput,int nbfirst,int nb,
						int NbStage, int NbOscill,double prec)
{
	
	return NormalHeston(rho,lambdaV,thetaV,kappaV, V0,   S0, k, T, callput, lambdaB, nbfirst, nb,NbStage,  NbOscill, prec);
}



void ARM_SABRToHestonSmileCalibration::buildStrikesAndVol(double alpha, double beta, double rho, double nu, int flag)
{
	strikes.resize(5,0.);
	mktvols.resize(5,0.);
	strikes[0] = fwdrate * 0.5;
	strikes[1] = fwdrate * 0.75;
	strikes[2] = fwdrate;
	strikes[3] = fwdrate * 1.25;
	strikes[4] = fwdrate * 1.5;
	
	for(int k = 0; k < 5; k++) mktvols[k] = CptSABR_implicit_vol_direct(fwdrate,strikes[k],resetTime,alpha,beta,rho,nu,flag);

	atmvol = mktvols[2];
}

void ARM_SABRToHestonSmileCalibration::calibrate()
{
	ARM_SmileCalibration_Params_Heston params(initialVar, Kappa, Theta, Rho, Vvol, Shift, 1., 
											  false, /* calibkappa */
											  CalibTheta, 
											  CalibRho, 
											  true, /* calibvvol */
											  CalibShift);

	ARM_SmileCalibration_Heston func;

	func.Init(resetTime, fwdrate, mktvols, strikes, true, fwdrate, atmvol, &params);
	func.Calibrate();

	caliberror = func.QuadraticErr();

	Theta = params.theta();
	Vvol = params.nu();
	Rho = params.rho();
	Shift = params.shifts()[0];
	Level = params.levels()[0];
}

string ARM_SABRToHestonSmileCalibration::toString(const string& indent, const string& nextIndent) const
{
	CC_Ostringstream os;

    os << "\n\n";
    os << indent << "shifted heston smile calibration on SABR params \n";
    os << indent << "----------------------------\n";
	os << "\n\n";

	return os.str();
}

CC_END_NAMESPACE()
 


/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
