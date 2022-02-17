/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 * Version initiale 01/15/2004
 *
 *  basic functions for the closed form framework 
 *
 *	\file heston.h
 *
 *  \brief
 *
 *	\author  O. Croissant
 *	\version 1.0
 *	\date January 2004
 */
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

 
#ifndef _GP_CF_HESTON_INTERFACE_H
#define _GP_CF_HESTON_INTERFACE_H


#include "firsttoinc.h"
#include "gpbase/port.h"
#include "gpbase/env.h"
#include "gpbase/gplinalgtypedef.h"
#include "gpbase/typedef.h"

CC_BEGIN_NAMESPACE(ARM)

/////////////////////////////////////////////////////////////////////////////////////////////////////
///  
///			Begining Exportable  Pricing Functions 
///
/////////////////////////////////////////////////////////////////////////////////////////////////////
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
									int callorput,
									int nb);

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
									int callorput,
									int nb);

///////////////////////////////////////////////////////////////////////
///  
///			2nd Derivatives  
///
///////////////////////////////////////////////////////////////////////
double Export_GHeston_VanillaOption(int i,int j,
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
									int nb);

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
										 std::vector<double>* Maturity_Vec,
										 std::vector<double>* F_Vec,
										 std::vector<double>* InitialVol_Vec,
										 std::vector<double>* longtermV_Vec,
										 std::vector<double>* theta_Vec,
										 std::vector<double>* ksi_Vec,
										 std::vector<double>* rho_Vec,
										 std::vector<double>* lambda_Vec,
										 std::vector<double>* muJ_Vec,
										 std::vector<double>* sigmaJ_Vec,
										 int nbsteps);

double Export_GHeston_Implicit_Volatility_ModelVector(
										 double K,
										 double T,
										 int Interpolation_Method,
										 std::vector<double>* Maturity_Vec,
										 std::vector<double>* F_Vec,
										 std::vector<double>* InitialVol_Vec,
										 std::vector<double>* longtermV_Vec,
										 std::vector<double>* theta_Vec,
										 std::vector<double>* ksi_Vec,
										 std::vector<double>* rho_Vec,
										 std::vector<double>* lambda_Vec,
										 std::vector<double>* muJ_Vec,
										 std::vector<double>* sigmaJ_Vec,
										 int nbsteps);



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
						int callorput,int nb1,int nb,int NbStage, int NbOscill,double prec);

///////////////////////////////////////////////////////////////////////
///  
///			Shifted Heston 
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
						int callorput,int nb1,int nb,int NbStage, int NbOscill,double prec);

///////////////////////////////////////////////////////////////////////
///  
///			Generalized Heston 
///
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
						int callorput,int nb1,int nb,int NbStage, int NbOscill,double prec);

double Export_Heston_OptionPrice(
		double resetTime, 
		double forward, 
		double strike, 
		int callorput, 
		double v0, 
		double kappa, 
		double theta,
		double vvol,
		double rho,
		double shift,
		const std::vector<double>* times, 
		const std::vector<double>* level);

double Export_Heston2b_OptionPrice(
		double resetTime, 
		double forward, 
		double strike, 
		int callorput, 
		double v01, 
		double kappa1, 
		double theta1,
		double vvol1,
		double rho1,
		double v02, 
		double kappa2, 
		double theta2,
		double vvol2,
		double rho2,
		double shift,
		const std::vector<double>* times, 
		const std::vector<double>* level);

double Export_MixteHeston_OptionPrice(
		double resetTime, 
		double forward, 
		double strike, 
		int callorput, 
		double sigma, 
		double v0, 
		double kappa, 
		double theta,
		double vvol,
		double rho,
		double shift,
		const std::vector<double>* times, 
		const std::vector<double>* level);

///////////////////////////////////////////////////////////////////////
///  
///			Normal Heston 
///
///////////////////////////////////////////////////////////////////////

double Export_Normal_Heston_VanillaOption(double rho,double lambdaV,double thetaV,
						double kappaV,double V0,  double S0,double k,double T,double lambdaB,int callput,int nbfirst,int nb,
						int NbStage, int NbOscill,double prec);

class ARM_SABRToHestonSmileCalibration : public ARM_RootObject
{
private:
	double	resetTime;
	double	fwdrate;
	double	atmvol;
	double	initialVar;
	double	Kappa;
	double	Theta;
	double	Rho;
	double	Shift;
	double	Vvol;
	double	Level;
	bool	CalibTheta;
	bool	CalibRho;
	bool	CalibShift;
	double	caliberror;

	std::vector<double>	mktvols;
	std::vector<double>	strikes;

public:
	ARM_SABRToHestonSmileCalibration(
		double reset, double fwd, double volatm, double alpha, double beta, double rho, double nu, int flag,
		double initvar, double meanrev, double vinfini, double correl, double shift, 
		bool calibtheta, bool calibrho, bool calibshift)
	{
		resetTime	= reset;
		fwdrate		= fwd;
		atmvol		= volatm;
		initialVar	= initvar;
		Kappa		= meanrev;
		Theta		= vinfini;
		Rho			= correl;
		Shift		= shift;
		CalibTheta	= calibtheta;
		CalibRho	= calibrho;
		CalibShift	= calibshift;

		// construction des strikes et des vols
		buildStrikesAndVol(alpha, beta, rho, nu, flag);

		// calibration
		calibrate();
	}

	ARM_SABRToHestonSmileCalibration(const ARM_SABRToHestonSmileCalibration& rhs) 
	{
	}

	~ARM_SABRToHestonSmileCalibration()
	{
	}

	virtual ARM_Object*		Clone() const { return new ARM_SABRToHestonSmileCalibration(*this); };
	virtual string			ExportShortName() const;
	virtual	string			toString(const string& indent="",const string& nextIndent="") const;

	double	GetKappa() const {return Kappa;};
	double	GetV0() const {return initialVar;};
	double	GetTheta() const {return Theta;};
	double	GetRho() const {return Rho;};
	double	GetLevel() const {return Level;};
	double	GetShift() const {return Shift;};
	double	GetVVol() const {return Vvol;};
	double	GetCalibError() const {return caliberror;};

private:
	void	buildStrikesAndVol(double alpha, double beta, double rho, double nu, int flag);
	void	calibrate();

};

inline string ARM_SABRToHestonSmileCalibration::ExportShortName() const
{
	return "LSHCA";
}

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
