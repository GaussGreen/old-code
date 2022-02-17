#ifndef _GP_CF_VANILLA_NORMAL_H
#define _GP_CF_VANILLA_NORMAL_H


#include <glob/firsttoinc.h>
#include "gpbase/port.h"
#include <stdlib.h>


CC_BEGIN_NAMESPACE(ARM)

///////////////////////////////////////////////////////////////////////
///  
///			 Fonctions de base 
///
///   because we are dealing with normal models, the volatility information needs
///   to be introduces as a standard deviation .
///   therfore to convert from a lognormal volatility information :
///   standardDev=vol*F
///
///////////////////////////////////////////////////////////////////////


/// callput =  K_CALL  for call
/// callput =  K_PUT  for put


double VanillaDigitalOption_N(double F,double stddev,double k,double t,int callput);
double VegaVanillaDigitalOption_N(double S1,double stddev,double k,double t,int callput);
double VanillaIndexPaying_DigitalOption_N(double F,double stddev,double k,double t,int callput );
double VanillaOption_N( double f, double stddev, double k, double t, int callput );
double DeltaVanillaOption_N( double f, double stddev, double k, double t, int callput );
double GammaVanillaOption_N( double F, double stdDev, double k, double t, int callOrPut);
double VegaVanillaOption_N( double F, double stdDev,double k, double t,  int callOrPut );
double VanillaImpliedVol_N( double f, double target, double k, double t, int callput, double* guess = NULL, bool* success = NULL );
double VanillaCall_N_BSphi_Inverse(double param,double objective, double accuracy);
double DigitalCall_N_ImpliedVol(double f,double k, double opt,int callput);

/// Analytics for Double Corridors pricing
double DoubleDigital_N (double maturity, 
						double Kx, double spread_x,
						double Ky, double spread_y,
						double X0, double volXplus, double volXminus,
						double Y0, double volYplus, double volYminus,
						double correl, int XCap, int YCap,
						double* volX=NULL, double* volY=NULL,
						double* KxShifted=NULL,double* KyShifted=NULL,
						double* probaXY=NULL,double* probaX=NULL,double* probaY=NULL);

CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/
