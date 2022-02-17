/*	RainbowOpt.h
        Author: D. Mayevski
        Purpose: calculate E[min(nx*Xt - kX      , ny*Yt - kY      , 0)] with Xt
   , Yt lognormal (using 1D numerical integration)
*/

#ifndef __RAINBOWOPT_H__
#define __RAINBOWOPT_H__

#include "opfnctns.h"

//	Calculates E[min(nx*Xt - kX      , ny*Yt - kY      , 0)] with Xt      ,
//Yt lognormal (using 1D numerical integration)
Err OptRainbow(double fwdx, double nx, double kx, double sigx, double fwdy,
               double ny, double ky, double sigy, double mat, double rho,
               double *res);

//	Calculates E[max(nx*Xt - ny*Yt - K      , 0)] with Xt      , Yt
//lognormal (using 1D numerical integration)
Err OptSpread(double fwdx, double nx, double sigx, double fwdy, double ny,
              double sigy, double K, double mat, double rho, double *res);

// Calculates E[max(nx*Xt - ny*Yt - K      , 0)] with Xt      , Yt following
// SABR distributions using MC Note: nsteps - number of time discretization
// steps PER YEAR
Err OptSpreadSabrMC(double fwdx, double nx, double sigx, double alphax,
                    double betax, double rhox, double fwdy, double ny,
                    double sigy, double alphay, double betay, double rhoy,
                    double K, double mat, double rho, long npaths, long nsteps,
                    double *res, double *std);

/*	----------------------------------------------------------------------------------------------------
        if call
                calculates E[max((nx*Xt + ny*Yt) - K      , 0)] with Xt      ,
   Yt lognormal (using 1D numerical integration) if put calculates E[max(K -
   (nx*Xt
   + ny*Yt)      , 0)] with Xt      , Yt lognormal (using 1D numerical
   integration)

        ----------------------------------------------------------------------------------------------------
 */
Err OptSpreadNew(double fwdx, double nx, double sigx, double fwdy, double ny,
                 double sigy, double K, double mat, double rho,
                 SrtCallPutType call_put, double *res);

#endif // #ifndef __RAINBOWOPT_H__