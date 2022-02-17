/*******************************************************************************
**                      Include Files
*******************************************************************************/

#include "utallhdr.h"
#include <OPFNCTNS.H>
#include <math.h"

/******************************************************************************/

/****************************************************************************
computes the following probability :

        Prob ( x(mux  ,sigx  ,tx) > k || min y (muy  ,sigy  ,ty) >l )
                     cp                          du
               x(mux  ,sigx  ,tx) < k || Max y (muy  ,sigy  ,ty) <l
****************************************************************************/

static double probascudout(double mux, double sigx, double matx, double k,
                           double muy, double sigy, double maty, double l,
                           double cpx, /*1 if call : x > k */
                           double duy, /*1 if down : min y > l */
                           double rho);

static double probascudout(double mux, double sigx, double matx, double k,
                           double muy, double sigy, double maty, double l,
                           double cpx, /*1 if call : x > k */
                           double duy, /*1 if down : min y > l */
                           double rho) {
  double result;
  double sigx_sqrt = sigx * sqrt(matx);
  double sigy_sqrt = sigy * sqrt(maty);

  result = bivar(cpx * (-k + mux * matx) / sigx_sqrt,
                 duy * (-l + muy * maty) / sigy_sqrt,
                 cpx * duy * rho * sqrt(maty / matx));
  result -=
      bivar(cpx * (-k + mux * matx + 2 * rho * sigx / sigy * l) / sigx_sqrt,
            duy * (l + muy * maty) / sigy_sqrt,
            duy * cpx * rho * sqrt(maty / matx)) *
      exp(2 * muy / (sigy * sigy) * l);

  return (result);
}

/*******************************************************************************
*
* FUNCTION     	: srt_f_optscdspr(14)
*
* PURPOSE      	: SCUD Option : X is the underlying
                                the scud is on the spread between Z and Y
                                (Z is the foreign asset; Y the domestic one)
*
* DESCRIPTION  	: Barrier on a spread between two underlyings
* 		  Premium is computed for OUT and then transformed for IN
*                 if necessary  , using the fact that OUT + IN = BS
*
* CALLS		: bivar()
*
* PARAMETERS   	: fwdx     	- forward price of 1st underlying x
*              	: fwdy          - forward price of 2nd underlying y
*              	: fwdz          - forward price of 3rd underlying z
*              	: spoty         - spot price of 2nd underlying y
*              	: spotz         - spot price of 3rd underlying z
*              	: strike        - strike price
*              	: barrier       - barrier level
*              	: sigx          - vol of 1st underlying x
*              	: sigy          - vol of 2nd underlying y
*              	: sigz          - vol of 3rd underlying z
*              	: rhoxy         - correlation between x and y
*              	: rhoyz         - correlation between y and z
*              	: rhoxz         - correlation between x and z
*              	: matx          - maturity of option  , in years
*              	: matyz         - maturity of scud spread in years
*              	: disc          - discount factor
*              	: call_put      - type of option: 0 call  , 1 put
*              	: scud_type     - type of scud: down - up
*		: ext_lit	- type of scud: out - in
*              	: greek        	- prem  , delta  , gamma  ,...
*
* RETURNS      	: double        - depends on greek
*
*******************************************************************************/

double srt_f_optscdspr(double fwdx, double fwdy, double fwdz, double spoty,
                       double spotz, double strike, double barrier, double sigx,
                       double sigy, double sigz, double rhoxy, double rhoyz,
                       double rhoxz, double matx, double matyz, double disc,
                       SrtCallPutType call_put, /* CALL or PUT */
                       int scud_type,           /* DOWN or UP*/
                       int ext_lit,             /* EXTING or LIGHT*/
                       SrtGreekType greek) {
  double voly, volz, volyz;

  double rhoxyz;

  double l;
  double k;

  double mux;
  double muyz;

  double mu1x;
  double mu2x;
  double mu1yz;
  double mu2yz;
  double premium;

  double deltax;
  double deltay;
  double deltaz;
  double gammax;
  double gammay;
  double gammaz;
  double gammaxy;
  double gammaxz;
  double gammayz;
  double vegax;
  double vegay;
  double vegaz;
  double theta;

  double shift;
  double shiftx;
  double shifty;
  double shiftz;

  int cpx;
  int duyz;

  cpx = (call_put == SRT_CALL) ? 1 : -1; /* 1 if call; -1 if put */
  duyz = (scud_type == 0) ? 1 : -1;      /* 1 if down; -1 if up  */

  if (matyz < 0) {
    premium =
        srt_f_optblksch(fwdx, strike, sigx, matx, disc, call_put, PREMIUM);
    if (((scud_type == 0) && (ext_lit == 0)) || /* Down and out*/
        ((scud_type == 1) && (ext_lit == 1)))   /* Up and in */
    {
      if ((spotz - spoty) < barrier)
        premium = 0;
    } else if (((scud_type == 0) && (ext_lit == 1)) || /* Down and in*/
               ((scud_type == 1) && (ext_lit == 0)))   /* Up and out */
    {
      if ((spotz - spoty) > barrier)
        premium = 0;
    }
    return (premium);
  } else if (((spotz - spoty) - barrier) * duyz < 0)
    premium = 0;
  else {

    voly = fwdy * sqrt((exp(sigy * sigy * matyz) - 1) / matyz);
    volz = fwdz * sqrt((exp(sigz * sigz * matyz) - 1) / matyz);

    volyz = sqrt(volz * volz + voly * voly - 2 * rhoyz * voly * volz);

    l = barrier - (spotz - spoty);
    k = log(strike / fwdx);

    /**** 1 corresponds to the domestic forward probability measure ****/
    /**** 2 corresponds to the foreign  forward probability measure ****/

    mux = 0.0;
    mu1x = mux - (0.5 * sigx * sigx);
    mu2x = mux + (0.5 * sigx * sigx);

    muyz = (fwdz - spotz) / matyz - (fwdy - spoty) / matyz;

    mu1yz = muyz - rhoxz * sigx * volz;
    mu2yz = muyz - rhoxy * sigx * voly;

    rhoxyz = (rhoxz * volz - rhoxy * voly) /
             sqrt(volz * volz + voly * voly - 2 * rhoyz * volz * voly);

    premium = fwdx * probascudout(mu2x, sigx, matx, k, mu2yz, volyz, matyz, l,
                                  cpx, duyz, rhoxyz);

    premium -= strike * probascudout(mu1x, sigx, matx, k, mu1yz, volyz, matyz,
                                     l, cpx, duyz, rhoxyz);

    premium *= disc;
    premium *= cpx;
  }

  switch (greek) {
    /****************************     PREMIUM 	******************************/
  case PREMIUM:
    if (ext_lit == 1) {
      premium = -premium;
      premium +=
          srt_f_optblksch(fwdx, strike, sigx, matx, disc, call_put, PREMIUM);
    }
    return (premium);
    break;

    /****************************     DELTAS 	******************************/

  case DELTAX:
    shift = fwdx / 10000;
    deltax =
        (srt_f_optscdspr(fwdx + shift, fwdy, fwdz, spoty, spotz, strike,
                         barrier, sigx, sigy, sigz, rhoxy, rhoyz, rhoxz, matx,
                         matyz, disc, call_put, scud_type, ext_lit, PREMIUM) -
         premium) /
        shift;
    return (deltax);
    break;

  case DELTAY:
    shift = spoty / 10000;
    deltay = (srt_f_optscdspr(fwdx, fwdy * (1 + shift / spoty), fwdz,
                              spoty + shift, spotz, strike, barrier, sigx, sigy,
                              sigz, rhoxy, rhoyz, rhoxz, matx, matyz, disc,
                              call_put, scud_type, ext_lit, PREMIUM) -
              premium) /
             shift;
    return (deltay);
    break;

  case DELTAZ:
    shift = spotz / 10000;
    deltaz = (srt_f_optscdspr(fwdx, fwdy, fwdz * (1 + shift / spotz), spoty,
                              spotz + shift, strike, barrier, sigx, sigy, sigz,
                              rhoxy, rhoyz, rhoxz, matx, matyz, disc, call_put,
                              scud_type, ext_lit, PREMIUM) -
              premium) /
             shift;
    return (deltaz);
    break;

    /****************************     GAMMAS 	******************************/
  case GAMMAX:
    shift = fwdx / 1000;
    gammax =
        srt_f_optscdspr(fwdx + shift, fwdy, fwdz, spoty, spotz, strike, barrier,
                        sigx, sigy, sigz, rhoxy, rhoyz, rhoxz, matx, matyz,
                        disc, call_put, scud_type, ext_lit, PREMIUM);
    gammax +=
        srt_f_optscdspr(fwdx - shift, fwdy, fwdz, spoty, spotz, strike, barrier,
                        sigx, sigy, sigz, rhoxy, rhoyz, rhoxz, matx, matyz,
                        disc, call_put, scud_type, ext_lit, PREMIUM);
    gammax -= 2 * premium;
    gammax /= shift * shift;
    return (gammax);
    break;

  case GAMMAY:
    shift = spoty / 1000;
    gammay = srt_f_optscdspr(fwdx, fwdy * (1 + shift / spoty), fwdz,
                             spoty + shift, spotz, strike, barrier, sigx, sigy,
                             sigz, rhoxy, rhoyz, rhoxz, matx, matyz, disc,
                             call_put, scud_type, ext_lit, PREMIUM);
    gammay += srt_f_optscdspr(fwdx, fwdy * (1 - shift / spoty), fwdz,
                              spoty - shift, spotz, strike, barrier, sigx, sigy,
                              sigz, rhoxy, rhoyz, rhoxz, matx, matyz, disc,
                              call_put, scud_type, ext_lit, PREMIUM);
    gammay -= 2 * premium;
    gammay /= shift * shift;
    return (gammay);
    break;

  case GAMMAZ:
    shift = spotz / 1000;
    gammaz = srt_f_optscdspr(fwdx * (1 + shift / spotz), fwdy, fwdz, spoty,
                             spotz + shift, strike, barrier, sigx, sigy, sigz,
                             rhoxy, rhoyz, rhoxz, matx, matyz, disc, call_put,
                             scud_type, ext_lit, PREMIUM);
    gammaz += srt_f_optscdspr(fwdx, fwdy, fwdz * (1 - shift / spotz), spoty,
                              spotz - shift, strike, barrier, sigx, sigy, sigz,
                              rhoxy, rhoyz, rhoxz, matx, matyz, disc, call_put,
                              scud_type, ext_lit, PREMIUM);
    gammaz -= 2 * premium;
    gammaz /= shift * shift;
    return (gammaz);
    break;

    /****************************     CROSS GAMMAS
     * ******************************/

  case GAMMAXY:
    shifty = spoty / 1000;
    shiftx = fwdx / 1000;
    gammaxy = srt_f_optscdspr(fwdx + shiftx, fwdy * (1 + shifty / spoty), fwdz,
                              spoty + shifty, spotz, strike, barrier, sigx,
                              sigy, sigz, rhoxy, rhoyz, rhoxz, matx, matyz,
                              disc, call_put, scud_type, ext_lit, PREMIUM);

    gammaxy -=
        srt_f_optscdspr(fwdx + shiftx, fwdy, fwdz, spoty, spotz, strike,
                        barrier, sigx, sigy, sigz, rhoxy, rhoyz, rhoxz, matx,
                        matyz, disc, call_put, scud_type, ext_lit, PREMIUM);

    gammaxy -= srt_f_optscdspr(fwdx, fwdy * (1 + shifty / spoty), fwdz,
                               spoty + shifty, spotz, strike, barrier, sigx,
                               sigy, sigz, rhoxy, rhoyz, rhoxz, matx, matyz,
                               disc, call_put, scud_type, ext_lit, PREMIUM);
    gammaxy += premium;
    gammaxy /= shiftx * shifty;
    return (gammaxy);
    break;

  case GAMMAXZ:
    shiftz = spotz / 1000;
    shiftx = fwdx / 1000;
    gammaxz = srt_f_optscdspr(fwdx + shiftx, fwdy, fwdz * (1 + shiftz / spotz),
                              spoty, spotz + shiftz, strike, barrier, sigx,
                              sigy, sigz, rhoxy, rhoyz, rhoxz, matx, matyz,
                              disc, call_put, scud_type, ext_lit, PREMIUM);

    gammaxz -=
        srt_f_optscdspr(fwdx + shiftx, fwdy, fwdz, spoty, spotz, strike,
                        barrier, sigx, sigy, sigz, rhoxy, rhoyz, rhoxz, matx,
                        matyz, disc, call_put, scud_type, ext_lit, PREMIUM);

    gammaxz -= srt_f_optscdspr(fwdx, fwdy, fwdz * (1 + shiftz / spotz), spoty,
                               spotz + shiftz, strike, barrier, sigx, sigy,
                               sigz, rhoxy, rhoyz, rhoxz, matx, matyz, disc,
                               call_put, scud_type, ext_lit, PREMIUM);
    gammaxz += premium;
    gammaxz /= shiftx * shiftz;
    return (gammaxz);
    break;

  case GAMMAYZ:
    shiftz = spotz / 1000;
    shifty = spoty / 1000;
    gammayz = srt_f_optscdspr(fwdx, fwdy * (1 + shifty / spoty),
                              fwdz * (1 + shiftz / spotz), spoty + shifty,
                              spotz + shiftz, strike, barrier, sigx, sigy, sigz,
                              rhoxy, rhoyz, rhoxz, matx, matyz, disc, call_put,
                              scud_type, ext_lit, PREMIUM);

    gammayz -= srt_f_optscdspr(fwdx, fwdy * (1 + shifty / spoty), fwdz,
                               spoty + shifty, spotz, strike, barrier, sigx,
                               sigy, sigz, rhoxy, rhoyz, rhoxz, matx, matyz,
                               disc, call_put, scud_type, ext_lit, PREMIUM);

    gammayz -= srt_f_optscdspr(fwdx, fwdy, fwdz * (1 + shiftz / spotz), spoty,
                               spotz + shiftz, strike, barrier, sigx, sigy,
                               sigz, rhoxy, rhoyz, rhoxz, matx, matyz, disc,
                               call_put, scud_type, ext_lit, PREMIUM);
    gammayz += premium;
    gammayz /= shifty * shiftz;
    return (gammayz);
    break;

    /****************************     VEGAS 	******************************/

  case VEGAX:
    vegax =
        (srt_f_optscdspr(fwdx, fwdy, fwdz, spoty, spotz, strike, barrier,
                         sigx + 0.01, sigy, sigz, rhoxy, rhoyz, rhoxz, matx,
                         matyz, disc, call_put, scud_type, ext_lit, PREMIUM) -
         premium);
    return (vegax);
    break;

  case VEGAY:
    vegay =
        (srt_f_optscdspr(fwdx, fwdy, fwdz, spoty, spotz, strike, barrier, sigx,
                         sigy + 0.01, sigz, rhoxy, rhoyz, rhoxz, matx, matyz,
                         disc, call_put, scud_type, ext_lit, PREMIUM) -
         premium);
    return (vegay);
    break;

  case VEGAZ:
    vegaz =
        (srt_f_optscdspr(fwdx, fwdy, fwdz, spoty, spotz, strike, barrier, sigx,
                         sigy, sigz + 0.01, rhoxy, rhoyz, rhoxz, matx, matyz,
                         disc, call_put, scud_type, ext_lit, PREMIUM) -
         premium);
    return (vegaz);
    break;

    /****************************     THETA 	******************************/

  case THETA:
    shift = 1.0 / 365.0;
    theta =
        srt_f_optscdspr(fwdx, fwdy, fwdz, spoty, spotz, strike, barrier, sigx,
                        sigy, sigz, rhoxy, rhoyz, rhoxz, matx - shift,
                        matyz - shift, disc * exp(-shift * log(disc) / matx),
                        call_put, scud_type, ext_lit, PREMIUM) -
        premium;
    return (theta);
    break;

  default:
    return premium;
    break;
  }

} /* END srt_f_optscdout() */

/******************************************************************************/
