/*
 * $Log: bssmiled.cpp,v $
 * Revision 1.46  2004/05/26 08:26:33  mab
 * New version ob Beta <> 1
 *
 * Revision 1.45  2004/04/26 12:39:54  emezzine
 * Store volatility in Storeinfocaplet in bp.
 *
 * Revision 1.44  2004/04/20 17:26:27  emezzine
 *  store the good volatility
 *
 * Revision 1.43  2004/04/07 07:18:42  rguillemot
 * Smile interpole
 *
 * Revision 1.42  2004/04/05 13:38:49  mab
 * double alpha0 = 0.0;
 * replaced by:
 * double alpha0 = sigmaATM*pow(f, 1.0-beta);
 *
 * Revision 1.41  2004/03/31 16:00:32  rguillemot
 * CMT Bug Fix
 *
 * Revision 1.40  2004/03/19 13:45:54  mab
 * Take in account SABR Sigma in Greeks calculation
 *
 * Revision 1.39  2004/03/05 17:28:16  mab
 * if  (beta) itsBeta   = (ARM_VolCurve *) beta->Clone();
 * replaced by : SetBeta(beta);
 *
 * Revision 1.38  2004/03/03 17:59:34  rguillemot
 * Caplet Bug Fix
 *
 * Revision 1.37  2004/03/03 14:51:45  emezzine
 * correct bug.
 *
 * Revision 1.36  2004/03/02 11:27:24  mab
 * Copy File Error: Corrected
 *
 * Revision 1.35  2004/03/02 10:50:08  mab
 * Improvements in computing Alpha (case Beta <> 1)
 *
 * Revision 1.33  2004/02/27 10:22:57  mab
 * in ComputeAlpha: Derivative test
 * if ((!finite(denomDeriv)) || IS_ZERO(denomDeriv))
 *
 * Revision 1.32  2004/02/26 13:46:55  mab
 * Improvements in : EuroCaplet
 *
 * Revision 1.31  2004/02/26 09:59:24  mab
 * Correction SABR formulaes (Beta <> 1)
 *
 * Revision 1.30  2004/02/20 12:54:12  rguillemot
 * Warning Fix
 *
 * Revision 1.29  2004/02/20 09:31:03  mab
 * Changes : itsAlphaOrSigmaInput in:
 * itsSigmaOrAlphaInput
 *
 * Revision 1.28  2004/02/17 16:02:48  mab
 * Added : Euro caplet
 *
 * Revision 1.27  2004/02/17 09:56:28  mab
 * Management of alpha and ATM Vol
 *
 * Revision 1.26  2004/02/04 11:27:19  mab
 *  itsSABRSig must be in 100 base
 *
 * Revision 1.25  2004/02/04 11:05:12  mab
 * Replace in ComputeVol sigma by: sigma*100
 *
 * Revision 1.24  2004/01/30 16:36:31  mab
 * volatility divided by 100
 * in case Beta <> 1
 *
 * Revision 1.23  2004/01/29 11:23:33  mab
 * Improcements
 *
 * Revision 1.22  2004/01/29 10:41:31  mab
 * Improvements for case Beta <> 1
 *
 * Revision 1.21  2004/01/29 10:23:19  mab
 * Management of case Beta <> 1
 *
 * Revision 1.20  2004/01/07 09:53:02  arm
 * ComputeFunc() : Transfered from .h
 *
 * Revision 1.19  2003/11/18 09:59:02  mab
 * Management of ARM_Digital(s)
 *
 * Revision 1.18  2003/10/01 07:51:30  sgasquet
 * Modif ComputeVol pour prendre en compte dimension de la vol ATM
 *
 * Revision 1.17  2003/09/23 18:01:47  ebenhamou
 * move computevol to cpp and corrected in the case of nu = 0
 *
 * Revision 1.16  2003/08/20 17:27:22  sgasquet
 * Modif pour prise en compte Cas Arithmetique
 *
 * Revision 1.15  2003/08/13 09:18:22  sgasquet
 * Ajout SABR avec approx arithmetique
 *
 * Revision 1.14  2003/07/01 19:09:51  mab
 * improvements
 *
 * Revision 1.13  2003/06/13 17:25:09  mab
 * improvements: more control tests
 *
 * Revision 1.12  2003/06/02 08:11:42  mab
 * Changes in the code in order to manage
 * correctly already fitted parameters
 *
 * Revision 1.11  2003/03/13 17:10:31  mab
 * int k = 0; replaced by : int k = 1; // We begin at 1
 *
 * Revision 1.10  2003/03/07 13:17:02  mab
 * Version with conditional Smoothing
 *
 * Revision 1.9  2003/02/14 15:48:01  mab
 * Improvements
 *
 * Revision 1.8  2003/02/14 09:58:14  mab
 * Integrating bootstrapping calibration
 *
 * Revision 1.7  2003/02/11 16:53:00  mab
 * Improvements
 *
 * Revision 1.6  2003/01/17 10:06:26  mab
 * Initial version
 *
 * Revision 1.5  2002/12/17 13:54:10  mab
 * This version integrate the calibration
 *
 * Revision 1.4  2002/12/11 11:00:57  mab
 * Generalisation by using : ARM_VolCurve
 *
 * Revision 1.3  2002/10/09 09:37:58  mab
 * Take in account SABR modelisation
 *
 * Revision 1.2  2002/03/28 13:20:19  mab
 * Alpha replaced by Rho's
 *
 * Revision 1.1  2002/02/25 13:53:53  mab
 * Initial revision
 *
 */




/*----------------------------------------------------------------------------*
    bssmiled.cpp
 
    This file implements the ARM_BSSmiledModel class

    (The shifted Black-Sholes Model)

*----------------------------------------------------------------------------*/


#include "firsttoinc.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "util.h"
#include "armglob.h"
#include "bssmiled.h"
#include "volflat.h"
#include "zeroflat.h"
#include "armdigital.h"

#include "fromto.h"

#include "cmsleg.h"

// SABR Formulaes

#include "gpclosedforms/sabrbdiff1.h"


/*---------------------------------------------------------------------------*/


#define PF_MATU_SIZE 3000






/*---------------------------------------------------------------------------*/
/*   Constructor : this is the first "intuitive" constructor corresponding   */
/*                       implicitly to ONE tenor                             */ 
/*---------------------------------------------------------------------------*/

ARM_BSSmiledModel::ARM_BSSmiledModel(ARM_Date& startDate,
                                     double spot,
                                     ARM_ZeroCurve* dividend, 
                                     ARM_ZeroCurve* discountRate, 
                                     ARM_VolCurve*  volatility,
                                     int type,
                                     ARM_Vector* matu,
                                     ARM_Vector* rho,
                                     ARM_Vector* nu,
                                     int SABR,
                                     ARM_CorrelManager* correlManager,
									 ARM_ZeroCurve* realDiscountRate,
									 ARM_VolCurve* AdjCvxVol,
									 int ConvToAdjVolWhenAlpha)
                                     :ARM_BSModel(startDate,
                                                  spot,
                                                  dividend,
                                                  discountRate,
                                                  volatility, type,
												  realDiscountRate)

{
    Init();

    // Initialize the Correlation Manager

    bool isShared = true;

    SetCorrelManager(correlManager, isShared);
                                  

    itsSmiledModelType = SABR;


    ARM_Vector* yearTermsRho = (ARM_Vector *) matu->Clone();
    ARM_Vector* yearTermsNu  = (ARM_Vector *) matu->Clone();

    ARM_Matrix* theRho       = new ARM_Matrix(*rho);

    ARM_Vector* undMatuRho   = new ARM_Vector(1);
    ARM_Vector* undMatuNu    = new ARM_Vector(1);

    ARM_Matrix* theNu        = new ARM_Matrix(*nu);

    itsRho = new ARM_VolLInterpol(startDate, 
                                  yearTermsRho,
                                  undMatuRho,
                                  theRho);

    itsNu  = new ARM_VolLInterpol(startDate, 
                                  yearTermsNu, 
                                  undMatuNu,
                                  theNu);


    itsLogRho = (ARM_VolCurve *) itsRho->Clone();
    itsLogNu  = (ARM_VolCurve *) itsNu->Clone();

	if (itsSmiledModelType == K_LD)
		UpdateLogRhoNu();

    SaveInitialCalibData();

	// Set the  volatility curve of convexity adjustments
	//---------------------------------------------------------------------

	if ( AdjCvxVol )
	{
		SetCvxAdjVolatility(AdjCvxVol);
	}
	else
	{
		ARM_VolCurve* tmpSigmaCurve = NULL;

		// 2 possible cases if AdjCvxVol == NULL:
		// -> Case Alpha: convert Alpha to Sigma and use the resulting curve 
		//			   for convexity adjustments
		// -> Case Sigma: use the vol ATM for convexity adjustments
		//---------------------------------------------------------------------
		if ( !itsSigmaOrAlphaInput )
		{
			if ( ConvToAdjVolWhenAlpha )
			{
				tmpSigmaCurve = FromAlphaToSigma(volatility);
			}
		}
		else
		{
			tmpSigmaCurve = (ARM_VolCurve*) volatility->Clone();
		}

		SetCvxAdjVolatility(tmpSigmaCurve);

		if ( tmpSigmaCurve )
		{
			delete tmpSigmaCurve;
			tmpSigmaCurve = NULL;
		}
	}
}



/*---------------------------------------------------------------------------*/
/*    This is a more general constructor without the preceding assumption    */ 
/*---------------------------------------------------------------------------*/

ARM_BSSmiledModel::ARM_BSSmiledModel(ARM_Date& startDate, double spot,
                                     ARM_ZeroCurve* dividend,
                                     ARM_ZeroCurve* discountRate,		// forecast curve
                                     ARM_VolCurve* volatility,
                                     int volType,
                                     ARM_VolCurve* rho,
                                     ARM_VolCurve* nu,
                                     int SABR,
                                     ARM_VolCurve* beta,
                                     double SabrWeight,
                                     int sigmaOrAlphaInput,
                                     ARM_CorrelManager* correlManager,
									 ARM_ZeroCurve* realDiscountRate,	// discount curve
									 ARM_VolCurve* AdjCvxVol,
									 int ConvToAdjVolWhenAlpha)	
                                     : ARM_BSModel(startDate,
                                                   spot,
                                                   dividend,
                                                   discountRate,
                                                   volatility, volType,
												   realDiscountRate)
{
    Init();

    // Initialize the Correlation Manager

    bool isShared = true;

    SetCorrelManager(correlManager, isShared);


    if (( SabrWeight > 1.0 )
        ||
        ( SabrWeight < 0.0 )
       )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
       "ARM_BSSmiledModel: SABR Weight must be in [0, 1]");
    } 
   
    itsSmiledModelType    = SABR;

    itsRho    = (ARM_VolCurve *) rho->Clone();

    itsNu     = (ARM_VolCurve *) nu->Clone();

    SetBeta(beta);

    itsLogRho = (ARM_VolCurve *) itsRho->Clone();
    itsLogNu  = (ARM_VolCurve *) itsRho->Clone();

    itsSABRWeight = SabrWeight;

    itsSigmaOrAlphaInput = sigmaOrAlphaInput;

	if (itsSmiledModelType == K_LD)
		UpdateLogRhoNu();

    SaveInitialCalibData();

	// Set the  volatility curve of convexity adjustments
	//---------------------------------------------------------------------

	if ( AdjCvxVol )
	{
		SetCvxAdjVolatility(AdjCvxVol);
	}
	else
	{
		ARM_VolCurve* tmpSigmaCurve = NULL;

		// 2 possible cases if AdjCvxVol == NULL:
		// -> Case Alpha: convert Alpha to Sigma and use the resulting curve 
		//			   for convexity adjustments
		// -> Case Sigma: use the vol ATM for convexity adjustments
		//---------------------------------------------------------------------
		if ( !itsSigmaOrAlphaInput )
		{
			if ( ConvToAdjVolWhenAlpha )
			{
				tmpSigmaCurve = FromAlphaToSigma(volatility);
			}
		}
		else
		{
			tmpSigmaCurve = (ARM_VolCurve*) volatility->Clone();
		}

		SetCvxAdjVolatility(tmpSigmaCurve);

		if ( tmpSigmaCurve )
		{
			delete tmpSigmaCurve;
			tmpSigmaCurve = NULL;
		}
	}
}


/*---------------------------------------------------------------------------*/
/*    General constructor with SABRVol structure						     */ 
/*---------------------------------------------------------------------------*/

ARM_BSSmiledModel::ARM_BSSmiledModel(ARM_Date& startDate, double spot,
                                     ARM_ZeroCurve* dividend,
                                     ARM_ZeroCurve* discountRate,		// forecast curve
									 ARM_SABRVol*	SABRVol,
									 int volType,
                                     ARM_CorrelManager* correlManager,
									 ARM_ZeroCurve* realDiscountRate,	// discount curve
									 ARM_VolCurve* AdjCvxVol,
									 int ConvToAdjVolWhenAlpha)
                                     : ARM_BSModel(startDate,
                                                   spot,
                                                   dividend,
                                                   discountRate,
                                                   (ARM_VolCurve *) SABRVol->GetSigmaOrAlpha(),
												   volType,
												   realDiscountRate)
{
    Init();

    // Initialize the Correlation Manager

    bool isShared = true;

	double SabrWeight = SABRVol->GetWeight();

    SetCorrelManager(correlManager, isShared);


    if (( SabrWeight > 1.0 )
        ||
        ( SabrWeight < 0.0 )
       )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
       "ARM_BSSmiledModel: SABR Weight must be in [0, 1]");
    } 
   
    itsSmiledModelType    = SABRVol->GetModelType();

    itsRho    = (ARM_VolCurve *) SABRVol->GetRho()->Clone();

    itsNu     = (ARM_VolCurve *) SABRVol->GetNu() ->Clone();

    SetBeta((ARM_VolCurve *) SABRVol->GetBeta());

    itsLogRho = (ARM_VolCurve *) itsRho->Clone();
    itsLogNu  = (ARM_VolCurve *) itsRho->Clone();

    itsSABRWeight = SabrWeight;

    itsSigmaOrAlphaInput = SABRVol->GetSigmaOrAlphaFlag();

	if (itsSmiledModelType == K_LD)
		UpdateLogRhoNu();

    SaveInitialCalibData();

	// Set the  volatility curve of convexity adjustments
	//---------------------------------------------------------------------

	if ( AdjCvxVol )
	{
		SetCvxAdjVolatility(AdjCvxVol);
	}
	else
	{
		ARM_VolCurve* tmpSigmaCurve = NULL;

		// 2 possible cases if AdjCvxVol == NULL:
		// -> Case Alpha: convert Alpha to Sigma and use the resulting curve 
		//				  for convexity adjustments
		// -> Case Sigma: use the Sigma curve for convexity adjustments
		//---------------------------------------------------------------------
		if ( !itsSigmaOrAlphaInput )
		{
			if ( ConvToAdjVolWhenAlpha )
			{
				tmpSigmaCurve = FromAlphaToSigma( (ARM_VolCurve *) SABRVol->GetSigmaOrAlpha() );
			}
		}
		else
		{
			tmpSigmaCurve = (ARM_VolCurve *) SABRVol->GetSigmaOrAlpha()->Clone();
		}

		SetCvxAdjVolatility(tmpSigmaCurve);

		if ( tmpSigmaCurve )
		{
			delete tmpSigmaCurve;
			tmpSigmaCurve = NULL;
		}
	}
}

/*----------------------------------------------------------------------------*
     Flat inputs constructor
*----------------------------------------------------------------------------*/

ARM_BSSmiledModel::ARM_BSSmiledModel(ARM_Date& startDate, double spot, 
                                     double flatdiv,
                                     double flatdiscrate,
                                     double flatvol, int type,
                                     ARM_Vector* matu,
                                     ARM_Vector* rho,
                                     ARM_Vector* nu,
                                     int SABR,
                                     ARM_CorrelManager* correlManager)
                                     :ARM_BSModel(startDate,
                                                  spot,
                                                  flatdiv,
                                                  flatdiscrate,
                                                  flatvol,
                                                  type) 
{
    Init();

    // Initialize the Correlation Manager

    bool isShared = true;

    SetCorrelManager(correlManager, isShared);


    // TMP SetConvAdjustManager(NULL, false);

    itsSmiledModelType = SABR;

    ARM_Vector* yearTermsRho = (ARM_Vector *) matu->Clone();
    ARM_Vector* yearTermsNu  = (ARM_Vector *) matu->Clone();
    ARM_Matrix* theRho       = new ARM_Matrix(*rho);

    ARM_Vector* undMatuRho   = new ARM_Vector(1);
    ARM_Vector* undMatuNu    = new ARM_Vector(1);
    ARM_Matrix* theNu        = new ARM_Matrix(*nu);

    itsRho = new ARM_VolLInterpol(startDate,
                                  yearTermsRho,
                                  undMatuRho,
                                  theRho);

    itsNu  = new ARM_VolLInterpol(startDate,
                                  yearTermsNu,
                                  undMatuNu,
                                  theNu);

    itsLogRho = (ARM_VolCurve *) itsRho->Clone();
    itsLogNu  = (ARM_VolCurve *) itsNu->Clone();

	if (itsSmiledModelType == K_LD)
		UpdateLogRhoNu();

    SaveInitialCalibData();
}



/*----------------------------------------------------------------------------*
    Constructor    (copy).
*----------------------------------------------------------------------------*/

ARM_BSSmiledModel::ARM_BSSmiledModel(const ARM_BSSmiledModel& bs) 
                  : ARM_BSModel(bs)
{
    Init();

    BitwiseCopy(&bs);
}



/*----------------------------------------------------------------------------*
    Assignment operator
*----------------------------------------------------------------------------*/

ARM_BSSmiledModel& ARM_BSSmiledModel::operator = (const ARM_BSSmiledModel& bs)
{
    (*this).ARM_BSModel::operator = (bs);

    BitwiseCopy(&bs);

    return(*this);
}



double z(double K, double f, double alpha, double beta)
{
    double _z;

    _z = 1.0/(alpha*(1.0-beta))*(pow(f, 1.0-beta)-pow(K, 1.0-beta));

    return(_z);
}



double x(double K, double f, double alpha, double beta, 
         double rho, double nu)
{
    double _z;

    _z = z(K, f, alpha, beta);


    double _x;

    _x = 1.0/nu*log((sqrt(1.0-2.0*rho*nu*_z+SQR(nu)*SQR(_z))-rho+nu*_z)
             /(1.0-rho));

    return(_x);
}


// Hagan derivative
double dsigBSappATM(double f, double K,
                    double mat,
                    double alpha, double rho, double nu, double beta)
{
   
    double dsig = 1.0/(pow(f*K, (1.0-beta)/2.0)*(1.0+SQR((1.0-beta)*log(f/K))/24.0
		          +SQR(SQR((1.0-beta)*log(f/K)))/1920.0))
                  *(1.0+(SQR((1.0-beta)*alpha)/(24*pow(f*K, 1.0-beta))+
                  rho*beta*nu*alpha/(4*pow(f*K,(1.0-beta)/2.0)) 
                  +(2.0-3.0*SQR(rho))*SQR(nu)/24.0)*mat)
				  +alpha/(pow(f*K, (1.0-beta)/2.0)*(1.0+SQR((1.0-beta)*log(f/K))/24.0
                  +SQR(SQR((1.0-beta)*log(f/K)))/1920.0))
                  *(SQR(1.0-beta)*alpha/(12.0*pow(f*K, 1.0-beta))+
                  rho*beta*nu/(4.0*pow(f*K, (1.0-beta)/2.0)))*mat;

    return(dsig);
}


// Impl LN 

double dsigBS2(double f, double K,
               double mat,
               double alpha, double rho, double nu, double beta)
{
    double dsig = 1.0/(pow(f*K,0.5*(1-beta))
			     *sqrt(1.0+2.0*mat*(-1/24*SQR(alpha
				 *(beta-1.0))/pow(f, 2.0*(1.0-beta))-1.0/4.0
				 *alpha*beta*rho*nu/pow(f, 1.0-beta)
				 +SQR(nu)*(-1.0/12.0+1.0/8.0*SQR(rho)))))
                 -0.5*alpha/(pow(f*K, 0.5*(1.0-beta))
				 *pow(1.0+2.0*mat*(-1/24*SQR(alpha
				 *(beta-1.0))/pow(f, 2.0*(1.0-beta))-1.0/4.0
				 *alpha*beta*rho*nu/pow(f, 1.0-beta)
				 +SQR(nu)*(-1.0/12.0+1.0/8.0*SQR(rho))), 1.5))
				 *(2.0*mat*(-1.0/12.0*alpha*SQR(beta-1.0)
				  /pow(f, 2.0*(1.0-beta))-1.0/4.0*
				  beta*rho*nu/pow(f, 1.0-beta)));

    return(dsig);
}



double ComputeSigATMBetaEqOne(double f, double K,
                              double maturity, 
                              double alpha,
                              double Rho_t, double Nu_t,
                              int type)
{
    double spot   = f;
    double strike = K;

    double resSigATM = 0.0;


    // Beta == 1

    if ( type != K_LD ) // SABR type
    {       
       if (IS_ONE(Rho_t))
       {
          Rho_t = double(0.99);
       }
       
       if (IS_ZERO(Rho_t))
       {
          resSigATM = alpha*(1.0+(SQR(Nu_t)*maturity/12.0));
         
          return(resSigATM);
       }
       else
       {
          double rho2_3_24;
            
          rho2_3_24 = ((2.0-3.0*SQR(Rho_t))/24.0)*SQR(Nu_t)
                       *maturity;

          resSigATM = (SQR(0.5*alpha*Rho_t*Nu_t*maturity+(1+rho2_3_24))-SQR(1+rho2_3_24))
                       /(Rho_t*Nu_t*maturity);
                      

          return(resSigATM);
       }  
    }
    
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                   "ComputeAlphaBetaEqOne() unavailable for LD type");
    
    return(0.0);
}



double ARM_BSSmiledModel::CalcSABRVolBetaDiffOne(double maturity, 
                                                 double volTenor,
                                                 double fwd, 
                                                 double strike)
{
    double SABRVol = 0.0;


    switch(itsSmiledModelType)
    {
        case K_SABR_ARITH :
        {
            SABRVol = CalcSABRVolBetaDiffOneArith(maturity,
                                                  volTenor,
                                                  fwd,
                                                  strike);

            return(SABRVol);
        };
        break;

        case K_SABR_IMPLNVOL :
        {
            SABRVol = CalcSABRVolBetaDiffOneImpLNVol(maturity,
                                                     volTenor,
                                                     fwd,
                                                     strike);

            return(SABRVol);

        };
        break;

        case K_SABR_WEIGHT :
        {
            SABRVol = CalcSABRVolBetaDiffOneWeighted(maturity,
                                                     volTenor,
                                                     fwd,
                                                     strike);

            return(SABRVol);

        };
        break;

        default :
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
       "ARM_BSSmiledModel::CalcSABRVolBetaDiffOne(): Unrecognized smile model type!?"); 
        }
        break;
    };

    return(SABRVol);
}



double ARM_BSSmiledModel::CalcSABRVolBetaDiffOneWeighted(double maturity, 
                                                         double volTenor,
                                                         double fwd, 
                                                         double strike)
{
    double sig1 = CalcSABRVolBetaDiffOneArith(maturity,
                                              volTenor,
                                              fwd,
                                              strike);

    double sig2 = CalcSABRVolBetaDiffOneImpLNVol(maturity,
                                                 volTenor,
                                                 fwd,
                                                 strike);

    double w = itsSABRWeight;

    double SABRVol = w*sig1+(1.0-w)*sig2;

    return(SABRVol);
}



double b0(double K, double f, double beta)
{
    double b0Res =  pow(0.5*(f+K), beta);

    return(b0Res);
}



double b1(double K, double f, double beta)
{
    double b1Res = beta*pow(0.5*(f+K), 2.0*beta-1.0);

    return(b1Res);
}



double b2(double K, double f, double beta)
{
    double b2Res = beta*(2.0*beta-1.0)*pow(0.5*(f+K), 3.0*beta-2.0);

    return(b2Res);
}



double theta(double K, double f, double alpha, double beta, 
             double rho, double nu)
{
    double _x;
    double _z;
    double _theta;

    _x = x(K, f, alpha, beta, rho, nu);

    _z = z(K, f, alpha, beta);

    _theta = log(alpha*_z*pow(K*f, beta/2.0)/(f-K))
              +log(_x*pow(1.0-2.0*rho*nu*_z+SQR(nu*_z), 0.25)/_z)
              +0.25*rho*nu*alpha*b1(K, f, beta)/b0(K, f, beta)*SQR(_z);

    return(_theta);
}



/*----------------------------------------------------------------------------*
    Implied Log Normal Vol case : Correction to the Hagan model
*----------------------------------------------------------------------------*/
double ARM_BSSmiledModel::CalcSABRVolBetaDiffOneImpLNVol(double maturity, 
                                                         double volTenor,
                                                         double fwd, 
                                                         double strike)
{
    double K = strike*0.01;

    double f = fwd*0.01;

    double SABRVol = 0.0;

    double volatility = GetVolatility()->ComputeVolatility(maturity,
                                                           volTenor);

    double alpha;


    // Retrieve the data corresponding to the tenor

    double rho = itsRho->ComputeVolatility(maturity, volTenor);

    double nu  = itsNu->ComputeVolatility(maturity, volTenor);

    double beta = itsBeta->ComputeVolatility(maturity, volTenor);


    try
    {
        if (itsSigmaOrAlphaInput) // The input is ATM Vol
        {
           alpha = ComputeAlpha(f, K, maturity, volatility/100.0, 
                                rho, nu, beta, itsSABRWeight, K_SABR_IMPLNVOL);
        }
        else
        {
           alpha = volatility/100.0;
        }
    }

    catch(Exception& anExpt)
    {
        throw anExpt;

        return(0.0);
    }

    double psi;
    double sig;
    double _x;
    double _theta;


    if ( fabs(K-f) < 0.000001*f )
    {
       sig = alpha/(pow(f*K, 0.5*(1-beta))
            *sqrt(1.0+2.0*maturity*(-1.0/24.0*pow(alpha*(beta-1.0), 2.0)
            /pow(f, 2.0*(1.0-beta))-1.0/4.0*alpha*beta*rho*nu/pow(f,1.0-beta)
             +nu*nu*(-1.0/12.0+1.0/8.0*rho*rho))));

    }
    else
    {
       psi = log(log(f/K)*sqrt(f*K)/(f-K));

       _x = x(K, f, alpha, beta, rho, nu);

       _theta = theta(K, f, alpha, beta, rho, nu);

       sig = log(f/K)/(_x*sqrt(1.0-2.0*_theta*maturity/SQR(_x)+2.0*psi
             *maturity/SQR(_x)));
    }

    SABRVol = sig;

    return(SABRVol);
}



double ARM_BSSmiledModel::ComputeAlpha(double f, double K,
                    double mat,
                    double sigmaATM,
                    double rho, double nu, double beta, 
                    double weight, int type)
{
	double alpha = ::ComputeAlpha(f, K, mat, sigmaATM, rho, nu, beta, weight, type);

	return alpha;
}

/*----------------------------------------------------------------------------*
    Hagan method case
*----------------------------------------------------------------------------*/
double ARM_BSSmiledModel::CalcSABRVolBetaDiffOneArith(double maturity, 
                                                      double volTenor,
                                                      double fwd, 
                                                      double strike)
{
    double volatility = GetVolatility()->ComputeVolatility(maturity,
                                                           volTenor);

    // Retrieve the data corresponding to the tenor

    double rho = itsRho->ComputeVolatility(maturity,
                                           volTenor);

    double nu  = itsNu->ComputeVolatility(maturity,
                                          volTenor);
  
    double beta = itsBeta->ComputeVolatility(maturity,
                                             volTenor);

    double alpha;

    double f = fwd*0.01;
    double K = strike*0.01;

    try
    {
        if (itsSigmaOrAlphaInput) // The input is ATM Vol
        {
           alpha = ComputeAlpha(f, K, maturity, volatility/100.0,
                                rho, nu, beta, itsSABRWeight, K_SABR_ARITH);
        }
        else
        {
           alpha = volatility/100.0;
        }
    }

    catch(Exception& theExpt)
    {
        throw theExpt;

        return(0.0);
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
             "ComputeAlpha() call in : ARM_BSSmiledModel::CalcSABRVolBetaDiffOneArith: \
              failed!?");  
    }

    double SABRVol = 0.0;

    double ratio;


    if ( fabs(f-K) < (0.001*f) )
    {
       ratio = 1.0;
    }
    else
    {
       ratio = z(K, f, alpha, beta)/x(K, f, alpha, beta, rho, nu);
    }

    double oneMoinsBeta = 1.0-beta;

    double SQROneMoinsBeta = oneMoinsBeta*oneMoinsBeta;

    double sig = alpha/(pow(f*K, oneMoinsBeta/2.0)*(1+SQROneMoinsBeta
                 *SQR(log(f/K))/24.0+SQR(SQROneMoinsBeta)
                 *SQR(log(f/K))*SQR(log(f/K))/1920.0))
                 *(ratio)
                 *(1.0+(SQR(oneMoinsBeta)*SQR(alpha)/(24*pow(f*K, oneMoinsBeta))
                 +rho*beta*nu*alpha/(4*pow(f*K, oneMoinsBeta/2.0))
                 +(2.0-3.0*SQR(rho))*SQR(nu)/24)*maturity);

    SABRVol = sig;

    return(SABRVol);
}



double ARM_BSSmiledModel::bsOptionBetaDiffOne(double spot,
                                              double strike,
                                              double volatility,
                                              double dividend,
                                              double discountRate,
                                              double maturity,
                                              double CallPut,
                                              double volTenor)
{
    double optionPrice = 0.0;

    double SABRSigma;

    SABRSigma = CalcSABRVolBetaDiffOne(maturity,
                                       volTenor,
                                       spot,
                                       strike);

    optionPrice = ::bsOption(spot, strike, SABRSigma,
                             dividend, discountRate,
                             maturity, CallPut);

    itsSABRSig = SABRSigma*100.0;

    return(optionPrice);
}



double ARM_BSSmiledModel::CptSABRVolBetaDiffOneDirectExact(double maturity,
                                                           double volTenor,
                                                           double fwd,
                                                           double strike)
{
    double K = strike*0.01;

    double f = fwd*0.01;

    double SABRVol = 0.0;

    double volatility = GetVolatility()->ComputeVolatility(maturity,
                                                           volTenor);

    double alpha;


    // Retrieve the data corresponding to the tenor

    double rho = itsRho->ComputeVolatility(maturity, volTenor);

    double nu  = itsNu->ComputeVolatility(maturity, volTenor);

    double beta = itsBeta->ComputeVolatility(maturity, volTenor);


    try
    {
        if (itsSigmaOrAlphaInput) // The input is ATM Vol
        {
           alpha = ComputeAlpha(f, K, maturity, volatility/100.0,
                                rho, nu, beta, itsSABRWeight, K_SABR_IMPLNVOL);
        }
        else
        {
           alpha = volatility/100.0;
        }
    }

    catch(Exception& except)
    {
        throw except;

        return(0.0);
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
             "ComputeAlpha() call in : ARM_BSSmiledModel::CptSABRVolBetaDiffOneDirectExact: \
              failed!?");  
    }

    double sig;

    sig = CptSABR_implicit_vol_direct(f, K, maturity, alpha, beta,
                                      rho, nu,
                                      K_SABR_IMPLNVOL); // TMP

    SABRVol = sig;

    return(SABRVol);
}



double ARM_BSSmiledModel::CptSABRVolBetaDiffOneNormalExact(double maturity,
                                                           double volTenor,
                                                           double fwd,
                                                           double strike)
{
    double K = strike*0.01;

    double f = fwd*0.01;

    double SABRVol = 0.0;

    double volatility = GetVolatility()->ComputeVolatility(maturity,
                                                           volTenor);

    double alpha;


    // Retrieve the data corresponding to the tenor

    double rho = itsRho->ComputeVolatility(maturity, volTenor);

    double nu  = itsNu->ComputeVolatility(maturity, volTenor);

    double beta = itsBeta->ComputeVolatility(maturity, volTenor);


    try
    {
        if (itsSigmaOrAlphaInput) // The input is ATM Vol
        {
           alpha = ComputeAlpha(f, K, maturity, volatility/100.0,
                                rho, nu, beta, itsSABRWeight, K_SABR_GEO);
        }
        else
        {
           alpha = volatility/100.0;
        }
    }

    catch(Exception& armExpt)
    {
        throw armExpt;

        return(0.0);
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
             "ComputeAlpha() call in : ARM_BSSmiledModel::CptSABRVolBetaDiffOneNormalExact: \
              failed!?");  
    }
    
    double sig;

    sig = CptSABR_implicit_vol_normal(f, K, maturity, alpha, beta,
                                      rho, nu,
                                      K_SABR_GEO); // TMP 

    SABRVol = sig;

    return(SABRVol);
}



double ARM_BSSmiledModel::CptSABRVolBetaDiffOneNormalOrDirectExact(double maturity,
                                                                   double volTenor,
                                                                   double fwd,
                                                                   double strike)
{
    double SABRVol = 0.0;


    switch(itsSmiledModelType)
    {
		case K_SABR_GEO : 
        {
            SABRVol = CptSABRVolBetaDiffOneNormalExact(maturity,
                                                       volTenor,
                                                       fwd,
                                                       strike);
            return(SABRVol);
        };
        break;

		case K_SABR_ARITH:
		{
			SABRVol = CalcSABRVolBetaDiffOneArith(maturity,
                                                  volTenor,
                                                  fwd,
                                                  strike);
            return(SABRVol);
        };
        break;

        case K_SABR_IMPLNVOL :
        {
            // case "Direct exact"

            SABRVol = CptSABRVolBetaDiffOneDirectExact(maturity,
                                                       volTenor,
                                                       fwd,
                                                       strike);


            return(SABRVol);

        };
        break;

        case K_SABR_WEIGHT :
        {
            double sig1  = CptSABRVolBetaDiffOneNormalExact(maturity,
                                                            volTenor,
                                                            fwd,
                                                            strike);

            double sig2 = CptSABRVolBetaDiffOneDirectExact(maturity,
                                                           volTenor,
                                                           fwd,
                                                           strike);

            double w = itsSABRWeight;

            SABRVol = w*sig1+(1.0-w)*sig2;

            return(SABRVol);

        };
        break;

        default :
        {
            throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
             "ARM_BSSmiledModel::CptSABRVolBetaDiffOneNormalOrDirectExact(): (Beta<>1) \
                          Unrecognized smile model type!?");  
        }
        break;
    };

    return(SABRVol);
}



double ARM_BSSmiledModel::bsOptionBetaDiffOneNormalOrDirectExact(double spot,
                                                                 double strike,
                                                                 double volatility,
                                                                 double dividend,
                                                                 double discountRate,
                                                                 double maturity,
                                                                 double CallPut,
                                                                 double volTenor)
{
    double optionPrice = 0.0;

    double SABRSigma;


    SABRSigma = CptSABRVolBetaDiffOneNormalOrDirectExact(maturity,
                                                         volTenor,
                                                         spot,
                                                         strike);


    optionPrice = ::bsOption(spot, strike, SABRSigma,
                             dividend, discountRate,
                             maturity, CallPut);

    itsSABRSig = SABRSigma*100.0;

    return(optionPrice);
}



void ARM_BSSmiledModel::ComputeFunc(double* params, double* FUNC,
                                    ARM_Model* eventualPricingModel)
{
    if (gTrace)
       fprintf(stdout, "\n ====> IN : ARM_BSSmiledModel::ComputeFunc \n");

    ARM_Model* pricingModel = ( eventualPricingModel ? eventualPricingModel :
                                itsFatherBSModel );

    ARM_Model::ComputeFunc(params, FUNC, pricingModel);


    // Smooth now if necessary

    if (( itsCalibMatus->GetSize() >= 3 ) && (itsSmoothFlag))
    {
       ARM_Vector* smoothFuncValue = SmoothFunction();

       double NP = 1.0;

       int sizePf = itsCalibrationPortfolio->GetSize();

       int index = 0;

       for (int i = 0; i < sizePf; i++)
       {
           FUNC[i] *= 1.0/(sqrt(NP));

           index++;
       }

       int sz = smoothFuncValue->GetSize();

       for (int j = 0; j < sz; j++)
       {
           FUNC[index] = smoothFuncValue->Elt(j);

           index++;
       }

       if (smoothFuncValue)
          delete smoothFuncValue;
    }

    if (gTrace)
       fprintf(stdout, "\n <==== OUT : ARM_BSSmiledModel::ComputeFunc \n");
}



int AssetMatuAlreadyExists(double matu, double* portfolioMatu, int pfSize)
{
    int i;

    int found = 0;

    i = 0;

    while (( i < pfSize ) && (!(found)))
    {
        if (fabs(portfolioMatu[i]-matu) <= 1e-4 )
        {
           found = 1; 
        }
        else
        {
           i++;
        }
    }

    if (( i == 0 ) && (!(found)))
    {
       return(0);
    }
    else
    {
       return(found);
    }
}



ARM_Vector* ARM_BSSmiledModel::GenCalibMatus(ARM_Portfolio* pf,
                                             double H_TimeStep)
{
    ARM_Date curAssetMatu;
    ARM_Date AsOfDate = GetStartDate();

    ARM_Vector* generatedMatus = NULL;

   

    double portfolioMatu[PF_MATU_SIZE];


    if ( H_TimeStep <= 0 )
    {
       int size = pf->GetSize();

       int pfSize = 0;
       int i;

       for (i = 0; i < size; i++)
       {
           ARM_Security* curAsset = pf->GetAsset(i);

           if (( curAsset->GetName() == ARM_CAPFLOOR ) 
               || 
               ( curAsset->GetName() == ARM_SWAPTION ) 
               ||
               ( curAsset->GetName() == ARM_DIGITAL )
              )
           {
              curAssetMatu = curAsset->GetExpiryDate();
     
              double matuTerm = (curAssetMatu-AsOfDate)/K_YEAR_LEN; 

              if (!(AssetMatuAlreadyExists(matuTerm,
                                           portfolioMatu, pfSize)))
              {
                 portfolioMatu[pfSize] = matuTerm;

                 pfSize++;
              }
           }
           else // Error
           {
              throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
      "Only CAPS/FLOORS, DIGITALS or SWAPTIONS expected for BS smiled calibration");
           }
       }

       generatedMatus = new ARM_Vector(pfSize, portfolioMatu);

       return(generatedMatus);
    }

    // H_TimeStep > 0
    // Generation of the lacking maturities
        
    int size = pf->GetSize();

    int i; 
    int pfSize = 0;

    for (i = 0; i < size; i++)
    {
        ARM_Security* curAsset = pf->GetAsset(i);

        if (( curAsset->GetName() == ARM_CAPFLOOR ) 
            || 
            ( curAsset->GetName() == ARM_SWAPTION )
            ||
            ( curAsset->GetName() == ARM_DIGITAL )
           )
        {
           curAssetMatu = curAsset->GetExpiryDate();

           double matuTerm = (curAssetMatu-AsOfDate)/K_YEAR_LEN;

           if (!(AssetMatuAlreadyExists(matuTerm,
                                        portfolioMatu, pfSize)))
           {
              portfolioMatu[pfSize] = matuTerm; 

              pfSize++;
           }
        }
        else // Error
        {
           throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
       "Only CAPS/FLOORS, DIGITALS or SWAPTIONS expected for BS smiled calibration");
        }
    }

    double tmpGenMatus[PF_MATU_SIZE];
    int    genMatusSize = 0;


    // First Tenor
    double firstGenMatu = portfolioMatu[0]-itsCalibrationTenor;

    double lowerBoundMatu;

    int k = 0;

    for (i = 0; i < pfSize; i++)
    {
        if ( i == 0 )
        {
           lowerBoundMatu = firstGenMatu;
        }
        else
        {
           lowerBoundMatu = portfolioMatu[i-1];
        }


        double pfMatu = portfolioMatu[i];

        double newGenMatu = lowerBoundMatu;

        int j = 1;

        while ( newGenMatu < pfMatu )
        {
            tmpGenMatus[k] = newGenMatu;

            newGenMatu = lowerBoundMatu+(H_TimeStep*j);

            j++;

            k++;                  
        }
    }

    // Manage last interval

    tmpGenMatus[k] = portfolioMatu[pfSize-1];
    genMatusSize = k+1;

    ARM_Vector* newMatus = new ARM_Vector(genMatusSize,
                                          tmpGenMatus);

    return(newMatus);
}



void ARM_BSSmiledModel::UpdateSigRhoNuByInterpol(double H_TimeStep,
                                                 int FitVolOrNot,
                                                 int FitRhoOrNot,
                                                 int FitNuOrNot,
                                                 ARM_Vector* generatedMatus,
                                                 int interpolMeth,
                                                 int FitBetaOrNot)
{
    int genMatusSize = generatedMatus->GetSize();


    // SIGMA

    if (FitVolOrNot)
    {
       ARM_Matrix* theATMVol    = new ARM_Matrix(genMatusSize, 1);
       ARM_Vector* yearTermsATM = (ARM_Vector *) generatedMatus->Clone();
       ARM_Vector* undMatuATM   = new ARM_Vector(1, itsCalibrationTenor);

       ARM_VolCurve* originalVol = GetVolatility();

       for (int i = 0; i < genMatusSize; i++)
       {
           if (( itsCalibMeth == ARM_CALIB_BOOTSTRAP )
               &&
               ( i > itsPreviousMaxCalibIndex )
              )
           {
              theATMVol->Elt(i, 0) =
                  itsInitialSig->ComputeVolatility(generatedMatus->Elt(i),
                                                   itsCalibrationTenor);
           }
           else
           {
              theATMVol->Elt(i, 0) =
                  originalVol->ComputeVolatility(generatedMatus->Elt(i),
                                                 itsCalibrationTenor);
           }
       }

       ARM_VolLInterpol* newVol = new ARM_VolLInterpol(GetStartDate(),
                                                       yearTermsATM,
                                                       undMatuATM,
                                                       theATMVol);
       SetVolatility(newVol);
       delete newVol;
    }

    // RHO

    if (FitRhoOrNot)
    {
       ARM_Matrix* theRho       = new ARM_Matrix(genMatusSize, 1);
       ARM_Vector* yearTermsRho = (ARM_Vector *) generatedMatus->Clone();
       ARM_Vector* undMatuRho   = new ARM_Vector(1, itsCalibrationTenor);

       for (int i = 0; i < genMatusSize; i++)
       {
           if (( itsCalibMeth == ARM_CALIB_BOOTSTRAP )
               &&
               ( i > itsPreviousMaxCalibIndex )
              )
           {
              theRho->Elt(i, 0) =
                    itsInitialRho->ComputeVolatility(generatedMatus->Elt(i),
                                                     itsCalibrationTenor);
           }
           else
           {
              theRho->Elt(i, 0) =
                itsRho->ComputeVolatility(generatedMatus->Elt(i),
                                          itsCalibrationTenor);
           }
       }

       if (itsRho)
          delete itsRho;

       itsRho = new ARM_VolLInterpol(GetStartDate(),
                                     yearTermsRho,
                                     undMatuRho,
                                     theRho);
    }

    // NU

    if (FitNuOrNot)
    {
       ARM_Matrix* theNu        = new ARM_Matrix(genMatusSize, 1);
       ARM_Vector* yearTermsNu  = (ARM_Vector *) generatedMatus->Clone();
       ARM_Vector* undMatuNu    = new ARM_Vector(1, itsCalibrationTenor);

       for (int i = 0; i < genMatusSize; i++)
       {
           if (( itsCalibMeth == ARM_CALIB_BOOTSTRAP )
               &&
               ( i > itsPreviousMaxCalibIndex )
              )
           {
              theNu->Elt(i, 0) =
                   itsInitialNu->ComputeVolatility(generatedMatus->Elt(i),
                                                   itsCalibrationTenor);
           }
           else
           {           
              theNu->Elt(i, 0) =
                 itsNu->ComputeVolatility(generatedMatus->Elt(i),
                                          itsCalibrationTenor);
           }
       }

       if (itsNu)
          delete itsNu;

       itsNu  = new ARM_VolLInterpol(GetStartDate(),
                                     yearTermsNu,
                                     undMatuNu,
                                     theNu);
    }

    // BETA

    if (FitBetaOrNot)
    {
       ARM_Matrix* theBeta       = new ARM_Matrix(genMatusSize, 1);
       ARM_Vector* yearTermsBeta = (ARM_Vector *) generatedMatus->Clone();
       ARM_Vector* undMatuBeta   = new ARM_Vector(1, itsCalibrationTenor);

       for (int i = 0; i < genMatusSize; i++)
       {
           if (( itsCalibMeth == ARM_CALIB_BOOTSTRAP )
               &&
               ( i > itsPreviousMaxCalibIndex )
              )
           {
              theBeta->Elt(i, 0) =
                   itsInitialBeta->ComputeVolatility(generatedMatus->Elt(i),
                                                   itsCalibrationTenor);
           }
           else
           {           
              theBeta->Elt(i, 0) =
                 itsBeta->ComputeVolatility(generatedMatus->Elt(i),
                                            itsCalibrationTenor);
           }
       }

       if (itsBeta)
          delete itsBeta;

       itsBeta = new ARM_VolLInterpol(GetStartDate(),
                                      yearTermsBeta,
                                      undMatuBeta,
                                      theBeta);
    }

    if (itsLogRho)
       delete itsLogRho;

    itsLogRho = (ARM_VolCurve *) itsRho->Clone();


    if (itsLogNu)
       delete itsLogNu;

    itsLogNu  = (ARM_VolCurve *) itsNu->Clone();

	if (itsSmiledModelType == K_LD)
	    UpdateLogRhoNu();

    // Initialize maximum index of maturities calib.

    if ( itsCalibMeth == ARM_CALIB_BOOTSTRAP )
       itsPreviousMaxCalibIndex = generatedMatus->GetSize()-1;
}



void ARM_BSSmiledModel::ThisIsWhatToFit(int FitVolOrNot,
                                        int FitRhoOrNot, 
                                        int FitNuOrNot,
                                        int FitBetaOrNot)
{
    int i;

    int nbParamsToFit = 0;
    int paramsToFit[MAX_PARAMS_FIT];

    int paramSizeSig   = GetVolatility()->GetExpiryTerms()->GetSize();
    int paramSizeRho   = itsRho->GetExpiryTerms()->GetSize();
    int paramSizeNu    = itsNu->GetExpiryTerms()->GetSize();

    int paramSizeBeta  = 0;

    if (itsBeta)
       paramSizeBeta  = itsBeta->GetExpiryTerms()->GetSize();


    if ( itsCalibMeth == ARM_CALIB_GLOB )
    {
       if (FitVolOrNot)
       {
          nbParamsToFit = paramSizeSig;
       }

       if (FitRhoOrNot)
       {
          nbParamsToFit += paramSizeRho;
       }

       if (FitNuOrNot)
       {
          nbParamsToFit += paramSizeNu;
       }

       if (itsBeta)
       {
          if (FitBetaOrNot)
          {
             nbParamsToFit += paramSizeBeta;
          }
       }
    }
    else
    {
       /*
        * N.B : In the case of bootstraping we calibrate at most
        *       3 parameters : Sigma, Rho, Nu (and may be Beta!)
        */

       if (FitVolOrNot)
       {
          nbParamsToFit++;
       }

       if (FitRhoOrNot)
       {
          nbParamsToFit++;
       }

       if (FitNuOrNot)
       {
          nbParamsToFit++;
       }
       
       if (itsBeta)
       {
          if (FitBetaOrNot)
          {
             nbParamsToFit++;
          }
       }
    }

    if ( nbParamsToFit == 0 )
    {
       throw Exception(__LINE__, __FILE__,
                      ERR_INVALID_ARGUMENT,
             "At least one param. must be fitted : Sig, Rho, Nu (or Beta)");
    }

    if ( itsCalibMeth == ARM_CALIB_GLOB )
    {
       for (i = 0; i < paramSizeSig; i++)
       {
           paramsToFit[i] = (FitVolOrNot? 1 : 0);
       }

       for (i = 0; i < paramSizeRho; i++)
       {
           paramsToFit[i+paramSizeSig] = (FitRhoOrNot? 1 : 0);
       }

       for (i = 0; i < paramSizeNu; i++) 
       {
           paramsToFit[i+paramSizeSig+paramSizeRho] = (FitNuOrNot? 1 : 0);
       }

       if (itsBeta)
       {
          for (i = 0; i < paramSizeBeta; i++) 
          {
              paramsToFit[i+paramSizeSig+paramSizeRho+paramSizeNu] = (FitBetaOrNot? 1 : 0);
          }
       }
    }
    else
    {
       for (i = 0; i < paramSizeSig; i++)
       {
           paramsToFit[i] = 0;
       }

       if (FitVolOrNot)
          paramsToFit[paramSizeSig-1] = 1;


       for (i = 0; i < paramSizeRho; i++)
       {
           paramsToFit[i+paramSizeSig] = 0;
       }

       if (FitRhoOrNot)
          paramsToFit[paramSizeSig+(paramSizeRho-1)] = 1;


       for (i = 0; i < paramSizeNu; i++)
       {
           paramsToFit[i+paramSizeSig+paramSizeRho] = 0;
       }

       if (FitNuOrNot)
          paramsToFit[paramSizeSig+paramSizeRho+(paramSizeNu-1)] = 1;

       // BETA

       if (itsBeta)
       {
          for (i = 0; i < paramSizeBeta; i++)
          {
              paramsToFit[i+paramSizeSig+paramSizeRho+paramSizeNu] = 0;
          }

          if (FitBetaOrNot)
             paramsToFit[paramSizeSig+paramSizeRho+paramSizeNu+(paramSizeBeta-1)] = 1;
       }
    }

    this->SetParametersToFit(nbParamsToFit, paramsToFit);
}



void ARM_BSSmiledModel::BootstrapCalibrate(ARM_Portfolio* pf,
                                           int calVolOrNot,
                                           int calRhoOrNot,
                                           int calNuOrNot,
                                           double volTenor,
                                           double H_TimeStep,
                                           double minSig,
                                           double maxSig,
                                           double minRho,
                                           double maxRho,
                                           double minNu,
                                           double maxNu,
                                           int interpolMeth,
                                           double tol,
                                           long maxIter,
                                           int gradCalc,
                                           double lambda,
                                           double beginSmoothMatu,
                                           int calBetaOrNot,
                                           double minBeta,
                                           double maxBeta)
{
    ARM_Date AsOfDate = GetStartDate();

    // Prepare Calibration

    double timeStep = -1;

    ARM_Vector* generatedMatus = GenCalibMatus(pf, timeStep);
    int genMatusSize           = generatedMatus->GetSize();

   
    if ( volTenor < 0 )
    {
       ARM_Security* curAsset = pf->GetAsset(0);
       
       if (( curAsset->GetName() == ARM_CAPFLOOR ) 
           || 
           ( curAsset->GetName() == ARM_SWAPTION )
           ||
           ( curAsset->GetName() == ARM_DIGITAL ) 
          )
       {       
           if ( curAsset->GetName() == ARM_CAPFLOOR )
           {
              ARM_CapFloor* capFloor = (ARM_CapFloor *) curAsset;

              if ( capFloor->GetSwapLeg()->GetName() == ARM_CMSLEG ) 
              { 
                 itsCalibrationTenor = ((ARM_CMSLeg *) 
                                        capFloor->GetSwapLeg())->GetSwapYearTerm();
              }
              else
              {
                 itsCalibrationTenor = capFloor->GetSwapLeg()->GetIRIndex()->GetYearTerm();
              }
           }
           else if ( curAsset->GetName() == ARM_DIGITAL ) 
           {
              ARM_Digital* digital = (ARM_Digital *) curAsset;

              itsCalibrationTenor  = digital->GetSwapLeg()->GetIRIndex()->GetYearTerm();
           }
           else
           {
              ARM_Swaption* swaption = (ARM_Swaption *) curAsset;

              itsCalibrationTenor = (swaption->GetEndDate().GetJulian()
                                     -swaption->GetStartDate().GetJulian())/K_YEAR_LEN;
           }
       }
       else // Error
       {       
          throw Exception(__LINE__, __FILE__,
                          ERR_INVALID_ARGUMENT,
          "Only CAPS/FLOORS, DIGITALS or SWAPTIONS expected for BS smiled calibration");
       }
    }
    else
    {
       itsCalibrationTenor = volTenor;
    }
    
    PrepareNormVectors(generatedMatus);

    ARM_Portfolio* curPF = NULL;

    for (int i = 0; i < genMatusSize; i++)
    {
        curPF = pf->FilterByMaturity(generatedMatus->Elt(i),
                                     AsOfDate); 

        try
        {
            this->Calibrate(curPF,
                            calVolOrNot,
                            calRhoOrNot,
                            calNuOrNot,
                            volTenor,
                            H_TimeStep,
                            minSig,
                            maxSig,
                            minRho,
                            maxRho,
                            minNu,
                            maxNu,
                            interpolMeth,
                            tol,
                            maxIter,
                            gradCalc,
                            lambda,
                            beginSmoothMatu,
                            calBetaOrNot,
                            minBeta,
                            maxBeta);

            if (curPF)
               delete curPF;
        }

        catch(Exception& anARMExpt)
        {
            if (generatedMatus)
               delete generatedMatus;

            if (curPF)
               delete curPF;

            throw anARMExpt;

            return;
        }

        catch(...)
        {
            if (generatedMatus)
               delete generatedMatus;

            if (curPF)
               delete curPF;

            throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Calibrate(): Unrecognized Failure?!");
        }
    }

    // Clean up
    delete generatedMatus;
}



void ARM_BSSmiledModel::Calibrate(ARM_Portfolio* pf,
                                  int calVolOrNot,
                                  int calRhoOrNot,
                                  int calNuOrNot,
                                  double volTenor,
                                  double H_TimeStep,
                                  double minSig,
                                  double maxSig,
                                  double minRho,
                                  double maxRho, 
                                  double minNu,
                                  double maxNu, 
                                  int interpolMeth,
                                  double tol, 
                                  long maxIter, 
                                  int gradCalc,
                                  double lambda,
                                  double beginSmoothMatu,
                                  int calBetaOrNot,
                                  double minBeta,
                                  double maxBeta)
{
    if ( lambda < 1e-15 )
    {
       itsSmoothFlag = 0;
    }
    else
    {
       itsSmoothFlag = 1;

       itsLambda     = lambda;
    }

    itsCalibrationPortfolio = pf;

    if ( volTenor < 0 )
    {
       ARM_Security* curAsset = pf->GetAsset(0);
      
       if (( curAsset->GetName() == ARM_CAPFLOOR ) 
           || 
           ( curAsset->GetName() == ARM_SWAPTION )
           ||
           ( curAsset->GetName() == ARM_DIGITAL )
          )
       {       
          if ( curAsset->GetName() == ARM_CAPFLOOR )
          {
             ARM_CapFloor* capFloor = (ARM_CapFloor *) curAsset;

             if ( capFloor->GetSwapLeg()->GetName() == ARM_CMSLEG ) 
             { 
                itsCalibrationTenor = ((ARM_CMSLeg *) 
                                        capFloor->GetSwapLeg())->GetSwapYearTerm();
             }
             else
             {
                itsCalibrationTenor = capFloor->GetSwapLeg()->GetIRIndex()->GetYearTerm();
             }
          }
          else if ( curAsset->GetName() == ARM_DIGITAL )
          {
             ARM_Digital* digital = (ARM_Digital *) curAsset;
             itsCalibrationTenor  = digital->GetSwapLeg()->GetIRIndex()->GetYearTerm();
          }
          else
          {
             ARM_Swaption* swaption = (ARM_Swaption *) curAsset;
             itsCalibrationTenor = (swaption->GetEndDate().GetJulian()
                                      -swaption->GetStartDate().GetJulian())/K_YEAR_LEN;
          }
       }
       else // Error
       {       
          throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
                "Only CAPS/FLOORS, DIGITALS or SWAPTIONS expected for BS smiled calibration");
       }
    }
    else
    {
       itsCalibrationTenor = volTenor;
    }

    // Prepare Calibration
    ARM_Vector* generatedMatus = GenCalibMatus(pf, H_TimeStep);
    int genMatusSize           = generatedMatus->GetSize();


    UpdateSigRhoNuByInterpol(H_TimeStep,
                             calVolOrNot, calRhoOrNot, calNuOrNot,
                             generatedMatus, interpolMeth,
                             calBetaOrNot);

    SetCalibMatus(generatedMatus);


    // Find the Index of the first Smooth Maturity
    if ((itsSmoothFlag) && ( beginSmoothMatu > 0 ))
    {
       // Update itsBeginSmoothMatuIndex
  
       int sz = generatedMatus->GetSize();

       int k = 1; // We begin at 1

       while (( k < sz ) 
              && 
              ( generatedMatus->Elt(k) < beginSmoothMatu ))
       {
           k++;
       }

       if ( k != sz ) // found
       {
          itsBeginSmoothMatuIndex = k; 
       } 
       else
       {
          itsSmoothFlag = 0;
       }
    }

    delete generatedMatus;


    ARM_Vector* params = BuildParameters();

    ARM_Model::SetParameters(NULL);
    ARM_Model::SetParameters(params);

    if (params)
       delete params;

    ThisIsWhatToFit(calVolOrNot, calRhoOrNot, calNuOrNot,
                    calBetaOrNot);

    int nbParamsToFit = GetNbParamsToFit();

    ARM_Vector LB(nbParamsToFit);
    ARM_Vector UB(nbParamsToFit);

    int index = 0;

    /*
     * N.B : In the case of bootstraping we calibrate at most
     *       3 parameters : Sigma, Rho, Nu
     */

    int idxBootStrap = 0;

    if (calVolOrNot)
    {
       if ( itsCalibMeth == ARM_CALIB_GLOB )
       {
          for (int i = 0; i < genMatusSize; i++)
          {
              // Sigma
              LB[index] = minSig;
              UB[index] = maxSig;
              index++;
          }
       }
       else
       {
          // Sigma
          LB[idxBootStrap] = minSig;
          UB[idxBootStrap] = maxSig;
          idxBootStrap++;
       }
    }

    if (calRhoOrNot)
    {
       if ( itsCalibMeth == ARM_CALIB_GLOB )
       {
          for (int i = 0; i < genMatusSize; i++)
          {
              // Rho
              LB[index] = minRho;
              UB[index] = maxRho;
              index++;
          }
       }
       else
       {
          // Rho
          LB[idxBootStrap] = minRho;
          UB[idxBootStrap] = maxRho;
          idxBootStrap++;
       }
    }

    if (calNuOrNot)
    {
       if ( itsCalibMeth == ARM_CALIB_GLOB )
       {
          for (int i = 0; i < genMatusSize; i++)
          {
              // Nu
              LB[index] = minNu;
              UB[index] = maxNu;
              index++;
          }
       }
       else
       {
          // Nu
          LB[idxBootStrap] = minNu;
          UB[idxBootStrap] = maxNu;
          idxBootStrap++;
       }
    }

    // BETA

    if (itsBeta && calBetaOrNot)
    {
       if ( itsCalibMeth == ARM_CALIB_GLOB )
       {
          for (int i = 0; i < genMatusSize; i++)
          {
              // Beta
              LB[index] = minBeta;
              UB[index] = maxBeta;
              index++;
          }
       }
       else
       {
          // Beta
          LB[idxBootStrap] = minBeta;
          UB[idxBootStrap] = maxBeta;
          idxBootStrap++;
       }
    }

    ARM_MinimizerOnPortfolio Calibrator(this,
                                        pf,
                                        tol,
                                        maxIter,
                                        &LB,
                                        &UB,
                                        gradCalc);

    try
    {
       if (itsSmoothFlag)
       {
          long nbSmoothFuncs = 3*(itsCalibMatus->GetSize()-2);

          if ( nbSmoothFuncs >= 3 )
          {
             Calibrator.SetNbSmoothingFuncs(nbSmoothFuncs);
          }
       }

       Calibrator.Calibrate();
    }

    catch(Exception& expt)
    {
        throw expt;
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
             "Unrecognized Failure in : ARM_BSSmiledModel::Calibrate()");
    }


    SetOptimizeDiag(Calibrator.GetOptimizeDiag());

	if (itsSmiledModelType == K_LD)
		this->UpdateLogRhoNu();
}



double ARM_BSSmiledModel::ComputeVol(double maturity, double volTenor, 
                                     double fwd, double strike, int mode )
{
    double Bt = 0.0;
    double At = 0.0;
    double Wt = 0.0;


	if ( maturity <= 0.0 )
    {
	   return(0.0);
    }

    double spot = fwd;
    
    double volatility = GetVolatility(mode)->ComputeVolatility(maturity,
                                                               volTenor);
    if ( volatility < 0.0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
            "ARM_BSSmiledModel::ComputeVol() : volatility is negative!?");
    }
    
	if (IS_ZERO(strike) || ( strike < 0 ))
    {
	   strike = 0.05; // strike = 0.00001;
    }

    if (itsBeta)
    {
       double beta = itsBeta->ComputeVolatility(maturity, volTenor);

       if (!(IS_ONE(beta)))
       {
#ifdef SABR_APPROX_NORMAL // JDA Approach

          double sigma = CalcSABRVolBetaDiffOne(maturity,
                                                volTenor,
                                                fwd,
                                                strike);
#else
          double sigma = CptSABRVolBetaDiffOneNormalOrDirectExact(maturity,
                                                                  volTenor,
                                                                  spot,
                                                                  strike);
#endif

          itsSABRSig = sigma*100.0;

          return(itsSABRSig);
       }
    }
 
    // Beta == 1

    if ( itsSmiledModelType != K_LD ) // SABR type
    {
       double sigSABR;
       double ATMVol;
       double Alpha_t;
        
       ATMVol = volatility*0.01;
        
       // Retrieve the data corresponding to the tenor
        
       double Rho_t = itsRho->ComputeVolatility(maturity,
                                                volTenor);
        
       double Nu_t  = itsNu->ComputeVolatility(maturity,
                                               volTenor);
        
        
       if (IS_ZERO(Nu_t))
       {
          sigSABR = ATMVol;
            
          itsSABRSig = sigSABR*100.0;

          return(itsSABRSig);
       }
       else
       {
          if (IS_ONE(Rho_t))
          {
             Rho_t = double(0.99);
          }
            
          if (IS_ZERO(Rho_t))
          {
             Alpha_t = ATMVol/(1.0+(SQR(Nu_t)*maturity/12.0));
          }
          else
          {
             double sqrtTerm;
                
             double rho2_3_24;
                
             rho2_3_24 = ((2.0-3.0*SQR(Rho_t))/24.0)*SQR(Nu_t)
                           *maturity;
                
             sqrtTerm = sqrt(SQR(1.0+rho2_3_24)
                          +(ATMVol*Rho_t*Nu_t*maturity));
                
             Alpha_t = (2.0/(Rho_t*Nu_t*maturity))*(sqrtTerm
                        -(1.0+rho2_3_24));
          }
            
          double KSI;
            
          if ( itsSmiledModelType == K_SABR_GEO ) // SABR geometric
          {
             KSI = (Nu_t/Alpha_t)*log(spot/strike);
          }
          else
          {
             KSI = (Nu_t/Alpha_t)*2.0*(spot-strike)/(spot+strike);
          }
            
          double M;
            
          if ( fabs(KSI) < SABR_PREC )
          {
             M = 1.0-(Rho_t*KSI/2.0)
                    +((1.0/6.0)-SQR(Rho_t)/4.0)*SQR(KSI);
          }
          else
          {
             double denom;
                
             denom = log((sqrt(1.0-2.0*Rho_t*KSI
                      +SQR(KSI))+KSI-Rho_t)
                      /(1.0-Rho_t));
                
             M = KSI/denom;
          }
            
          sigSABR = ATMVol*M*100.0;
            
          itsSABRSig = sigSABR;

          return(sigSABR);
       }
    }
    
    throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                   "ARM_BSSmiledModel::ComputeVol() unavailable for LD type");
}



double ARM_BSSmiledModel::EuroCaplet(double Settlement,
                                     double ResetDate,
                                     double StartDate, 
                                     double EndDate,
                                     double PayDate,
                                     int OptionType, 
                                     double Fwd, 
                                     double OptionStrike, 
                                     int DomOrFrg,
                                     int CompMeth, 
                                     int DayCount,
                                     bool IsCMS,
                                     bool isTreasury,
                                     int FixFrequency,
                                     double UnderlyingTenor,
                                     int YieldDecomp,
                                     double Margin,
                                     ARM_Currency* ccy,
                                     StoreFwdRateAndCapletInfo* StoreInfo)
{
    double price = 0.0;

    double asofDate = GetStartDate().GetJulian();

    ARM_Date startDate(asofDate+StartDate*K_YEAR_LEN);
    ARM_Date endDate;

    if (IsCMS)
    {
       endDate = startDate;
       endDate.AddYears(UnderlyingTenor);
    }
    else
    {
       endDate = ARM_Date(asofDate+EndDate*K_YEAR_LEN);
    }
    
    double fwdAdj = 0.0;
    double fwdAdjNoDecap = 0.0;
    double strikeAdj = 0.0;
    double volAdj = 0.0;

    ComputeCapletBSParameters(ResetDate,
                              StartDate,
                              EndDate,
                              PayDate,
                              OptionType,
                              Fwd, 
                              OptionStrike,
                              DomOrFrg,
                              CompMeth,
                              DayCount,
                              IsCMS,
                              isTreasury,
                              FixFrequency,
                              UnderlyingTenor,
                              YieldDecomp,
                              Margin,
                              ccy,
                              fwdAdj,
                              fwdAdjNoDecap,
                              strikeAdj,
                              volAdj,
                              StoreInfo);

    // If the option didn't reach its maturity date
    if ( ResetDate > 0 )
    {
       // We interpolate volatility
       double vol = 0.0;

       double rawFwdRate = 0.0;

       if (IsCMS)
       {
          rawFwdRate = SwapRate(startDate, endDate,ccy->GetFixedPayFreq(), 
                                ccy->GetFwdRule(), ccy)*100.0;
       }
       else
       {
          rawFwdRate = ExpectedFwdYield(startDate,
                                        endDate,
                                        endDate,
                                        CompMeth,
                                        DayCount,
                                        1,
                                        1,
                                        1,
                                        0.0);
        }

        double strikeAdjRecap = strikeAdj;

        if ( YieldDecomp > 0 )
        {
           strikeAdjRecap = FromRateToRate(strikeAdjRecap, 1.0, YieldDecomp, 
                                           K_COMP_ANNUAL);
        }

        // vol=Could be ATM vol or Alpha

        vol = GetVolatility(K_IRG)->ComputeVolatility(ResetDate, 
                                                 UnderlyingTenor)/100.0;

        if (itsSigmaOrAlphaInput)
        {
           vol *= volAdj;
        }

        // We finally apply the Black&Scholes formulae
        double delta = 0.0;
        double gamma = 0.0;
        double kappa = 0.0;
        double vega = 0.0;
        double theta = 0.0;
        double rho = 0.0;

        if ( strikeAdj >= 0 )
        {
           price = bsOption(fwdAdj, 
                            strikeAdj, 
                            vol, 
                            0.0, 
                            0.0,
                            ResetDate, 
                            OptionType,
							UnderlyingTenor);


           // Compute Greeks
            
           delta = bsDelta(fwdAdj,
                           strikeAdj,
                           itsSABRSig/100.0,
                           0.0,
                           0.0,
                           ResetDate,
                           OptionType);

           gamma = bsGamma(fwdAdj,
                           strikeAdj,
                           itsSABRSig/100.0,
                           0.0,
                           0.0,
                           ResetDate,
                           OptionType);

           kappa = bsKappa(fwdAdj,
                           strikeAdj,
                           itsSABRSig/100.0,
                           0.0,
                           0.0,
                           ResetDate,
                           OptionType);

           theta = bsTheta(fwdAdj,
                           strikeAdj,
                           itsSABRSig/100.0,
                           0.0,
                           0.0,
                           ResetDate,
                           OptionType);

           vega = bsVega(fwdAdj,
                         strikeAdj,
                         itsSABRSig/100.0,
                         0.0,
                         0.0,
                         ResetDate,
                         OptionType);

           rho = bsRho(fwdAdj,
                       strikeAdj,
                       itsSABRSig/100.0,
                       0.0,
                       0.0,
                       ResetDate,
                       OptionType);
        }
        else
        {
           if ( OptionType == K_CALL )
           {
              price = fwdAdj-strikeAdj;
           }
           else
           {
              price = 0.0;
           }
        }

        if (StoreInfo)
        {
           StoreInfo->StoreCaplet(fwdAdj,
                                  strikeAdj,
                                  itsSABRSig/100/volAdj,
                                  itsSABRSig/100,
                                  delta,
                                  gamma,
                                  kappa,
                                  vega,
                                  theta,
                                  rho);
        }
    }
    else // If the option did reach its maturity date
    {
       // We apply the payoff formulae
       if ( OptionType*Fwd > strikeAdj*OptionType )
       {
          price = OptionType*Fwd-strikeAdj*OptionType;
       }
       else
       {
          price = 0.0;
       }

	   if (StoreInfo)
	   {
	       StoreInfo->StoreCaplet(Fwd,
								  strikeAdj,
								  0.0,
								  0.0,
								  0.0,
								  0.0,
								  0.0,
								  0.0,
								  0.0,
								  0.0);
	    }
    }

    return(price);
}



ARM_VolCurve* ARM_BSSmiledModel::FromSigmaToAlpha(ARM_VolCurve* sigma,
                                                  double strike)
{
    ARM_VolCurve* alphaRes = NULL;

    ARM_VolCurve* originalVol = NULL;



    if (sigma)
    {
       alphaRes = (ARM_VolCurve *) sigma->Clone();

       originalVol = sigma;
    }
    else
    {
       originalVol = (ARM_VolCurve *) GetVolatility();
       
       alphaRes = (ARM_VolCurve *) GetVolatility()->Clone();
    }

    
    ARM_Vector* tenors = ((ARM_VolLInterpol *) originalVol)->GetStrikes();
    ARM_Vector* matuOp = originalVol->GetExpiryTerms();
    ARM_Matrix* matVol = alphaRes->GetVolatilities();

    double F;
    double K;
    double mat;
    double sigmaATM;
    double rho;
    double nu;
    double beta; 
    double weight = itsSABRWeight;
    int    type   = itsSmiledModelType;

    double tenor;

    int i, j, matuSize, sizeTenor;

    matuSize  = matuOp->GetSize();
    sizeTenor = tenors->GetSize();

    for (i = 0; i < matuSize; i++)
    {
        for (j = 0; j < sizeTenor; j++)
        {
            mat   = matuOp->Elt(i);
            tenor = tenors->Elt(j);

            if ( strike < 0.0 )
            {
               double curFwd = CalcFwd(mat, tenor)/100.0;

               K = F = curFwd;
            }
            else
            {
               K = F = strike;
            }

            sigmaATM = originalVol->ComputeVolatility(mat, tenor)/100.0;
            
            rho = itsRho->ComputeVolatility(mat, tenor);
            nu  = itsNu->ComputeVolatility(mat, tenor);

            if (itsBeta) // Beta <> 1
            {
               beta = itsBeta->ComputeVolatility(mat, tenor);
            }
            else
            {
               beta = 1.0;
            }

            double curAlpha;

            try
            {
                if (!(IS_ONE(beta)))
                {
                   curAlpha = ComputeAlpha(F, K,
                                           mat,
                                           sigmaATM,
                                           rho, nu, beta, 
                                           weight, type);
                }
                else
                {                
                   curAlpha = ComputeAlphaBetaEqOne(F, K,
                                                    mat, 
                                                    sigmaATM,
                                                    rho, nu,
                                                    type);
                }
            }

            catch(Exception& anExpt)
            {
                throw anExpt;

                return(NULL);
            }

            catch(...)
            {
                throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                "ARM_BSSmiledModel::FromSigmaToAlpha : Unrecognized Failure!?");

                return(NULL);
            }

            matVol->Elt(i, j) = curAlpha*100.0;
        }
    }
    
    return(alphaRes);
}
 


ARM_VolCurve* ARM_BSSmiledModel::FromAlphaToSigma(ARM_VolCurve* alpha,
                                                  double strike)
{
    ARM_VolCurve* sigmaRes = NULL;

    ARM_VolCurve* originalVol = NULL;


    if (alpha)
    {
       sigmaRes = (ARM_VolCurve *) alpha->Clone();

       originalVol = alpha;
    }
    else
    {
       originalVol = (ARM_VolCurve *) GetVolatility();
       
       sigmaRes = (ARM_VolCurve *) GetVolatility()->Clone();
    }

    
    ARM_Vector* tenors = ((ARM_VolLInterpol *) originalVol)->GetStrikes();
    ARM_Vector* matuOp = originalVol->GetExpiryTerms();
    ARM_Matrix* matVol = sigmaRes->GetVolatilities();

    double F;
    double K;
    double mat;
  
    double tenor;

    int i, j, matuSize, sizeTenor;

    matuSize  = matuOp->GetSize();
    sizeTenor = tenors->GetSize();

    for (i = 0; i < matuSize; i++)
    {
        for (j = 0; j < sizeTenor; j++)
        {
            mat   = matuOp->Elt(i);
            tenor = tenors->Elt(j);

            if ( strike < 0.0 )
            {
               double curFwd = CalcFwd(mat, tenor);

               K = F = curFwd;
            }
            else
            {
               K = F = strike;
            }

            double beta;

            if (itsBeta) // Beta <> 1
            {
               beta = itsBeta->ComputeVolatility(mat, tenor);
            }
            else
            {
               beta = 1.0;
            }

            double rho     = itsRho->ComputeVolatility(mat, tenor);
            double nu      = itsNu->ComputeVolatility(mat, tenor);

            double alpha_t = originalVol->ComputeVolatility(mat, tenor)*0.01;
                
            double curSigma;

            try
            {
                if (!(IS_ONE(beta)))
                {
                   curSigma = CptSABRVolBetaDiffOneNormalOrDirectExact(mat,
                                                                       tenor,
                                                                       K, F);
                }
                else
                {
                   curSigma = ComputeSigATMBetaEqOne(F*0.01, K*0.01,
                                                     mat, 
                                                     alpha_t,
                                                     rho, nu,
                                                     itsSmiledModelType);
                }
            }

            catch(Exception& anExpt)
            {
                throw anExpt;

                return(NULL);
            }

            catch(...)
            {
                throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                "ARM_BSSmiledModel::FromAlphaToSigma : Unrecognized Failure!?");

                return(NULL);
            }

            matVol->Elt(i, j) = curSigma*100.0;
        }
    }
    
    return(sigmaRes);
}



/*---------------------------------------------------------------------------*/
/*---- End Of File ----*/
