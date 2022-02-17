/*
 * $Log: bssmiled.h,v $
 * Revision 1.43  2004/05/25 15:26:29  mab
 * Correction : free objects in BuildParameters
 *
 * Revision 1.42  2004/03/31 16:00:41  rguillemot
 * CMT Bug Fix
 *
 * Revision 1.41  2004/03/11 15:19:51  mab
 * Improvement in View
 *
 * Revision 1.40  2004/03/05 17:27:13  mab
 * Modification in : View
 *
 * Revision 1.39  2004/02/20 09:30:21  mab
 * Change : itsAlphaOrSigmaInput to:
 * itsSigmaOrAlphaInput
 *
 * Revision 1.38  2004/02/17 16:02:26  mab
 * Added : Euro caplet
 *
 * Revision 1.37  2004/02/17 09:56:02  mab
 * Management of alpha and ATM Vol
 *
 * Revision 1.36  2004/02/04 11:26:48  mab
 * itsSABRSig must be in 100 base
 *
 * Revision 1.35  2004/02/03 15:55:39  mab
 * correction of SetBeta
 *
 * Revision 1.34  2004/01/29 10:57:25  mab
 * Improvement
 *
 * Revision 1.33  2004/01/29 10:47:12  mab
 * Comment added on SABR flag
 *
 * Revision 1.32  2004/01/29 10:41:04  mab
 * Improvelments for case Beta <> 1
 *
 * Revision 1.31  2004/01/29 10:22:56  mab
 * Management of case Beta <> 1
 *
 * Revision 1.30  2004/01/19 09:49:52  mab
 * IS_ZERO modified
 *
 * Revision 1.29  2004/01/09 10:35:04  arm
 * IS_ONE corrected
 *
 * Revision 1.28  2004/01/07 09:51:51  arm
 * ComputeFunc in .cpp
 *
 * Revision 1.27  2003/09/23 18:01:26  ebenhamou
 * move computevol to cpp and corrected in the case of nu = 0
 *
 * Revision 1.26  2003/09/23 06:21:15  ebenhamou
 * remove unused var
 *
 * Revision 1.25  2003/09/22 15:33:21  mab
 * Added : double ComputeVol()
 *
 * Revision 1.24  2003/09/11 13:40:20  mab
 * bssmiled.h : change in : UpdateLogRhoNu()
 *
 * Revision 1.23  2003/08/20 17:34:56  sgasquet
 * Ajout include de swaption.h
 *
 * Revision 1.22  2003/08/19 12:34:30  jpriaudel
 * use SABR flag
 *
 * Revision 1.21  2003/08/13 09:17:57  sgasquet
 * Ajout SABR avec approx arithmetique
 *
 * Revision 1.20  2003/07/01 19:09:28  mab
 * improvements
 *
 * Revision 1.19  2003/06/13 17:23:24  mab
 * Added : if ( N_SIG == 0.0 ) N_SIG = 1.0;
 * same thing for N_RHO , N_NU
 *
 * Revision 1.18  2003/06/02 08:10:11  mab
 * Changes in the code in order to manage
 * correctly already fitted parameters
 *
 * Revision 1.17  2003/05/05 13:00:32  mab
 * Added in bs Option : if ( volatility < 0.0 ) ...
 *
 * Revision 1.16  2003/03/31 11:12:29  mab
 * generatedPf->SortByExpiry(); ADDED
 *
 * Revision 1.15  2003/03/13 17:07:28  mab
 * Handling of conditional smoothing (see beginSmoothMatuIndex)
 *
 * Revision 1.14  2003/03/07 13:13:39  mab
 * Version with conditional Smoothing
 *
 * Revision 1.13  2003/02/14 17:45:48  mab
 * Added : if (smoothFuncValue) delete smoothFuncValue;
 *
 * Revision 1.12  2003/02/14 15:47:39  mab
 * Improvements
 *
 * Revision 1.11  2003/02/14 09:57:45  mab
 * Integrating bootstrapping calibration
 *
 * Revision 1.10  2003/02/11 16:51:32  mab
 * New version of bsOption
 *
 * Revision 1.9  2003/01/17 10:07:10  mab
 *  Initial version
 *
 * Revision 1.8  2002/12/19 10:01:59  mab
 * Version with calibration and bsOption old specifications
 *
 * Revision 1.7  2002/12/17 13:53:40  mab
 * This version integrate the calibration
 *
 * Revision 1.6  2002/12/11 10:59:25  mab
 * Generalisation by using : ARM_VolCurve
 *
 * Revision 1.5  2002/11/15 10:21:44  mab
 * added double GetSABRSigma(void),
 * if ( fabs(spot-strike) < 1.0e-3 )
 * replaced by : if ( fabs(spot-strike) < 1.0e-16 )
 *
 * Revision 1.4  2002/10/11 08:14:37  mab
 * Correction in Alpha_t Formulae
 *
 * Revision 1.3  2002/10/09 09:37:30  mab
 * Take in account SABR modelisation
 *
 * Revision 1.2  2002/03/28 13:19:43  mab
 * Alpha replaced by Rho's
 *
 * Revision 1.1  2002/02/25 13:53:35  mab
 * Initial revision
 *
 */


/*----------------------------------------------------------------------*

   bssmiled.h
 
   The smiled Black & Scholes framework(SABR)
 
*-----------------------------------------------------------------------*/
#ifndef _BSSMILEDMODEL_H
#define _BSSMILEDMODEL_H




#include "armglob.h"
#include "linalg.h"
#include "zerocurv.h"
#include "volcurv.h"
#include "volint.h"
#include "ycmodel.h"
#include "bsmodel.h"
#include "bsxtic.h"
#include "capfloor.h"
#include "swaption.h"
#include "merge.h"
#include "genminimizer.h"
#include "fromto.h"

// ARM_SABRVol structure
#include "SABRVol.h"



#define SABR_PREC 1e-7 



inline int IS_ZERO(double x) // at 1e-8
{
    if ( fabs(x) < 1e-8 )
    {
       return(1);
    }
    else
    {
       return(0);
    }
}



inline int IS_ONE(double x) // at 1e-6
{
    long truncated_x = long(floor(x*100000.0));

    return( truncated_x == 100000L );
}


extern double ComputeSigATMBetaEqOne(double f, double K,
                                     double maturity, 
                                     double alpha,
                                     double Rho_t, double Nu_t,
                                     int type);

class ARM_Date;
class ARM_Vector;


class ARM_BSSmiledModel : public ARM_BSModel
{
    private:

        int itsSmiledModelType;

        // 0: LD type (no SABR)
        // 1: SABR geometric
        // 2: SABR arithmetic
       
        /* case Beta < 1 */
        // 3: K_SABR_IMPLNVOL
        // 4: K_SABR_WEIGHT
 
        ARM_VolCurve* itsInitialSig;
        ARM_VolCurve* itsInitialRho;
        ARM_VolCurve* itsInitialNu;
        ARM_VolCurve* itsInitialBeta;

        ARM_VolCurve* itsRho;
        ARM_VolCurve* itsLogRho;
       
        ARM_VolCurve* itsNu;
        ARM_VolCurve* itsLogNu;

        ARM_VolCurve* itsBeta;

        // the SABR Weight between arith smile and
        // implied LN vol smile 

        double itsSABRWeight;
 
        double itsCalibrationTenor;

        ARM_Portfolio* itsCalibrationPortfolio;

        ARM_Vector*    itsCalibMatus;


        double itsSABRSig;

        int itsSmoothFlag;

        double itsLambda;

        double itsNorm_SIG;
        double itsNorm_RHO;
        double itsNorm_NU;

        // The following vectors are for optimizations needs
        // when Normalizing

        ARM_Vector*  itsNormCoefSIG;
        ARM_Vector*  itsNormCoefRHO;
        ARM_Vector*  itsNormCoefNU;


        int          itsBeginSmoothMatuIndex;

        int          itsCalibMeth; // ARM_CALIB_GLOB , ARM_CALIB_BOOTSTRAP

        int          itsPreviousMaxCalibIndex;

        char         itsOptimizationDiag[ARM_OPTIMIZ_BUF_DIAG_SIZE];

        int          itsSigmaOrAlphaInput;			// 1: Sigma ;  0 : Alpha

		int			 itsConvToAdjVolInCaseOfAlpha;	// 0: no conversion ; 1 : conversion
													// This flag is used to launch or not the conversion function:
													// FromAlphaToSigma(...)

        ARM_BSModel* itsFatherBSModel; // The eventual BS father model
                                       // A SHARED object so don't delete



        void Init(void)
        {
            SetName(ARM_BSSMILEDMODEL);

            SetBSCategory(1);

            itsSmiledModelType = K_LD;

            itsInitialSig  = NULL;
            itsInitialRho  = NULL;
            itsInitialNu   = NULL;
            itsInitialBeta = NULL;

            itsRho    = NULL;
            itsLogRho = NULL;

            itsNu    = NULL;
            itsLogNu = NULL;

            itsBeta  = NULL;

            itsSABRWeight = 0.5;

            itsCalibrationTenor = 0.0;

            itsCalibrationPortfolio = NULL;

            itsCalibMatus = NULL;

            itsSABRSig    = -10000;

            itsSmoothFlag = 0;

            itsLambda     = 0.0;

            itsNorm_SIG   = -1.0;
            itsNorm_RHO   = -1.0;
            itsNorm_NU    = -1.0;

            itsNormCoefSIG = NULL;
            itsNormCoefRHO = NULL;
            itsNormCoefNU  = NULL;

            itsBeginSmoothMatuIndex = 1;

            itsCalibMeth = ARM_CALIB_GLOB;

            itsPreviousMaxCalibIndex = 0;

            strcpy(itsOptimizationDiag, 
                   "===> OK: An optimal Solution has been found");

            itsSigmaOrAlphaInput = 1; // the input is Sigma0

			itsConvToAdjVolInCaseOfAlpha = 1; // default: convert from alpha to sigma

            itsFatherBSModel     = NULL;
        }

        void SaveInitialCalibData(void)
        {
            ARM_VolCurve* vol = GetVolatility();


            if (itsInitialSig)
            {
               delete itsInitialSig;

               itsInitialSig = NULL;
            }

            if (vol)
               itsInitialSig = (ARM_VolCurve *) vol->Clone();


            if (itsInitialRho)
            {
               delete itsInitialRho;

               itsInitialRho = NULL;
            }
    
            if (itsRho)
               itsInitialRho = (ARM_VolCurve *) itsRho->Clone();

            if (itsInitialNu)
            {
               delete itsInitialNu;

               itsInitialNu = NULL;
            }
          
            if (itsNu)
               itsInitialNu = (ARM_VolCurve *) itsNu->Clone();

            if (itsInitialBeta)
            {
               delete itsInitialBeta;

               itsInitialBeta = NULL;
            }

            if (itsBeta)
               itsInitialBeta = (ARM_VolCurve *) itsBeta->Clone();
        }

    public:

        ARM_BSSmiledModel(void)
        {
            Init();
        }

        ARM_BSSmiledModel(ARM_Date& startDate, double spot, 
                          ARM_ZeroCurve* dividend,
                          ARM_ZeroCurve* discountRate, 
                          ARM_VolCurve* volatility,
                          int volType = K_PRICE,
                          ARM_Vector* matu  = NULL,
                          ARM_Vector* rho  = NULL,
                          ARM_Vector* nu = NULL,
                          int SABR = 0,
                          ARM_CorrelManager* correlManager = NULL,
						  ARM_ZeroCurve* realDiscountRate = NULL,
						  ARM_VolCurve* AdjCvxVol = NULL,
						  int ConvToAdjVolWhenAlpha = 0);

     
        ARM_BSSmiledModel(ARM_Date& startDate, double spot,
                          double flatdiv,
                          double flatdiscrate, double flatvol,
                          int volType = K_PRICE,
                          ARM_Vector* matu  = NULL,
                          ARM_Vector* rho  = NULL,
                          ARM_Vector* nu = NULL,
                          int SABR = 0,
                          ARM_CorrelManager* correlManager = NULL);

        // A more general constructor : in fact accepting 
        // multi-tenors data when in the context of the Pricing
        ARM_BSSmiledModel(ARM_Date& startDate, double spot,
                          ARM_ZeroCurve* dividend,
                          ARM_ZeroCurve* discountRate,	// forecast curve
                          ARM_VolCurve* volatility,
                          int volType = K_PRICE,
                          ARM_VolCurve* rho  = NULL,
                          ARM_VolCurve* nu = NULL,
                          int SABR = 0,
                          ARM_VolCurve* beta = NULL,
                          double SabrWeight = 0.5,
                          int alphaOrSigmaInput = 1,
                          ARM_CorrelManager* correlManager = NULL,
						  ARM_ZeroCurve* realDiscountRate = NULL,	// discount curve
						  ARM_VolCurve* AdjCvxVol = NULL,
						  int ConvToAdjVolWhenAlpha = 0);		

		// General constructor with SABR volatility structure
		ARM_BSSmiledModel(ARM_Date& startDate, double spot,
						  ARM_ZeroCurve* dividend,
						  ARM_ZeroCurve* discountRate,		// forecast curve
						  ARM_SABRVol*	SABRVol,			// SABR volatility structure
						  int volType,
						  ARM_CorrelManager* correlManager,
						  ARM_ZeroCurve* realDiscountRate,
						  ARM_VolCurve* AdjCvxVol = NULL,
						  int ConvToAdjVolWhenAlpha = 0);

        ARM_BSSmiledModel(const ARM_BSSmiledModel& bs);
    
        ARM_BSSmiledModel& operator = (const ARM_BSSmiledModel& bs);
    
       ~ARM_BSSmiledModel(void) 
        {
            if (itsInitialSig)
               delete itsInitialSig;

            if (itsInitialRho)
               delete itsInitialRho;

            if (itsInitialNu)
               delete itsInitialNu;

            if (itsInitialBeta)
               delete itsInitialBeta;

            if (itsRho)
               delete itsRho;

            if (itsLogRho)
               delete itsLogRho;

            if (itsNu)
               delete itsNu;

            if (itsBeta)
               delete itsBeta;

            if (itsLogNu)
               delete itsLogNu;

            if (itsCalibMatus)
               delete itsCalibMatus;

            if (itsNormCoefSIG)
               delete itsNormCoefSIG;

            if (itsNormCoefRHO)
               delete itsNormCoefRHO;

            if (itsNormCoefNU)
               delete itsNormCoefNU;
        }

	    int GetConversionFlag(void)
		{
			return itsConvToAdjVolInCaseOfAlpha;
		}

		void SetConversionFlag(int flag)
		{
			itsConvToAdjVolInCaseOfAlpha = flag;
		}
        
        int GetSABRFlag(void)
        {
            return(itsSmiledModelType);
        }

        void SetSABRFlag(int SABRFlag)
        {
            itsSmiledModelType = SABRFlag;
        }

        void SetSigmaOrAlphaInput(int flag)
        {
            itsSigmaOrAlphaInput = flag;
        }

        int GetSigmaOrAlphaInput(void)
        {
           return(itsSigmaOrAlphaInput);
        }

        double GetWeight(void)
        {
            return(itsSABRWeight);
        }

        void SetWeight(double weight)
        {
            itsSABRWeight = weight;
        }

        void SetCalibMeth(int meth)
        {
            itsCalibMeth = meth;
        }

        void SetCalibMatus(ARM_Vector* matus)
        {
            if (itsCalibMatus)
               delete itsCalibMatus;

            if (matus)    
               itsCalibMatus = (ARM_Vector *) matus->Clone();
        }

        ARM_BSModel* GetBSFather(void)
        {
            return(itsFatherBSModel);
        }

        void SetBSFather(ARM_BSModel* bsModel)
        {
            itsFatherBSModel = bsModel;
        }

        void BitwiseCopy(const ARM_Object* srcBSModel)
        {
            ARM_BSSmiledModel* bsmodel = 
                              (ARM_BSSmiledModel *) srcBSModel;

            
            itsSmiledModelType = bsmodel->itsSmiledModelType;

 
            if (itsInitialSig)
            {
               delete itsInitialSig;

               itsInitialSig = NULL;
            }

            if (bsmodel->itsInitialSig)
            {
               itsInitialSig = (ARM_VolCurve *)
                               bsmodel->itsInitialSig->Clone();
            }

            if (itsInitialRho)
            {
               delete itsInitialRho;

               itsInitialRho = NULL;
            }
         
            if (bsmodel->itsInitialRho)
            {
               itsInitialRho = (ARM_VolCurve *)
                               bsmodel->itsInitialRho->Clone();
            }

            if (itsInitialNu)
            {
               delete itsInitialNu;

               itsInitialNu = NULL;
            }

            if (bsmodel->itsInitialNu)
            {
               itsInitialNu = (ARM_VolCurve *)
                               bsmodel->itsInitialNu->Clone();
            }

            if (itsInitialBeta)
            {
               delete itsInitialBeta;

               itsInitialBeta = NULL;
            }

            if (bsmodel->itsInitialBeta)
            {
               itsInitialBeta = (ARM_VolCurve *)
                                bsmodel->itsInitialBeta->Clone();
            }

            if (itsRho)
            {
               delete itsRho;

               itsRho = NULL;
            }

            if (bsmodel->itsRho)
            {
               itsRho = (ARM_VolCurve *) bsmodel->itsRho->Clone();
            }

            if (itsLogRho)
            {
               delete itsLogRho;

               itsRho = NULL;
            }

            if (bsmodel->itsLogRho)
            {
               itsLogRho = (ARM_VolCurve *) bsmodel->itsLogRho->Clone();
            }

            if (itsNu)
            {
               delete itsNu;

               itsNu = NULL;
            }            

            if (bsmodel->itsNu)
            {
               itsNu = (ARM_VolCurve *) bsmodel->itsNu->Clone();
            }

            if (itsLogNu)
            {
               delete itsLogNu;

               itsLogNu = NULL;
            }

            if (bsmodel->itsLogNu)
            {
               itsLogNu = (ARM_VolCurve *) bsmodel->itsLogNu->Clone();
            }

            if (itsBeta)
            {
               delete itsBeta;

               itsBeta = NULL;
            }

            if (bsmodel->itsBeta)
            {
               itsBeta = (ARM_VolCurve *) bsmodel->itsBeta->Clone();
            }

            itsSABRWeight = bsmodel->itsSABRWeight;

            itsCalibrationTenor = bsmodel->itsCalibrationTenor;

            if (bsmodel->itsCalibrationPortfolio)
               itsCalibrationPortfolio = bsmodel->itsCalibrationPortfolio;

            itsSABRSig          = bsmodel->itsSABRSig;

            if (itsCalibMatus)
            {
               delete itsCalibMatus;

               itsCalibMatus = NULL;
            }

            if (bsmodel->itsCalibMatus)
               itsCalibMatus = (ARM_Vector *) bsmodel->itsCalibMatus->Clone();

            itsSmoothFlag = bsmodel->itsSmoothFlag;

            itsLambda     = bsmodel->itsLambda;

            itsNorm_SIG   = bsmodel->itsNorm_SIG;
            itsNorm_RHO   = bsmodel->itsNorm_RHO;
            itsNorm_NU    = bsmodel->itsNorm_NU;

            if (itsNormCoefSIG)
            {
               delete itsNormCoefSIG;

               itsNormCoefSIG = NULL;
            }

            if (bsmodel->itsNormCoefSIG)
               itsNormCoefSIG = (ARM_Vector *) bsmodel->itsNormCoefSIG->Clone();

            if (itsNormCoefRHO)
            {
               delete itsNormCoefRHO;

               itsNormCoefRHO = NULL;
            }

            if (bsmodel->itsNormCoefRHO)
               itsNormCoefRHO = (ARM_Vector *) bsmodel->itsNormCoefRHO->Clone();

            if (itsNormCoefNU)
            {
               delete itsNormCoefNU;

               itsNormCoefNU  = NULL;
            }

            if (bsmodel->itsNormCoefNU)
               itsNormCoefNU = (ARM_Vector *) bsmodel->itsNormCoefNU->Clone();


            itsBeginSmoothMatuIndex = bsmodel->itsBeginSmoothMatuIndex;

            itsCalibMeth = bsmodel->itsCalibMeth;

            itsPreviousMaxCalibIndex = bsmodel->itsPreviousMaxCalibIndex;

            strcpy(itsOptimizationDiag, 
                   bsmodel->itsOptimizationDiag);

            itsSigmaOrAlphaInput = bsmodel->itsSigmaOrAlphaInput;

			itsConvToAdjVolInCaseOfAlpha = bsmodel->itsConvToAdjVolInCaseOfAlpha;

            // The BS father model is a shared object so we don't clonne it

            itsFatherBSModel     = bsmodel->itsFatherBSModel;
        }
 
        void Copy(const ARM_Object* srcBSModel)
        {
            ARM_BSModel::Copy(srcBSModel);
 
            BitwiseCopy(srcBSModel);
        }
 
        ARM_Object* Clone(void)
        {
            ARM_BSSmiledModel* theClone = new ARM_BSSmiledModel();
 
 
            theClone->Copy(this);
 
            return(theClone);
        }

        void SetOptimizeDiag(char* diagnoMsg)
        {
            strcpy(itsOptimizationDiag, diagnoMsg);
        }

        ARM_Portfolio* GetCalibrationPortfolio(void)
        {
            return(itsCalibrationPortfolio);
        }

        ARM_VolCurve* GetRho(void)
        {
            return(itsRho);
        }

        ARM_VolCurve* GetNu(void)
        {
            return(itsNu);
        }

        ARM_VolCurve* GetBeta(void)
        {
            return(itsBeta);
        }

		void Set_XXX_Rho(ARM_VolCurve* rho)
        {
            * itsRho->GetVolatilities() = * rho->GetVolatilities();
        }

        void Set_XXX_Nu(ARM_VolCurve* nu)
        {
			* itsNu->GetVolatilities() = * nu->GetVolatilities();
        }


        void Set_XXX_Beta(ARM_VolCurve* beta)
        {
			* itsBeta->GetVolatilities() = * beta->GetVolatilities();
        }

        void SetRho(ARM_VolCurve* rho)
        {
            if (itsRho)
            {
               delete itsRho;

               itsRho = NULL;
            }

            if (rho)
            {
               itsRho = (ARM_VolCurve *) rho->Clone();
            }
        }

        void SetNu(ARM_VolCurve* nu)
        {
            if (itsNu)
            {
               delete itsNu;

               itsNu = NULL;
            }

            if (nu)
            {
               itsNu = (ARM_VolCurve* ) nu->Clone();
            }

        }


        void SetBeta(ARM_VolCurve* beta)
        {
            if (itsBeta)
            {
               delete itsBeta;

               itsBeta = NULL;
            }

            if (beta)
            {
               if ( beta->GetName() != ARM_VOL_FLAT )
               {
                  itsBeta = (ARM_VolCurve *) beta->Clone();
               }
               else
               {
                  ARM_Vector* expiries = NULL;


                  ARM_Vector* yearTermsBeta = NULL;
                  
                  if (GetVolatility()->GetExpiryTerms())
                  {
                     expiries = GetVolatility()->GetExpiryTerms();
                  }
                  else if (itsRho->GetExpiryTerms())
                  {
                     expiries = itsRho->GetExpiryTerms();  
                  }
                  else
                  {
                     expiries = itsNu->GetExpiryTerms();
                  }
                  
                  yearTermsBeta = (ARM_Vector *) expiries->Clone();

                  ARM_Vector* undMatus  = NULL;

                  if (((ARM_VolLInterpol *) GetVolatility())->GetStrikes())
                  {
                     undMatus = ((ARM_VolLInterpol *) GetVolatility())->GetStrikes();
                  }
                  else if (((ARM_VolLInterpol *) itsRho)->GetStrikes())
                  {
                     undMatus = ((ARM_VolLInterpol *) itsRho)->GetStrikes();
                  }
                  else
                  {
                     undMatus = ((ARM_VolLInterpol *) itsNu)->GetStrikes();
                  }

                  ARM_Vector* undMatuBeta = (ARM_Vector *) undMatus->Clone();


                  itsBeta = new ARM_VolLInterpol((ARM_VolFlat *) beta, 
                                                 yearTermsBeta,
                                                 undMatuBeta);
 
               }
            }

            if (itsInitialBeta)
            {
               delete itsInitialBeta;

               itsInitialBeta = NULL;
            }

            if (itsBeta)
               itsInitialBeta = (ARM_VolCurve *) itsBeta->Clone();
        }

        void PrepareNormVectors(ARM_Vector* calibMatus)
        {
            int i;

            double Fwd_Mi;
            double Sig_Mi;
            double tmpVal;

            int sz = calibMatus->GetSize();
 
            if (itsNormCoefSIG)
               delete itsNormCoefSIG;

            if (itsNormCoefSIG)
               delete itsNormCoefSIG;

            if (itsNormCoefNU)
               delete itsNormCoefNU;

            itsNormCoefSIG = new ARM_Vector(sz);
            itsNormCoefRHO = new ARM_Vector(sz);
            itsNormCoefNU  = new ARM_Vector(sz);
   
            for (i = itsBeginSmoothMatuIndex; i < (sz-1); i++)
            {
                Fwd_Mi = CalcFwd(calibMatus->Elt(i))/100.0;
                Sig_Mi = CalcVol(calibMatus->Elt(i))/100.0;

                tmpVal = SmoothF(Fwd_Mi*Sig_Mi);
                itsNormCoefSIG->Elt(i) = SQR(tmpVal);
            }

            double Rho_Mi;

            for (i = itsBeginSmoothMatuIndex; i < (sz-1); i++)
            {
                Rho_Mi = itsRho->ComputeVolatility(calibMatus->Elt(i),
                                             itsCalibrationTenor);

                tmpVal = SmoothG(Rho_Mi);
                itsNormCoefRHO->Elt(i) = SQR(tmpVal);
            }

            double Nu_Mi;

            for (i = itsBeginSmoothMatuIndex; i < (sz-1); i++)
            {
                Nu_Mi = itsNu->ComputeVolatility(calibMatus->Elt(i),
                                                 itsCalibrationTenor);

                tmpVal = SmoothH(Nu_Mi);
                itsNormCoefNU->Elt(i) = SQR(tmpVal);
            }
        }

		double xVal(double z, double Rho_t)
        {
            double res;

            res = log((z-Rho_t+sqrt(1.0-2.0*z*Rho_t+(z*z)))/(1.0-Rho_t));

            return(res);
        }

        virtual double ComputeVol(double maturity, double volTenor, 
                          double fwd, double strike, int mode = -1);

        double CalcSABRVolBetaDiffOneWeighted(double maturity, 
                                              double volTenor,
                                              double fwd, 
                                              double strike);

        double CalcSABRVolBetaDiffOneImpLNVol(double maturity, 
                                              double volTenor,
                                              double fwd, 
                                              double strike);

        double CalcSABRVolBetaDiffOneArith(double maturity, 
                                           double volTenor,
                                           double fwd, 
                                           double strike);

        
        double CalcSABRVolBetaDiffOne(double maturity, 
                                      double volTenor,
                                      double fwd, 
                                      double strike);

        double bsOptionBetaDiffOne(double spot,
                                   double strike,
                                   double volatility,
                                   double dividend,
                                   double discountRate,
                                   double maturity,
                                   double CallPut,
                                   double volTenor = 0.0);


        // The latest approach 

        double CptSABRVolBetaDiffOneNormalExact(double maturity,
                                                double volTenor,
                                                double fwd,
                                                double strike);

        double CptSABRVolBetaDiffOneDirectExact(double maturity,
                                                double volTenor,
                                                double fwd,
                                                double strike);

		
		double ComputeAlpha(double f, double K,
                    double mat,
                    double sigmaATM,
                    double rho, double nu, double beta, 
                    double weight, int type);


        double CptSABRVolBetaDiffOneNormalOrDirectExact(double maturity,
                                                        double volTenor,
                                                        double fwd,
                                                        double strike);

        double bsOptionBetaDiffOneNormalOrDirectExact(double spot,
                                                      double strike,
                                                      double volatility,
                                                      double dividend,
                                                      double discountRate,
                                                      double maturity,
                                                      double CallPut,
                                                      double volTenor = 0.0);
        double bsOption(double spot,
                        double strike,
                        double volatility,
                        double dividend,
                        double discountRate,
                        double maturity,
                        double CallPut,
                        double volTenor = 0.0)
        {
            double optionPrice;

            if ( volatility < 0.0 )
            {
               throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                "ARM_BSSmiledModel::bsOption() : volatility is negative!?");
            }
    
 	        if (IS_ZERO(strike) || ( strike < 0 ))
            {
	           strike = 0.05; // strike = 0.00001;
            }
            
            if (itsBeta)
            {
               double beta = itsBeta->ComputeVolatility(maturity,
                                                        volTenor);
               if (!(IS_ONE(beta)))
               {
#ifdef SABR_APPROX_NORMAL // JDA Approach

                  optionPrice = bsOptionBetaDiffOne(spot,
                                                    strike,
                                                    volatility,
                                                    dividend,
                                                    discountRate,
                                                    maturity,
                                                    CallPut,
                                                    volTenor);
#else

                  optionPrice = bsOptionBetaDiffOneNormalOrDirectExact(spot,
                                                                       strike,
                                                                       volatility,
                                                                       dividend,
                                                                       discountRate,
                                                                       maturity,
                                                                       CallPut,
                                                                       volTenor);

#endif

                  return(optionPrice);
               }
            }

            // Beta == 1
            if ( itsSmiledModelType == K_SABR_IMPLNVOL )
            {
               throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
    "ARM_BSSmiledModel::bsOption() : Implied LN Vol flag not valid in case Beta==1");
            }
            else if ( itsSmiledModelType == K_SABR_WEIGHT )
            {
               throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
    "ARM_BSSmiledModel::bsOption() : Weighted method flag not valid in case Beta==1");
            } 
     
            double sigmaPlus  = 0.0;
            double sigmaMinus = 0.0;


            ARM_Date asOfDate = GetStartDate();

            double delta;
            double Bt = 0.0;
            double At = 0.0;
            double Wt = 0.0;

            if ( volatility < 0.0 )
            {
               throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
              "ARM_BSSmiledModel::bsOption() : volatility is negative!?");
            }

            if ( itsSmiledModelType != K_LD ) // SABR type
            {
               double sigSABR;
               double ATMVol;
               double Alpha_t;


               ATMVol = volatility;

               // Retrieve the data corresponding to the tenor

               double Rho_t = itsRho->ComputeVolatility(maturity, 
                                                        volTenor);

               double Nu_t  = itsNu->ComputeVolatility(maturity,
                                                       volTenor);

               if ( itsSigmaOrAlphaInput == 0 )
               {
                  double alpha_t = volatility;

                  double ATM_VAL = strike*0.01; // spot*0.01, strike*0.01,

                  ATMVol = ComputeSigATMBetaEqOne(ATM_VAL, ATM_VAL,
                                                  maturity, 
                                                  alpha_t,
                                                  Rho_t, Nu_t,
                                                  itsSmiledModelType);
               }

               if (IS_ZERO(Nu_t))
               {
                  sigSABR = ATMVol;

                  itsSABRSig = sigSABR*100.0;
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
                     if ( strike <= 0.0 )
                     {
						if ( CallPut == K_CALL )
						{
                           return(spot);
						}
						else
						{
						   return(0.0);
						}
                     }

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

                     if (( Rho_t == 1.0 ) // Just in case ...
                         ||
                         (!finite(denom)) // The relevant test
                         ||
                         (!finite(KSI))
                        )
                     {
                        throw Exception(__LINE__, __FILE__, ERR_INVALID_DATA,
                        "bsOption() : Numerical problem an NaN or INFINITY encoutered");
                     }
 
                     M = KSI/denom;
                  }

                  sigSABR = ATMVol*M;

                  itsSABRSig = sigSABR*100.0;
               }

               optionPrice = ::bsOption(spot, strike, sigSABR,
                                        dividend, discountRate,
                                        maturity, CallPut);

               return(optionPrice);
            }

            Bt = itsLogRho->ComputeVolatility(maturity, volTenor);
            Bt = exp(Bt);

            At = spot*((1-Bt)/Bt);
            
            Wt = itsLogNu->ComputeVolatility(maturity, volTenor); 
            Wt = exp(Wt);
 
            delta = volatility*sqrt(maturity)*Wt; 


            // Compute sigmaPlus

            sigmaPlus = volatility*(1.0+log(spot/strike)
                                        *(At/(2.0*(spot+At)))) 
                                  *(1.0+(maturity*delta*delta)/8.0)
                                  +delta;

            // Compute sigmaMinus

            sigmaMinus = sigmaPlus-(2.0*delta); 


            optionPrice = 0.5*(::bsOption(spot, strike, sigmaPlus, 
                                          dividend, discountRate, 
                                          maturity, CallPut)
                              +::bsOption(spot, strike, sigmaMinus, 
                                          dividend, discountRate,
                                          maturity, CallPut)
                              );

            return(optionPrice);
        }

        double GetSABRSigma(void)
        {
            return(itsSABRSig);
        }

        ARM_Vector* GenCalibMatus(ARM_Portfolio* pf,
                                  double H_TimeStep);

        void UpdateSigRhoNuByInterpol(double H_TimeStep,
                                      int FitVolOrNot,
                                      int FitRhoOrNot,
                                      int FitNuOrNot,
                                      ARM_Vector* generatedTenors, 
                                      int interpolMeth,
                                      int FitBetaOrNot = 0);

        void ThisIsWhatToFit(int FitVolOrNot,
                             int FitRhoOrNot, 
                             int FitNuOrNot,
                             int FitBetaOrNot = 0);

        // Calibration utilities

        void Calibrate(ARM_Portfolio* pf,
                       int calVolOrNot   = 1, // YES
                       int calRhoOrNot   = 1, // YES
                       int calNuOrNot    = 1, // YES
                       double volTenor   = -1,
                       double H_TimeStep = -1,
                       double minSig     = LOW_INFINITE_BOUND,
                       double maxSig     = UPPER_INFINITE_BOUND,
                       double minRho     = LOW_INFINITE_BOUND,
                       double maxRho     = UPPER_INFINITE_BOUND,
                       double minNu      = LOW_INFINITE_BOUND, 
                       double maxNu      = UPPER_INFINITE_BOUND,
                       int interpolMeth  = K_LINEAR,
                       double tol        = -1,  // Retrieve default tol
                       long maxIter      = ARM_DEF_MAX_ITER,
                       int gradCalc      = 1,  // Yes
                       double lambda     = 0.0,
                       double beginSmoothMatu = -1.0,
                       int FitBetaOrNot  = 0, // NO
                       double minBeta    = LOW_INFINITE_BOUND,
                       double maxBeta    = UPPER_INFINITE_BOUND);

        void BootstrapCalibrate(ARM_Portfolio* pf,
                                int calVolOrNot   = 1, // YES
                                int calRhoOrNot   = 1, // YES
                                int calNuOrNot    = 1, // YES
                                double volTenor   = -1,
                                double H_TimeStep = -1,
                                double minSig     = LOW_INFINITE_BOUND,
                                double maxSig     = UPPER_INFINITE_BOUND,
                                double minRho     = LOW_INFINITE_BOUND,
                                double maxRho     = UPPER_INFINITE_BOUND,
                                double minNu      = LOW_INFINITE_BOUND,
                                double maxNu      = UPPER_INFINITE_BOUND,
                                int interpolMeth  = K_LINEAR,
                                double tol        = -1,  // Retrieve default tol
                                long maxIter      = ARM_DEF_MAX_ITER,
                                int gradCalc      = 1,  // Yes
                                double lambda     = 0.0,
                                double beginSmoothMatu = -1.0,
                                int calBetaOrNot  = 0, // NO,
                                double minBeta    = LOW_INFINITE_BOUND,
                                double maxBeta    = UPPER_INFINITE_BOUND);

        ARM_Vector* BuildParameters(void)
        {
            int paramSizeSig = GetVolatility()->GetExpiryTerms()->GetSize();
            int paramSizeRho = itsRho->GetExpiryTerms()->GetSize();
            int paramSizeNu  = itsNu->GetExpiryTerms()->GetSize();

            int paramSizeBeta = 0;

            if (itsBeta)
               paramSizeBeta = itsBeta->GetExpiryTerms()->GetSize();

            ARM_Vector* builtParams = new ARM_Vector(paramSizeSig
                                                     +paramSizeRho
                                                     +paramSizeNu
                                                     +paramSizeBeta);

            // First ATM Vol, Rho's , next Nu's

            ARM_Vector* ATMTenorVol = GetVolatility()->GetCol(itsCalibrationTenor);

            ARM_Vector* rhoTenor    = itsRho->GetCol(itsCalibrationTenor);

            ARM_Vector* nuTenor     = itsNu->GetCol(itsCalibrationTenor);

            ARM_Vector* betaTenor   = NULL;

            if (itsBeta)
               betaTenor = itsBeta->GetCol(itsCalibrationTenor);

            int i;

            for (i = 0; i < paramSizeSig; i++)
            {
                builtParams->Elt(i) = ATMTenorVol->Elt(i);
            }

            for (i = 0; i < paramSizeRho; i++)
            {
                builtParams->Elt(i+paramSizeSig) = rhoTenor->Elt(i);
            }

            for (i = 0; i < paramSizeNu; i++)
            {
                builtParams->Elt(i+paramSizeSig
                                 +paramSizeRho) = nuTenor->Elt(i);
            }

            if (itsBeta)
            {
               for (i = 0; i < paramSizeBeta; i++)
               {
                   builtParams->Elt(i+paramSizeSig
                                     +paramSizeRho
                                     +paramSizeNu) = betaTenor->Elt(i);
               }
            }

            delete ATMTenorVol;
            delete rhoTenor;
            delete nuTenor;

            if (betaTenor)
               delete betaTenor;

            return(builtParams);
        }
 
        void SetParameters(ARM_Vector* params)
        {
            ARM_Model::SetParameters(params);


            int paramSizeSig = GetVolatility()->GetExpiryTerms()->GetSize();
            int paramSizeRho = itsRho->GetExpiryTerms()->GetSize();
            int paramSizeNu  = itsNu->GetExpiryTerms()->GetSize();

            int paramSizeBeta = 0;

            if (itsBeta)
               paramSizeBeta = itsBeta->GetExpiryTerms()->GetSize();

            
            // Set the model, Real parameters

            ARM_Vector* newParams = GetParameters();

            ARM_Vector volATM(paramSizeSig);
            ARM_Vector rho(paramSizeRho);
            ARM_Vector nu(paramSizeNu);
 
            ARM_Vector beta(paramSizeBeta);


            int i;

            for (i = 0; i < paramSizeSig; i++)
            {
                volATM.Elt(i) = newParams->Elt(i);
            }

            for (i = 0; i < paramSizeRho; i++)
            {
                rho.Elt(i)    = newParams->Elt(i+paramSizeSig);
            }

            for (i = 0; i < paramSizeNu; i++)
            {
                nu.Elt(i)     = newParams->Elt(i+paramSizeSig
                                               +paramSizeRho);
            }

            if (itsBeta)
            {
               for (i = 0; i < paramSizeBeta; i++)
               {
                   beta.Elt(i) = newParams->Elt(i+paramSizeSig
                                                 +paramSizeRho
                                                 +paramSizeNu);
               }
            }

            GetVolatility()->UpdateCol(&volATM, itsCalibrationTenor);
            itsRho->UpdateCol(&rho, itsCalibrationTenor);
            itsNu->UpdateCol(&nu, itsCalibrationTenor);

            if (itsBeta)
            {
               itsBeta->UpdateCol(&beta, itsCalibrationTenor);
            }
        }

        void UpdateLogRhoNu(void)
        {
            ARM_Matrix* logRhoMat = itsLogRho->GetVolatilities();
            ARM_Matrix* logNuMat  = itsLogNu->GetVolatilities();

            int nbLin = logRhoMat->GetNumLines();
            int nbCol = logRhoMat->GetNumCols();
            int i;

            for (i = 0; i < nbLin; i++)
            {
                for (int j = 0; j < nbCol; j++)
                {
                    logRhoMat->Elt(i, j) =
                             log(itsRho->GetVolatilities()->Elt(i, j));
                }
            }

            nbLin = logNuMat->GetNumLines();
            nbCol = logNuMat->GetNumCols();

            for (i = 0; i < nbLin; i++)
            {
                for (int j = 0; j < nbCol; j++)
                {
                    logNuMat->Elt(i, j)  =
                             log(itsNu->GetVolatilities()->Elt(i, j));
                }
            }
        }

        double SmoothF(double x)
        {
            return(x);
        }

        double SmoothG(double x)
        {
            return(x);
        }

        double SmoothH(double x)
        {
            return(x);
        }

        double CalcFwd(double matu, double tenor = -1)
        {
            ARM_ZeroCurve* zc = GetZeroCurve();

            double fwd;
            double fwdDate = matu; 


            double newTenor;

            if ( tenor < 0.0 )
			{
               newTenor = itsCalibrationTenor;
            }
            else
            {
               newTenor = tenor;
            }

            int IsSwapFwd = ( newTenor < 1.0 ) ? 0 : 1;

            if (IsSwapFwd) // A Swap fwd Rate
            {
               fwd = ApproximatedForwardRate(fwdDate,
                                             newTenor,
                                             zc,
                                             IsSwapFwd);
               return(fwd);

            }

			if ( tenor < 0.0 )
			{
               fwd = zc->ForwardYield(fwdDate, 
                                   fwdDate+itsCalibrationTenor, -1);
            }
			else
			{
			   fwd = zc->ForwardYield(fwdDate, 
                                   fwdDate+tenor, -1);
			}

            char* ISOccy = zc->GetCurrencyUnit()->GetCcyName();

            if (( strcmp(ISOccy, "ZAR") != 0 )
                && ( strcmp(ISOccy, "PLN") != 0 )
                && ( strcmp(ISOccy, "BEF") != 0 )
                && ( strcmp(ISOccy, "GBP") != 0 )
                && ( strcmp(ISOccy, "PTE") != 0 )
                && ( strcmp(ISOccy, "AUD") != 0 )
               )
            {
               fwd = fwd*360.0/365.0;
            }

            return(fwd);
        }

        double CalcVol(double matu)
        {
            double vol;


            vol = GetVolatility()->ComputeVolatility(matu, 
                                             itsCalibrationTenor);

            return(vol);
        }
 
        ARM_Vector* SmoothFunction(void)
        {
            int beginSmoothMatuIndex = itsBeginSmoothMatuIndex;

            double LAMBDA = itsLambda;
            double N_SIG  = 1.0;
            double N_RHO  = 1.0;
            double N_NU   = 1.0;

            double F_term = 0.0;
            double G_term = 0.0;
            double H_term = 0.0;

            int sz = itsCalibMatus->GetSize();

            ARM_Vector* smoothVector = new ARM_Vector(3*(sz-2));

            if ( beginSmoothMatuIndex == (sz-1) )
            {
               return(smoothVector);
            }

            double Fwd_Miplus1;
            double Sig_Miplus1;

            double Fwd_Mimoins1;
            double Sig_Mimoins1;

            double Fwd_Mi;
            double Sig_Mi;

            double tmpVal;

            int index = 0;
            int i;

            // Calculate the Sigma Normalization term

            if ( itsCalibMeth == ARM_CALIB_GLOB )
            {
               if ( itsNorm_SIG < 0.0 ) // in fact the first time
               {
                  N_SIG = 0.0;
 
                  for (i = beginSmoothMatuIndex; i < (sz-1); i++)
                  {
                      Fwd_Mi = CalcFwd(itsCalibMatus->Elt(i))/100.0;
                      Sig_Mi = CalcVol(itsCalibMatus->Elt(i))/100.0; 
   
                      tmpVal = SmoothF(Fwd_Mi*Sig_Mi);            
                      N_SIG += SQR(tmpVal);
                  }

                  itsNorm_SIG = N_SIG;
               }
               else
               {
                  N_SIG = itsNorm_SIG;
               }
            }
            else
            {
               N_SIG = 0.0;

               for (i = beginSmoothMatuIndex; i < (sz-1); i++)
               {
                   N_SIG += itsNormCoefSIG->Elt(i);
               }
            }

            if ( N_SIG == 0.0 )
            {
               N_SIG = 1.0;
            }

            // Integrate F
            for (i = beginSmoothMatuIndex; i < (sz-1); i++)
            {
                Fwd_Miplus1  = CalcFwd(itsCalibMatus->Elt(i+1))/100.0;
                Sig_Miplus1  = CalcVol(itsCalibMatus->Elt(i+1))/100.0;

                Fwd_Mimoins1 = CalcFwd(itsCalibMatus->Elt(i-1))/100.0; 
                Sig_Mimoins1 = CalcVol(itsCalibMatus->Elt(i-1))/100.0;

                Fwd_Mi       = CalcFwd(itsCalibMatus->Elt(i))/100.0;
                Sig_Mi       = CalcVol(itsCalibMatus->Elt(i))/100.0;


                double curTerm = SmoothF(Fwd_Miplus1*Sig_Miplus1)
                                 +SmoothF(Fwd_Mimoins1*Sig_Mimoins1)
                                 -2.0*SmoothF(Fwd_Mi*Sig_Mi);

                double denom = itsCalibMatus->Elt(i+1)
                               -itsCalibMatus->Elt(i-1);

                denom = SQR(denom);

                curTerm = curTerm*sqrt(itsCalibMatus->Elt(i+1)
                                           -itsCalibMatus->Elt(i)); 

                curTerm = 4.0*curTerm/denom;

                curTerm = sqrt(LAMBDA/N_SIG)*curTerm;

                smoothVector->Elt(index) = curTerm;

                index++;
            }


            // Integrate G

            double Rho_Miplus1;

            double Rho_Mimoins1;

            double Rho_Mi;


            // Calculate the Rho Normalization term
            if ( itsCalibMeth == ARM_CALIB_GLOB )
            {
               if ( itsNorm_RHO < 0.0 )
               {
                  N_RHO = 0.0;

                  for (i = beginSmoothMatuIndex; i < (sz-1); i++)
                  {
                      Rho_Mi = 
                        itsRho->ComputeVolatility(
                                    itsCalibMatus->Elt(i),
                                    itsCalibrationTenor);

                      tmpVal = SmoothG(Rho_Mi);
                      N_RHO += SQR(tmpVal); 
                  }

                  itsNorm_RHO = N_RHO;
               }
               else
               {
                  N_RHO = itsNorm_RHO;
               }
            }
            else
            {
               N_RHO = 0.0;

               for (i = beginSmoothMatuIndex; i < (sz-1); i++)
               {
                   N_RHO += itsNormCoefSIG->Elt(i); 
               }
            }

            if ( N_RHO == 0.0 )
            {
               N_RHO = 1.0;
            }

            for (i = beginSmoothMatuIndex; i < (sz-1); i++)
            {
                Rho_Miplus1  = 
                  itsRho->ComputeVolatility(itsCalibMatus->Elt(i+1),
                                            itsCalibrationTenor);

                Rho_Mimoins1 = 
                  itsRho->ComputeVolatility(itsCalibMatus->Elt(i-1),
                                            itsCalibrationTenor);

                Rho_Mi       = 
                  itsRho->ComputeVolatility(itsCalibMatus->Elt(i),
                                            itsCalibrationTenor);


                double curTerm = SmoothG(Rho_Miplus1)
                                 +SmoothG(Rho_Mimoins1)
                                 -2.0*SmoothG(Rho_Mi);

                double denom = itsCalibMatus->Elt(i+1)
                               -itsCalibMatus->Elt(i-1);

                denom = SQR(denom);

                curTerm = curTerm*sqrt(itsCalibMatus->Elt(i+1)
                                           -itsCalibMatus->Elt(i));

                curTerm = 4.0*curTerm/denom;

                curTerm = sqrt(LAMBDA/N_RHO)*curTerm;

                smoothVector->Elt(index) = curTerm;

                index++;
            }

            // Integrate H

            double Nu_Miplus1;

            double Nu_Mimoins1;

            double Nu_Mi;


            // Calculate the Nu Normalization term
            if ( itsCalibMeth == ARM_CALIB_GLOB )
            {
               if ( itsNorm_NU < 0.0 )
               {
                  N_NU = 0.0;

                  for (i = beginSmoothMatuIndex; i < (sz-1); i++)
                  {
                      Nu_Mi = 
                      itsNu->ComputeVolatility(itsCalibMatus->Elt(i),
                                               itsCalibrationTenor);

                      tmpVal = SmoothH(Nu_Mi);
                      N_NU  += SQR(tmpVal); 
                  }

                  itsNorm_NU = N_NU;
               }
               else
               {
                  N_NU = itsNorm_NU;
               }
            }
            else
            {
               N_NU = 0.0;

               for (i = beginSmoothMatuIndex; i < (sz-1); i++)
               {
                   N_NU += itsNormCoefNU->Elt(i);
               }
            }

            if ( N_NU == 0.0 )
            {
               N_NU = 1.0;
            }

            for (i = beginSmoothMatuIndex; i < (sz-1); i++)
            {
                Nu_Miplus1  =
                  itsNu->ComputeVolatility(itsCalibMatus->Elt(i+1),
                                            itsCalibrationTenor);

                Nu_Mimoins1 =
                  itsNu->ComputeVolatility(itsCalibMatus->Elt(i-1),
                                            itsCalibrationTenor);

                Nu_Mi       =
                  itsNu->ComputeVolatility(itsCalibMatus->Elt(i),
                                            itsCalibrationTenor);


                double curTerm = SmoothH(Nu_Miplus1)
                                 +SmoothH(Nu_Mimoins1)
                                 -2.0*SmoothH(Nu_Mi);

                double denom = itsCalibMatus->Elt(i+1)
                               -itsCalibMatus->Elt(i-1);

                denom = SQR(denom);

                curTerm = curTerm*sqrt(itsCalibMatus->Elt(i+1)
                                           -itsCalibMatus->Elt(i));

                curTerm = 4.0*curTerm/denom;

                curTerm = sqrt(LAMBDA/N_NU)*curTerm;

                smoothVector->Elt(index) = curTerm;

                index++;
            }

            return(smoothVector);
        }

        void ComputeFunc(double* params, double* FUNC,
                         ARM_Model* eventualPricingModel = NULL);

        void View(char* id = NULL, FILE* ficOut = NULL)
        {
            FILE* fOut;
            char fOutName[200];


            if ( ficOut == NULL )
            {
               ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

               fOut = fopen(fOutName, "w");
            }
            else
            {
               fOut = ficOut;
            }

            fprintf(fOut, "\n ====> BS Smiled Model :\n\n");

            if (itsSmiledModelType == K_LD) 
               fprintf(fOut, "\n ====> Smiled Model Type: LD +  2 vol states \n\n");
            else
            {
               if (itsSmiledModelType == K_SABR_GEO)
                  fprintf(fOut, "\n ====> Smiled Model Type: SABR geometric \n\n");
               else if (itsSmiledModelType == K_SABR_ARITH)
                  fprintf(fOut, "\n ====> Smiled Model Type: SABR arithmetic \n\n");
               else if (itsSmiledModelType == K_SABR_WEIGHT)
                  fprintf(fOut, "\n ====> Smiled Model Type: SABR weight \n\n");
			   else
				   fprintf(fOut, "\n ====> Smiled Model Type: SABR Implied LogNormal Vol \n\n");
            }

            fprintf(fOut, "\n\n>>>>>> CALIBRATION DIAGNOSTIC <<<<<<\n\n");
            fprintf(fOut, "%s \n\n", itsOptimizationDiag);
            fprintf(fOut, ">>>>>> END CALIBRATION DIAGNOSTIC <<<<<<\n\n");

            fprintf(fOut, "\n ====> SIGMA :\n");
            GetVolatility()->View(id, fOut);

            fprintf(fOut, "\n ====> RHO :\n");
            itsRho->View(id, fOut);

            fprintf(fOut, "\n ====> NU :\n");
            itsNu->View(id, fOut);

            fprintf(fOut, "\n ====> BETA :\n");
            if (itsBeta) 
               itsBeta->View(id, fOut);
            else
               fprintf(fOut, "\n ====> NO BETA <> 1");

            if (itsSigmaOrAlphaInput && itsBeta)
            {
               fprintf(fOut, "\n ====> SigmaOrAlphaInput : SIGMA\n");
            }
            else
            {
               if (itsBeta)
               {
                  fprintf(fOut, "\n ====> SigmaOrAlphaInput : ALPHA considered\n");
               }
               else
               {
                  fprintf(fOut, "\n ====> SigmaOrAlphaInput : SIGMA considered\n");
               }
            }

            if ( ficOut == NULL )
            {
               fclose(fOut);
            }
        }

        double EuroCaplet(double Settlement,
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
                          StoreFwdRateAndCapletInfo* StoreInfo);

        ARM_VolCurve* FromSigmaToAlpha(ARM_VolCurve* sigma = NULL, 
                                       double strike = -1.0);
        
        ARM_VolCurve* FromAlphaToSigma(ARM_VolCurve* alpha = NULL,
                                       double strike = -1.0);
};



inline ARM_Portfolio* GenRealPortfolio(ARM_Portfolio* pf)
{
    ARM_Security* curAsset = NULL;

    int size = pf->GetSize();

    int newSize = 0;
    int i;


    for (i = 0; i < size; i++)
    {
        curAsset = pf->GetAsset(i);

        double Wi = pf->GetWeights()->Elt(i);

        if (!( Wi <= 1e-15 )) // Wi == 0
        {
           newSize++;
        }
    }

    if ( newSize == 0 )
    {
       throw Exception(__LINE__, __FILE__, ERR_INVALID_ARGUMENT,
            "All Assets Weights are NULL");

       return(NULL);
    }

    ARM_Vector*    Weights = pf->GetWeights();
    ARM_Vector*    Prices  = pf->GetMktPrices();

    ARM_Portfolio* generatedPf  = new ARM_Portfolio(newSize);

    ARM_Vector*    genWeights   = generatedPf->GetWeights();
    ARM_Vector*    genMktPrices = generatedPf->GetMktPrices();

    int j = 0;

    for (i = 0; i < size; i++)
    {
        curAsset = pf->GetAsset(i);

        double Wi = pf->GetWeights()->Elt(i);

        if (!( Wi <= 1e-15 )) // Wi == 0
        {
           generatedPf->SetAsset(curAsset, j);

           genWeights->Elt(j)   = Weights->Elt(i);

           genMktPrices->Elt(j) = Prices->Elt(i);  

           j++;
        }
    }

    generatedPf->SortByExpiry();

    return(generatedPf);
}



inline ARM_BSSmiledModel* Calibrate(ARM_BSSmiledModel* model, 
                                    ARM_Portfolio* pf,
                                    int calVolOrNot   = 1, // YES
                                    int calRhoOrNot   = 1, // YES
                                    int calNuOrNot    = 1, // YES
                                    double volTenor   = -1,
                                    double H_TimeStep = -1,
                                    double minSig     = LOW_INFINITE_BOUND,
                                    double maxSig     = UPPER_INFINITE_BOUND,
                                    double minRho     = LOW_INFINITE_BOUND,
                                    double maxRho     = UPPER_INFINITE_BOUND,
                                    double minNu      = LOW_INFINITE_BOUND,
                                    double maxNu      = UPPER_INFINITE_BOUND,
                                    int interpolMeth  = K_LINEAR,
                                    double tol = -1,  // Retrieve default tol
                                    long maxIter      = ARM_DEF_MAX_ITER,
                                    int gradCalc      = 1,   // Yes
                                    double lambda     = 0.0,
                                    int globOrBootStrap = 0, // Yes : Glob
                                    double beginSmoothMatu = -1,
                                    int FitBetaOrNot = 0, // NO
                                    double minBeta   = LOW_INFINITE_BOUND,
                                    double maxBeta   = UPPER_INFINITE_BOUND)
{
    ARM_BSSmiledModel* modelClone = (ARM_BSSmiledModel *) model->Clone();

    ARM_Portfolio* realPortfolio = NULL;

    try
    {
        realPortfolio = GenRealPortfolio(pf);

        if (globOrBootStrap) // Global calibration
        { 
           modelClone->SetCalibMeth(ARM_CALIB_GLOB);

           modelClone->Calibrate(realPortfolio, 
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
                                 FitBetaOrNot,
                                 minBeta,
                                 maxBeta);
        }
        else // Calibrate by bootstrapping
        {
           modelClone->SetCalibMeth(ARM_CALIB_BOOTSTRAP);

           modelClone->BootstrapCalibrate(realPortfolio,
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
                                          FitBetaOrNot,
                                          minBeta,
                                          maxBeta);
        }
    }

    catch(Exception& x)
    {
        if (modelClone)
           delete modelClone;

        if (realPortfolio)
           delete realPortfolio;

        // x.DebugPrint();

        throw x;

        return((ARM_BSSmiledModel *) NULL);
    }

    if (realPortfolio)
       delete realPortfolio;

    return(modelClone);
}


/*
ARM_BSSmiledModel* ChangeSABRModelSABRFlagWithBeta(ARM_BSSmiledModel* inSABRMod, 
                                                   int newSABRFlag,
                                                   ARM_VolCurve* newBeta = NULL, // NULL if Beta == 1
                                                   int SigmaOrALphaInput = 1)
{
    ARM_BSSmiledModel* newMod = NULL;


    newMod = new ARM_BSSmiledModel(inSABRMod->GetAsOfDate(),
                                   inSABRMod->GetSpot(),
                                   inSABRMod->GetDividend(),
                                   inSABRMod->GetDiscountCurve(),
                                   inSABRMod->GetVolatility(),
                                   inSABRMod->GetType(),
                                   inSABRMod->GetRho(),
                                   inSABRMod->GetNu(),
                                   newSABRFlag,
                                   newBeta, // NULL if Beta == 1
                                   0.5, // The SABR default Weight
                                   SigmaOrALphaInput,
                                   inSABRMod->GetCorrelManager());


    return(newMod);
}

*/


#endif
/*----------------------------------------------------------------------*/
/*---- End of file ----*/