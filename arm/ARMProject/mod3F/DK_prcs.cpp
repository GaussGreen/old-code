/*
 * $Log: DK_prcs.cpp,v $
 * Revision 1.26  2004/05/18 15:52:37  mab
 * Correction In CALIBRATION
 *
 * Revision 1.25  2004/05/11 14:52:56  mab
 * Just Formatting
 *
 * Revision 1.24  2004/05/10 12:42:06  mcampet
 * MC from SP add PRCS3F_ConvertObjFwdVolToSpotVol
 *
 * Revision 1.23  2004/05/04 13:30:45  mcampet
 * MC from SP modif prcs3f_ConvertObjSpotVolToFwdVol
 *
 * Revision 1.22  2004/04/29 07:26:01  mcampet
 *  MC from SP add HW3F_Calibration
 *
 * Revision 1.21  2004/04/22 11:59:51  mcampet
 * MC for SP updated
 *
 * Revision 1.20  2004/04/14 13:30:59  mcampet
 *  MC update SP release
 * add ImpliedFwdCorrelation_VFDK_HW1To3F
 *
 * Revision 1.19  2004/03/25 16:53:05  mcampet
 * MC integration SwaptionPrice_VFDK_HW1To3F
 *
 * Revision 1.17  2004/02/26 09:48:14  mab
 * Version of : 26 fev 2004
 *
 * Revision 1.16  2004/01/22 11:07:22  mab
 * Correction: ARM_Vector* matuUnder = (ARM_Vector *) volFx->GetStrikes()->Clone();
 * replaced by : ARM_Vector* matuUnder = new ARM_Vector(1, 100.0);
 *
 * Revision 1.15  2004/01/13 18:23:11  mab
 * the FX vol is now in base 100 so we have to divide it by 100!
 *
 * Revision 1.14  2004/01/05 15:29:46  arm
 * In vol conversion replacing GetCol in the vol object by:
 * GetVolatilities->GetCol(0)
 *
 * Revision 1.13  2003/11/13 17:02:48  mab
 * Retrieve the ATM VOl From Matrix struct (GetColumn(0))
 *
 * Revision 1.12  2003/10/21 10:52:40  mab
 * last release from Dimitri
 *
 * Revision 1.11  2003/10/17 12:54:01  mab
 * Take reset dates in account in vol conversion
 *
 * Revision 1.9  2003/10/13 13:55:36  mab
 * new Version
 *
 * Revision 1.8  2003/10/08 13:38:21  mab
 * Systematic conversion of Vol Fx from spot to Fwd (lattice also)
 *
 * Revision 1.7  2003/10/06 09:33:24  mab
 * Better management of Exceptions
 *
 * Revision 1.6  2003/09/23 10:04:38  mab
 * One correction in returning an Exception
 *
 * Revision 1.5  2003/09/11 10:43:45  mab
 * New Version
 *
 * Revision 1.4  2003/08/21 09:21:09  mab
 * last version
 *
 * Revision 1.2  2003/08/01 13:01:36  mab
 * Improvements
 *
 * Revision 1.1  2003/06/30 16:23:44  mab
 * Initial revision
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/*    DK_prcs.cpp : Interface to PRCS DK pricing using ARM objects            */
/*                                                                            */
/*----------------------------------------------------------------------------*/

 

#include "DK_prcs.h"

#include "hw_vfdk_LDHD_lattice.h"

#include "hw_vfdk_analytics.h"

#include "dk_utils.h"

#include "HW3F_Calibration.h"

#include "volflat.h"

#include "zeroint.h"

#include "volint.h"




/*----------------------------------------------------------------------------*/



int Get3FCurrencyId(ARM_Currency* ccy)
{
    char* ccyName = ccy->GetCcyName();

    if ( strcmp(ccyName, "JPY") == 0 )
    {
       return(0);
    }
    else if ( strcmp(ccyName, "USD") == 0 )
    {
       return(1);
    }
    else if ( strcmp(ccyName, "AUD") == 0 )
    {
       return(2);
    }
    else if ( strcmp(ccyName, "EUR") == 0 )
    {
       return(3);
    }
    else
    {
       throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
        "Unexpected currency");       
    }
}



// Function to convert year terms to Julian dates

void ConvertYearTermsToJulianDates(double julianAsOf, ARM_Vector* yearTerms)
{
    int size = yearTerms->GetSize();


    for (int i = 0; i < size; i++)
    {
        (*yearTerms)[i] = (*yearTerms)[i]*365.0+julianAsOf;
    }
}



// Function to convert Julian dates to year terms

void ConvertJulianDatesToYearTerms(double julianAsOf, ARM_Vector* julianDates)
{
    int size = julianDates->GetSize();


    for (int i = 0; i < size; i++)
    {
        (*julianDates)[i] = ((*julianDates)[i]-julianAsOf)/365.0;
    }
}



ARM_Vector* ConvertToARMVector(DKMaille<double>& inVect)
{
    ARM_Vector* newVect = NULL;

    int size = inVect.entries();

    if ( size == 0 )
    {
       return(NULL);
    }
    else
    {
       newVect = new ARM_Vector(size);
       
       for (int i = 0; i < size; i++)
       {
           (*newVect)[i] = inVect[i];
       }
    }

    return(newVect);
}



void FromARMVectorToDKMaille(ARM_Vector* vect, DKMaille<double>& maille)
{
    int sz = vect->GetSize();


    maille.resize(sz);

    for (int i = 0; i < sz; i++)
    {
        maille[i] = vect->Elt(i);
    }
}



void FromARMVolToDKMaille2D(double julianAsOf, 
                            ARM_VolLInterpol* volCrv,
                            DKMaille2D<double>& vect2D, 
                            DKMaille<double>& X, 
                            DKMaille<double>& Y)
{
    vect2D.resize(volCrv->GetExpiryTerms()->GetSize(), 
                  volCrv->GetStrikes()->GetSize());

    ARM_Matrix* matVol = volCrv->GetVolatilities();

    
    for (int i = 0; i < matVol->GetNumLines(); i++)
    {
        for (int j = 0; j < matVol->GetNumCols(); j++)
        {
            vect2D(i, j) = matVol->Elt(i, j);
        }
    }
   

    ARM_Vector matusOp  = *volCrv->GetExpiryTerms();
    ARM_Vector matusUnd = *volCrv->GetStrikes();

    // ConvertYearTermsToJulianDates(julianAsOf, &matusOp);
    FromARMVectorToDKMaille(&matusOp, X);

    // ConvertYearTermsToJulianDates(julianAsOf, &matusUnd);
    FromARMVectorToDKMaille(&matusUnd, Y);
}





/*----------------------------------------------------------------------------*/


ARM_Vector* PRCS3F_Lattice_HWVFDK_Pricing(ARM_Vector* dLatticeGeometryData,
                                          double dNumTimeLinesBeforeFirstNotice,
                    double dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                                          double dNumTimeLines, 
                                          ARM_Date& AsOfDate, 
                                          ARM_ZeroCurve* dBaseYieldCurve,
                                          ARM_ZeroCurve* dForeignYieldCurve,
                                          ARM_ZeroCurve* dBaseRatesNoBasisCurve,
                                          ARM_ZeroCurve* dForeignRatesNoBasisCurve,
                                          ARM_VolCurve* volSwopBase,
                                          ARM_VolCurve* volSwopForeign,
                                          ARM_VolCurve* volFx,
                                          ARM_Vector* dRedemptionData,
                                          double dStrike,
                                          double dType,
                                          ARM_Date& dOptionExpiry, 
                                          double dMeanReversionBase, 
                                          double dMeanReversionForeign,  
                                          double dSpotFX,
                                          ARM_VolCurve* dBaseForeignCorrelation,
                                          ARM_VolCurve* dBaseSpotFXCorrelation,
                                          ARM_VolCurve* dForeignSpotFXCorrelation,
                                          double dProductModelCode,
                                          double dX1Limit,
                                          double dX2Limit,
                                          double dX3Limit,
                                          double dI1Limit,
                                          double dI2Limit,
                                          double dI3Limit,
                                          double dADPLimit,
                                          double dOptimal,
                                          double dTimeBoost,
                                          double dDeltaFlag,
                                          double smoothingValue,
                                          long calcProbaSurvOrNot,
                                          double QBaseSmile,
                                          double QForeignSmile,
                                          double CutOff,
                                          double LongDatedSpotFxVol,
                                          double dStringModel,
                                          ARM_Matrix* dBoosterData,
                                          int CalibSwoptBasis)
{
     ARM_Vector* res = NULL;


     DKMaille<double>   dStdDevBaseX;
     DKMaille<double>   dStdDevBaseY;
     DKMaille2D<double> dStdDevBaseZ;

     DKMaille<double>   dStdDevForeignX;
     DKMaille<double>   dStdDevForeignY;
     DKMaille2D<double> dStdDevForeignZ;

   
     double dSmoothing = smoothingValue;


     // Interface DK function

     double dSpotDate = AsOfDate.GetJulian();

     DKMaille<double> pricingRes;

     unsigned vectSize = 0;
     double*  vectElements = NULL;
   
     
     vectSize     = dLatticeGeometryData->GetSize();
     vectElements = dLatticeGeometryData->GetElt();
     DKMaille<double> dLatticeGeometryDataConv(vectSize, vectElements);

     ARM_Vector vectCopy = *(dBaseYieldCurve->GetYearTerms());
     ConvertYearTermsToJulianDates(dSpotDate, &vectCopy);
     vectSize     = vectCopy.GetSize();
     vectElements = vectCopy.GetElt();
     DKMaille<double> dBaseDates(vectSize, vectElements);

     vectSize     = dBaseYieldCurve->GetDiscountFactors()->GetSize();
     vectElements = dBaseYieldCurve->GetDiscountFactors()->GetElt();
     ARM_Vector* baseYearTerms = dBaseYieldCurve->GetYearTerms();
     DKMaille<double> dBaseRates(vectSize, vectElements);


     vectCopy = *(dForeignYieldCurve->GetYearTerms());
     ARM_Vector* forYearTerms = dForeignYieldCurve->GetYearTerms();
     ConvertYearTermsToJulianDates(dSpotDate, &vectCopy);
     vectSize     = vectCopy.GetSize();
     vectElements = vectCopy.GetElt();
     DKMaille<double> dForeignDates(vectSize, vectElements);

     vectSize     = dForeignYieldCurve->GetDiscountFactors()->GetSize();
     vectElements = dForeignYieldCurve->GetDiscountFactors()->GetElt();
     DKMaille<double> dForeignRates(vectSize, vectElements);

     // Curves without Basis
    
/* OLD CODE
     vectSize     = dBaseRatesNoBasisCurve->GetDiscountFactors()->GetSize();
     vectElements = dBaseRatesNoBasisCurve->GetDiscountFactors()->GetElt();
*/
     ARM_Vector newBaseDFNoBasis(baseYearTerms->GetSize());
     int k;
     for (k = 0; k < baseYearTerms->GetSize(); k++)
     {
         newBaseDFNoBasis.Elt(k) = 
                dBaseRatesNoBasisCurve->DiscountPrice(baseYearTerms->Elt(k));
     }
     
     vectSize = baseYearTerms->GetSize();
     vectElements = newBaseDFNoBasis.GetElt();
     DKMaille<double> dBaseRatesNoBasis(vectSize, vectElements);


     ARM_Vector newForDFNoBasis(forYearTerms->GetSize());

     for (k = 0; k < forYearTerms->GetSize(); k++)
     {
         newForDFNoBasis.Elt(k) = 
              dForeignRatesNoBasisCurve->DiscountPrice(forYearTerms->Elt(k));
     }
/* OLD CODE
     vectSize     = dForeignRatesNoBasisCurve->GetDiscountFactors()->GetSize();
     vectElements = dForeignRatesNoBasisCurve->GetDiscountFactors()->GetElt();
*/

     
     vectSize = forYearTerms->GetSize();
     vectElements = newForDFNoBasis.GetElt();
     DKMaille<double> dForeignRatesNoBasis(vectSize, vectElements);

     // End Curves without Basis


     vectSize     = dRedemptionData->GetSize();
     vectElements = dRedemptionData->GetElt();


     DKMaille<double> dRedemptionDataConv(vectSize, vectElements);

     DKMaille2D<double> dBoosterDataConv(dBoosterData->GetNumLines(), 
                                         dBoosterData->GetNumCols());

     for (int i = 0; i < dBoosterData->GetNumLines(); i++)
     {
         for (int j = 0; j < dBoosterData->GetNumCols(); j++)
         {
             dBoosterDataConv(i, j) = dBoosterData->Elt(i, j);
         }
     }

    
     // The correls are objects but similar to flat vol object

     double dBaseForeignCorrelationConv   = 0.0;
     if ( dBaseForeignCorrelation->GetName() != ARM_VOL_FLAT )
     {
        throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
        "Only Flat objects are handled so far for BaseForeignCorrelation");
     }
     else
     {
        dBaseForeignCorrelationConv = ((ARM_VolFlat *) 
                                dBaseForeignCorrelation)->GetVolatility();
     }

     
     double dBaseSpotFXCorrelationConv    = 0.0;
     if ( dBaseSpotFXCorrelation->GetName() != ARM_VOL_FLAT )
     {
        throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                 "Only Flat objects are handled so far for BaseSpotFXCorrelation");
     }
     else
     {
        dBaseSpotFXCorrelationConv = 
                          ((ARM_VolFlat *) dBaseSpotFXCorrelation)->GetVolatility();
     }

     double dForeignSpotFXCorrelationConv = 0.0;
     if ( dForeignSpotFXCorrelation->GetName() != ARM_VOL_FLAT )
     {
        throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
              "Only Flat objects are handled so far for ForeignSpotFXCorrelation");
     }
     else
     {
        dForeignSpotFXCorrelationConv = 
               ((ARM_VolFlat *) dForeignSpotFXCorrelation)->GetVolatility();
     }
     

     FromARMVolToDKMaille2D(dSpotDate, (ARM_VolLInterpol *) volSwopBase, 
                            dStdDevBaseZ, dStdDevBaseX, dStdDevBaseY);
                            

     FromARMVolToDKMaille2D(dSpotDate, (ARM_VolLInterpol *) volSwopForeign, 
                            dStdDevForeignZ, dStdDevForeignX, dStdDevForeignY);

     DKMaille<double> dSpotFXVolDatesTD;
     DKMaille<double> dSpotFXVolTD;

     ARM_Vector vectDates = *volFx->GetExpiryTerms();
     // ConvertYearTermsToJulianDates(dSpotDate, &vectDates);
     FromARMVectorToDKMaille(&vectDates, dSpotFXVolDatesTD);

     ARM_Vector* vectValues = volFx->GetVolatilities()->GetColumn(0);

     // The inputs are in base 100
     ARM_Vector vectValFlagAndVol(vectValues->GetSize()+2);

     for (int h = 0; h < vectValFlagAndVol.GetSize()-2; h++)
     {
         vectValFlagAndVol.Elt(h) = vectValues->Elt(h)/100.0;
     }

     vectValFlagAndVol.Elt(vectValFlagAndVol.GetSize()-2) = (double) CutOff;
     vectValFlagAndVol.Elt(vectValFlagAndVol.GetSize()-1) = LongDatedSpotFxVol;

     FromARMVectorToDKMaille(&vectValFlagAndVol, dSpotFXVolTD);

     if (vectValues)
        delete vectValues;

     try
     {
         pricingRes = Lattice_HWVFDK_3F_(dLatticeGeometryDataConv,
                                         dNumTimeLinesBeforeFirstNotice,
                        dNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                                         dNumTimeLines, 
                                         dSpotDate, 
                                         dBaseDates,
                                         dBaseRates,
                                         dForeignDates,
                                         dForeignRates,
                                         dRedemptionDataConv,
                                         dStrike,
                                         dType,
                                         dOptionExpiry.GetJulian(),
                                         dStdDevBaseX, 
                                         dStdDevBaseY,
                                         dStdDevBaseZ,
                                         dMeanReversionBase,  
                                         dStdDevForeignX, 
                                         dStdDevForeignY,
                                         dStdDevForeignZ,
                                         dMeanReversionForeign,  
                                         dSpotFX,  
                                         dSpotFXVolDatesTD,
                                         dSpotFXVolTD,
                                         dBaseForeignCorrelationConv,
                                         dBaseSpotFXCorrelationConv,
                                         dForeignSpotFXCorrelationConv,
                                         dProductModelCode,
                                         dX1Limit,
                                         dX2Limit,
                                         dX3Limit,
                                         dI1Limit,
                                         dI2Limit,
                                         dI3Limit,
                                         dADPLimit,
                                         dOptimal,
                                         dTimeBoost,
                                         dDeltaFlag,
                                         dStringModel,
                                         dSmoothing,
                                         (double) calcProbaSurvOrNot,
                                         dBoosterDataConv,
                                         dBaseRatesNoBasis,
                                         dForeignRatesNoBasis,
                                         QBaseSmile,
                                         QForeignSmile,
                                         CalibSwoptBasis);
     }

     catch(char* a3FException)
     {
         char buf[255];

         strcpy(buf, "3F Model: ");

         strcat(buf, a3FException);

         throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL, buf);
     }

     catch(Exception& theExpt)
     {
         throw theExpt;
     }

     catch(...)
     {
         throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                           "3Factors model pricing Failed!");
     }

     // Convert pricingRes from DKMaille<double> to an ARM_Vector
      
     res = ConvertToARMVector(pricingRes);

     if ( res == NULL )
     {
        throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                         "PRCS pricing has failed!?");
     }

     return(res);
}



/*----------------------------------------------------------------------------*/
/*                                                                            */
/*   Convert a Spot volatility to a Forward Volatility                        */
/*                                                                            */
/*----------------------------------------------------------------------------*/



ARM_Vector* PRCS3F_FromSpotVolToFwdVol(ARM_Date& AsOfDate,
                                       ARM_ZeroCurve* dBaseRatesNoBasisCurve,
                                       ARM_ZeroCurve* dBaseYieldCurve,
                                       ARM_ZeroCurve* dForeignRatesNoBasisCurve,
                                       ARM_ZeroCurve* dForeignYieldCurve,
                                       ARM_VolCurve*  volSwopBase,
                                       ARM_VolCurve*  volSwopForeign,
                                       ARM_VolCurve*  volFx, // This is a Spot Vol!.
                                       double dMeanReversionBase, 
                                       double dMeanReversionForeign,
                                       ARM_VolCurve* dBaseSpotFXCorrelation,
                                       ARM_VolCurve* dForeignSpotFXCorrelation,
                                       ARM_VolCurve* dBaseForeignCorrelation,
                                       ARM_Vector* NoticeDates,
                                       ARM_Vector* FXCouponResetDates,
                                       ARM_Vector* FXCouponPaymentDates,
                                       int    CutOff,
                                       double LongDatedSpotFxVol,
                                       int CalibSwoptBasis,
                                       ARM_Vector* dForwardVolDates)
{
     ARM_Vector* res = NULL;

     //double dSwaptionCalibrationDate = 30.0; // For Now

     // Dates Vectors
     DKMaille<double>   dStartDates;
     DKMaille<double>   dEndDates;
     DKMaille<double>   dForwardMaturityDates;

     DKMaille<double>   dStdDevBaseX;
     DKMaille<double>   dStdDevBaseY;
     DKMaille2D<double> dStdDevBaseZ;

     DKMaille<double>   dStdDevForeignX;
     DKMaille<double>   dStdDevForeignY;
     DKMaille2D<double> dStdDevForeignZ;


     // Interface DK Fwd Volatility function

     double dSpotDate = AsOfDate.GetJulian();


     unsigned vectSize = 0;
     double*  vectElements = NULL;
 

     // Curves with basis adjustment
    
     ARM_Vector vectCopy = *(dBaseYieldCurve->GetYearTerms());
     ConvertYearTermsToJulianDates(dSpotDate, &vectCopy);
     vectSize     = vectCopy.GetSize();
     vectElements = vectCopy.GetElt();
     DKMaille<double> dBaseDates(vectSize, vectElements);

     vectSize     = dBaseYieldCurve->GetDiscountFactors()->GetSize();
     vectElements = dBaseYieldCurve->GetDiscountFactors()->GetElt();
     ARM_Vector* baseYearTerms = dBaseYieldCurve->GetYearTerms();
     DKMaille<double> dBaseRates(vectSize, vectElements);


     vectCopy = *(dForeignYieldCurve->GetYearTerms());
     ARM_Vector* forYearTerms = dForeignYieldCurve->GetYearTerms();
     ConvertYearTermsToJulianDates(dSpotDate, &vectCopy);
     vectSize     = vectCopy.GetSize();
     vectElements = vectCopy.GetElt();
     DKMaille<double> dForeignDates(vectSize, vectElements);

     vectSize     = dForeignYieldCurve->GetDiscountFactors()->GetSize();
     vectElements = dForeignYieldCurve->GetDiscountFactors()->GetElt();
     DKMaille<double> dForeignRates(vectSize, vectElements);


     // Curves without Basis

     ARM_Vector newBaseDFNoBasis(baseYearTerms->GetSize());
     int k;
     for (k = 0; k < baseYearTerms->GetSize(); k++)
     {
         newBaseDFNoBasis.Elt(k) = 
              dBaseRatesNoBasisCurve->DiscountPrice(baseYearTerms->Elt(k));
     }
     
     vectSize     = baseYearTerms->GetSize();
     vectElements = newBaseDFNoBasis.GetElt();

     DKMaille<double> dBaseRatesNoBasis(vectSize, vectElements);


     ARM_Vector newForDFNoBasis(forYearTerms->GetSize());

     for (k = 0; k < forYearTerms->GetSize(); k++)
     {
         newForDFNoBasis.Elt(k) = 
                 dForeignRatesNoBasisCurve->DiscountPrice(forYearTerms->Elt(k));
     }
     
     vectSize = forYearTerms->GetSize();
     vectElements = newForDFNoBasis.GetElt();
     DKMaille<double> dForeignRatesNoBasis(vectSize, vectElements);

     // End Curves without Basis


    
     // The correls are objects but similar to flat vol object

     double dBaseForeignCorrelationConv   = 0.0;
     if ( dBaseForeignCorrelation->GetName() != ARM_VOL_FLAT )
     {
        throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
        "Only Flat objects are handled so far for BaseForeignCorrelation");
     }
     else
     {
        dBaseForeignCorrelationConv = 
                 ((ARM_VolFlat *) dBaseForeignCorrelation)->GetVolatility();
     }

     
     double dBaseSpotFXCorrelationConv    = 0.0;
     if ( dBaseSpotFXCorrelation->GetName() != ARM_VOL_FLAT )
     {
        throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
                 "Only Flat objects are handled so far for BaseSpotFXCorrelation");
     }
     else
     {
        dBaseSpotFXCorrelationConv = 
                   ((ARM_VolFlat *) dBaseSpotFXCorrelation)->GetVolatility();
     }

     double dForeignSpotFXCorrelationConv = 0.0;
     if ( dForeignSpotFXCorrelation->GetName() != ARM_VOL_FLAT )
     {
        throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
           "Only Flat objects are handled so far for ForeignSpotFXCorrelation");
     }
     else
     {
        dForeignSpotFXCorrelationConv = 
                   ((ARM_VolFlat *) dForeignSpotFXCorrelation)->GetVolatility();
     }
     

     FromARMVolToDKMaille2D(dSpotDate, (ARM_VolLInterpol *) volSwopBase, 
                            dStdDevBaseZ, dStdDevBaseX, dStdDevBaseY);
                            

     FromARMVolToDKMaille2D(dSpotDate, (ARM_VolLInterpol *) volSwopForeign, 
                            dStdDevForeignZ, dStdDevForeignX, dStdDevForeignY);

     DKMaille<double> dSpotFXVolDatesTD;
     DKMaille<double> dSpotFXVolTD;

     ARM_Vector vectDates = *volFx->GetExpiryTerms();

     ARM_Vector zeroVect(dForwardVolDates->GetSize());

     FromARMVectorToDKMaille(&vectDates, dSpotFXVolDatesTD);

     // Obtain FX year terms vectors
     
     FromARMVectorToDKMaille(&zeroVect, dStartDates);

     vectDates = *dForwardVolDates;

     ConvertJulianDatesToYearTerms(dSpotDate, &vectDates);

     FromARMVectorToDKMaille(&vectDates, dEndDates);
     FromARMVectorToDKMaille(&vectDates, dForwardMaturityDates);


     // OLD CODE : ARM_Vector* vectValues = volFx->GetCol(1.0);

     ARM_Vector* vectValues = volFx->GetVolatilities()->GetColumn(0);

     ARM_Vector vectValFlagAndVol(vectValues->GetSize()+2);

     for (int h = 0; h < vectValFlagAndVol.GetSize()-2; h++)
     {
         vectValFlagAndVol.Elt(h) = vectValues->Elt(h)/100.0;
     }

     vectValFlagAndVol.Elt(vectValFlagAndVol.GetSize()-2) = (double) CutOff;
     vectValFlagAndVol.Elt(vectValFlagAndVol.GetSize()-1) = LongDatedSpotFxVol;

     FromARMVectorToDKMaille(&vectValFlagAndVol, dSpotFXVolTD);

     if (vectValues)
        delete vectValues;

     DKMaille<double> dNoticeDates;
     FromARMVectorToDKMaille(NoticeDates, dNoticeDates);
      
     DKMaille<double> dFXCouponResetDates;
     FromARMVectorToDKMaille(FXCouponResetDates, dFXCouponResetDates);

     DKMaille<double> dFXCouponPaymentDates;
     FromARMVectorToDKMaille(FXCouponPaymentDates, dFXCouponPaymentDates);

     DKMaille<double> fwdVol; 
      
     double dCurrencyPair1;
     double dCurrencyPair2;

     try 
     {
         dCurrencyPair1 = Get3FCurrencyId(dBaseRatesNoBasisCurve->GetCurrencyUnit());
         dCurrencyPair2 = Get3FCurrencyId(dForeignRatesNoBasisCurve->GetCurrencyUnit());

         fwdVol = FromSpotToForwardVol3F(dSpotDate, 
                                         dStartDates,
                                         dEndDates,
                                         dForwardMaturityDates,
                                         dBaseDates,
                                         dBaseRatesNoBasis,
                                         dBaseRates,
                                         dForeignDates,
                                         dForeignRatesNoBasis,
                                         dForeignRates,
                                         dStdDevBaseX, 
                                         dStdDevBaseY,
                                         dStdDevBaseZ,
                                         dStdDevForeignX, 
                                         dStdDevForeignY,
                                         dStdDevForeignZ,
                                         dSpotFXVolDatesTD,
                                         dSpotFXVolTD,
                                         dMeanReversionBase,
                                         dMeanReversionForeign,
                                         dBaseSpotFXCorrelationConv,
                                         dForeignSpotFXCorrelationConv,
                                         dBaseForeignCorrelationConv,
                                         dCurrencyPair1,
                                         dCurrencyPair2,
                                         dNoticeDates,
                                         dFXCouponResetDates,
                                         dFXCouponPaymentDates,
                                         CalibSwoptBasis);
     }

     catch(char* a3FException)
     {
         char buf[200];

         strcpy(buf, "3F Model Vol Conv: ");

         strcat(buf, a3FException);

         throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL, buf);
     }

     catch(Exception& anExpt)
     {
         throw(anExpt);
     }

     catch(...)
     {
         throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                         "Converting Spot Vol To Fwd Vol Failed!");
     }

     // Convert fwdVol to ARM_Vector

     ARM_Vector* resVect = ConvertToARMVector(fwdVol);

     int sz = resVect->GetSize();

     for (int i = 0; i < sz; i++)
     {
         resVect->Elt(i) *= 100.0;
     }

     return(resVect);
}



/*----------------------------------------------------------------------------*/


ARM_VolLInterpol* PRCS3F_ConvertObjSpotVolToFwdVol(ARM_Date& AsOfDate,
                                              ARM_ZeroCurve* dBaseRatesNoBasisCurve,
                                              ARM_ZeroCurve* dBaseYieldCurve,
                                              ARM_ZeroCurve* dForeignRatesNoBasisCurve,
                                              ARM_ZeroCurve* dForeignYieldCurve,
                                              ARM_VolCurve*  volSwopBase,
                                              ARM_VolCurve*  volSwopForeign,
                      ARM_VolLInterpol*  volFx, // This is the inputed Spot Vol!.
                                              double dMeanReversionBase, 
                                              double dMeanReversionForeign,
                                              ARM_VolCurve* dBaseSpotFXCorrelation,
                                              ARM_VolCurve* dForeignSpotFXCorrelation,
                                              ARM_VolCurve* dBaseForeignCorrelation,
                                              ARM_Vector* NoticeDates,
                                              ARM_Vector* FXCouponResetDates,
                                              ARM_Vector* FXCouponPaymentDate,
                                              double CutOff,
                                              double LongDatedSpotFxVol,
                                              int CalibSwoptBasis,
                                              ARM_Vector* ForwardVolDates)
{
    ARM_VolLInterpol* fwdVol = NULL;

    int sz = volFx->GetExpiryTerms()->GetSize();


    double julAsOf = AsOfDate.GetJulian();
    int i;

    if ( ForwardVolDates == NULL )
        ForwardVolDates = FXCouponResetDates;


    ARM_Vector* fwdVolVect = NULL;
    
    try
    {
         fwdVolVect = PRCS3F_FromSpotVolToFwdVol(AsOfDate,
                                                 dBaseRatesNoBasisCurve,
                                                 dBaseYieldCurve,
                                                 dForeignRatesNoBasisCurve,
                                                 dForeignYieldCurve,
                                                 volSwopBase,
                                                 volSwopForeign,
                                                 volFx, // This is a Spot Vol!.
                                                 dMeanReversionBase, 
                                                 dMeanReversionForeign,
                                                 dBaseSpotFXCorrelation,
                                                 dForeignSpotFXCorrelation,
                                                 dBaseForeignCorrelation,
                                                 NoticeDates,
                                                 FXCouponResetDates,
                                                 FXCouponPaymentDate,
                                                 CutOff,
                                                 LongDatedSpotFxVol,
                                                 CalibSwoptBasis,
                                                 ForwardVolDates);
    }
    
    catch(Exception& anExpt)
    {
        throw(anExpt);
    }
        
    ARM_Vector* tmpVolMatu = new ARM_Vector(ForwardVolDates->GetSize());

    for (i = 0; i < ForwardVolDates->GetSize(); i++)
    {
        tmpVolMatu->Elt(i) = (ForwardVolDates->Elt(i)-julAsOf)/365.0;
    }

    ARM_Vector* matuOp    = tmpVolMatu;

    ARM_Matrix* resFwdVol  = new ARM_Matrix(ForwardVolDates->GetSize(), 1);

    for (i = 0; i < ForwardVolDates->GetSize(); i++)
    {
        for (int j = 0; j < 1; j++)
        {
            resFwdVol->Elt(i, j) =  fwdVolVect->Elt(i);               
        }
    }


    ARM_Vector* matuUnder = new ARM_Vector(1, 100.0);

    fwdVol = new ARM_VolLInterpol(AsOfDate, matuOp, matuUnder, resFwdVol);    
    
    if (fwdVolVect)
       delete fwdVolVect;

    return(fwdVol);
}





ARM_Vector* PRCS3F_Bootstrapping_VFDK_HW1To3F(ARM_VolCurve* volCurve,
		 									 ARM_ZeroCurve* zeroCurve,
											 ARM_Vector* NoticeDates,
                                             ARM_Vector* StartDates,
                                             ARM_Vector* EndDates,
                                             ARM_Vector* HW3FParams,
                                             ARM_Date& AsOfDate)
{
    DKMaille<double>   dStdDevX;
    DKMaille<double>   dStdDevY;
    DKMaille2D<double> dStdDevZ;

    // Dates Vectors
    DKMaille<double>   dNoticeDates;
    DKMaille<double>   dStartDates;
    DKMaille<double>   dEndDates;
     
    DKMaille<double>   dHW3FParams;

    double dSpotDate = AsOfDate.GetJulian();


    FromARMVolToDKMaille2D(dSpotDate, (ARM_VolLInterpol *) volCurve, 
                           dStdDevZ, dStdDevX, dStdDevY);


    unsigned vectSize = 0;
    double*  vectElements = NULL;
 

    // Curves with basis adjustment
    
    ARM_Vector vectCopy = *(zeroCurve->GetYearTerms());
    ConvertYearTermsToJulianDates(dSpotDate, &vectCopy);
    vectSize     = vectCopy.GetSize();
    vectElements = vectCopy.GetElt();
    DKMaille<double> dDates(vectSize, vectElements);

    vectSize     = zeroCurve->GetDiscountFactors()->GetSize();
    vectElements = zeroCurve->GetDiscountFactors()->GetElt();
    ARM_Vector* baseYearTerms = zeroCurve->GetYearTerms();
    DKMaille<double> dRates(vectSize, vectElements);


    FromARMVectorToDKMaille(NoticeDates, dNoticeDates);
    FromARMVectorToDKMaille(StartDates, dStartDates);
    FromARMVectorToDKMaille(EndDates, dEndDates);
    FromARMVectorToDKMaille(HW3FParams, dHW3FParams);

    DKMaille<double> dSigma;

    try 
    {
        dSigma = Bootstrapping_VFDK_HW1To3F(dStdDevZ,
                                            dStdDevX,
                                            dStdDevY,
                                            dDates,
                                            dRates,
                                            dRates,
                                            dNoticeDates,
                                            dStartDates,
                                            dEndDates,
                                            dHW3FParams,
                                            dSpotDate);
        int dSize = dSigma.entries();

         double lastvalue = dSigma.at(dSize-1);
         dSigma.insert(lastvalue);
         dSigma.insert(lastvalue);
     }

     catch(char* a3FException)
     {
         char buf[200];

         strcpy(buf, "3F Model Bootstrapping Vol Swaption: ");

         strcat(buf, a3FException);

         throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL, buf);
     }

     catch(Exception& anExpt)
     {
         throw(anExpt);
     }

     catch(...)
     {
         throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                         "Bootstrapping Vol swaption failed!");
     }

     // Convert fwdVol to ARM_Vector

     ARM_Vector* resVect = ConvertToARMVector(dSigma);

     return(resVect);    
}



double PRCS3F_SwaptionPrice_VFDK_HW1To3F(double dSwaptionExpiryInYears,
                                         double dSwaptionTenorInYears,
                                         double dNoticePeriodInDays,
                                         double dStrike,
                                         double dCallPut,
                                         ARM_ZeroCurve* zeroCurve,
                                         ARM_Vector* NoticeDates,
                                         ARM_Vector* Sigma,
                                         ARM_Vector* HW3FParams,
                                         ARM_Date& ObservationDate)
{
    // Dates Vectors
    DKMaille<double>   dNoticeDates;
    DKMaille<double>   dSigma;
    DKMaille<double>   dHW3FParams;

    double dObservationDate = ObservationDate.GetJulian();


    unsigned vectSize = 0;
    double*  vectElements = NULL;
 
    // Curves with basis adjustment
    
    ARM_Vector vectCopy = *(zeroCurve->GetYearTerms());
    ConvertYearTermsToJulianDates(dObservationDate, &vectCopy);
    vectSize     = vectCopy.GetSize();
    vectElements = vectCopy.GetElt();
    DKMaille<double> dDates(vectSize, vectElements);

    vectSize     = zeroCurve->GetDiscountFactors()->GetSize();
    vectElements = zeroCurve->GetDiscountFactors()->GetElt();
    ARM_Vector* baseYearTerms = zeroCurve->GetYearTerms();
    DKMaille<double> dRates(vectSize, vectElements);

    FromARMVectorToDKMaille(NoticeDates, dNoticeDates);
    FromARMVectorToDKMaille(Sigma, dSigma);
    FromARMVectorToDKMaille(HW3FParams, dHW3FParams);

    double Price;
        
    try 
    {
        Price = SwaptionPrice_VFDK_HW1To3F(dSwaptionExpiryInYears,
                                           dSwaptionTenorInYears,
                                           dNoticePeriodInDays,
                                           dStrike,
                                           dCallPut,
                                           dDates,
                                           dRates,
                                           dRates,
                                           dNoticeDates,
                                           dSigma,
                                           dHW3FParams,
                                           dObservationDate);
    }

    catch(char* a3FException)
    {
        char buf[200];

        strcpy(buf, "3F Model Swaption Price: ");

        strcat(buf, a3FException);

        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL, buf);
    }

    catch(Exception& anExpt)
    {
        throw(anExpt);
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                        "Swaption Price with HW3F failed!");
    }

    return(Price);    
}



double PRCS3F_ImpliedFwdCorrelation_VFDK_HW1To3F(double dSwaptionExpiryInYears,
                                                 double dSwaptionTenorInYears,
                                                 double dSwaptionTenor2InYears,
                                                 double dNoticePeriodInDays,
                                                 ARM_ZeroCurve* zeroCurve,
                                                 ARM_Vector* NoticeDates,
                                                 ARM_Vector* Sigma,
                                                 ARM_Vector* HW3FParams,
                                                 ARM_Date& ObservationDate)
{
    // Dates Vectors
    DKMaille<double>   dNoticeDates;
    DKMaille<double>   dSigma;
    DKMaille<double>   dHW3FParams;

    double dObservationDate = ObservationDate.GetJulian();


    unsigned vectSize = 0;
    double*  vectElements = NULL;
 

    // Curves with basis adjustment
    
    ARM_Vector vectCopy = *(zeroCurve->GetYearTerms());
    ConvertYearTermsToJulianDates(dObservationDate, &vectCopy);
    vectSize     = vectCopy.GetSize();
    vectElements = vectCopy.GetElt();
    DKMaille<double> dDates(vectSize, vectElements);

    vectSize     = zeroCurve->GetDiscountFactors()->GetSize();
    vectElements = zeroCurve->GetDiscountFactors()->GetElt();
    ARM_Vector* baseYearTerms = zeroCurve->GetYearTerms();
    DKMaille<double> dRates(vectSize, vectElements);


    FromARMVectorToDKMaille(NoticeDates, dNoticeDates);
    FromARMVectorToDKMaille(Sigma, dSigma);
    FromARMVectorToDKMaille(HW3FParams, dHW3FParams);

    double Price;
        
    try 
    {
       Price = ImpliedFwdCorrelation_VFDK_HW1To3F(dSwaptionExpiryInYears,
                                                  dSwaptionTenorInYears,
                                                  dSwaptionTenor2InYears,
                                                  dNoticePeriodInDays,
                                                  dDates,
                                                  dRates,
                                                  dRates,
                                                  dNoticeDates,
                                                  dSigma,
                                                  dHW3FParams,
                                                  dObservationDate);
    }

    catch(char* a3FException)
    {
        char buf[200];

        strcpy(buf, "3F Model Implied Fwd Correlation: ");

        strcat(buf, a3FException);

        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL, buf);
    }

    catch(Exception& anExpt)
    {
        throw(anExpt);
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                        "Implied Fwd Correlation  failed!");
    }

    return(Price);    
}



ARM_Vector* HW3F_CALIBRATION(ARM_Date& AsOfDate,
                             ARM_ZeroCurve* zeroCurve,
                             ARM_Vector* HW3FParamsIn,
                             ARM_VolCurve* volCurve,
                             ARM_VolCurve* correlationCurve,
                             ARM_VolCurve* volWeightCurve,
                             ARM_VolCurve* correlationWeightCurve,
                             ARM_Vector* NoticeDates,
                             ARM_Vector* SwapStartDates,
                             ARM_Vector* SwapEndDates)
{
    DKMaille<double>   dStdDevX;
    DKMaille<double>   dStdDevY;
    DKMaille2D<double> dStdDevZ;
    DKMaille2D<double> dStdDevWeight;

    DKMaille<double>   dCorrelationX;
    DKMaille<double>   dCorrelationZ;
    DKMaille<double>   dCorrelationWeight;

    // Dates Vectors
    DKMaille<double>   dNoticeDates;
    DKMaille<double>   dSwapStartDates;
    DKMaille<double>   dSwapEndDates;
     

    DKMaille<double>   dHW3FParamsIn;

    double dSpotDate = AsOfDate.GetJulian();


    FromARMVolToDKMaille2D(dSpotDate, (ARM_VolLInterpol *) volCurve, 
                           dStdDevZ, dStdDevX, dStdDevY);

    FromARMVolToDKMaille2D(dSpotDate, (ARM_VolLInterpol *) volWeightCurve, 
                            dStdDevWeight, dStdDevX, dStdDevY);

    if (( dStdDevZ.rows() != dStdDevWeight.rows() ) 
        || 
        ( dStdDevZ.columns() != dStdDevWeight.columns() )
       )
    {
       throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
        "Vol Swopt Curve and Vol Swopt Weight Curve must have the same size");
    }

    ARM_Vector vectDates = *correlationCurve->GetExpiryTerms();
    // ConvertYearTermsToJulianDates(dSpotDate, &vectDates);
    FromARMVectorToDKMaille(&vectDates, dCorrelationX);

    ARM_Vector* vectValues = correlationCurve->GetVolatilities()->GetColumn(0);
    FromARMVectorToDKMaille(vectValues, dCorrelationZ);

    if (vectValues)
        delete vectValues;

    vectValues = correlationWeightCurve->GetVolatilities()->GetColumn(0);
    FromARMVectorToDKMaille(vectValues, dCorrelationWeight);

    if (vectValues)
        delete vectValues;

    if ( dCorrelationZ.entries() != dCorrelationWeight.entries() )
    {
       throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
      "Fwd Correlation Curve and Fwd Correlation Weight Curve must have the same size");
    }

    unsigned vectSize = 0;
    double*  vectElements = NULL;
   
    ARM_Vector vectCopy = *(zeroCurve->GetYearTerms());
    ConvertYearTermsToJulianDates(dSpotDate, &vectCopy);
    vectSize     = vectCopy.GetSize();
    vectElements = vectCopy.GetElt();
    DKMaille<double> dDates(vectSize, vectElements);

    vectSize     = zeroCurve->GetDiscountFactors()->GetSize();
    vectElements = zeroCurve->GetDiscountFactors()->GetElt();
    ARM_Vector* baseYearTerms = zeroCurve->GetYearTerms();
    DKMaille<double> dRates(vectSize, vectElements);


    FromARMVectorToDKMaille(NoticeDates, dNoticeDates);
    FromARMVectorToDKMaille(SwapStartDates, dSwapStartDates);
    FromARMVectorToDKMaille(SwapEndDates, dSwapEndDates);
    FromARMVectorToDKMaille(HW3FParamsIn, dHW3FParamsIn);

    DKMaille<double> dHW3FParamsOut;

    try 
    {
        dHW3FParamsOut = HW3F_CALIBRATION(dSpotDate,
                                          dDates,
                                          dRates,
                                          dHW3FParamsIn,
                                          dStdDevX,
                                          dStdDevY,
                                          dStdDevZ,
                                          dCorrelationX,
                                          dCorrelationZ,
                                          dStdDevWeight,
                                          dCorrelationWeight,
                                          dNoticeDates,
                                          dSwapStartDates,
                                          dSwapEndDates);
    }

    catch(char* a3FException)
    {
        char buf[200];

        strcpy(buf, "HW3F Calibration: ");

        strcat(buf, a3FException);

        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL, buf);
    }

    catch(Exception& anExpt)
    {
        throw(anExpt);
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                        "HW3F Calibration failed!");
    }

    // Convert fwdVol to ARM_Vector

    ARM_Vector* resVect = ConvertToARMVector(dHW3FParamsOut);

    return(resVect);    
}


ARM_VolLInterpol* PRCS3F_ConvertObjFwdVolToSpotVol(ARM_Date& AsOfDate,
                                             ARM_ZeroCurve* dBaseRatesNoBasisCurve,
                                             ARM_ZeroCurve* dBaseYieldCurve,
                                             ARM_ZeroCurve* dForeignRatesNoBasisCurve,
                                             ARM_ZeroCurve* dForeignYieldCurve,
                                             ARM_VolCurve*  volSwopBase,
                                             ARM_VolCurve*  volSwopForeign,
                           ARM_VolLInterpol*  volFx, // This is the inputed Fwd Vol!.
                                             double dMeanReversionBase, 
                                             double dMeanReversionForeign,
                                             ARM_VolCurve* dBaseSpotFXCorrelation,
                                             ARM_VolCurve* dForeignSpotFXCorrelation,
                                             ARM_VolCurve* dBaseForeignCorrelation,
                                             ARM_Vector* NoticeDates,
                                             ARM_Vector* FXCouponResetDates,
                                             ARM_Vector* FXCouponPaymentDate,
                                             double CutOff,
                                             double LongDatedSpotFxVol,
                                             int CalibSwoptBasis)
{
    ARM_VolLInterpol* SpotVol = NULL;

    int sz = volFx->GetExpiryTerms()->GetSize();


    double julAsOf = AsOfDate.GetJulian();
    int i;

    ARM_Vector* dSpotVolDates;
    ARM_Vector* dSpotVol;
    
    try
    {
         PRCS3F_FromFwdVolToSpotVol(AsOfDate,
                                    dBaseRatesNoBasisCurve,
                                    dBaseYieldCurve,
                                    dForeignRatesNoBasisCurve,
                                    dForeignYieldCurve,
                                    volSwopBase,
                                    volSwopForeign,
                                    volFx, // This is a Spot Vol!.
                                    dMeanReversionBase, 
                                    dMeanReversionForeign,
                                    dBaseSpotFXCorrelation,
                                    dForeignSpotFXCorrelation,
                                    dBaseForeignCorrelation,
                                    NoticeDates,
                                    FXCouponResetDates,
                                    FXCouponPaymentDate,
                                    CutOff,
                                    LongDatedSpotFxVol,
                                    CalibSwoptBasis,
                                    dSpotVolDates,
                                    dSpotVol);
    }
    
    catch(Exception& anExpt)
    {
        throw(anExpt);
    }
        
    ARM_Vector* tmpVolMatu = new ARM_Vector(dSpotVolDates->GetSize());

    for (i = 0; i < dSpotVolDates->GetSize(); i++)
    {
        tmpVolMatu->Elt(i) = (dSpotVolDates->Elt(i)-julAsOf)/365.0;
    }

    ARM_Vector* matuOp    = tmpVolMatu;

    ARM_Matrix* resFwdVol  = new ARM_Matrix(dSpotVolDates->GetSize(), 1);

    for (i = 0; i < dSpotVolDates->GetSize(); i++)
    {
        for (int j = 0; j < 1; j++)
        {
            resFwdVol->Elt(i, j) =  dSpotVol->Elt(i);               
        }
    }


    ARM_Vector* matuUnder = new ARM_Vector(1, 100.0);

    SpotVol = new ARM_VolLInterpol(AsOfDate, matuOp, matuUnder, resFwdVol);    
    
    if (dSpotVolDates)
       delete dSpotVolDates;

    if (dSpotVol)
       delete dSpotVol;

    return(SpotVol);
}



void PRCS3F_FromFwdVolToSpotVol(ARM_Date& AsOfDate,
                                ARM_ZeroCurve* dBaseRatesNoBasisCurve,
                                ARM_ZeroCurve* dBaseYieldCurve,
                                ARM_ZeroCurve* dForeignRatesNoBasisCurve,
                                ARM_ZeroCurve* dForeignYieldCurve,
                                ARM_VolCurve*  volSwopBase,
                                ARM_VolCurve*  volSwopForeign,
                                ARM_VolCurve*  volFx, // This is a Spot Vol!.
                                double dMeanReversionBase, 
                                double dMeanReversionForeign,
                                ARM_VolCurve* dBaseSpotFXCorrelation,
                                ARM_VolCurve* dForeignSpotFXCorrelation,
                                ARM_VolCurve* dBaseForeignCorrelation,
                                ARM_Vector* NoticeDates,
                                ARM_Vector* FXCouponResetDates,
                                ARM_Vector* FXCouponPaymentDates,
                                int    CutOff,
                                double LongDatedSpotFxVol,
                                int CalibSwoptBasis,
                                ARM_Vector*& dSpotVolDates,
                                ARM_Vector*& dSpotVol)
{
    // Dates Vectors
    DKMaille<double>   dStdDevBaseX;
    DKMaille<double>   dStdDevBaseY;
    DKMaille2D<double> dStdDevBaseZ;

    DKMaille<double>   dStdDevForeignX;
    DKMaille<double>   dStdDevForeignY;
    DKMaille2D<double> dStdDevForeignZ;


    // Interface DK Fwd Volatility function

    double dSpotDate = AsOfDate.GetJulian();


    unsigned vectSize = 0;
    double*  vectElements = NULL;

    // Curves with basis adjustment
    
    ARM_Vector vectCopy = *(dBaseYieldCurve->GetYearTerms());
    ConvertYearTermsToJulianDates(dSpotDate, &vectCopy);
    vectSize     = vectCopy.GetSize();
    vectElements = vectCopy.GetElt();
    DKMaille<double> dBaseDates(vectSize, vectElements);

    vectSize     = dBaseYieldCurve->GetDiscountFactors()->GetSize();
    vectElements = dBaseYieldCurve->GetDiscountFactors()->GetElt();
    ARM_Vector* baseYearTerms = dBaseYieldCurve->GetYearTerms();
    DKMaille<double> dBaseRates(vectSize, vectElements);


    vectCopy = *(dForeignYieldCurve->GetYearTerms());
    ARM_Vector* forYearTerms = dForeignYieldCurve->GetYearTerms();
    ConvertYearTermsToJulianDates(dSpotDate, &vectCopy);
    vectSize     = vectCopy.GetSize();
    vectElements = vectCopy.GetElt();
    DKMaille<double> dForeignDates(vectSize, vectElements);

    vectSize     = dForeignYieldCurve->GetDiscountFactors()->GetSize();
    vectElements = dForeignYieldCurve->GetDiscountFactors()->GetElt();
    DKMaille<double> dForeignRates(vectSize, vectElements);


    // Curves without Basis

    ARM_Vector newBaseDFNoBasis(baseYearTerms->GetSize());
    int k;
    for (k = 0; k < baseYearTerms->GetSize(); k++)
    {
        newBaseDFNoBasis.Elt(k) = 
              dBaseRatesNoBasisCurve->DiscountPrice(baseYearTerms->Elt(k));
    }
     
    vectSize     = baseYearTerms->GetSize();
    vectElements = newBaseDFNoBasis.GetElt();

    DKMaille<double> dBaseRatesNoBasis(vectSize, vectElements);


    ARM_Vector newForDFNoBasis(forYearTerms->GetSize());

    for (k = 0; k < forYearTerms->GetSize(); k++)
    {
        newForDFNoBasis.Elt(k) = 
               dForeignRatesNoBasisCurve->DiscountPrice(forYearTerms->Elt(k));
    }
     
    vectSize = forYearTerms->GetSize();
    vectElements = newForDFNoBasis.GetElt();
    DKMaille<double> dForeignRatesNoBasis(vectSize, vectElements);

    // End Curves without Basis

    
    // The correls are objects but similar to flat vol object

    double dBaseForeignCorrelationConv   = 0.0;
    if ( dBaseForeignCorrelation->GetName() != ARM_VOL_FLAT )
    {
       throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
        "Only Flat objects are handled so far for BaseForeignCorrelation");
    }
    else
    {
       dBaseForeignCorrelationConv = 
            ((ARM_VolFlat *) dBaseForeignCorrelation)->GetVolatility();
    }

    double dBaseSpotFXCorrelationConv    = 0.0;
    if ( dBaseSpotFXCorrelation->GetName() != ARM_VOL_FLAT )
    {
       throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
         "Only Flat objects are handled so far for BaseSpotFXCorrelation");
    }
    else
    {
       dBaseSpotFXCorrelationConv = 
              ((ARM_VolFlat *) dBaseSpotFXCorrelation)->GetVolatility();
    }

    double dForeignSpotFXCorrelationConv = 0.0;
    if ( dForeignSpotFXCorrelation->GetName() != ARM_VOL_FLAT )
    {
       throw Exception(__LINE__, __FILE__, ERR_UNIMP_METHOD_CALL,
        "Only Flat objects are handled so far for ForeignSpotFXCorrelation");
    }
    else
    {
       dForeignSpotFXCorrelationConv = 
                ((ARM_VolFlat *) dForeignSpotFXCorrelation)->GetVolatility();
    }

    FromARMVolToDKMaille2D(dSpotDate, (ARM_VolLInterpol *) volSwopBase, 
                           dStdDevBaseZ, dStdDevBaseX, dStdDevBaseY);
                            

    FromARMVolToDKMaille2D(dSpotDate, (ARM_VolLInterpol *) volSwopForeign, 
                           dStdDevForeignZ, dStdDevForeignX, dStdDevForeignY);

    DKMaille<double> dFwdFXVolDatesTD;
    DKMaille<double> dFwdFXVolTD;

    ARM_Vector vectDates = *volFx->GetExpiryTerms();

    ARM_Vector zeroVect(FXCouponResetDates->GetSize());

    FromARMVectorToDKMaille(&vectDates, dFwdFXVolDatesTD);

    // Obtain FX year terms vectors
     
    ARM_Vector* vectValues = volFx->GetVolatilities()->GetColumn(0);

    ARM_Vector vectValFlagAndVol(vectValues->GetSize()+2);

    for (int h = 0; h < vectValFlagAndVol.GetSize()-2; h++)
    {
        vectValFlagAndVol.Elt(h) = vectValues->Elt(h)/100.0;
    }

    vectValFlagAndVol.Elt(vectValFlagAndVol.GetSize()-2) = (double) CutOff;
    vectValFlagAndVol.Elt(vectValFlagAndVol.GetSize()-1) = LongDatedSpotFxVol;

    FromARMVectorToDKMaille(&vectValFlagAndVol, dFwdFXVolTD);

    if (vectValues)
       delete vectValues;

    DKMaille<double> dNoticeDates;
    FromARMVectorToDKMaille(NoticeDates, dNoticeDates);
      
    DKMaille<double> dFXCouponResetDates;
    FromARMVectorToDKMaille(FXCouponResetDates, dFXCouponResetDates);

    DKMaille<double> dFXCouponPaymentDates;
    FromARMVectorToDKMaille(FXCouponPaymentDates, dFXCouponPaymentDates);

    DKMaille<double> SpotVolDates;
    DKMaille<double> SpotVol;
      
    double dCurrencyPair1;
    double dCurrencyPair2;

    try 
    {
        dCurrencyPair1 = Get3FCurrencyId(dBaseRatesNoBasisCurve->GetCurrencyUnit());
        dCurrencyPair2 = Get3FCurrencyId(dForeignRatesNoBasisCurve->GetCurrencyUnit());

        FromForwardToSpotVol3F(dSpotDate, 
                               dBaseDates,
                               dBaseRatesNoBasis,
                               dBaseRates,
                               dForeignDates,
                               dForeignRatesNoBasis,
                               dForeignRates,
                               dStdDevBaseX, 
                               dStdDevBaseY,
                               dStdDevBaseZ,
                               dStdDevForeignX, 
                               dStdDevForeignY,
                               dStdDevForeignZ,
                               dFwdFXVolDatesTD,
                               dFwdFXVolTD,
                               dMeanReversionBase,
                               dMeanReversionForeign,
                               dBaseSpotFXCorrelationConv,
                               dForeignSpotFXCorrelationConv,
                               dBaseForeignCorrelationConv,
                               dCurrencyPair1,
                               dCurrencyPair2,
                               dNoticeDates,
                               dFXCouponResetDates,
                               dFXCouponPaymentDates,
                               CalibSwoptBasis,
                               SpotVolDates,
                               SpotVol);
    }

    catch(char* a3FException)
    {
        char buf[200];

        strcpy(buf, "3F Model Vol Conv: ");

        strcat(buf, a3FException);

        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL, buf);
    }

    catch(Exception& anExpt)
    {
        throw(anExpt);
    }

    catch(...)
    {
        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
                        "Converting Fwd Vol To Spot Vol Failed!");
    }

    // Convert fwdVol to ARM_Vector

    ARM_Vector* dSpotVolDatesLocal = ConvertToARMVector(SpotVolDates);
    ARM_Vector* dSpotVolLocal = ConvertToARMVector(SpotVol);

    int sz = dSpotVolLocal->GetSize();

    for (int i = 0; i < sz; i++)
    {
        dSpotVolLocal->Elt(i) *= 100.0;
    }

    dSpotVolDates = dSpotVolDatesLocal;
    dSpotVol      = dSpotVolLocal;
}



/*----------------------------------------------------------------------------------*/
/*---- End Of File ----*/
