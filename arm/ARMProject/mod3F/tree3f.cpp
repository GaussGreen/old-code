/*
 * $Log: tree3f.cpp,v $
 * Revision 1.11  2004/06/14 15:15:57  mab
 * Added: ARM_VolLInterpol* ARM_Tree3F::ConvertSpotVolTo3FFwdVol
 *
 * Revision 1.10  2003/10/21 11:03:24  mab
 * last release from Dimitri
 *
 * Revision 1.9  2003/10/08 13:32:22  mab
 *  Systematic conversion of Vol Fx From spot To Fwd
 *
 * Revision 1.8  2003/10/06 09:35:49  mab
 * Added: void ARM_Tree3F::View(char* id, FILE* ficOut)
 *
 * Revision 1.7  2003/09/12 16:55:18  mab
 * Supression of itsNoticeDates
 *
 * Revision 1.6  2003/09/02 16:00:49  mab
 * Added : SetStartDate(itsAsOfDate) in constructor
 *
 * Revision 1.5  2003/08/14 09:53:37  mab
 * Improvements
 *
 * Revision 1.4  2003/08/01 13:05:08  mab
 * Improvements
 *
 * Revision 1.2  2003/06/30 16:48:26  mab
 * modif for RCS
 *
 * Revision 1.1  2003/06/30 16:33:33  mab
 * Initial revision
 *
*/


/*-----------------------------------------------------------------------------*
   tree3f.cpp

   Header of the 3 factors tree model class. this model
   allows (for now) the pricing of products such as PRCS

*-----------------------------------------------------------------------------*/

#include "firsttoinc.h"

#include "DK_prcs.h"
#include "tree3f.h"



ARM_Tree3F::ARM_Tree3F(const ARM_Tree3F& tree3F)
{
    Init();
}



ARM_Tree3F& ARM_Tree3F::operator = (const ARM_Tree3F& tree3F)
{
   (*this).ARM_Model::operator = (tree3F);

    BitwiseCopy(&tree3F);

    return(*this);
}



ARM_Tree3F::ARM_Tree3F(ARM_Date& asOfDate,
                       ARM_DFBSModel* dFBSModel,
                       ARM_VolCurve*  volSwopBase,
                       ARM_VolCurve*  volSwopForeign,
                       ARM_VolCurve*  baseForeignCorrelation,
                       ARM_Vector*    latticeGeometryData,
                       double         numTimeLinesBeforeFirstNotice,
                       double numTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal,
                       double         numTimeLines,
                       double         meanReversionBase,
                       double         meanReversionForeign,
                       double         x1Limit,
                       double         x2Limit,
                       double         x3Limit,
                       double         i1Limit,
                       double         i2Limit,
                       double         i3Limit,
                       double         aDPLimit,
                       double         optimal,
                       double         timeBoost,
                       double         deltaFlag,
                       double         smoothingValue,
                       double         calcProbaSurvOrNot,
					   double         QBaseSmile,
					   double         QForeignSmile,
                       double         CutOff,
                       double         LongDatedSpotFxVol,
                       int            FwdVolConvFlag,
                       int            SwopConvFlag,
                       int            CalibSwoptBasis)
{
    Init();

    itsAsOfDate = asOfDate;

    SetStartDate(itsAsOfDate);

    itsDFBSModel = (ARM_DFBSModel *) dFBSModel->Clone();

    if (SwopConvFlag)
    {
       // Conversion of the lognormal SUMMIT vol, as found in the BS model,
       // into absolute unitary volatility surfaces,
       // as awaited by Dimitri's PRCS3F_Lattice_HWVFDK_Pricing function.

       ARM_VolLInterpol* DomAbsVol = volSwopBase->ConvertToNormalVol(
                                     dFBSModel->GetDBSModel()->GetDiscountCurve());

       itsVolSwopBase = (ARM_VolCurve *) DomAbsVol;

       ARM_VolLInterpol* ForAbsVol = volSwopForeign->ConvertToNormalVol(
                                     dFBSModel->GetFBSModel()->GetDiscountCurve());

       itsVolSwopForeign = (ARM_VolCurve *) ForAbsVol;
    }
    else
    {
       itsVolSwopBase    = (ARM_VolCurve *) volSwopBase->Clone();
       itsVolSwopForeign = (ARM_VolCurve *) volSwopForeign->Clone(); 
    } 

    itsBaseForeignCorrelation = (ARM_VolCurve *) baseForeignCorrelation->Clone();

    itsLatticeGeometryData = (ARM_Vector *) latticeGeometryData->Clone();

    itsNumTimeLinesBeforeFirstNotice = numTimeLinesBeforeFirstNotice;
    itsNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal =
           numTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal;
    itsNumTimeLines = numTimeLines;
    itsMeanReversionBase = meanReversionBase;
    itsMeanReversionForeign = meanReversionForeign;
    itsX1Limit = x1Limit;
    itsX2Limit = x2Limit;
    itsX3Limit = x3Limit;
    itsI1Limit = i1Limit;
    itsI2Limit = i2Limit;
    itsI3Limit = i3Limit;
    itsADPLimit = aDPLimit;
    itsOptimal  = optimal;
    itsTimeBoost = timeBoost;
    itsDeltaFlag  = deltaFlag;
    itsSmoothingValue = smoothingValue;

    itsCalcProbaSurvOrNot = calcProbaSurvOrNot;

    itsCutOff = CutOff;

    itsFwdVolConvFlag = FwdVolConvFlag;

    itsLongDatedSpotFxVol = LongDatedSpotFxVol;

    itsCalibSwoptBasis = CalibSwoptBasis;

    itsQBaseSmile = QBaseSmile;
				   
    itsQForeignSmile = QForeignSmile;
}



ARM_VolLInterpol* ARM_Tree3F::ConvertSpotVolTo3FFwdVol(ARM_Vector* noticeDates,
                                                       ARM_Vector* FXCouponResetDates,
                                                       ARM_Vector* FXCouponPaymentDates,
                                                       ARM_Vector* FwdVolDates)
{
    ARM_VolLInterpol* fwdVol = NULL;


    try
    {
        fwdVol = PRCS3F_ConvertObjSpotVolToFwdVol(this->itsAsOfDate,
                         this->itsDFBSModel->GetDomYieldCurve(),
                         this->itsDFBSModel->GetDBsCrv(),
                         this->itsDFBSModel->GetForeignYieldCurve(),
                         this->itsDFBSModel->GetFBsCrv(),
                         this->itsVolSwopBase,
                         this->itsVolSwopForeign,
                 (ARM_VolLInterpol *) this->itsDFBSModel->GetFxVol(),
                         this->itsMeanReversionBase,
                         this->itsMeanReversionForeign,
                         this->itsDFBSModel->GetdFxCorr(),// BaseSpotFXCorrel
                         this->itsDFBSModel->GetfFxCorr(),// ForgnSpotFXCorrel
                         this->itsBaseForeignCorrelation,
                         noticeDates,
                         FXCouponResetDates,
                         FXCouponPaymentDates,
                         this->itsCutOff,
                         this->itsLongDatedSpotFxVol,
                         this->itsCalibSwoptBasis,
                         FwdVolDates);
    }

    catch(char* a3FException)
    {
        if (fwdVol)
           delete fwdVol;

        char buf[200];

        strcpy(buf, "ARM_Tree3F::ConvertSpotVolTo3FFwdVol: 3F Model: ");

        strcat(buf, a3FException);

        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL, buf);

        return(NULL);
    }

    catch(Exception& theExpt)
    {
        if (fwdVol)
           delete fwdVol;

        throw theExpt;

        return(NULL);
    }

    catch(...)
    {
        if (fwdVol)
           delete fwdVol;

        throw Exception(__LINE__, __FILE__, ERR_INVALID_MODEL,
           "ARM_Tree3F::ConvertSpotVolTo3FFwdVol: 3F Model: Failed!");

        return(NULL);
    }

    return(fwdVol);
}



void ARM_Tree3F::View(char* id, FILE* ficOut)
{
    FILE* fOut;

    char fOutName[200];


    if ( ficOut == NULL )
    {
       ARM_GetViewFile(ARM_VIEW_FILE, id, fOutName);

       (void) unlink(fOutName);

       fOut = fopen(fOutName, "w");
    }
    else
    {
       fOut = ficOut;
    }

    fprintf(fOut, "\n ====> Converted FX Vol: \n");
    if (itsCalcFxFwdVol)
       itsCalcFxFwdVol->View(id, fOut); 

    fprintf(fOut, "\n\n ====> VolSwopBase: \n");

    itsVolSwopBase->View(id, fOut);

    fprintf(fOut, "\n\n ====> VolSwopForeign: \n");

    itsVolSwopForeign->View(id, fOut);


    if ( ficOut == NULL )
    {
       fclose(fOut);
    }
}













/*-----------------------------------------------------------------------------*/
/*---- End Of File ----*/
