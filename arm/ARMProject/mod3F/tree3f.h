/*
 * $Log: tree3f.h,v $
 * Revision 1.10  2004/06/14 15:14:43  mab
 * Added :  ARM_VolLInterpol* ConvertSpotVolTo3FFwdVol(..)
 *
 * Revision 1.9  2004/03/05 17:29:03  mab
 * Added : void SetDiscPricingMode(int discPricingMode)
 *
 * Revision 1.8  2003/10/21 11:03:36  mab
 * last release from Dimitri
 *
 * Revision 1.7  2003/10/08 13:31:41  mab
 * Systematic conversion of Vol Fx From spot To Fwd
 *
 * Revision 1.6  2003/10/06 09:35:06  mab
 * Added : View(..)
 *
 * Revision 1.5  2003/09/12 16:54:56  mab
 * Supression of itsNoticeDates
 *
 * Revision 1.4  2003/08/14 09:51:48  mab
 * Added: ARM_VolCurve*  itsCalcFxFwdVol;
 *
 * Revision 1.3  2003/08/01 13:05:33  mab
 * Improvements
 *
 * Revision 1.1  2003/06/30 16:49:16  mab
 * Initial revision
 *
 */


/*-----------------------------------------------------------------------------*
   tree3f.h

   Header of the 3 factors tree model class. this model
   allow (for now) the pricing of products such as PRCS 
   
*-----------------------------------------------------------------------------*/
#ifndef _TREE3F_H
#define _TREE3F_H


#include "xbsfx.h"



class ARM_Tree3F : public ARM_Model
{
//    private:

      public : // For now!

        ARM_Date itsAsOfDate;

        ARM_DFBSModel* itsDFBSModel;

        ARM_VolCurve*  itsVolSwopBase;
        ARM_VolCurve*  itsVolSwopForeign;
        ARM_VolCurve*  itsBaseForeignCorrelation;

        ARM_Vector*    itsLatticeGeometryData;

        double         itsNumTimeLinesBeforeFirstNotice;
        double         itsNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal;
        double         itsNumTimeLines;
        double         itsMeanReversionBase;
        double         itsMeanReversionForeign;
        double         itsX1Limit;
        double         itsX2Limit;
        double         itsX3Limit;
        double         itsI1Limit;
        double         itsI2Limit;
        double         itsI3Limit;
        double         itsADPLimit;
        double         itsOptimal;
        double         itsTimeBoost;
        double         itsDeltaFlag;
        double         itsSmoothingValue;

        int            itsCalcProbaSurvOrNot;

		double         itsQBaseSmile;
				   
		double         itsQForeignSmile;

        double         itsCutOff;

        int            itsFwdVolConvFlag;

        int            itsSwopConvFlag;

        double         itsLongDatedSpotFxVol;

        ARM_VolCurve*  itsCalcFxFwdVol;

        int            itsCalibSwoptBasis;


        void Init(void)
        {
            SetName(ARM_TREE3F);

            itsDFBSModel              = NULL;

            itsVolSwopBase            = NULL;
            itsVolSwopForeign         = NULL;
            itsBaseForeignCorrelation = NULL;

            itsLatticeGeometryData    = NULL;

            itsNumTimeLinesBeforeFirstNotice = 0;
            itsNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal = 0;
            itsNumTimeLines = 0;
            itsMeanReversionBase    = 0.0;
            itsMeanReversionForeign = 0.0;

            itsX1Limit  = 0.0;
            itsX2Limit  = 0.0;
            itsX3Limit  = 0.0;
            itsI1Limit  = 0.0;
            itsI2Limit  = 0.0;
            itsI3Limit  = 0.0;
            itsADPLimit = 0.0;
            itsOptimal  = 0.0;
            itsTimeBoost = 0.0;

            itsDeltaFlag          = 0;

            itsSmoothingValue     = 0.0;

            itsCalcProbaSurvOrNot = 0;

			itsQBaseSmile         = 0.0;
				   
		    itsQForeignSmile      = 0.0;
       
            itsCutOff             = 0.0;
 
            itsFwdVolConvFlag     = 0; // The default is
                                       // no conversion to fwd Vol 

            itsSwopConvFlag       = 0; // The default is
                                       // no conversion to fwd Vol

            itsLongDatedSpotFxVol = 0.0;

            itsCalcFxFwdVol       = NULL;

            itsCalibSwoptBasis    = 1;
        }

    public :

        ARM_Tree3F(void)
        {
            Init();
        }

        ARM_Tree3F(const ARM_Tree3F& tree3F);

        ARM_Tree3F& operator = (const ARM_Tree3F& tree3F);

        ARM_Tree3F(ARM_DFBSModel* bsFxMod);

       ~ARM_Tree3F(void)
        {
            if (itsDFBSModel)
               delete itsDFBSModel;

            if (itsVolSwopBase)
               delete itsVolSwopBase;

            if (itsVolSwopForeign)
               delete itsVolSwopForeign;

            if (itsBaseForeignCorrelation)
               delete itsBaseForeignCorrelation;

            if (itsLatticeGeometryData)
               delete itsLatticeGeometryData;

            if (itsCalcFxFwdVol)
               delete itsCalcFxFwdVol;
        }
   
        ARM_Tree3F(ARM_Date& asOfDate,
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
				   double         QBaseSmile       = 0.0,
				   double         QForeignSmile    = 0.0,
                   double         CutOff           = 0.0,
                   double         LongDatedSpotFxVol = 0.0,
                   int            FwdVolConvFlag   = 0,
                   int            SwopConvFlag  = 0,
                   int            CalibSwoptBasis = 1);

        void BitwiseCopy(const ARM_Object* src3FModel)
        {
            ARM_Tree3F* tree3F = (ARM_Tree3F *) src3FModel;


            if (itsDFBSModel)
               delete itsDFBSModel;
            itsDFBSModel              = NULL;

            if (tree3F->itsDFBSModel)
               itsDFBSModel = (ARM_DFBSModel *) tree3F->itsDFBSModel->Clone();

            if (itsVolSwopBase)
               delete itsVolSwopBase;
            itsVolSwopBase            = NULL;

            if (tree3F->itsVolSwopBase)
               itsVolSwopBase = (ARM_VolCurve *) tree3F->itsVolSwopBase->Clone();

            if (itsVolSwopForeign)
               delete itsVolSwopForeign;
            itsVolSwopForeign         = NULL;

            if (tree3F->itsVolSwopForeign)
               itsVolSwopForeign = 
                    (ARM_VolCurve *) tree3F->itsVolSwopForeign->Clone();

            if (itsBaseForeignCorrelation)
               delete itsBaseForeignCorrelation;
            itsBaseForeignCorrelation = NULL;

            if (tree3F->itsBaseForeignCorrelation)
               itsBaseForeignCorrelation = 
                      (ARM_VolCurve *) tree3F->itsBaseForeignCorrelation->Clone();

            if (itsLatticeGeometryData)
               delete itsLatticeGeometryData;
            itsLatticeGeometryData = NULL;

            if (tree3F->itsLatticeGeometryData)
               itsLatticeGeometryData = 
                          (ARM_Vector *) tree3F->itsLatticeGeometryData->Clone();
   
            itsNumTimeLinesBeforeFirstNotice = 
                  tree3F->itsNumTimeLinesBeforeFirstNotice;
            itsNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal = 
                  tree3F->itsNumTimeLinesPerYearAfterFirstNoticeAndBeforeOptimal;

            itsNumTimeLines = tree3F->itsNumTimeLines;
            itsMeanReversionBase    = tree3F->itsMeanReversionBase;
            itsMeanReversionForeign = tree3F->itsMeanReversionForeign;

            itsX1Limit  = tree3F->itsX1Limit;
            itsX2Limit  = tree3F->itsX2Limit;
            itsX3Limit  = tree3F->itsX3Limit;
            itsI1Limit  = tree3F->itsI1Limit;
            itsI2Limit  = tree3F->itsI2Limit;
            itsI3Limit  = tree3F->itsI3Limit;
            itsADPLimit = tree3F->itsADPLimit;
            itsOptimal  = tree3F->itsOptimal;
            itsTimeBoost = tree3F->itsTimeBoost;

            itsDeltaFlag = tree3F->itsDeltaFlag;

            itsSmoothingValue = tree3F->itsSmoothingValue;

            itsCalcProbaSurvOrNot = tree3F->itsCalcProbaSurvOrNot;

            itsCutOff    = tree3F->itsCutOff;
 
            itsFwdVolConvFlag     = tree3F->itsFwdVolConvFlag;

            itsSwopConvFlag       = tree3F->itsSwopConvFlag;

            itsLongDatedSpotFxVol   = tree3F->itsLongDatedSpotFxVol;
           
            if (itsCalcFxFwdVol)
            {
               delete itsCalcFxFwdVol;

               itsCalcFxFwdVol = NULL;
            }
            
            if (tree3F->itsCalcFxFwdVol)
               itsCalcFxFwdVol = 
                    (ARM_VolCurve *) tree3F->itsCalcFxFwdVol->Clone();

            itsCalibSwoptBasis = tree3F->itsCalibSwoptBasis;
        }

        void Copy(const ARM_Object* src3FModel)
        {
            ARM_Model::Copy(src3FModel);

            BitwiseCopy(src3FModel);
        }

        ARM_Object* Clone(void)
        {
            ARM_Tree3F* theClone = new ARM_Tree3F();


            theClone->Copy(this);

            return(theClone);
        }


        void SetDiscPricingMode(int discPricingMode)
        {
            ARM_Model::SetDiscPricingMode(discPricingMode);

            if (itsDFBSModel)
               itsDFBSModel->SetDiscPricingMode(discPricingMode);
        }

        ARM_VolCurve* GetCalcFxFwdVol(void)
        {
            return(itsCalcFxFwdVol);
        }

        void SetCalcFxFwdVol(ARM_VolCurve* fwdVol)
        {
            if (itsCalcFxFwdVol)
            {
               delete itsCalcFxFwdVol;

               itsCalcFxFwdVol = NULL;
            }
            
            if (fwdVol)
               itsCalcFxFwdVol = (ARM_VolCurve *) fwdVol->Clone();
        }

        ARM_VolLInterpol* ConvertSpotVolTo3FFwdVol(ARM_Vector* noticeDates,
                                                   ARM_Vector* FXCouponResetDates,
                                                   ARM_Vector* FXCouponPaymentDates,
                                                   ARM_Vector* FwdVolDates);


        void View(char* id = NULL, FILE* fOut = NULL);
};





























#endif
/*-----------------------------------------------------------------------------*/
/*---- End Of File ----*/
