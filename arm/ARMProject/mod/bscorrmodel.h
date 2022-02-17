/*
 * $Log: bscorrmodel.h,v $
 * Revision 1.4  2004/01/21 11:39:43  emezzine
 * added #include"volcurve.h"
 *
 * Revision 1.3  2004/01/13 19:43:32  ebenhamou
 * version to have a spread option working in ARM
 *
 * Revision 1.2  2002/11/13 10:35:33  mab
 * Formatting
 *
 * Revision 1.1  2002/05/23 15:52:59  sgasquet
 * Initial revision
 *
 */


/*----------------------------------------------------------------------------*

    bscorrmodel.h
    
    This is a Black model enabling to price spread options
    Have been added to standard Black model:    
    - a Cash Vol Matrix (added to an IRG Vol Matrix)
    - a Correlation Matrix 
    - a SpreadOption Vol Matrix
    - Flags for Pricing Type and for SpreadVol Type                                                   
*----------------------------------------------------------------------------*/

#ifndef _BSCORRMODEL_H
#define _BSCORRMODEL_H



#include "bsmodel.h"
#include "foncstd.h"
#include "volcurv.h"


class ARM_Date;
class ARM_Vector;



class ARM_BSCorrModel : public ARM_BSModel
{
    private:
/*
        ARM_VolCurve*  itsCorrelationMatrix; 
        ARM_VolCurve*  itsCashVolMatrix; 
        ARM_VolCurve*  itsSpreadVolCurve; 

        int itsModelType;
        int itsVolType;*/


        void Init(void)
        {
            SetBSCategory(1);

 /*           itsCorrelationMatrix  = NULL;  
            itsCashVolMatrix = NULL;
            itsSpreadVolCurve = NULL;*/
        }


    public:

        ARM_BSCorrModel(void)
        {
            Init();
            //SetName(ARM_BSCorrModel);
        }
        
        
        ARM_BSCorrModel(ARM_Date& startDate,                         
                        ARM_ZeroCurve* zeroCurve,        // ZC
                        ARM_VolCurve* spreadLock,        // SpreadLock Vol
                        ARM_VolCurve* capIRGVol,         // IRG Vol
                        ARM_VolCurve* capCashVol = NULL, // CASH Vol
                        ARM_VolCurve* indexVAdjol = NULL,// Adjustment Vol
                        ARM_VolCurve* SpreadVol = NULL, 
                        ARM_VolCurve* Correlations = NULL, // Correlation Structure
                        int ModelType = K_2LOG,
                        int VolType = K_INPUTED); 
                         // ModelType 2LOG: 2 Underlyings LOGNOR 
                         // NOR:    Spread gaussien
                         // LOG:    Spread lognormal             
                         // VolType K_INPUTED: we directly use 
                         // the inputed spread vol
                         // K_COMPUTED:we compute the spread 
                         // vol implied by the correl and the 2 underl vols                                    

        ARM_BSCorrModel& operator = (ARM_BSCorrModel& bs);

        
       ~ARM_BSCorrModel(void) 
        {
/*            if (itsCorrelationMatrix)
               delete itsCorrelationMatrix;

            itsCorrelationMatrix = NULL;

            if (itsCashVolMatrix)
               delete itsCashVolMatrix;

            itsCashVolMatrix = NULL;

            if (itsSpreadVolCurve)
               delete itsSpreadVolCurve;

            itsSpreadVolCurve = NULL;*/
        }


       
        void Set(ARM_VolCurve* capCashVol,           // CASH Vol 
                 ARM_VolCurve* SpreadVol = NULL,     // Spread Vol
                 ARM_VolCurve* Correlations = NULL,  // Correlation Structure
                 int ModelType = 0 /* K_2LOG    */,
                 int VolType = 0   /* K_INPUTED */);


        
        void BitwiseCopy(const ARM_Object* srcBSModel)
        {
            ARM_BSCorrModel* bsmodel = (ARM_BSCorrModel *) srcBSModel;
            ARM_BSModel::BitwiseCopy(srcBSModel);

//            itsModelType = bsmodel->GetModelType();
//            itsVolType = bsmodel->GetVolType();

/*            if (itsCashVolMatrix)
               delete itsCashVolMatrix;*/
            if (bsmodel->GetCashVolMatrix())
            {
               /*itsCashVolMatrix = (ARM_VolCurve *)
                           bsmodel->GetCashVolMatrix()->Clone();*/
				SetCashVolMatrix(bsmodel->GetCashVolMatrix());
            }
            else
            {
                SetCashVolMatrix(NULL);   
            }


/*            if (itsCorrelationMatrix)
               delete itsCorrelationMatrix;*/

            if (bsmodel->GetCorrelationMatrix())
            {
/*               itsCorrelationMatrix  = (ARM_VolCurve *)
                     bsmodel->GetCorrelationMatrix()->Clone();*/
				SetCorrelationMatrix(bsmodel->GetCorrelationMatrix());
            }
            else
            {
                SetCorrelationMatrix(NULL);
            }

/*            if (itsSpreadVolCurve)
                delete itsSpreadVolCurve;*/


            if (bsmodel->GetSpreadVolCurve())
            {
/*               itsSpreadVolCurve = (ARM_VolCurve *)
                           bsmodel->GetSpreadVolCurve()->Clone();*/
				SetSpreadVolCurve(bsmodel->GetSpreadVolCurve());
            }
            else
            {
                SetSpreadVolCurve(NULL);   
            }
        }
 

        void Copy(const ARM_Object* srcBSModel)
        {
            ARM_Model::Copy(srcBSModel);
 
            BitwiseCopy(srcBSModel);
        }
 


        ARM_Object* Clone(void)
        {
            ARM_BSCorrModel* theClone = new ARM_BSCorrModel();
  
            theClone->Copy(this);
 
            return(theClone);
        }



        double SpreadOptionPrice(double fwd1, double fwd2, 
             double vol1, 
             double vol2, double Correl, 
             double strike, 
             double optMat, int optType, int ComputedFormula = 1/* 0 formule d'Olivier, 1 Formule Classique*/); 

        /// static method to allow call directly from outside
        static double SpreadOptionFormula(double fwd1, double fwd2, 
             double vol1, 
             double vol2, double Correl, 
             double strike, 
             double optMat, int optType, int modelType, double SpreadVol=-1, int ComputedFormula = 1); 
};


#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/
