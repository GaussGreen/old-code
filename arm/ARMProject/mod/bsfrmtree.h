/*
 * $Log: bsfrmtree.h,v $
 * Revision 1.3  2003/07/10 06:16:36  ebenhamou
 * pragma for bool conversion
 *
 * Revision 1.2  2002/07/02 12:56:44  mab
 * Yan LI add-on
 *
 */


#ifndef _BSFRMTREE_H
#define _BSFRMTREE_H


/*----------------------------------------------------------*/



#include "frmtree.h"






class ARM_BSFrmTree : public virtual ARM_FrmTree1
{
    private :

       int                itsBasisDiscountFlag;
       ARM_ZeroCurve*     itsUnAdjustedCurve;
       ARM_ZeroCurve*     itsAdjustedCurve;
       ARM_ZeroLInterpol* itsAdjWorkingCurve;

       ARM_Forwards*      itsAdjustedFwds;

       ARM_Vector*        itsBasisSpreads; // spread = basis 
                                           // adjusted - plain libor


    public : 

        ARM_BSFrmTree(void)
        {
            Init();
        }

        double GetFirstPeriodBSAdj(void)
        {
            return exp(-itsBasisSpreads->Elt(0));
        }

        // Construction par calibrateur    

        ARM_BSFrmTree(ARM_ZeroCurve* zc, ARM_ZeroCurve* zd,
                      ARM_VolCurve* vol, ARM_VolCurve* smile,
                      int mfine,
                      int autoMode, ARM_Date& endDate,  
                      double decay = 0, double slope = 0, 
                      double asymptote = 0, int NbFactor = 1, 
                      ARM_Vector* CorrelatedIndexes = NULL,
                      ARM_Vector* indexes = NULL,
                      ARM_Matrix* correlations = NULL);

        ARM_BSFrmTree(ARM_ZeroCurve* zc, ARM_ZeroCurve* zd, 
                      ARM_VolCurve* vol, ARM_VolCurve* smile,
                      ARM_VolCurve* irgvol, 
                      ARM_VolCurve* irgsmile, int mfine,
                      int autoMode, ARM_Date& endDate,  
                      double decay = 0, double slope = 0,
                      double asymptote = 0, int NbFactor = 1, 
                      ARM_Vector* CorrelatedIndexes = NULL,
                      ARM_Vector* indexes = NULL,
                      ARM_Matrix* correlations = NULL);

        ARM_BSFrmTree(ARM_FrmAna* AnaModel, ARM_ZeroCurve* baZc,
                      ARM_Date& endDate, int fineSteps  = K_L_DEFAULT,
                      ARM_Vector* CorrelatedIndexes = NULL, 
                      ARM_Vector* indexes = NULL,
                      ARM_Matrix* correlations = NULL);

/*
    // Construction pour autocalibration (1 classe de risque en vol)
    ARM_BSFrmTree(ARM_ZeroCurve* zc, ARM_ZeroCurve* baZc,
                  ARM_VolCurve* vol, ARM_VolCurve* smile, 
                  int productType, ARM_Date &endDate,  
                  int Ntraj = 1000, int  MCGeneratorType = K_MC_FAURE,
                  double decay = 0, double slope = 0,
                  double asymptote = 0, int NbFactor = 1, 
                  ARM_Vector* CorrelatedIndexes = NULL,
                  ARM_Vector* indexes  = NULL,
                  ARM_Matrix* correlations = NULL,
                  int NoControl = 0, long seed = 10000);
     
    // Construction pour autocalibration (2 classes de risque en vol)    
    ARM_BSFrmTree(ARM_ZeroCurve* zc, ARM_ZeroCurve* baZc,
                  ARM_VolCurve* swoptVol, ARM_VolCurve *swoptSmile, 
                  ARM_VolCurve* irgVol, ARM_VolCurve *irgSmile, 
                  int productType, ARM_Date &endDate,  
                  int Ntraj = 1000, int  MCGeneratorType = K_MC_FAURE,
                  double decay = 0, double slope = 0,
                  double asymptote = 0, int NbFactor = 1, 
                  ARM_Vector* CorrelatedIndexes = NULL,
                  ARM_Vector* indexes  = NULL,
                  ARM_Matrix* correlations = NULL, 
                  int NoControl = 0, long seed = 10000);
*/

     ARM_BSFrmTree(ARM_BSFrmTree &inModel) : ARM_FrmTree1(inModel)
     {
         Init();

         BitwiseCopy((ARM_BSFrmTree *) &inModel);
     };


    ~ARM_BSFrmTree(void)
     {
          cleanit(itsBasisSpreads);

          cleanit(itsAdjustedFwds);

          cleanit(itsAdjWorkingCurve);
     };

/*
    ARM_BSFrmTree operator = (ARM_BSFrmTree &inModel)
    {
        (*this).ARM_FrmTree1::operator =(inModel);

        BitwiseCopy((ARM_BSFrmTree *) &inModel);
        
        return (*this);
    }
*/
    void Init(void)
    {
        SetName(ARM_FRM_TREE);

        itsBasisDiscountFlag = 1;
        itsUnAdjustedCurve = NULL;
        itsBasisSpreads    = NULL;
        itsAdjustedCurve   = NULL;
        itsAdjWorkingCurve = NULL;
        itsAdjustedFwds    = NULL;
    };

    ARM_Object *Clone(void)
    {
        ARM_BSFrmTree *theClone = new ARM_BSFrmTree(*this);
        return theClone;
    }


    ARM_ZeroCurve* GetDiscountCurve(void)
    {
        if (itsBasisDiscountFlag)
           return(itsAdjWorkingCurve);
        else
           return(GetWorkingCurve());
    }

    double FirstPeriodDiscountFactor(int nPath);

    double ZeroPrice(double calcDate, double zMaturity, int DomOrFrg = 1);

    double ZeroPriceOnDiscountCurve(double calcDate,
                                    double zMaturity, int DomOrFrg = 1);

    void BeFittedTo(ARM_Security* sec, char* ccyName, int mfine);

    void SetWorkingFor_n_Curve(ARM_Forwards *inFwds);

    bool IsBasisDiscountModel(void)
    {
		/// to remove the performance warning
		#pragma warning(disable : 4800)
        return (bool)itsBasisDiscountFlag;
    }

    void SetBasisDiscountFlag(int flag=1)
    {
        itsBasisDiscountFlag = flag;
    }

    void Propagate(ARM_Object* OutCurve,
                   ARM_Vector* dr1,
                   int time1,
                   int time2,
                   ARM_Object* InpCurve, 
                   bool propFullCurve,
                   ARM_Vector* driversTarget,
                   ARM_Vector* driversCurrent);

};



#endif

/*-----------------------------------------------------------------*/
/*---- End Of File ----*/
