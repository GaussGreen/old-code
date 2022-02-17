/*
 * $Log: bmc_frm.h,v $
 * Revision 1.14  2004/03/31 15:27:41  rguillemot
 * CMT Bug Fix
 *
 * Revision 1.13  2004/03/02 13:56:16  rguillemot
 * Spread Bug Fix
 *
 * Revision 1.12  2004/02/16 13:58:46  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.11  2004/02/09 08:54:23  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.10  2004/02/04 15:18:37  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.9  2004/01/26 13:39:20  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.8  2004/01/12 07:15:42  rguillemot
 * Replication Convexity Adjustment
 *
 * Revision 1.7  2003/07/10 06:15:16  ebenhamou
 * remove performance warning
 *
 * Revision 1.6  2003/06/26 18:09:34  ebenhamou
 * explicitly say which function to use
 *
 * Revision 1.5  2002/03/01 16:29:22  mab
 * idem
 *
 * Revision 1.3  2001/10/03 16:27:10  smysona
 * doc
 *
 * Revision 1.2  2001/09/25 13:40:30  mab
 * Dans Init(): Correction SetName :  SetName(ARM_BMCFRM);
 *
 */


/*****************************************************/
/*                                                   */
/*    This class is first attempt to deal with       */
/*    basis discounting.                             */
/*    We use the basis adjusted curve in the dynamic */
/*    of the model (the instruments must be modified */
/*    condsequently). For the dynamic, we adjust the */
/*    vol at each time step on each trajectory.      */
/*                                                   */
/*****************************************************/

#ifndef _BMC_FRM_H
#define _BMC_FRM_H


#include "smc_frm.h"


class ARM_BMCFrm : public virtual ARM_SMCFrm
{
private :

// YAN 02/2002
	int itsBasisDiscountFlag;
    ARM_ZeroCurve*   itsUnAdjustedCurve;
    ARM_ZeroCurve*   itsAdjustedCurve;
// YAN 02/2002
	ARM_ZeroLInterpol*   itsAdjWorkingCurve;

    ARM_Forwards*    itsAdjustedFwds;

    ARM_Vector*     itsBasisSpreads;        // spread = basis adjusted - plain libor


public : 

    ARM_BMCFrm(void)
    {
        Init();
    }

// YAN 02/2002
	double GetFirstPeriodBSAdj()
    {return exp( -itsBasisSpreads->Elt(0) );}
    // Construction par calibrateur    
    ARM_BMCFrm(ARM_ZeroCurve* zc, ARM_ZeroCurve* baZc, ARM_FRM_Calibrator* aCalibrator, 
               ARM_Date& endDate, int Ntraj,  
               int MCGeneratorType =   K_MC_SIMPLE, int NumFact = 1,
               bool decal = false, long seed = 10000);



    // Construction par modele anlytique    
     ARM_BMCFrm(ARM_FrmAna* AnaModel, ARM_ZeroCurve* baZc, ARM_Date& endDate,
                int Ntraj, int MCGeneratorType = K_MC_SIMPLE,
                ARM_Vector *CorrelatedIndexes = NULL, 
                ARM_Vector *indexes = NULL, ARM_Matrix *correlations = NULL,
                long seed = 10000);


    // Construction pour autocalibration (1 classe de risque en vol)
	ARM_BMCFrm(ARM_ZeroCurve *zc, ARM_ZeroCurve* baZc, ARM_VolCurve *vol, ARM_VolCurve *smile, 
                            int productType, ARM_Date &endDate,  
							int Ntraj = 1000, int  MCGeneratorType = K_MC_FAURE,
            		       double decay = 0, double slope = 0, double asymptote = 0, int NbFactor = 1, 
            		       ARM_Vector *CorrelatedIndexes = NULL, ARM_Vector *indexes  = NULL,
            			   ARM_Matrix *correlations = NULL, int NoControl = 0, long seed = 10000);
	 
    // Construction pour autocalibration (2 classes de risque en vol)    
	ARM_BMCFrm(ARM_ZeroCurve *zc, ARM_ZeroCurve* baZc, ARM_VolCurve *swoptVol, ARM_VolCurve *swoptSmile, 
               ARM_VolCurve *irgVol, ARM_VolCurve *irgSmile, 
               int productType, ARM_Date &endDate,  
               int Ntraj = 1000, int  MCGeneratorType = K_MC_FAURE,
               double decay = 0, double slope = 0, double asymptote = 0, int NbFactor = 1, 
               ARM_Vector *CorrelatedIndexes = NULL, ARM_Vector *indexes  = NULL,
               ARM_Matrix *correlations = NULL, int NoControl = 0, long seed = 10000);





     ARM_BMCFrm(ARM_BMCFrm &inModel) : ARM_SMCFrm(inModel)
     {
         Init();

         BitwiseCopy((ARM_BMCFrm *) &inModel);
     };

    ~ARM_BMCFrm(void)
     {
        cleanit(itsBasisSpreads);

        cleanit(itsAdjustedFwds);

// YAN 02/2002
		cleanit(itsAdjWorkingCurve);
     };


    ARM_BMCFrm operator = (ARM_BMCFrm &inModel)
    {
        (*this).ARM_SMCFrm::operator =(inModel);

        BitwiseCopy((ARM_BMCFrm *) &inModel);
        
        return (*this);
    }

    void Init(void)
    {
        SetName(ARM_BMCFRM);

// YAN 02/2002
		itsBasisDiscountFlag = 1;
        itsUnAdjustedCurve = NULL;
        itsBasisSpreads    = NULL;
        itsAdjustedCurve   = NULL;
// YAN 02/2002
		itsAdjWorkingCurve = NULL;
        itsAdjustedFwds    = NULL;
// YAN 02/2002
//        itsSpreadCurve     = NULL;
    };

    ARM_Object *Clone(void)
    {
        ARM_BMCFrm *theClone = new ARM_BMCFrm(*this);
        return theClone;
    }


    ARM_ZeroCurve* GetDiscountCurve(void)
    {
// YAN 02/2002
		if(itsBasisDiscountFlag)
			return(itsAdjWorkingCurve);
		else
			return(GetWorkingCurve());
// Fin YAN
    }

    double FirstPeriodDiscountFactor(int nPath);

	double ZeroPrice(double calcDate, double zMaturity, int DomOrFrg);

    double ZeroPriceOnDiscountCurve(double calcDate, double zMaturity, int DomOrFrg=1);

    void BeFittedTo(ARM_Security *sec);

// YAN 02/2002
	void PropagateCurve(int nPath, int nTimeStep);

	bool IsBasisDiscountModel(void)
	{
// YAN 02/2002
		/// to remove the performance warning
		#pragma warning(disable : 4800)
		return (bool)itsBasisDiscountFlag;
	}

// YAN 02/2002
	void SetBasisDiscountFlag(int flag=1)
	{
		itsBasisDiscountFlag = flag;
	}

	void InitTrajectory(int nPath);

	ARM_Vector* GetDiscSprds() { return itsBasisSpreads;}
// Fin YAN 02/2002


	
	
	///////////////////////////////////////////////
	/// to remove ambiguity of multiple inheritance
	////////////////////////////////////////////////

    double ExpectedFwdYield(ARM_Date& fwdDate, 
		ARM_Date& maturity, ARM_Date& payDate, 
		int compMeth = 0, int dayCount = KACTUAL_365, 
		int DomOrFrgRate = 1, int discYC=1,
		int YieldDecomp = K_COMP_PROP,
		double Margin = 0.0)
	{
		return ARM_FrmAna::ExpectedFwdYield(fwdDate, 
			maturity, payDate, compMeth, dayCount, DomOrFrgRate, discYC, YieldDecomp, Margin);
	}
			

    double ExpectedFwdYield(ARM_Date& calcDate, ARM_Date& resetDate, 
		ARM_Date& maturity, ARM_Date& payDate,
		int compMeth = 0, int dayCount = KACTUAL_365, 
		int DomOrFrgRate = 1, int discYC=1,
		int YieldDecomp = K_COMP_PROP,
		double Margin = 0.0)
	{		
		return ARM_FrmAna::ExpectedFwdYield(calcDate, resetDate, maturity, payDate,
			compMeth, dayCount, DomOrFrgRate, discYC, YieldDecomp, Margin);
	}

    double ExpectedFwdYield(double fwdDate, 
		double maturity, double payDate, 
		int compMeth = 0, int dayCount = KACTUAL_365, 
		int DomOrFrgRate = 1, int discYC = 1,
		int YieldDecomp = K_COMP_PROP,
		double Margin = 0.0)
	{		
		return ARM_FrmAna::ExpectedFwdYield(fwdDate, 
			maturity, payDate, compMeth, dayCount, DomOrFrgRate, discYC, YieldDecomp, Margin);
	}

	double ExpectedFwdYield(ARM_Date& currDate, double fwdDate, 
		double maturity, double payDate, 
		int compMeth = 0, int dayCount = KACTUAL_365, 
		int DomOrFrgRate = 1, int discYC = 1,
		int YieldDecomp = K_COMP_PROP,
		double Margin = 0.0)
	{		
		return ARM_FrmAna::ExpectedFwdYield(currDate, fwdDate, 
			maturity, payDate, compMeth, dayCount, DomOrFrgRate,discYC, YieldDecomp, Margin);
	}
	
    double EuroCaplet(ARM_Security *sec, double settlement,
		double startMaturity, double endMaturity, 
		int optionType, double Forwd, double optionStrike, 
		int DomOrFrg=1, int compMeth = 0,
		int dayCount=KACTUAL_365)
	{
		return ARM_FrmAna::EuroCaplet(sec, settlement,
			startMaturity, endMaturity, optionType, Forwd, optionStrike, 
			DomOrFrg, compMeth, dayCount);
	}


	double EuroCaplet(double VolStartMaturity, double resetMaturity,
		double startMaturity, double endMaturity,
		double payMaturity,
		int optionType, double Forwd, double Strike,
		int DomOrFrg, int compMeth, int dayCount,
		bool IsCMS = false,
        bool isTreasury = false,
		int FixFrequency = 1,
		double UnderlyingTenor = 0.0,
		int YieldDecomp = K_COMP_PROP,
		double Margin = 0.0,
		ARM_Currency* ccy = NULL,
		StoreFwdRateAndCapletInfo* StoreInfo = NULL)
	{
		return ARM_FrmAna::EuroCaplet(VolStartMaturity, resetMaturity,
		startMaturity, endMaturity, payMaturity, 
		optionType, Forwd, Strike,
		DomOrFrg, compMeth, dayCount,
		IsCMS, isTreasury, FixFrequency, UnderlyingTenor,
		YieldDecomp, Margin, ccy, StoreInfo);
	}

    double EuroSwaption(double startswap, double matswaption,
		int optionType, double swapfwd,
		double Strike, ARM_Vector* Forwds,
		ARM_Vector* poids, int DomOrFrg=1)
	{
		return ARM_FrmAna::EuroSwaption( startswap, matswaption,
			optionType, swapfwd, Strike, Forwds, poids, DomOrFrg );
	}

    double  EuroSwaption(ARM_Security *sec,double startswap,
		double matswaption,int optionType,
		double swapfwd, double Strike,
		ARM_Vector* Forwds,ARM_Vector* poids,
		int DomOrFrg=1)
	{	
		return ARM_FrmAna::EuroSwaption(sec, startswap, matswaption, optionType,
			swapfwd, Strike, Forwds, poids, DomOrFrg);
	}


    void SetParams(double* para)
	{	
		ARM_FrmAna::SetParams( para );
	}


	void UpdatesModelDatas(ARM_Date* currentDate)
	{
		ARM_FrmAna::UpdatesModelDatas( currentDate );	
	}


	void View(char* id, FILE* ficOut)
	{
		ARM_FrmAna::View(id, ficOut);
	}


    /// explicitly specify which function to use
	double CFPrice(double calcDate, double* cfTerms, 
		double* cfValues, int size, int DomOrFrg)
	{
		return ARM_LogDecal::CFPrice(calcDate, cfTerms, cfValues, size, DomOrFrg);
	}
	
	
	
};

#endif
