/*
* $Log: armfrmmcmodel.h,v $
* Revision 1.9  2004/03/16 13:56:19  emezzine
* Use GetFwdResetdates to localise the terminal probabilty.
*
* Revision 1.8  2003/11/25 14:57:58  emezzine
* initialize correctly guassian number
*
* Revision 1.7  2003/11/13 18:16:57  emezzine
* modif constructor
*
* Revision 1.6  2003/10/24 16:29:12  rguillemot
* Prise en compte du spread de la mixture dans MC
*
* Revision 1.5  2003/06/17 14:07:23  emezzine
*  Add ExpectedFwdYield()
*
* Revision 1.4  2003/06/11 13:33:27  emezzine
* New version using ARM_FRMModel only
*
* Revision 1.3  2003/05/26 09:40:04  emezzine
*   Update MC in Multifactors
*
* Revision 1.2  2003/04/30 11:36:39  emezzine
* Pb de version
*
* Revision 1.1  2003/04/28 14:57:19  emezzine
* Initial revision
*
*/
/*--------------------------------------------------------------------------*/
/*                                                                          */
/* armfrmmcmodel.h: interface for the ARM_FRMMCModel class.             */
/*                                                                          */
/*--------------------------------------------------------------------------*/
#ifndef _ARMFRMMCMODEL_H
#define _ARMFRMMCMODEL_H

#include "mcmodel.h"
#include "armfrmvol.h"

class ARM_FRMMCModel : public ARM_MCModel  
{
	// Data :
	private : 
        ARM_PathGenerator *itsRandGen;
        int  itsRandGenType ;
        
        ARM_Matrix* itsCorrelCoef;
        ARM_Vector* itsVarianceForMC;
        ARM_Vector* itsShift;
		ARM_Vector* itsSpread;
        ARM_Vector* itsInterestPeriod;

        int itsStepsIn;
        
        ///Security parameters
        ARM_Vector* itsResetDates;
        ARM_Vector* itsDiffusionDates;
        ARM_Vector* itsEvenementDates;

        /// Fwds Parameters
        ARM_Vector* itsFwdResetDates;
        ARM_Vector* itsFwdStartDates;
        ARM_Vector* itsFwdEndDates;

        ARM_FRMVol* itsFRMVol;
        ARM_Security* itsSecurity;  
        
        ARM_Matrix* itsFwdRates;
        ARM_Matrix* itsZCFwds;

        ARM_PRICER_TYPE itsPricerType;
        int itsdayCount;	

	public : 
		
		ARM_FRMMCModel(void);
	    ARM_FRMMCModel(const ARM_FRMMCModel& ARM_FRMMCModel);

	    ARM_FRMMCModel(ARM_FRMModel* FRMModel,long nbTraj, 
                        long nbStepIn,ARM_PRICER_TYPE PricerType, long mcMethod = K_MC_SIMPLE);

	   ~ARM_FRMMCModel(void);
        void FreeMemory(void);
	 
	    void Init(void);

	    void BitwiseCopy(const ARM_Object* FRMMCModel);
	    void Copy(const ARM_Object* FRMMCModel);
        ARM_Object* Clone(void);

        double ExpectedFwdYield(double startDate,double endDate,
                                double payDate,int compMeth,int dayCount,
                                int domOrFrg,int discYC);

        double TimeStep(int current = 0);
        double GetSecExpiry(void);
        void BeFittedTo(ARM_Security* sec);
        void ComputeInterestPeriod(void);
        ARM_Vector* GetInterestPeriod(void);
        
        void ComputeCovariance(void);   
        void CptOnePath(double * bmInnov);
        void PropagateFwds(double * bmInnov);

        ARM_PathGenerator* Generator(void);
        int IsEvntDate(double& date);       

        double ZeroPrice(double settlement, double maturity, 
                                            int DomOrFrg);

        double ForwardYield(double calcDate, double resetDate, 
                            double maturityDate = 0.0, double yieldMaturity = 0.0,
                            int compMeth = -1, int dayCount= KACTUAL_365,
                            int DomOrFrg = 1); 

        inline ARM_PRICER_TYPE PricerType(void)        
        {
            return(itsPricerType);
        }
	    
};



#endif
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/

