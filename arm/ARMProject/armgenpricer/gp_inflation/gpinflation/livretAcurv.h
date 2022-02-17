/*
 * $Log: livretAcurv.h,v $
 * Revision 1.3  2004/06/09 13:56:36  mab
 * Correction in : int ARM_LivretACurve::GetRateDateNb(ARM_Date& dateIn)
 *
 * Revision 1.2  2004/06/04 08:30:41  mab
 * Passing dates by ref.
 *
 * Revision 1.1  2004/05/11 10:22:54  jpriaudel
 * Initial revision
 *
 */

/*----------------------------------------------------------------------------*/
/*                                                                            */
/* FILE        : livretAcurv.h                                                */
/*                                                                            */
/* DESCRIPTION : Header for the ARM_LivretACurve class, a class for dealing   */
/*               with  livretA rate curves.                                   */
/*                                                                            */
/* DATE        : Tue May 4 2004                                               */
/*                                                                            */
/*----------------------------------------------------------------------------*/

#ifndef _INGPINFLATION_LIVRETACURVE_H
#define _INGPINFLATION_LIVRETACURVE_H

#include "infcurv.h"
#include <util/refvalue.h>
#include <ccy/currency.h>
#include "gpbase/port.h"

CC_BEGIN_NAMESPACE( ARM )
        
class ARM_LivretACurve : public ARM_ZeroCurve
{

    private:
     
        long itsFlagRounding;
    
        void Init(void);
    
        void BuildCurve(ARM_InfCurv* infCurv, ARM_ZeroCurve* euriborCurv,
                        ARM_ReferenceValue* fixingLivretA, 
                        ARM_ReferenceValue* fixingEuribor,
						long monthForAugustId = K_MONTH_DEFAULT,
						long monthForFebruaryId = K_MONTH_DEFAULT);

        void BuildCurve(ARM_InfCurv* infCurv, ARM_ZeroCurve* euriborCurv,
                        ARM_ResetManager* LivretAResetManager, 
                        ARM_ResetManager* EuriborResetManager,
						long monthForAugustId = K_MONTH_DEFAULT,
						long monthForFebruaryId = K_MONTH_DEFAULT);
		
        double GetBeginCalculSchedule(void);

        double Rounding(double valeurIn);

        double RateLivretA(ARM_InfCurv* infCurv, 
                           ARM_ZeroCurve* euriborCurv,
                           ARM_Date& dateIn, 
                           ARM_ReferenceValue* fixingEuribor);

        double RateLivretA(ARM_InfCurv* infCurv, 
                           ARM_ZeroCurve* euriborCurv,
                           ARM_Date& dateIn, 
                           ARM_ResetManager* EuriborResetManager);

        double RateLivretA(ARM_InfCurv* infCurv, 
                           ARM_ZeroCurve* euriborCurv, 
                           double dateJulIn,
                           ARM_ReferenceValue* fixingEuribor);
						   
        double LivretACPIRatio(ARM_Date& numDate, ARM_InfCurv* infCurv);
                                
        // get the ratio of the CPI
        double LivretACPIRatio(double numJulDate, ARM_InfCurv* infCurv);

        // get value of CPI from a inflation curve
        double ValueCPI(const ARM_Date& resetDate, ARM_InfCurv* infCurv);
    
        ////    functions for txEuribor
        double AverageZCRate(ARM_ZeroCurve* euriborCurv, 
                             const ARM_Date& dateEuribor,
                             ARM_ReferenceValue* fixingEuribor);

        ////    functions for txEuribor
        double AverageZCRate(ARM_ZeroCurve* euriborCurv, 
                             const ARM_Date& dateEuribor,
                             ARM_ResetManager* EuriborResetManager);

        double LivretAInterpolate(double yearTerm);
        double LivretAInterpolate(double yearTerm, int& nb);

        double LivretAInterpolate(const ARM_Date& rateDate);
        double LivretAInterpolate(const ARM_Date& rateDate, int& nb);

        void SetFlagRounding(long flag);
        long GetFlagRounding()const;


    public:
    
        // constructor
        ARM_LivretACurve(void);
    
        ARM_LivretACurve(ARM_Date& asOfDate, ARM_InfCurv* infCurv, 
                         ARM_ZeroCurve* euriborCurv,
                         long flagRounding, 
  					     ARM_ResetManager* ResetManager=NULL,
						 ARM_ReferenceValue* fixingLivretA=NULL, 
                         ARM_ReferenceValue* fixingEuribor=NULL,
						 long monthForAugustId = K_MONTH_DEFAULT,
						 long monthForFebruaryId = K_MONTH_DEFAULT);

        ARM_LivretACurve(ARM_Date& asOfDate, ARM_InfCurv* infCurv, 
                         ARM_ZeroCurve* euriborCurv,
                         long flagRounding, 
  					     ARM_ResetManager* InflResetManager=NULL,
						 ARM_ResetManager* LivretAResetManager=NULL, 
                         ARM_ResetManager* EuriborResetManager=NULL,
						 long monthForAugustId = K_MONTH_DEFAULT,
						 long monthForFebruaryId = K_MONTH_DEFAULT);
		
        // Copy constructor
        ARM_LivretACurve(const ARM_LivretACurve& srcLivreACurv);

        // destructor
        ~ARM_LivretACurve(void);
    
        // operator = 
        ARM_LivretACurve& operator = (const ARM_LivretACurve& srcLivreACurv);
    
        // function for ARM_Objet compatibility
        void BitwiseCopy(const ARM_Object* src);
        
        virtual void Copy(const ARM_Object* src);    

        virtual ARM_Object* Clone(void);

        virtual void View(char* id = NULL, FILE* ficOut = NULL);

        int GetRateDateNb(ARM_Date& dateIn);

        double DiscountFunction(double yearTerm);

        double DiscountFunctionLivretA(double yearTerm, int compMeth);
        
        double DiscountFunctionLivretA(const ARM_Date& dateIn, int compMeth);

        virtual double DiscountYield(double yearTerms, int compMeth = 0);

        double DiscountYield(ARM_Date& maturity, int compMeth = 0);
        
        virtual double ForwardYield(double yearTerm1, double yearTerm2,
                                    int compMeth = 0);
        
        double ForwardYield(ARM_Date& maturity1, ARM_Date& maturity2,
                            int compMeth = 0);
};

CC_END_NAMESPACE()

#endif
/*----------------------------------------------------------------------------*/
/*---- End Of File ----*/
