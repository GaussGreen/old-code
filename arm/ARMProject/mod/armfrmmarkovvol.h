/*
 * $Log: armfrmmarkovvol.h,v $
 * Revision 1.7  2003/08/06 10:47:29  emezzine
 *  Modif de constructeur
 *
 * Revision 1.6  2003/06/10 17:13:54  emezzine
 * update new version
 *
 * Revision 1.5  2003/05/26 09:39:42  emezzine
 * Update all method of calibration
 *
 * Revision 1.4  2003/04/15 15:50:08  jpriaudel
 * modif dans le commentaire RCS
 * pour compil
 *
 * Revision 1.3  2003/04/11 07:54:39  emezzine
 * Enlever qui bloquait la compilation
 *
 *Revision 1.2  2003/04/04 15:41:30  emezzine
 *ajout 
 *
 */

/*--------------------------------------------------------------------------*/
/*                                                                          */
/* armfrmmarkovvol.h interface for the ARM_FRMMarkovVol class.             */
/*                                                                          */
/*--------------------------------------------------------------------------*/
#ifndef _ARMFRMMARKOVVOL_H
#define _ARMFRMMARKOVVOL_H


#include "armfrmvol.h"
#include "refvalue.h"
#include "zerocurv.h"
#include "volcurv.h"

#include "genminimizer.h"



double SmoothFunction(void);
double ObjectFunction(void);

class ARM_FRMMarkovVol : public ARM_FRMVol  
{
    private :

		 ARM_Date itsAsOfDate;
	
		 ARM_ReferenceValue* itsACurvePrice;
		 ARM_ReferenceValue* itsKCurvePrice;

         ARM_ReferenceValue* itsACurve;
		 ARM_ReferenceValue* itsKCurve;

		 ARM_ReferenceValue* itsACurveToCalibrate;
		 ARM_ReferenceValue* itsKCurveToCalibrate;
	
		 ARM_Vector* itsMainRevA;
		 ARM_Vector* itsMainRevK;

		 ARM_Matrix* itsBounds;

	     // Parameters to minimize

		 double itsSmoothACurve;
		 double itsSmoothKCurve;
		 
		 size_t itsNbParms;

         
    public :
        
		 ARM_FRMMarkovVol(void);
	
		 ARM_FRMMarkovVol(const ARM_FRMMarkovVol& FRMVol);
	
	     ARM_FRMMarkovVol(double NbFactors,ARM_Date AsOfDate,
						  ARM_Vector* ACalibrationSchedule,
						  ARM_Vector* KCalibrationSchedule,
						  ARM_ReferenceValue* AInitialise,
						  ARM_ReferenceValue* XlKInitialise,
						  ARM_Vector* MeanReversionA,
						  ARM_Vector* MeanReversionK,
                          ARM_Vector* SchedToPrice,
                          ARM_Vector* Mdec,
						  ARM_Matrix* Bounds,
						  ARM_Vector* PowerData,
						  ARM_Vector* SmoothData,
                          ARM_IRIndex* Index);

		 ~ARM_FRMMarkovVol(void);

	      // Standard functions
	      void Init(void);
	      void BitwiseCopy(const ARM_Object* ARM_FRMMarkovVol);
	      void Copy(const ARM_Object* ARM_FRMMarkovVol);
          ARM_Object* Clone(void);
    
	      // FRMVol 
		  void SetSmoothParameters(double SmoothACurve,double SmoothKCurve);
		  void GetSmoothParameters(double& SmoothACurve, double& SmoothKCurve);
		  
		  void SetCurvesPrice(ARM_ReferenceValue* ACurve,ARM_ReferenceValue* KCurve);
		  void GetCurvesPrice(ARM_ReferenceValue*& ACurve,ARM_ReferenceValue*& KCurve);

		  void SetCurvesToCalibrate(ARM_ReferenceValue* ACurve,ARM_ReferenceValue* KCurve);
		  void GetCurvesToCalibrate(ARM_ReferenceValue*& ACurve,ARM_ReferenceValue*& KCurve);

          void SetACurve(ARM_ReferenceValue* ACurve);
          void SetKCurve(ARM_ReferenceValue* KCurve);
          ARM_ReferenceValue* GetACurve(void);
          ARM_ReferenceValue* GetKCurve(void);

		  void SetMeanRevParams(ARM_Vector* MainRevA,ARM_Vector* MainRevK);
		  void GetMeanRevParams(ARM_Vector*& MainRevA,ARM_Vector*& MainRevK);


		  void SetCurvesToPrice(ARM_Vector* vect);
		  ARM_Vector* SetCurvesToOptimize(ARM_Vector*& UBounds,ARM_Vector*& LBounds,size_t);

		  size_t GetNbParams(void);
		  void   SetNbParms(size_t nbParams);

		  size_t GetNbParamstoNag(void);

		 ARM_Vector* VolatilityVector(double t,double T,
									  ARM_CapFloor* Data);

		 ARM_Vector* VolatilityVector(double t,double T,
									  ARM_Vector* Mu,
									  ARM_Swaption* swaption);

		 
		 ARM_Vector* CalculateSmoothError(void);

		 void UpdateCurves(size_t Idx =0);

         ARM_Vector* StdDevVect(double s,double t,double T);
         ARM_Vector* CorrelCoefVect(double T);
         ARM_Vector* VarianceForMC(double s, double t);

     // Functions for Tree method one Factor
        double VarianceForTree(double s, double t);
        double CorrelCoef(double T);
        double TimeCoef(double s);
        double StdDev(double s,double t,double T);
		double VarianceToTime(double var);
     
		 ARM_Vector* InterpolateCurves(double t0, double T, double u);
		 	    		 
		/* double IntegrateCurves(double tmin,double tmax,
											 double T1,double T2);*/

	     void Integratepower2Curve(ARM_ReferenceValue* Curve,ARM_Vector* MeanRev);

		 size_t GetNbFactors(void);
	 
};






#endif 
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/
	
