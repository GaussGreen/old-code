/*
$Log: armfrmhwvol.h,v $
Revision 1.21  2004/02/24 12:28:25  emezzine
 Clean and modif  constructor

Revision 1.20  2004/02/09 18:12:02  emezzine
modif constructor

Revision 1.19  2003/11/25 14:54:40  emezzine
global factor normalization

Revision 1.18  2003/10/31 08:12:20  emezzine
New version to calibrate LogProba Curve.

Revision 1.17  2003/10/27 09:27:42  emezzine
 Correct a bug and improvment

Revision 1.16  2003/10/07 13:24:47  emezzine
Added a new functions Set/GetInitSigmaCurve.

Revision 1.15  2003/10/07 10:33:03  rguillemot
Spread parametre par defaut

Revision 1.13  2003/07/07 09:22:27  emezzine
add Log

*/
/*--------------------------------------------------------------------------*/
/*                                                                          */
/* armfrmmarkovvol.h interface for the ARM_FRMHWVol class.                  */
/*                                                                          */
/*--------------------------------------------------------------------------*/
#ifndef _ARMFRMHWVOL_H
#define _ARMFRMHWVOL_H


#pragma warning(disable :4786)



#include "armfrmvol.h"

class ARM_ReferenceValue;

class ARM_FRMHWVol : public ARM_FRMVol  
{
    private :

		 ARM_Date itsAsOfDate;
  		 ARM_ReferenceValue* itsInitSigmaCurve;
	
		 ARM_ReferenceValue* itsCurvePrice;
		 ARM_ReferenceValue* itsCurveCurrent;
         ARM_ReferenceValue** itsCurveFactors;
		 double itsMeanRev;

         size_t itsVolType;
         
    public :
			ARM_FRMHWVol(void);
			ARM_FRMHWVol(const ARM_FRMHWVol& FRMHWVol);            
			ARM_FRMHWVol(ARM_Date AsOfDate,
                         ARM_INDEX_TYPE liborType,
                         ARM_Portfolio* PfSigma,
                         ARM_Portfolio* PfLogProba,
					     ARM_Vector* CurveInit,
						 ARM_Vector* LogProbaCurveInit, 
                         int flagBetaOrShift,
					     double MeanRev,
                         long nbfactor,
                         ARM_Matrix* correlMatrix,
                         size_t VolType = 0,
						 ARM_Vector* SpreadCurveInit = NULL);

			ARM_FRMHWVol(ARM_Date AsOfDate,
                         ARM_INDEX_TYPE liborType,
                         ARM_Portfolio* pf1,
						 int nbProductsPerMaturity,
					     ARM_Vector* CurveInit,
						 ARM_Vector* MCurveInit,                        
						 ARM_Vector* SpreadCurveInit,
					     double MeanRev,
                         long nbfactor,
                         ARM_Matrix* correlMatrix,
                         size_t VolType = 0);

			~ARM_FRMHWVol(void);

            void FreeRefValues(void);

			// Standard functions
			void Init(void);
			void BitwiseCopy(const ARM_Object* ARM_FRMHWVol);
			void Copy(const ARM_Object* ARM_FRMHWVol);
			ARM_Object* Clone(void);

			// FRMHWVol 
            void SetInitSigmaCurve(ARM_ReferenceValue* Curve);
            ARM_ReferenceValue* GetInitSigmaCurve(void);
			
			void SetCurvesPrice(ARM_ReferenceValue* ACurve);
			ARM_ReferenceValue* GetCurvesPrice(void);

			void SetCurrentCurve(ARM_ReferenceValue* Curve);
			ARM_ReferenceValue* GetCurrentCurve(void);

			void SetMeanRevParams(double MainRev);
			double GetMeanRevParams(void);

			size_t GetVolType();

			void SetEltToPrice(double x,size_t Idx);
            
			ARM_Vector* SetCurvesToOptimize(size_t Idx);

          	double buildVolatility(ARM_Security* security,
                                   double startdate);

			ARM_Vector* VolatilityVector(double t,
                                         double T,
                                         ARM_CapFloor* Data);

			ARM_Vector* VolatilityVector(double t,
                                         double T,
                                         ARM_Vector* Mu,
                                         ARM_Swaption* swaption);

			void UpdateCurves(size_t Idx);
			
			//ARM_Vector* Interpolate(double t0, double T, double U);
            ARM_Vector* InterpolateCurves(double t0,
                                          double T,
                                          double U);

            // Functions for Tree method
            double VarianceForTree(double s, double t);
            double CorrelCoef(double T);
            double TimeCoef(double s);
            double StdDev(double s,double t,double T);
			double VarianceToTime(double var);

            // Functions for MC method
            ARM_Vector* StdDevVect(double s,double t,double T);
            ARM_Vector* CorrelCoefVect(double T);
            void PrepareFactorsFromACP(ARM_Matrix* correlMat,
                                       double asofdate);

			void Integratepower2Curve(ARM_ReferenceValue* Curve,
                                      double MeanRev,
                                      size_t Idx);

            void View(char* id, FILE* ficOut);

};






#endif 
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/
	