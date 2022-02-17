/*
$Log: armfrmvol.h,v $
Revision 1.20  2004/05/17 13:14:53  emezzine
formatting.

Revision 1.19  2004/02/24 12:34:46  emezzine
change the interpolation methods and function's names

Revision 1.18  2004/02/09 18:13:27  emezzine
added SetBetaCurve().

Revision 1.17  2003/10/31 08:12:50  emezzine
New version to calibrate LogProba Curve.

Revision 1.16  2003/10/07 13:25:13  emezzine
Added a new functions Set/GetInitSigmaCurve.

Revision 1.15  2003/10/06 15:52:07  rguillemot
Ajout d'un spread dans le modele FRM (Modele Mixture)

Revision 1.14  2003/09/17 13:11:19  emezzine
Added itsBeta and its functions.

Revision 1.13  2003/07/11 11:45:04  emezzine
Cut SearchIndex()

Revision 1.12  2003/06/17 14:15:46  emezzine
 Add view virtual function.

Revision 1.11  2003/06/10 17:21:14  emezzine
define  a new virtaul method

Revision 1.10  2003/05/26 09:43:40  emezzine
added  virtual functions

Revision 1.9  2003/04/25 13:59:28  emezzine
Ajout de Log pour commentaire

*/

/*--------------------------------------------------------------------------*/
/*                                                                          */
/* armfrmvol.h: interface for the ARM_FRMVol class.                         */
/*                                                                          */
/*--------------------------------------------------------------------------*/
#ifndef _ARMFRMVOL_H
#define _ARMFRMVOL_H

#ifdef unix
#include <ieeefp.h>
#else
#define finite _finite
#endif


#include "armglob.h"

class ARM_ReferenceValue;
class ARM_Swaption;
class ARM_CapFloor;
class ARM_IRIndex;
class ARM_Matrix;
class ARM_Vector;

ARM_Vector* FuncSmoothOneCurve(const ARM_Vector* Curve,
							   double SmoothCoef,
							   double SmoothPower);


class ARM_FRMVol : public ARM_Object  
{
	private	: 
      
	    // Function parameters to minimize        
        size_t  itsSmoothPower;
	    size_t  itsMarketPower;

        ARM_INDEX_TYPE itsIndexType;
		ARM_ReferenceValue* itsInitShift;// Shift Curve
        ARM_ReferenceValue* itsLogProbaCurve;// Probability of lognormality   

		ARM_ReferenceValue* itsSpreadCurve;

        int itsNbFactors;           
        
    protected:
        /// ITSCeoff is the forward value of the rate  
        ARM_Vector* ITSCeoff;
        
    public :

		ARM_FRMVol(void);
		ARM_FRMVol(const ARM_FRMVol& FRMVol);
		ARM_FRMVol(ARM_ReferenceValue* BetaCurve);
	   ~ARM_FRMVol(void);
	
	    // Standard functions
	    void Init(void);
	    void BitwiseCopy(const ARM_Object* FRMVol);
	    void Copy(const ARM_Object* FRMVol);
        ARM_Object* Clone(void);
		
		ARM_ReferenceValue* GetSpreadCurve(void);
		void SetSpreadCurve(ARM_ReferenceValue* SpreadCurve);

        int GetNbFactors(void);
        void SetNbFactors(int Nbfactors);

		double Mdec(double date);
        double LogProbaParam(double date);
  		double Spread(double date);

        ARM_ReferenceValue* GetLogProbaCurve(void);
        void SetLogProbaCurve (ARM_ReferenceValue* LogProbaCurve);        

        ARM_ReferenceValue* GetMCurve(void);              
        void SetMCurve(ARM_ReferenceValue* MCurve);

        void ComputeShiftCurve(void); 
        void ComputeBetaCurve(void); 
        
        inline ARM_INDEX_TYPE GetIndexType() const {return itsIndexType;}
        void SetIndexType(ARM_INDEX_TYPE indextype) {itsIndexType = indextype;}
            
        // For FRMHWVol
        virtual void SetMeanRevParams(double MainRev);
        virtual	double GetMeanRevParams(void);

        virtual void SetInitSigmaCurve(ARM_ReferenceValue* Curve);
        virtual ARM_ReferenceValue* GetInitSigmaCurve(void);

        virtual void SetCurrentCurve(ARM_ReferenceValue* Curve);
        virtual ARM_ReferenceValue* GetCurrentCurve(void);

        virtual void SetEltToPrice(double x,size_t Idx);         

        virtual ARM_Vector*  VolatilityVector(double t, double T, ARM_CapFloor* Data); 
        virtual ARM_Vector*  VolatilityVector(double t, double T, ARM_Vector* Mu, ARM_Swaption* Data);            
        
        virtual double VarianceForTree(double s,double t);
        virtual double StdDev(double s,double t,double T);
        virtual double CorrelCoef(double T);
        virtual double TimeCoef(double s);
        virtual double VarianceToTime(double var);


        //// For FRMMarkovVol        
        void SetPowers( size_t MarketPower, size_t SmoothPower);
		size_t GetMarketPower(void);
		size_t GetSmoothPower(void);

        virtual    void CurvesToOptimizerParam(ARM_Vector& Data);
        virtual   void UpdateCurvesToCalibrate(ARM_Vector& Try);    
        virtual size_t GetNbParams(void);
        virtual void SetNbParams(size_t nbParmas);
        virtual void SetCurvesToPrice(ARM_Vector* Try); 
        virtual void SetCurvesPrice(ARM_ReferenceValue* ACurve,ARM_ReferenceValue* KCurve);
        virtual void GetCurvesPrice(ARM_ReferenceValue*& ACurve,ARM_ReferenceValue*& KCurve);
        virtual ARM_Vector* SetCurvesToOptimize(ARM_Vector*& UBounds,ARM_Vector*& LBounds, size_t Idx = 0);

        virtual ARM_Vector*  CalculateSmoothError(void); 
        
        virtual void UpdateCurves(size_t Idx =0);
        virtual void OutPutCurves(void);
        virtual ARM_Vector* StdDevVect(double s,double t,double T);
        virtual ARM_Vector* CorrelCoefVect(double T);
        virtual ARM_Vector* VarianceForMC(double s, double t);

        virtual ARM_Matrix* MatrixVolatility(void);
        virtual	ARM_Vector* InterpolateCurves(double t0,double T,double U);

        virtual void View(char* id, FILE* fOut);


};




#endif 
/*--------------------------------------------------------------------------*/
/*---- End of file ----*/
	
