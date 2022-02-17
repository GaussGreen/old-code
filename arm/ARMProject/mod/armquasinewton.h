/*
*$Log: armquasinewton.h,v $
*Revision 1.10  2004/05/12 14:03:18  emezzine
*Added a new parameter to control vega perturbation
*
*Revision 1.9  2004/02/17 10:57:08  emezzine
* new version for optimiser
*
*Revision 1.8  2003/10/07 13:26:16  emezzine
* Added a new argument for constructor to initilise epsilon to calculate gradient
*
*Revision 1.7  2003/07/31 08:09:57  emezzine
**** empty log message ***
*
*Revision 1.6  2003/07/25 15:17:43  emezzine
* Modif of constrector of ARMQuasinewton
*
*Revision 1.5  2003/06/10 17:20:44  emezzine
*new version
*
* Revision 1.1  2003/02/25 16:34:29  Ezzine
* Initial revision
*
*/
/*--------------------------------------------------------------------------*/
/*                                                                          */
/* armquasinewton.h: interface for the ARM_QuasiNewton class.               */
/*                                                                          */
/*--------------------------------------------------------------------------*/
#ifndef _ARMQUASINEWTON_H
#define _ARMQUASINEWTON_H



#include "calibration.h"

class  ARM_FRMModel;

class ARM_QuasiNewton : public ARM_Calibration
{
    private :

        ARM_Vector*    itsLowerBound;
        ARM_Vector*    itsUpperBound;

        double         itsPrecision;
        double         itsVegaLevel;
        double         itsInitGuess;
        size_t         itsNbMaxIter;
        double         itsEpsilon;
 
    public :

		void Init(void); 

        ARM_QuasiNewton(void)
		{
			Init();
		}

        ARM_QuasiNewton(double tol,
                        long maxIter   = ARM_DEF_MAX_ITER,
                        ARM_Vector* LB = NULL,
                        ARM_Vector* UB = NULL); 

        ARM_QuasiNewton(const ARM_QuasiNewton& calib); 

       ~ARM_QuasiNewton(void);
	    
		void BitwiseCopy(const ARM_Object* calib);
		void Copy(const ARM_Object* calib);
		ARM_Object* Clone(void);

        void SetPrecision(double Precision) {itsPrecision = Precision;};
        void SetFistDerivativeLevel(double VegaLevel) {itsVegaLevel = VegaLevel;};
        void SetInitGuess(double Initguess){itsInitGuess = Initguess;};
        void SetEpsilon(double epsilon){itsEpsilon = epsilon;};
        inline double GetEpsilon() const {return itsEpsilon;};

		double Funcd(double x,ARM_FRMModel* model);
		double Func(double x, ARM_FRMModel* model);
        double f(double x) { return 0.0;};

		ARM_Matrix* Fun2mjac(ARM_Vector* x, size_t Idx);
		ARM_Vector* Fun2m(ARM_Vector* x, size_t Idx);

		//Optimizer (Brackting, Newton Raphson)

		void Zbrak(double x1, double x2,size_t n, 
		     	   ARM_Vector* xb1,ARM_Vector* xb2, 
				   size_t& nroot,size_t Idx);

        ARM_Vector* MultiNewton(ARM_Vector* X0, size_t Idx);
		
		double BNROptimizer(ARM_FRMModel* model);
		ARM_Vector* NROptimizer2m(ARM_Vector* X0,size_t Idx);

        double BrentMinimizer(const double ax, const double bx, const double cx,
                             double tol,double& xmin);
		
        };


#endif 
/*--------------------------------------------------------------------------*/
/*---- End Of File ----*/
