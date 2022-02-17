/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *  function solution of the Riccati Equation with Piece-Wise functions A(t), B(t) and C(t)
 *	df/dt = A(t) * f(t)² + B(t) * f(t) + C(t)
 *  
 *	\file riccati.h
 *
 *  \brief
 *
 *	\author  A. Triki
 *	\version 1.0
 *	\date October 2005
 */
 
#ifndef _GP_CF_RICCATI_H
#define _GP_CF_RICCATI_H

#include <glob/firsttoinc.h>

#include "gpbase/countedptr.h"
#include "gpbase/functor.h"
#include "gpbase/typedef.h"
#include "gpbase/gpmatrix.h"
#include "gpbase/curve.h"

#include "gpnumlib/odefunctions.h"


#include <vector>
CC_USING_NS(std,vector)

#include <map>
CC_USING_NS(std,map)



CC_BEGIN_NAMESPACE(ARM)

class ARM_Riccati:public ARM_ODEFunc
{
	public :

    /// Internal structures for function description
    struct ARM_FunctData
    {
        double              itsTime;
        double              itsAt;
        double              itsBt;
		double              itsCt;
        double              itsDelta;
        int                 itsNbRoots;
        double              itsRootInf;
        double              itsRootSup;
		// To avoid calculation repetition
		double				itsSumRoots; //(x1+x2)/2
		double				itsDiffRoots; //(x1-x2)/2
    };

    struct ARM_FunctType
    {
        ARM_FunctType(const ARM_FunctData* data=NULL) : itsData(data){}
        const ARM_FunctData*  itsData;
		virtual ~ARM_FunctType () 
		{
		}
    };

	struct ARM_FunctTypeRiccati : public ARM_FunctType
    {
        ARM_FunctTypeRiccati(const ARM_FunctData* data=NULL) : ARM_FunctType(data) , itsFlag(0) {}

        virtual void InitCst(double t, double a){};
        virtual double ComputeA(double t) const{ return 0.;};
		virtual double ComputeB(double t1, double t2);
        double          itsCst;
        short           itsFlag;
        double          itsCACoef;
    };

    struct ARM_Funct0 : public ARM_FunctTypeRiccati
    {
        /// Discriminant < 0 : no solution
        ARM_Funct0(const ARM_FunctData* data=NULL) : ARM_FunctTypeRiccati(data) {}
        virtual void InitCst(double t, double a);
        virtual double ComputeA(double t) const;
    };

    struct ARM_Funct1 : public ARM_FunctTypeRiccati
    {
        /// Discriminant = 0 : one solution
        ARM_Funct1(const ARM_FunctData* data=NULL) : ARM_FunctTypeRiccati(data) {}
        virtual void InitCst(double t, double a);
        virtual double ComputeA(double t) const;
    };

    struct ARM_Funct2 : public ARM_FunctTypeRiccati
    {
        /// Discriminant > 0 : two solutions
        ARM_Funct2(const ARM_FunctData* data=NULL) : ARM_FunctTypeRiccati(data) {}
    };

    struct ARM_Funct2Out : public ARM_Funct2
    {
        /// Solution outside roots
        ARM_Funct2Out(const ARM_FunctData* data=NULL) : ARM_Funct2(data) {}
        virtual void InitCst(double t, double a);
        virtual double ComputeA(double t) const;
    };

    struct ARM_Funct2In : public ARM_Funct2
    {
        /// Solution inside roots
        ARM_Funct2In(const ARM_FunctData* data=NULL) : ARM_Funct2(data) {}
        virtual void InitCst(double t, double a);
        virtual double ComputeA(double t) const;
    };

    struct ARM_Funct2Sup : public ARM_Funct2
    {
        /// Solution equal to the greater root
        ARM_Funct2Sup(const ARM_FunctData* data=NULL) : ARM_Funct2(data) {}

        virtual void InitCst(double t, double a);

        virtual double ComputeA(double t) const;
    };

    struct ARM_Funct2Inf : public ARM_Funct2
    {
        /// Solution equal to the lower root
        ARM_Funct2Inf(const ARM_FunctData* data=NULL) : ARM_Funct2(data) {}
        virtual void InitCst(double t, double a);
        virtual double ComputeA(double t) const;
    };

    typedef ARM_CountedPtr< ARM_FunctTypeRiccati >      ARM_FunctTypeRiccatiPtr;
    typedef ARM_CountedPtr< ARM_FunctType >				ARM_FunctTypePtr;
    typedef map< double,vector< ARM_FunctTypePtr > >	ARM_FunctionsMap;
    typedef ARM_FunctionsMap::iterator                  ARM_FunctionsMapIter;

    /// Function parameters for all sampling intervals
    vector< ARM_FunctData > itsFunctionDatas;

    /// For each maturity, function description over all sampling intervals
    /// Mutable because GenerateFunction() will update this map
    CC_IS_MUTABLE ARM_FunctionsMap itsFunctions;

	public : 
		
	ARM_Curve*		   itsAtCurve;
	ARM_Curve*		   itsBtCurve;
	ARM_Curve*		   itsCtCurve;
	ARM_GP_Vector	   itsSchedule;

	/// Variables for time-dependent Riccati Resolution
	double					itsTempAt;
	double					itsTempBt;
	CC_IS_MUTABLE double	itsTempCt;
	
	/// Accessors
	const ARM_GP_Vector& GetSchedule() { return itsSchedule;}


	/// Standard ARM object support
	ARM_Riccati(const ARM_Curve& At,const ARM_Curve& Bt,const ARM_Curve& Ct,bool IsInitFunctionData = true);
	ARM_Riccati(const ARM_Riccati& rhs);
	virtual ~ARM_Riccati();
    virtual ARM_Riccati& operator = (const ARM_Riccati& rhs);
	virtual void CopyNoCleanUp(const ARM_Riccati& rhs);
	virtual ARM_Object* Clone() const;
	virtual string toString(const string& indent="",const string& nextIndent="") const;

	virtual void GenerateSchedule();

	/// General function for Runge Kutta Resolution
	virtual void derivs(double x, ARM_GP_Vector* yt, ARM_GP_Vector* dyt) const;
	/// Update the values of Temp Parameters TempAt, TempBt, TempCt at each time Step
	virtual void UpdateTempParams(double x) const {};

    /// Collect & save model parameters for function definition
    virtual void InitFunctionDatas();

	/// Accessors
	virtual void SetTempAt (double at ) {itsTempAt = at;}
    /// Build the function solution of the Riccati equation
    virtual ARM_FunctionsMapIter GenerateFunction(double maturity);
	virtual double ComputeA(double time,const vector< ARM_FunctTypePtr >& function) const;

    /// Computation of A(t,T), Solution of the Riccati Equation
    virtual double A(double t, double T);
};


CC_END_NAMESPACE()



#endif

/*---------------------------------------------------------------------------*/
/*---- End of file ----*/