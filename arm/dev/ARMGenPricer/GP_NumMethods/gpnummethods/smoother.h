/*!
 *
 * Copyright (c) CDC IXIS CM July 2003 Paris
 *
 *	\file smoother.h
 *
 *  \brief
 *
 *	\author  J-M Prié
 *	\version 1.0
 *	\date February 2005
 */


#ifndef _INGPNUMMETHODS_SMOOTHER_H
#define _INGPNUMMETHODS_SMOOTHER_H

/// gpbase
#include "gpbase/port.h"
#include "gpbase/rootobject.h"
#include "gpbase/gplinalgtypedef.h"


CC_BEGIN_NAMESPACE( ARM )


struct ARM_SmootherBase : public ARM_RootObject
{
    enum SchedulerType
	{
		DoNothing = 0,
		Linear,
        Quadratic,
        Cubic
    };
    virtual double Compute(const ARM_GP_Vector& exerFct, int exerIdx, double coef, ARM_GP_Vector& smooth) const = 0;
};

struct ARM_SmootherDoNothing : public ARM_SmootherBase
{
    virtual double Compute(const ARM_GP_Vector& exerFct, int exerIdx, double coef, ARM_GP_Vector& smooth) const { return 0.5; }

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_SmootherDoNothing(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_SmootherDoNothing"; }
};
struct ARM_SmootherLinear : public ARM_SmootherBase
{
    virtual double Compute(const ARM_GP_Vector& exerFct, int exerIdx, double coef, ARM_GP_Vector& smooth) const;

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_SmootherLinear(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_SmootherLinear"; }
};
struct ARM_SmootherQuadratic : public ARM_SmootherLinear
{
    virtual double Compute(const ARM_GP_Vector& exerFct, int exerIdx, double coef, ARM_GP_Vector& smooth) const;

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_SmootherQuadratic(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_SmootherQuadratic"; }
};
struct ARM_SmootherCubic : public ARM_SmootherQuadratic
{
    /// Nested derivative function class
    class ExerciseD1Function
    {
        double itsA,itsB,itsC;

    public:
        ExerciseD1Function(double a, double b, double c) : itsA(a),itsB(b),itsC(c) {}
        double operator () ( double x ) const { double x2=x*x; return itsA*x2+itsB*x+itsC; }
    };

    /// Nested function class
    class ExerciseFunction
    {
        double itsA,itsB,itsC,itsD;
        ExerciseD1Function *itsD1Function;

    public:
        ExerciseFunction(double a, double b, double c, double d) : itsA(a),itsB(b),itsC(c),itsD(d),itsD1Function(new ExerciseD1Function(3*a,2*b,c)) {}
        virtual ~ExerciseFunction() { delete itsD1Function;}

        double operator () ( double x ) const { double x2=x*x,x3=x2*x; return itsA*x3+itsB*x2+itsC*x+itsD; }
        ExerciseD1Function* Derivative() const { return itsD1Function; }
    };

    virtual double Compute(const ARM_GP_Vector& exerFct, int exerIdx, double coef, ARM_GP_Vector& smooth) const;

	/// Standard ARM support
	virtual ARM_Object* Clone() const { return new ARM_SmootherCubic(*this); }
    virtual string toString(const string& indent="",const string& nextIndent="") const { return "ARM_SmootherCubic"; }
};


CC_END_NAMESPACE()

#endif

/*----------------------------------------------------------------------------*/
/*---- End of file ----*/

