//////////////////////////////////////////////////////////////////////////////////----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Interpolator.hpp
//
//   Description : 
//
//   Date        : 06 June 2002
//
//
//   $Log: Interpolator.hpp,v $
//   Revision 1.7  2004/06/02 17:24:08  mvenardos
//   Initialized n to 0 in the default constructors
//
//   Revision 1.6  2004/06/02 14:35:57  mvenardos
//   x and y are now DoubleArray and not DoubleArraySP
//
//   Revision 1.5  2004/06/01 16:47:57  mvenardos
//   Redesign interpolant so that it can have either
//   a virtual base or non virtual. This allows for
//   faster interpolants when their type is known at
//   compile time
//
//   Revision 1.4  2003/05/21 13:35:41  rguichar
//   Made some methods const, as they should
//
//   Revision 1.3  2002/10/08 14:41:23  evenardo
//   Added constructor to Interpolant
//
//   Revision 1.2  2002/10/08 08:46:03  evenardo
//   Interpolant is a CObject
//
//   Revision 1.1  2002/06/07 14:12:15  evenardo
//   Interpolator and Interpolant interfaces
//
//
//----------------------------------------------------------------------------
#include "edginc/Object.hpp"
#include "edginc/Format.hpp"
#include "edginc/smartPtr.hpp"
#include <math.h>
#include "edginc/AtomicArray.hpp"

#ifndef EDR_INTERPOLATOR_HPP
#define EDR_INTERPOLATOR_HPP

DRLIB_BEGIN_NAMESPACE


/** Interface for classes that return an interpolant */
class UTIL_DLL Interpolator: virtual public IObject{
public:
    static CClassConstSP const TYPE;

    /** Interface for interpolants */
    class UTIL_DLL Interpolant: virtual public IObject{
    public:
        static CClassConstSP const TYPE;
        
        /** Returns interpolated value for xx */
        virtual double value(double xx) const = 0;

        /** Returns nth derivative evaluated at xx */
        virtual double value(double xx,
                             int    deriv) const = 0;

        /** Returns nth derivative evaluated at each
            value of xvec */
        virtual void value(const CDoubleArray& xvec,           
                           int                 deriv,
                           CDoubleArray&       valuevec) const = 0;


		virtual void editY(int, double) {};

    private:
        /** For Reflection */
        static void load(CClassSP& clazz);
    };

    typedef smartPtr<Interpolant> InterpolantSP;
    typedef smartConstPtr<Interpolant> InterpolantConstSP;

    /** Produces an interpolant from the arrays of x and y */
    virtual InterpolantConstSP computeInterp(const CDoubleArray& xdata,
                                             const CDoubleArray& fdata) const = 0;

    /** Class used to specify non-virtual interpolation methods. */
    class UTIL_DLL NullBase: virtual public IObject {
    public:
        static CClassConstSP const TYPE;
        
    private:
        /** For Reflection */
        static void load(CClassSP& clazz);
    };
    
private:
    /** For reflection */
    static void load(CClassSP& clazz);
};

////////////////////////////////////////////////////////////////////////////////

/** Concrete class that contains data relevant to interpolation.
    Template parameter Base controls whether methods will be
    virtual or not. Setting _TBase = Interpolant (default) will
    produce virtual methods. Setting _TBase = InterpolantNonVirtual
    will result in non virtual methods. */
template<class _TBase = Interpolator::Interpolant>
class InterpolantBase: public CObject, virtual public _TBase {
public:
    static CClassConstSP const TYPE;
    
    /** Default constructor */
    InterpolantBase(): CObject(TYPE), n(0) {}
    
    /** Constructor for reflection */
    InterpolantBase(const CClassConstSP& clazz): CObject(clazz), n(0) {}
    
    /** Full constructor from data */
    InterpolantBase(const CClassConstSP& clazz,
                    const CDoubleArray&  xx,
                    const CDoubleArray&  yy): CObject(clazz), x(xx), y(yy) {

        static string routine = "InterpolantBase::InterpolantBase";
        n = xx.size();
        if(n != yy.size()) {
            throw ModelException(routine, "Sizes of x and y arrays do not match");
        }
    }

    /** Allows access to x array */
    const CDoubleArray& getXarray() const {
        return x;
    }
    
    /** Allows access to y array */
    const CDoubleArray& getYarray() const {
        return y;
    }
    
    /** Returns size of arrays */
    int size() const {
        return n;
    }

protected:
    mutable CDoubleArray x; //!<    Array of x values
    mutable CDoubleArray y; //!<    Array of y values
    int                  n; //!<    Size of arrays
    
private:
    /** For reflection */
    static void load(CClassSP& clazz) {
        REGISTER(InterpolantBase<_TBase>, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(_TBase);
        FIELD(x, "X array");
        FIELD(y, "Y array");
        FIELD(n, "size");
    };
};

typedef InterpolantBase<Interpolator::Interpolant> InterpolantVirtual;
typedef InterpolantBase<Interpolator::NullBase>    InterpolantNonVirtual;

typedef smartPtr<Interpolator> InterpolatorSP;
typedef smartConstPtr<Interpolator> InterpolatorConstSP;

typedef array<Interpolator::InterpolantSP, Interpolator::Interpolant> InterpolantArray;
typedef smartPtr<InterpolantArray> InterpolantArraySP;

DRLIB_END_NAMESPACE

#endif
