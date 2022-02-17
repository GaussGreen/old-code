//----------------------------------------------------------------------------
//
//   File        : CompositeForwardRatePricer.hpp
//
//   Description : A container for models capable
//                 of generating and pricing deterministic forward rates
//                 Currently has slots for an IR and a CR form
//                 with ClosedForm being available by default
//
//----------------------------------------------------------------------------

#ifndef QLIB_COMPOSITEFORWARDRATEPRICER_HPP
#define QLIB_COMPOSITEFORWARDRATEPRICER_HPP

#include "edginc/IForwardRatePricer.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL CompositeForwardRatePricer : public CObject, // CModel ??
                                               virtual public IForwardRatePricer
{
public:

    static CClassConstSP const TYPE;

    //allow public construction
    CompositeForwardRatePricer();

    virtual void Forward(CObject *instrument, double *result);

/*
    class IProduct
    {
        //the method to return the forward rate
        virtual double fwd() const = 0;
    };

    //Interface that target objects must implement
    class IIntoProduct
    {
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedFormForwardRatePricer* pricer) = 0;
    };
*/
private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();

    //-------
    // Fields
    //-------

    IForwardRatePricerSP irModel; //for interest rate forwards
    IForwardRatePricerSP crModel; //for credit forwards
};

DECLARE(CompositeForwardRatePricer)

DRLIB_END_NAMESPACE

#endif
