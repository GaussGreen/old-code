//----------------------------------------------------------------------------
//
//   File        : ClosedFormForwardRatePricer.hpp
//
//   Description : A closed form model capable
//                 of generating and pricing deterministic forward rates
//
//----------------------------------------------------------------------------

#ifndef QLIB_CLOSEDFORMFORWARDRATEPRICER_HPP
#define QLIB_CLOSEDFORMFORWARDRATEPRICER_HPP

#include "edginc/IForwardRatePricer.hpp"
#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL ClosedFormForwardRatePricer : public CObject, // CModel ??
                                                public IForwardRatePricer
{
public:

    static CClassConstSP const TYPE;

    //allow public construction
    ClosedFormForwardRatePricer();

    //Forward method
    virtual void Forward(CObject *instrument, double *result);

    class IProduct
    {
    public:
        //the method to return the forward rate
        virtual double fwd(ClosedFormForwardRatePricer* model) const = 0;
    };

    //Interface that target objects must implement
    class IIntoProduct
    {
    public:
        static CClassConstSP const TYPE;
        virtual IProduct* createProduct(ClosedFormForwardRatePricer* pricer) const = 0;
    };

private:
    static void load(CClassSP& clazz);
    static IObject* defaultConstructor();
};

DECLARE(ClosedFormForwardRatePricer)

DRLIB_END_NAMESPACE

#endif
