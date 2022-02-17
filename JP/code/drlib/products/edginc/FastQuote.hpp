//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FastQuote.cpp
//
//   Description : model tuned for real time vanilla options pricing.
//
//   Author      : Ning Shen
//
//   Date        : 24 June 2002
//
//
//----------------------------------------------------------------------------

#ifndef FAST_QUOTE_H
#define FAST_QUOTE_H

#include "edginc/config.hpp"
#include "edginc/FastQuoteEnv.hpp"
#include "edginc/ClientRunnable.hpp"
#include "edginc/RegressionTest.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

/** fast quote class computes array of options */
class PRODUCTS_DLL FastQuote: public CObject, 
                 public ClientRunnable, 
                 public IRegressionTest {
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);

    static IObject* defaultFastQuote(){
        return new FastQuote();
    }

    // EdrAction version of addin
    virtual IObjectSP run();

    /** calculate options */
    CResultsArraySP ComputeQuote() const;

    /** addin for calculate options */
    static IObjectSP ComputeQuoteAddin(FastQuote*);

    /** for regression run */
    virtual IObjectSP runTest() const;

private:
    FastQuote():  CObject(TYPE){}

    double                  spot;
    FastQuoteInputArraySP   input;
    FastQuoteEnvSP          quoteEnv;
};


DRLIB_END_NAMESPACE

#endif
