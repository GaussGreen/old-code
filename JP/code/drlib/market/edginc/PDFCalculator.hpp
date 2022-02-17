//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PDFCalculator.hpp
//
//   Description : 
//
//   Author      : Andrew J Swain
//
//   Date        : 15 March 2002
//
//
//----------------------------------------------------------------------------

#ifndef PDFCALCULATOR_HPP
#define PDFCALCULATOR_HPP

#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Lattice.hpp"

DRLIB_BEGIN_NAMESPACE

class Correlation;
class PDFCalculator;
class VolProcessed;
class CAsset;
class PDFRequest;
class FXAsset;
class Correlation;

class MARKET_DLL PDFCalculatorException : public ModelException {
public:
    PDFCalculatorException(const string& routine, const string& msg) : 
        ModelException(routine, msg) {}
    PDFCalculatorException(PDFCalculatorException& e, const string &routine):
        ModelException(e, routine) {}
};

// defines an interface to be implemented by concrete vols
class MARKET_DLL IPDFCalculator: virtual public IObject {
public:        
    static CClassConstSP const TYPE;

    IPDFCalculator();
    virtual ~IPDFCalculator();

    virtual PDFCalculator* getPDFCalculator(const PDFRequest* request,
                                            const CAsset*     asset) const = 0;
        
#if 0
    // don't think I need this now
    /* ccy struck version */
    virtual PDFCalculator* getPDFCalculator(const PDFRequest*  request,
                                            const CAsset*      eqAsset,
                                            const FXAsset*     fxAsset,
                                            const Correlation* eqFXCorr) const = 0;
#endif

private:
    static void load(CClassSP& clazz);
};

class MARKET_DLL PDFCalculator : public CObject {
public:
    static CClassConstSP const TYPE;


    /** Calculate the cumulative probability at each strike.
        Each prob is estimated as a "small" call spread (i.e., by tweaking the strike)
        Returns false if if the difference between 2 consecutive probs is negative.
    */
    virtual void probabilities(const DoubleArray& strikes,
                               const DateTime&    maturity,
                               DoubleArray&       probs) const = 0;

    virtual void probabilities(const CLatticeDouble&   strikes,
                               const DateTimeArray&    maturities,
                               CLatticeDouble&         probs) const = 0;


    /** Calculate the "local" density at each strike.
        Each density is estimated as a "small" butterfly (using same strikes 
        as above together with the center strike).
        Returns false if any of the densities is negative
    */
    virtual void localDensity(const DoubleArray& strikes,
                              const DateTime&    maturity,
                              DoubleArray&       density) const = 0;

    /** Calculate the "integrated" density at each strike.
        Each density is computed by taking the difference between 2 consecutive
        probabilities divided by the difference between the 2 consecutive strikes.
        false should be returned if any of the so-computed densities is negative;
    */        
    virtual void integratedDensity(const DoubleArray& strikes,
                                   const DateTime&    maturity,
                                   DoubleArray&       density) const = 0;

    virtual ~PDFCalculator();

protected:
    PDFCalculator(const CClassConstSP& clazz);
private:
    friend class PDFCalculatorHelper;
    PDFCalculator(const PDFCalculator &rhs);
    PDFCalculator& operator=(const PDFCalculator &rhs);
};

typedef smartConstPtr<PDFCalculator> PDFCalculatorConstSP;
typedef smartPtr<PDFCalculator> PDFCalculatorSP;

DRLIB_END_NAMESPACE
#endif
