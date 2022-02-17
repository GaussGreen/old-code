//----------------------------------------------------------------------------
//
//   Group       : QR Equities
//
//   Filename    : VegaMatrixLite.hpp
//
//   Description : Class to drive computation of 'Lite' Vega Matrix 
//                 sensitivity
//
//   Author      : Jon Dee
//
//   Date        : 6 October 2006
//
//
//----------------------------------------------------------------------------

#ifndef EDR_VEGAMATRIXLITE_HPP
#define EDR_VEGAMATRIXLITE_HPP

#include "edginc/Sensitivity.hpp"
#include "edginc/TweakQualifierID.hpp"
#include "edginc/Expiry.hpp"
#include "edginc/Untweakable.hpp"

DRLIB_BEGIN_NAMESPACE

/** Instruments supporting the VegaMatrixLite computation 
    should implement this interface. */
class RISKMGR_DLL ISupportVegaMatrixLite: virtual public IObject {
public:
    static CClassConstSP const TYPE;
    
    /** Gives the instrument a chance to throw a ModelException explaining
        why VegaMatrixLite cannot be computed e.g. ccy struck assets, fwdStarting
        options etc. */
    virtual void avoidVegaMatrixLite(const IModel* model) = 0;

    static void load(CClassSP &clazz);
};

/** Lite Vega Matrix controller: bypasses sensitivity framework */
class RISKMGR_DLL VegaMatrixLite :  public Sensitivity
{
    static void load(CClassSP& clazz);

    //Same name as ScalarShift::shiftSize, and plays the same role.
    //Will be used differently: instrument queries shiftSize directly
    //when computing the lite vega matrix
    double shiftSize;

public:

    double getShiftSize() const;

    static CClassConstSP const TYPE;
    const static string NAME;

    /** Implement Sensitivity interface */
    virtual void calculate(TweakGroup* tweakGroup, Results* results);

    /** identifies the name used for storing associated results in the output*/
    virtual const string& getSensOutputName() const;

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
    an approximately continuous one. Returns false */
    virtual bool discreteShift() const;

    VegaMatrixLite();
    virtual ~VegaMatrixLite(void) {}
};

DRLIB_END_NAMESPACE

#endif
