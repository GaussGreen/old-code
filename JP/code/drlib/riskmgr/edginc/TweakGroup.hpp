//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TweakGroup.hpp
//
//   Description : Holds objects to be tweaked for greeks - basically
//                 instrument and model
//
//   Author      : Mark A Robson
//
//   Date        : 9 May 2002
//
//
//----------------------------------------------------------------------------

#ifndef EDR_TWEAKGROUP_HPP
#define EDR_TWEAKGROUP_HPP

#include "edginc/Instrument.hpp"
#include "edginc/Model.hpp"
DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL TweakGroup: public CObject{
public:
    static CClassConstSP const TYPE;

    /** Constructor: takes references to inputs, not copies */
    TweakGroup(const InstrumentSP& inst,
               const IModelSP&     model);

    /** Returns the instrument */
    Instrument* getInstrument() const;

    /** Returns the instrument as SP - todo: merge this with getInstrument */
    InstrumentSP getInstrumentSP() const;

    /** Returns the IModel */
    IModel* getModel() const;
    /** Returns the Model as SP */
    IModelSP getModelSP() const;

private:
    // fields
    InstrumentSP inst;
    IModelSP     model;

    // methods
    static IObject* defaultTweakGroup();
    TweakGroup();
    static void load(CClassSP& clazz);
};

typedef smartPtr<TweakGroup> TweakGroupSP;
#ifndef QLIB_TWEAKGROUP_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<TweakGroup>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<TweakGroup>);
#endif

DRLIB_END_NAMESPACE

#endif

