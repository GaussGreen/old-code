/**
 * @file MultiTweakGroup.hpp
 */

#ifndef EDR_MultiTweakGroup_H
#define EDR_MultiTweakGroup_H

#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Model.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IInstrumentCollection)
FORWARD_DECLARE(MultiTweakGroup)

/**
 * Multi-instrument counterpart to TweakGroup
 *
 * An object from which everything in the risk manager's "world" (instruments,
 * model, market data) is reachable by following references.
 */

class RISKMGR_DLL MultiTweakGroup: public CObject {
public:
    static CClassConstSP const TYPE;

    virtual ~MultiTweakGroup();

    /** Constructor: takes references to inputs, not copies */
    MultiTweakGroup(IInstrumentCollectionSP inst,
                    IModelSP model);

    /** Constructor: takes references to inputs, not copies */
    static MultiTweakGroupSP SP(IInstrumentCollectionSP inst,
                                IModelSP model);

    /** Returns the instruments as SP - todo: merge this with getInstrument */
    IInstrumentCollectionSP getInstruments() const;

    /** Returns the IModel */
    IModel* getModel() const;
    /** Returns the Model as SP */
    IModelSP getModelSP() const;

private:
    MultiTweakGroup(const MultiTweakGroup& rhs);
    MultiTweakGroup& operator=(const MultiTweakGroup& rhs);
    // fields
    IInstrumentCollectionSP insts;
    IModelSP model;

    // methods
    static IObject* defaultMultiTweakGroup();
    MultiTweakGroup();
    static void load(CClassSP& clazz);
};

DRLIB_END_NAMESPACE

#endif

