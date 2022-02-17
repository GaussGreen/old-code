//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AggregateModel.hpp
//
//   Description : Way of pricing an array of instruments in one go
//
//   Author      : Mark A Robson
//
//   Date        : 30 Oct 2003
//
//
//----------------------------------------------------------------------------

#ifndef EDR_COMPOSITE_MODEL_HPP
#define EDR_COMPOSITE_MODEL_HPP
#include "edginc/MarketObject.hpp"
#include "edginc/Results.hpp"
#include "edginc/Control.hpp"
#include "edginc/Model.hpp"
#include "edginc/Instrument.hpp"

DRLIB_BEGIN_NAMESPACE
class MarketData;
class SensControl;

/** An AggregateModel is an object which supports the pricing of an
    array of [certain types of] instruments in one go. */
class RISKMGR_DLL IAggregateModel: public virtual IModel {
public:
    static CClassConstSP const TYPE;

    virtual ~IAggregateModel();

    /** Create an Instrument which can be passed to the price() method
        below. Note that the instruments (nor this model) will not
        have been populated with their market data at this point. See
        comments on price() regarding use of weights */
    virtual CInstrument* createAggregateInstrument(
        const CInstrumentArray& instruments,
        const DoubleArray&      weights) = 0;

    /** calculate n prices and store result in CResult. The 
        compositeInstrument passed in is that returned from 
        createAggregateInstrument. Note that results for each instrument
        (in the array) must be returned as if the weight was 1.0 ie do
        not scale the price by weight. */
    virtual void price(CInstrument*         compositeInstrument, 
                       Control*             control, 
                       CResultsArray&       results) = 0;

private:
    class MyAggregateModel;
    class MyAggInst;
    static void load(CClassSP& clazz);
};


// typedef for smart pointers to Model
typedef smartConstPtr<IAggregateModel> IAggregateModelConstSP;
typedef smartPtr<IAggregateModel> IAggregateModelSP;

DRLIB_END_NAMESPACE
#endif
