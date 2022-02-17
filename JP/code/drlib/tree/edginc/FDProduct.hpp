//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : FDProduct.hpp
//
//   Description : FD/tree product base class.
//
//   Author      : Ning Shen
//
//   Date        : February 8, 2005
//
//----------------------------------------------------------------------------
#ifndef FDPRODUCT_HPP
#define FDPRODUCT_HPP

#include "edginc/TreeSlice.hpp"
#include "edginc/IRVolBase.hpp"
#include <typeinfo>

//!!! this is to make old tree work, to be reviewed/removed
#include "edginc/VolRequestLN.hpp"
#include "edginc/YieldCurve.hpp"

DRLIB_BEGIN_NAMESPACE

class FDModel;
class FDProduct;
typedef refCountPtr< FDProduct > FDProductSP;
typedef vector< FDProductSP > FDProductArray;
typedef refCountPtr< FDProductArray > FDProductArraySP;

/** product base class for FD and tree
    support component instrument (state variable like and allow component valuation **/
class TREE_DLL FDProduct
{
public:
    // product 
    FDProduct(FDModel* model) : model(model), calcOff(false) {}
    virtual ~FDProduct() {}

    /** get model */
    FDModel * getModel() const { return model; }

    /** if this is elementary product - model supplied state variables */
    virtual bool isElementary() const { return false; }

    /** vanilla, quanto or struck */
    virtual string getCcyTreatment() const { return "V"; }

    /** returns all slices to DEV, method is used by the model */
    virtual const vector< TreeSliceSP > & getSlicesToDEV() const;

    /** called before each DEV and update() */
    virtual void preCalc( int step ) {}

    /** switches of the computation of this product for performance reasons */
    virtual void setCalcOff() { calcOff = true; }
    bool isCalcOff() { return calcOff; }

    ////////////////////  required interface methods //////////////////
    /** return product start date. i.e. first fixing date */
    virtual DateTime getStartDate() const {
        throw ModelException("FDProduct::getStartDate()", 
            "Virtual function not implemented by product " 
            + string(typeid(*this).name()));
    }

    // Specifies that a value (getValue) will be needed at modelResetDates[i] and 
    // an approximate values may be needed in range [modelResetDates[i]; lateResetDates[i]]
    virtual void addModelResetDates(
        const DateTimeArray &modelResetDates, 
        const DateTimeArray &lateResetDates
    ) {}
    inline void addModelResetDates(const DateTimeArray &modelResetDates) {
        addModelResetDates(modelResetDates, modelResetDates);
    }
    void addModelResetDate(DateTime modelResetDate, DateTime lateResetDate);
    inline void addModelResetDate(DateTime modelResetDate) {
        addModelResetDate(modelResetDate, modelResetDate);
    }

    /** initialization, called ONCE only for each new model instance 
    this gives product chance to init model params. 
    must not init product variables here, use initProd() instead */
    virtual void init (Control * control) const = 0;

    /** initializing and setting product variables */
    // this is called per pricing call before each pricing 
    virtual void initProd() = 0;

    virtual void printInfo(ostream& outputStream) const;

    /** backward induction or fwd induction, at time point = 0, t or T */
    enum UpdateType {BWD, BWD_T, FWD_0, FWD, BWD_NODE_INSERTION};
    /** update value array by product */
    virtual void update( int & step, UpdateType type ) = 0;

    /** get product value */
    virtual const TreeSlice & getValue( int step ) const;
    virtual const TreeSlice & getValue( int step, DateTime eventDate ) const;

    /** get product cash flow at current time, default empty implementation provided */
    virtual const TreeSlice & getCashFlow( int step ) const;

    /** output prices and any additional results */
    virtual void recordOutput(Control* control, YieldCurveConstSP disc, Results* results) = 0;

    /** function to assist runtime debugging by returning an outputName/instance name if setup */
    virtual string getOutputName() const;

    /** Product can decide if it can be priced forward or backward inductively 
     *  default implementation is backward induction                           */
    enum FDDirection { FWD_INDUCTION, BACK_INDUCTION, FWD_BACK_INDUCTION };
    virtual FDDirection getDirectionPref() const { return BACK_INDUCTION; }

    /** return array of points where we could use state prices to express DEV *
    *  should be for example European Exercise dates (default empty array)    */
    virtual DateTimeArray getStatePriceDates( ) const { return DateTimeArray(); };

    /** helper function to store the result of a slice into a requestResult */
    void recordSliceToOutputName(Control* ctrl, Results* results, FDModel* model,
        bool isMainPrice, const string &outputName, string ccy, const TreeSlice & price );
    ////////////////////  required interface methods //////////////////

protected:
    /** convenience function that rethrows an exception after 
     *  having added the function name to the stack */
    ModelException makeException(exception &e, const string &functionName) const;

    /** Tests if an historical value is available from "inst", put it in "slice" and return true
     *  If not, return false. */
    bool historicalValueAvailable(IProdCreator const &inst, TreeSlice &slice, DateTime eventDate) const;

    /** start doing DEV on a slice,
    returns whether slice was affected **/
    bool startDEV( const TreeSliceSP & slice );

    /** stop doing DEV on a slice,
    returns whether slice was affected **/
    bool stopDEV( const TreeSliceSP & slice );

    FDModel* model;
    bool calcOff;
    mutable vector< TreeSliceSP > slicesToDEV; // list of slices on which model performs DEV
};

DRLIB_END_NAMESPACE
#endif


