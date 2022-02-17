//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Sensitivity.hpp
//
//   Description : Base class for sensitivities
//
//   Author      : Mark A Robson
//
//   Date        : 11 Jun 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_SENSITIVITY_HPP
#define EDR_SENSITIVITY_HPP

#include "edginc/FORWARD_DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE
class CInstrument; // to avoid circular header dependencies
class IModel;
class Control;
class Results;
class MarketData;
class TweakGroup;
FORWARD_DECLARE(OutputName);
// fix for MSVC6
//#ifndef QLIB_OUTPUTNAME_CPP
//EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<OutputNameArray>);
//EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<OutputNameArray>);
//#endif

/** An abstract base class. It identifies a particular type of tweak
    to make. It is also responsible for driving how that tweak is
    calculated as well as how the result is returned.  */
class RISKMGR_DLL Sensitivity: public CObject{
public:
    static CClassConstSP const TYPE;

    virtual ~Sensitivity();

    /** identifies the name used for storing associated results in the output*/
    virtual const string& getSensOutputName() const = 0;

    /** Is this sensitivity made using a discrete shift (ie a jump) or a
        an approximately continuous one (return value: false) */
    virtual bool discreteShift() const = 0;

    /** identifies the packet in which the results are stored. By default this
        returns the same as getSensOutputName() */
    virtual const string& getPacketName() const;

    /** scale the results in the Results Object for this sensitivity
        by supplied factor. Default implementation is (provided instance
        of Sensitivity implements the Additive interface) to scale all
        the results in the packet provided by getPacketName() if this is
        the same as getSensOutputName() otherwise the result with name
        getSensOutputName() in packet getPacketName() is scaled */
    virtual void scaleResult(Results*     results,     // (M)
                             double       scaleFactor) const;

    /** Modify results in the Results Object for this sensitivity by
        adding all results in resultsToAdd as indicated by control */
    virtual void addResult(Results*           results,     // (M)
                           const Results*     resultsToAdd,
                           double             scaleFactor) const;

    /** Extracts 'price' for this sensitivity  from results set. Default
        uses results->retrievePrice() */
    virtual double getSensPrice(Results*     results,
                                CInstrument* inst,
                                IModel*      model,
                                Control*     control);
    
    /** calculates given sensitivity. Default implementation checks 
     validitity of pointers and calls calculate */
    virtual void calculateSens(IModel*          algorithm,
                               CInstrument*     instrument,
                               Control*         control,
                               Results*         results);

    /** returns a reference to the algorithm used for pricing. Only valid
        when calculating a sensitivity */
    IModel* getModel() const;

    /** returns a reference to the control. Only valid
        when calculating a sensitivity  */
    Control* getControl() const;

    /** populate from market cache - default provided */
    virtual void getMarket(const IModel* model, const MarketData* market);

    /** does this sensitivity have its own list of underlyings ? */
    virtual bool hasOverrideNames() const;

    /** return its own list of underlyings to tweak */
    virtual OutputNameArrayConstSP overrideNames() const;

    /** store list of underlyings to tweak */
    virtual void storeOverrideNames(OutputNameArraySP names);

    /** specialised copy that allows over-ride of algorithm */
    virtual Sensitivity * spawn(IModel* model) const;

    /** Removes any names in the array returned by overrideNames() that have
        already been calculated for the specified packet. Note that this
        object is modified accordingly */
    virtual void removeOverrideNames(const string& packetName,
                                     const Results* results);

    virtual void setAlgorithm(IModel*);
    virtual void setControl(Control*);

protected:
    /** Note Sensitivity is abstract */
    Sensitivity(const CClassConstSP& clazz);

    /** calculates given sensitivity - invoked by calculateSens */
    virtual void calculate(TweakGroup*      tweakGroup,
                           Results*         results) = 0;

    // protected fields
    /** This is set by calculateSens at the start of a calculation for
       a given senstivitiy. It is then cleared */
    IModel*      algorithm;  // $unregistered
    Control*     control;    // $unregistered

    /** to allow only certain underlyings to be tweaked */
    OutputNameArraySP toTweak;

private:
    friend class SensitivityHelper;
    friend class ScalarShiftTwoSided; // FIXME why do I need this?!?!??!
    Sensitivity();
    Sensitivity(const Sensitivity &rhs);
    Sensitivity& operator=(const Sensitivity& rhs);
};

DECLARE(Sensitivity)
#ifndef QLIB_SENSITIVITY_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<Sensitivity>);
EXTERN_TEMPLATE(class RISKMGR_DLL array<SensitivitySP _COMMA_ Sensitivity>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<SensitivityArray>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<SensitivityArray>);
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<Sensitivity>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL array<SensitivitySP _COMMA_ Sensitivity>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<SensitivityArray>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<SensitivityArray>);
#endif

DRLIB_END_NAMESPACE

#endif

