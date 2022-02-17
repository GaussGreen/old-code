
#ifndef EDG_CONTROL_H
#define EDG_CONTROL_H
#include "edginc/Object.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/Sensitivity.hpp"
#include "edginc/FORWARD_DECLARE.hpp"
#include "edginc/Results_forward.hpp"
#include "edginc/Control_forward.hpp"
#include <time.h>
#include <map>

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(IModel)
class CInstrument;
FORWARD_DECLARE(MarketData)
FORWARD_DECLARE(IInstrumentCollection)
FORWARD_DECLARE(RiskQuantityEvaluator)

/** General risk manager control */
class RISKMGR_DLL Control: public CObject{
public:
    friend class ControlHelper;
    static CClassConstSP const TYPE;

    static const string ASSET_THEO; // use theortical value of assets
    static const string ASSET_MTM;  // use mark-to-market value of assets
    static const string ASSET_THEO_MTM; // do both

    virtual ~Control();

    Control();

    /* Takes copy of sens array. sens and outReqs can be null */
    Control(const SensitivityArrayConstSP&   sens,
            const OutputRequestArrayConstSP& outReqs,
            bool                             writeToFile,
            const string&                    fileName);

    /* As above but takes value for skipShifts() return value. This is an
       obscure constructor - generally don't use this one.*/
    Control(const SensitivityArrayConstSP&   sens,
            const OutputRequestArrayConstSP& outReqs,
            bool                             writeToFile,
            const string&                    fileName,
            bool                             skipShifts); // skip shifting

    /* Takes copy of sens array. sens and outReqs can be null */
    Control(const SensitivityArrayConstSP&   sens,
            const OutputRequestArrayConstSP& outReqs,
            bool                             writeToFile,
            const string&                    fileName,
            const string&                    assetPriceSource); 

    /** calculates price and sensitivities */
    void calculate(IModel*      model,
                   CInstrument* instrument, 
                   Results*     results);

    /** calculates prices and sensitivities for multiple instruments */
    void calculateMulti(IModel* model,
                        IInstrumentCollectionSP instrument, 
                        CResultsArraySP results);
    
    /** calculates the supplied sensivity (the instrument must have already
        been priced). Calls to getCurrentSensitivity will return sens for the
        duration of the call */
    void calculateSens(const SensitivitySP&  sens,
                       IModel*               model,
                       CInstrument*          instrument, 
                       Results*              results);

    /** calculates exposure characteristics for multiple instruments */
    void calculateMultiExposures(IModel* model,
        IInstrumentCollectionSP instrument, 
        CResultsArraySP results);

    /** returns reference to SensControlArray */
    SensitivityArrayConstSP getSens();

    /** shuffle's the sensitivity array. Used for testing caching. */
    void shuffleSens(int seed);


    /** returns reference to OutputRequestArray */
    OutputRequestArrayConstSP getOutputRequests() const;

    /** checks whether a specific output has been requested. Returns null
        or the class representing the request if it has been requested.
        Do NOT free the pointer returned. (awkward to make it a const since
        we need to alter the object to flag that the result has been
        calculated) */
    OutputRequest* requestsOutput(const string& output);

    /** Deprecated. checks whether a specific output has been requested. 
        Do NOT free the resulting pointer */
    bool requestsOutput(const string& output, OutputRequest*& request);

    /** Returns flag indicating whether a regression test file should 
        be created. */
    bool getWriteToFile() const;

    /** Returns filename for regression file. In particular, if
        validate is true an exception will be thrown if writeToFile is
        true but the file should not be created (eg because a file
        exists already [in a spreadsheet environment]). Also, if
        validate is true throws an exception if the fileName is empty
        and writeToFile is true */
    const string& getFileName(bool validate) const;

    /** Switches off the write to file flag */
    void switchOffWriteToFile();

    /** Processes any command line options which alter the control. 
        The returned ControlSP is either this (and suitably modified) or is 
        a new instance of a Control. */
    CControlSP applyCommandLineOptions();

    /** object validation after construction */
    void validatePop2Object();

    /** Are we currently doing the initial pricing run or are we doing a 
        tweak */
    bool isPricing() const;

    /** Returns the first instance of the specified sensitivity in the
        control. The comparision is an exact comparison of types
        (rather than using isAssignableFrom) */
    SensitivitySP sensitivityRequested(const CClassConstSP& clazz) const;
    //// Same as above but with clazz = sensToDo->getClass()
    SensitivitySP sensitivityRequested(SensitivitySP sensToDo) const;

    void handleUnfulfilledRequests(IModel* model,
                                   IInstrumentCollectionSP instrument,
                                   CResultsArraySP results);

    /** Override clone method to copy non registered fields over */
    virtual IObject* clone() const;

    /** get the sensitivity which is currently being calculated. Returns
        a null SP if no sensitivity is currently being calculated (eg initial
        pricing pass) */
    SensitivitySP getCurrentSensitivity() const;

    /** populate from market cache */
    void getMarket(IModelConstSP model, MarketDataConstSP market,
                   IInstrumentCollectionConstSP instruments);

    /** get delta shift size of a control.
        returns 0 if delta tweak not requested */
    double getDeltaShiftSize() const;
      
    /** get shift size of a scalar shift */
    double scalarShiftSize(CClassConstSP shiftType) const;

    /** add a new sensitivity */
    void addSensitivity(SensitivitySP newSens);

    /** remove sensitivity */
    void removeSensitivity(CClassConstSP sensType);
        
    // add a new request */
    void addRequest(OutputRequestSP request);

    //// remove a request */
    void  removeRequest(const string& requestName);

    /** build up a control from single character flags e.g. D = Delta */
    static Control* makeFromFlags(const string& flags, double impVolTarget);

    /** Returns whether outputs should be scaled or not */
    bool scaleOutputs() const;

    /** Sets whether outputs should be scaled or not */
    void setScaleOutputs(bool scale);

    /** Sets whether assets should use theoretical price or mtm */
    void setUseTheoAssetPrice(bool useTheoAssetPrice);

    /** Indicates whether assets should use theoretical price or mtm for
        the current pricing */
    bool getUseTheoAssetPrice() const;

    /** Returns the source we should use of asset prices as requested by the
        client. Note the string has multiple values including 'both'. Its
        value should not be used to indicate the rule for the current
        pricing. Instead use getUseTheoAssetPrice */
    const string& getAssetPriceSource() const;

    /** Resets the control for a new pricing - clears current sensitivity and
        marks all output requests as not calculated */
    void reset();
    /** Requests that the specified packet is removed from the results at the
        end of the calculate() call */
    void removePacketAfterCalc(const string& packetName);
    
    /** Records the price time in the results (if requested) assuming
        that the pricing call has now finished. Useful for models that
        calculate their own greeks in the initial pricing
        call. Returns the stored value of priceTime (in seconds) */
    double recordPriceTime(CResultsArraySP results);

    //// work out the compute index for Pyramid - a series of lookups based on
    //// how large the compute estimate is and how many (equity) underlings
    int computeIndex(IModel* model, 
                     IInstrumentCollectionSP insts, 
                     double priceTime,
                     Results* results);

    void startTiming();
    void endTiming(IModel* model, IInstrumentCollectionSP inst,
                   CResultsArraySP results);

    /** If skipShifts() = true then when a shift to the market is requested
        it can be skipped. This is for a performance gain when assessing
        the total number of calculations an instrument requires */
    bool skipShifts() const;

private:
    class ComputeEstimator;
    friend class ComputeEstimator;
    Control(const Control &rhs);
    Control& operator=(const Control& rhs);

    SensitivityArraySP   sens;
    RiskQuantityEvaluatorSP riskQuantityEvaluator;
    OutputRequestArraySP outputRequests;
    string               fileName;
    string               assetPriceSource; // see values above
    // derived map used for faster access to output requests
    map<string, OutputRequest*>  outputRequestSet; // $unregistered
    SensitivitySP                currentSens;      // transient
    int currentInstrumentIndexInCollection; // transient: index into tempPackets
    StringArrayArray             tempPackets; /* transient: names of packets
                                                 to remove after pricing */
    clock_t                      startTime;   // transient: set in calculate $unregistered
#ifdef UNIX
    time_t                       startSec;   // transient: set in calculate $unregistered
#endif
    /* true: the instrument contains a Derivative Asset, and therefore
       requires the control object to be passed through before pricing
       and every tweak */
    bool                 useTheoAssetPrice;    // transient
    bool                 isPricingRun;     // transient
    bool                 scaleResults;
    bool                 writeToFile;
    /** If skipAllShifts = true then when a shift to the market is requested
        it can be skipped. This is for a performance gain when assessing
        the total number of calculations an instrument requires */
    bool                 skipAllShifts; // transient, default false

    friend class HypothesisTree; // FIXME temporary---so that we can set currentSens
};

DRLIB_END_NAMESPACE
#endif
