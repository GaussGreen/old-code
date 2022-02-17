
#ifndef EDG_RESULTS_H
#define EDG_RESULTS_H

#include "edginc/Results_forward.hpp"
#include "edginc/OutputName.hpp"
#include "edginc/OutputRequest.hpp"
#include "edginc/Control.hpp"
#include "edginc/NotApplicable.hpp"
#include "edginc/Atomic.hpp" // not strictly needed here

DRLIB_BEGIN_NAMESPACE

class SensControl;
class ResultsHelper;

/** Stores price + sensitivities etc. Results are stored within
    packets with each packet typically holding the results for a
    different sensitivity. As such results are identified by which
    packet they belong to and their OutputName which typically holds
    the name(s) of what was tweaked to get this result. Any oject
    derived from IObject can be stored as a result. However, unless
    this object is stored within a DecoratedResult it is assumed that
    the result is in the same units as the price. For example,
    indicative vol should be stored using storeCcylessResult()
    method. It is also possible to store results in another ccy from
    the price using the storeCcyResult() method */
class RISKMGR_DLL Results: public CObject{
public:
    friend class ResultsHelper;
    static CClassConstSP const TYPE;
    static const string INSTRUMENT_PACKET;
    static const string FWD_AT_MAT_PACKET;
    static const string DEBUG_PACKET;
    static const string DEBUG_PACKETS_PREFIX;
    static const string VALUE;
    static const string ESW_LEG_PRICE;
    static const string KNOWN_CASHFLOWS_PACKET;
    static const string CLEAN_DEFAULT_SPREAD_CURVE_PACKET;
    static const string SPI_PACKET;
    static const string TRANCHE_CONTINGENT_LEG_PACKET;
    static const string TRANCHE_FEE_LEG_PACKET;
    //// used for storing shift sizes when calculating two sided
    //// derivatives - appended to sensOutputName
    static const string SHIFT_SIZE_POSTFIX;

    static OutputNameSP  emptyName;        // empty output name to avoid repeated
                                           // memory allocation

    /** Creates an empty results set */
    Results();

    /** Returns the name of currency in which the value is in */
    const string& getCcyName() const;

    /** Stores the fair value of the instrument (aka price) */
    void storePrice(double         price,
                    const string&  ccyName);

    /** Deprecated - do not use. To be retired. Use storeGreek or
        storeRequestResult */
    void storeCcylessResult(IObjectSP                result,
                            const string&            packetName,
                            const OutputNameConstSP& outputName);

    /** Stores a generic greek without any decoration */ 
    void storeGreek(IObjectSP                result,
                    const string&            packetName,
                    const OutputNameConstSP& outputName);

    /** Removes a generic greek if it exists. The packet is also removed
        if it becomes empty */ 
    void removeGreek(const string&            packetName,
                     const OutputNameConstSP& outputName);

    /** Stores a 'scalar' greek. Equivalent to  
        storeGreek(CDouble::create(result), sens) */
    void storeScalarGreek(double         result,
                          SensControl*   sens);  // (I) identifies greek

    /** Stores a 'scalar' greek. Equivalent to  
        storeGreek(resultObj,
                   packetName, 
                   ouputName)
    where resultObj is CDouble containing result*/ 
    void storeScalarGreek(double                   result,
                          const string&            packetName,
                          const OutputNameConstSP& outputName);

    /** Stores a generic greek. */
    void storeGreek(IObjectSP      result,
                    SensControl*   sens);  // (I) identifies greek

    /** If a greek doesn't make sense for a product (e.g. delta for
        cashflows) store a NotApplicable object as the result    */
    void storeNotApplicable(const string&       packetName,
                            const OutputNameSP& outputName,
                            const IObjectSP&    na);
    /** If a greek doesn't make sense for a product (e.g. delta for
        cashflows) store a NotApplicable object as the result    */
    void storeNotApplicable(const Sensitivity* sens); // (I) identifies greek
    //// As above but can specify cause
    void storeNotApplicable(const Sensitivity*     sens,
                            NotApplicable::Cause   cause);

    /** If an output request doesn't make sense for a product 
        (e.g. effective strike for a vanilla) store a NotApplicable 
        object as the result */
    void storeNotApplicable(const OutputRequest*    request);

    //// As above but can specify cause
    void storeNotApplicable(const OutputRequest*    request,
                            NotApplicable::Cause    cause);

    /** Store not applicable for the given packet name */
    void storeNotApplicable(const string&  packetName);

    /** Retrieves the instrument price */
    double retrievePrice() const;

    /** Whether a price has been set using storePrice() */
    bool priceExists() const;

    /** Returns NotApplicable object if a result is stored for the
        specified request and it is of type NotApplicable else returns
        null */
    const NotApplicable* isNotApplicable(const OutputRequest* request) const;

    // has a greek been stored as NotApplicable?
    bool isNotApplicable(SensControl* sens) const;

    /** checks whether results hold a valid results for a given greek */
    bool isValidScalarGreek(
        const string&            packetName,
        const OutputNameConstSP& outputName) const; // (I) identifies greek

    bool   isValidScalarGreek(
        SensControl*   sens) const;  // (I) identifies greek

    /** Retrieves all scalar greeks (incl price) plus packet names, output names.
        Returns num of results retrieved */
    int retrieveAllScalarGreeks(StringArray&          packetNames,
                                 OutputNameArray&    outNames,
                                 DoubleArray&          results) const;

    /** Retrieves a 'scalar' greek */
    double retrieveScalarGreek(
        SensControl*   sens) const;  // (I) identifies greek

    /** Retrieves a 'scalar' greek. */
    double retrieveScalarGreek(
        const string&            packetName,
        const OutputNameConstSP& outputName) const;

    /** Retrieves a generic greek. Note returns reference to result */
    IObjectConstSP retrieveGreek(
        const string&            packetName,
        const OutputNameConstSP& outputName) const;

    /** write object out in XML format */
    virtual void write(const string& tag, Writer* writer) const;

    /** populate an empty object from XML description */
    virtual void import(Reader::Node* elem, Reader* reader);

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** Returns true if given result exists. Equivalent to 
        exists(sens->getSensOutputName(), sens->getMarketDataName()) */
    virtual bool exists(SensControl*   sens) const; 

    /** Returns true if given result exists.  */
    virtual bool exists(const string&            packetName,
                        const OutputNameConstSP& outputName) const;

    /** Returns true if a result exists for the specified output
        request */
    bool exists(const OutputRequest* outputRequest) const;

    
    /** Returns true if given packet exists.  */
    virtual bool packetExists(const string&  packetName) const;
    
    /** Removes the specified packet if it exists */
    void removePacket(const string&       packetName);

    virtual ~Results();

    /** Returns the list of packets within the Results */
    StringArraySP listAllPackets() const;

    /** Returns a list of output names within a packet */
    OutputNameArraySP packetContents(const string& packetName) const;

    /** Returns the list of packets within the Results for use by
        Pyramid (hides certain packets eg DEBUG) */
    vector<const string*> listPyramidPackets() const;

    /** Returns an array of all the results in the given packet */
    vector<pair<OutputNameConstSP, IObjectConstSP> > listPacketResults(
        const string&  packetName) const;

    void storeRequestResult(OutputRequest*           request,
                            const IObjectSP&         resultObj,
                            const OutputNameConstSP& outputName);

    // lazy version of the above
    void storeRequestResult(OutputRequest* request,
                            double         result,
                            const string&  outputName);
    
    /** stores double together with ccy ISO code */
    void storeRequestResult(OutputRequest*           request,
                            double                   result,
                            const string&            ccyISOCode,
                            const OutputNameConstSP& outputName);

#if 0
    /** stores object together with ccy ISO code */
    void storeRequestResult(OutputRequest*           request,
                            const IObjectSP&         result,
                            const string&            ccyISOCode,
                            const OutputNameConstSP& outputName);
#endif

    /** for requests that don't need further qualification e.g. IND_VOL */
    void storeRequestResult(OutputRequest*   request,
                            const IObjectSP& resultObj);

    /** for requests that don't need further qualification e.g. IND_VOL */
    void storeRequestResult(OutputRequest* request,
                            double         result);

    /** Retrieves a generic output request result. Note returns reference to result */
    IObjectConstSP retrieveRequestResult(
                const string& outputRequestName) const;
    
    /** scale all the Additive results (as determined by the control)
        inside a Results object by supplied factor. The singleInstStatistics
        indicates whether we're trying to caclulate the mean value of
        a result by pricing the same instrument multiple times */
    void scale(const CControlSP&  control,
               double             scaleFactor,
               bool               singleInstStatistics);

    /** scale all the CombinableResult results in the given packet
        by supplied factor */
    void scale(const string&       packetName,
               double              scaleFactor);

    /** scale the result with given sensName in the given packet
        by the supplied factor */
    void scale(const string&       packetName,
               const string&       sensName,
               double              scaleFactor);

    /** scale the outputs of a result set by the supplied factor -
        this is mainly used to convert the output conventions of
        convertible bonds and is driven by a different flag than the
        scale methods above, thus a different method */
    void scalePostProcess(const CControlSP&   control,
                          double              scaleFactor);

    /** Modify 'this' Results set by adding all results in resultsToAdd
        as indicated by control */
    void add(const Results*     resultsToAdd,
             const CControlSP&  control,  // used to obtain resultsToAdd
             double             scaleFactor,
             bool               sameInstrument); /* are the 2 results for
                                                    the same instrument */

    /** Modify 'this' Results set by adding all CombinableResult
        results in the packet of given name to this Results */
    void add(const string&      packetName,
             const Results*     resultsToAdd,
             double             scaleFactor);

    /** Modify 'this' Results set by adding the result with given
        sensName in the given packet by the supplied factor */
    void add(const string&      packetName,
             const string&      sensName,
             const Results*     resultsToAdd,
             double             scaleFactor);

    /** Modify 'this' Results set by merging all objects in the packet of
        given name with those in this Results object */
    void merge(const string&      packetName,
               const Results*     resultsToAdd);

    /** Merge all packets from resultsToAdd into this results object */
    void merge(const Results* resultsToAdd);

    /** Modify this Results by overwriting all results contained in
        resultsToCopyFrom */
    void overwrite(const Results* resultsToCopyFrom);

    //// Combines results from pricing an array of instruments
    static ResultsSP combineResults(
        const CControlArray&   ctrls,
        const DoubleArray&     weights,
        CResultsArray&         results);

private:
    ResultsHelper *data;
    void storeNotApplicable(const OutputRequest*    request, 
                            const IObjectSP&        na);
    void storeResult(const string&            packet,
                     const OutputNameConstSP& resultName,
                     const IObjectSP&         result);
    IObjectSP retrieveResult(
        const string&            packet,
        const OutputNameConstSP& resultName) const;
    //// helper method. Do we have results for when DerivativeAsset's use mtm
    bool mtmResultsExist() const;
    //// helper method. Get results for when DerivativeAsset's use mtm
    Results* getMtmResults();
    //// helper method. Get results for when DerivativeAsset's use mtm
    const Results* getMtmResults() const;

    void addSingleResultsSet(
        const Results*     resultsToAdd,
        const CControlSP&  control,  // used to obtain resultsToAdd
        double             scaleFactor,
        bool               sameInstrument); /* are the 2 results for
                                               the same instrument */
    

};


DRLIB_END_NAMESPACE

#endif
