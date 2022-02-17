//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MarketData.hpp
//
//   Description : Class controlling access to market data
//
//   Author      : Mark A Robson
//
//   Date        : 16 Jan 2001
//
//
//----------------------------------------------------------------------------

#ifndef MARKETDATA_HPP
#define MARKETDATA_HPP
#include "edginc/DateTime.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/WriteableMap.hpp"
#include "edginc/ReadableMap.hpp"
#include "edginc/Model.hpp"
#include "edginc/MarketData_forward.hpp"

DRLIB_BEGIN_NAMESPACE

class IModel;

/** Class controlling access to market data */
class RISKMGR_DLL MarketData: public CObject,
                  virtual public IReadableMap,
                  virtual public IWriteableMap {
public:
    static CClassConstSP const TYPE;
    friend class MarketDataHelper;
    friend class GetMarketDataAddin;
    friend class MarketDataProxy;

    ~MarketData();

    MarketData();

    MarketData(const DateTime& today);

    /** creates deep copy of MarketData */
    virtual IObject* clone() const;

    /** write object out to writer */
    virtual void write(const string& tag, Writer* writer) const;

    /** populate an empty object from XML description */
    virtual void import(Reader::Node* elem, Reader* reader);
    
    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    void outputWrite(const string& linePrefix,
                     const string& prefix, ostream& stream) const;

    /** Stores (a reference) to the supplied data in this market data
     cache. The getName method() on data is used for the object's id.
     Should consider making data a MarketObjectConstSP */
    void AddData(const MarketObjectSP& data); 

    /** Asks the supplied model to return a [deep] copy of the
        requested MarketObject object from this cache. The getMarket method
        is then invoked on the returned MarketObject */
    MarketObjectSP GetData(const IModel*         model,
                           const string&         name,
                           const CClassConstSP&  type) const;

    /** Returns a [deep] copy of the requested MarketObject object from
        this cache. If a unique object with the given name and of the
        given type (or derived from the given type) then that object is
        returned. Otherwise an exception is thrown. 
        'model' here is optional : if provided it is used to specialise 
        the data extraction. */
    MarketObjectSP GetData(const string&         name,
                           const CClassConstSP&  type) const;

    /** As above but can use model to specialise 
     */
    MarketObjectSP GetObjectData(const string&         name,
                                 const CClassConstSP&  type,
                                 const IModel*         model) const;

    /** Returns a [deep] copy of all objects in the cache with the 
        the supplied name. */
    MarketObjectArraySP GetAllDataWithName(const string& name)const;

    /** Return the reference date (aka value date, today) associated
        with this data */
    DateTime GetReferenceDate() const;

    /** Populates the supplied dateTime with the value date from the
        market data cache. Fails if no date is in the cache and the
        supplied date is empty or if the supplied date is not empty
        and does not match the one in the cache */
    void GetReferenceDate(DateTime& valueDate) const;

    /** Returns the [concrete] type of the market data requested - fails
        if the data does not exist */
    CClassConstSP getDataType(const string&         name,
                              const CClassConstSP&  type) const;

    /** Returns true if there is a item of market data with the supplied
        name and derived from the supplied type */
    bool hasData(const string&         name,
                 const CClassConstSP&  type) const;

    /** Merge supplied market with this one. Data in market overwrites data in
        this if conflicts exist */
    void addData(const MarketData& market);

    // may want to combine these get/set name functions into two functions 

    /** Returns the name identifying a correlation given the names of the two
        objects whose correlation is required */
    string getCorrelationName(const string& id1, 
                              const string& id2) const;

    bool   hasCorrelationData(const string& id1,
                              const string& id2) const;

    /** Sets the name identifying a correlation given the names of the two
        objects identifying the correlation */
    void setCorrelationName(const string& id1,
                            const string& id2,
                            const string& correlationName);

    /** Returns the name identifying a correlation given the names of the two
        objects whose correlation is required */
    string getCorrelationTermName(  const string& id1, 
                                    const string& id2) const;

    bool   hasCorrelationTermData(  const string& id1,
                                    const string& id2) const;

    /** Sets the name identifying a correlation given the names of the two
        objects identifying the correlation */
    void setCorrelationTermName(const string& id1,
                                const string& id2,
                                const string& correlationTermName);

	/** Returns the name identifying a corrSwapBasisAdj given the names of the two
		regions for which the sampling adj is required */
	string getCorrSwapSamplingAdjName(const string& id1, 
								      const string& id2) const;

	bool hasCorrSwapSamplingAdjData(const string& id1,
								    const string& id2) const;

	/** Sets the name identifying a corrSwapSamplingAdj given the names of the two
		regions for which the sampling adj is required  */
	void setCorrSwapSamplingAdjName(const string& id1,
								    const string& id2,
									const string& corrSwapBasisAdjName);


    /** get all equity/fx correlations for a given equity */
    MarketObjectArraySP getEquityFXCorrelations(
                                         const string& equityName,
                                         const CClassConstSP&  fxAssetType,
                                         const CClassConstSP&  correlationType);

    /** Returns the name identifying an fx asset given the names of the two
        yield curves for the two currencies making up the fx rate */
    string getFXName(const string& riskCcy, 
                     const string& baseCcy) const;

    bool   hasFXData(const string& riskCcy,
                     const string& baseCcy) const;

    /** Sets the name identifying a correlation given the names of the two
        objects identifying the correlation */
    void setFXName(const string& riskCcy,
                   const string& baseCcy,
                   const string& fxName);

    /** Record what the iso code is for a specific yield curve */
    void setYieldCurveISOCode(const string& ycName,
                              const string& ycISOCode);

    /** Get what the iso code is for a specific yield curve. This is useful
        if you want to choose what type of data to get inside a yield curve
        dependent upon the iso code */
    const string& getYieldCurveISOCode(const string& ycName) const;

    /** Returns all ObservableHistorys in the cache for a given
        market observable name and type*/
    MarketObjectArraySP getAllObservableHistories(const IModel* model,
                                                  const string& obsName,
                                                  const CClassConstSP&  type) const;

    /** Returns the name identifying the ObservableHistory in the cache for a given
        Observable name, type and source*/
    string getObservableHistoryName(const string& obsName,
                                   const string& source,
                                   const CClassConstSP& type) const;

    /** Checks whether we already have an ObservableHistory for the given
        Observable name, type and source - note this also throws an Exception
        if we try to add a history of a different type*/
    bool hasObservableHistoryData(const string& obsName,
                                  const string& source,
                                  const CClassConstSP& type) const;

    /** Sets the name identifying a correlation given the Observable name */
    void setObservableHistoryName(const string& obsName,
                                  const string& source,
                                  const CClassConstSP& type,
                                  const string& cacheName);

    /** Fetches any market data required by eg. Business252 */
    void fetchForNonMarketObjects(IObjectSP     obj,
                                  IModelConstSP model,
                                  const string& className) const;

private:

    IObjectConstSP _plugin(const type_info& which,
                           IObjectConstSP (*factory)(const MarketData&)) const;

public:

    /**
     * Return (creating if not already present) a domain-specific "extension
     * module"
     *
     * This is a way of storing information related to RiskMapping, potentially
     * also correlations, asset history etc., without having to hardwire
     * anything into the general-purpose MarketData class.
     *
     * For example, RiskMapping wants to index the available
     * RiskMappingMatrix's so it doesn't have to search through them all at the
     * start of each pricing --- obviously the index needs to be stored "in"
     * with the MarketData object so that it persists across RiskMgr calls, but
     * we want to keep the code privately in RiskMapping.cpp because it really
     * has nothing to do with market data in general.  So here's how it works:
     *
     *    -  RiskMapping.cpp defines a class RiskMappingMatrixCache, which
     *       has a static "constructor" method from MarketData, i.e.
     *       RiskMappingMatrixCache::SP(const MarketData&).  The constructor
     *       pulls out the RiskMappingMatrix's out of the market and indexes
     *       them.  A method RiskMappingMatrixCache::rmmsForParameters(...)
     *       is provided for retrieval.
     *
     *    -  During the getMarket() phase of a pricing, RiskMapping::fromMarket()
     *       calls
     *       <PRE>
     *       MarketData::plugin<RiskMappingMatrixCache>()->rmmsForParameters(...)
     *       </PRE>
     *       to retrieve the RiskMappingMatrix's it requires.  The first time
     *       it does that, a RiskMappingMatrixCache is automatically created
     *       an initialised from the MarketData --- i.e. the index is built.
     *       Subsequently the same RiskMappingMatrixCache is returned.
     *
     *    -  Whenever something new is added to the MarketData, all the cached
     *       plugins are discarded, so the next time they're requested, they
     *       are created afresh.
     *
     * The type PLUGIN must have the following signature:
     *
     * <PRE>
     *    // Must be a CObject directly or indirectly
     *
     *    class PLUGIN: public CObject {
     *    public:
     *
     *        // Usual type descriptor
     *        static const CClassConstSP TYPE;
     *
     *        // Static "constructor" or "factory" method
     *        static IObjectConstSP SP(const MarketData&);
     *    };
     * </PRE>
     *
     * Obviously the idea is that you add whatever domain-specific methods you
     * need for your application.
     *
     * <H3>Discussion</H3>
     *
     * This is really a first draft.  I've only briefly looked at how easy it
     * would be to migrate the domain-specific methods above (for correlations
     * etc.) to this scheme --- the main choice point in the interface design
     * is, do you go for a "throw away and build from scratch" model, as here,
     * or do you instead "build up incrementally", as the asset history stuff
     * seems to.  For me the latter is going to be a bit complicated as long
     * as it's possible to _overwrite_ market data as well as just _add_ it.
     */

    template <class PLUGIN>
    smartConstPtr<PLUGIN> plugin() const {
        return smartConstPtr<PLUGIN>::dynamicCast(_plugin(typeid(PLUGIN),
                                                          &PLUGIN::SP));
    }

    /** Returns a [shallow] copy of all objects in the cache with the 
        the supplied type. */
    MarketObjectArraySP GetAllDataWithType(const CClassConstSP& type) const;

    // Map interfaces
    // Builds an iterator
    virtual IMap::IIteratorSP createIterator();
    /** Is this object truly a map ie does toObject()/toMap() return this */
    virtual bool isTrueMap() const;    

    /** Add data to the market */
    virtual void put(const string& key, const IObjectSP& value);

    /** Get data from the market - returns an array */
    virtual IObjectSP get(const string& key) const;
private:
    friend class MarketIterator;
    // hide implementation in separate class
    class Imp;
    auto_ptr<Imp> my; // $unregistered

    MarketData(const MarketData& rhs);
    MarketData& operator=(const MarketData& rhs);
    //// methods ///
    static void load(CClassSP& clazz);
    static IObject* defaultMarketData();
    MarketObjectConstSP getData(const string&         name,
                                const CClassConstSP&  type,
                                const IModel*         model = 0) const;
    
    /// static fields ////
    static const string REF_DATE;
    static const string DATA_CACHE;
};

DRLIB_END_NAMESPACE
#endif
