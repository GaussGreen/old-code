//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AssetHistory.hpp
//
//   Description : Holds all the samples for an asset/exchange
//
//   Author      : Andrew McCleery
//
//   Date        : August 12, 2004
//
//
//----------------------------------------------------------------------------

#ifndef ASSETHISTORY_HPP
#define ASSETHISTORY_HPP

#include "edginc/config.hpp"
#include "edginc/MarketObject.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/WrapperNameCollector.hpp"
#include "edginc/MarketData.hpp"
#include "edginc/ObservableHistory.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL AssetHistory : public MarketObject,
                     virtual public IObservableHistory {
public:
    static CClassConstSP const TYPE;
    friend class AssetHistoryHelper;
    
    void validatePop2Object();
    
    ~AssetHistory();

    /** AssetHistory are essentially immutable objects so clone just returns this
        - more or less */
    virtual IObject* clone() const;
    
    /** Pull out the holidays from market data */
    void getMarket(const IModel* model, const MarketData* market);

    virtual string getName() const;

    virtual string getMarketObservableName() const;

    virtual string getSource() const;

    /** Initialises this piece of market data - records the source, asset name and
    history type identifying the history*/
    virtual void initialise(MarketData* market);

    // Populates CashFlowArray with samples
    void getSamples(CashFlowArraySP cfls, const DateTime& valueDate) const;

    // retrieve a single sample
    virtual double pastValue(const DateTime&           sampleDate,
                     const ObservationType*    obsType,
                     const FixingType*         fixType,
                     const SamplingConvention* sampleRule) const;

    // retrieve a single observation date - returns false if date omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const;

    // is the given date a holiday for the relevant source
    virtual bool isHoliday(const DateTime& sampleDate) const;

    AssetHistory(string name, int size);

private:
    AssetHistory();
    AssetHistory(const AssetHistory &rhs);
    AssetHistory& operator=(const AssetHistory& rhs);

    static void acceptWrapperNameCollector(AssetHistory* assetHist,
                                           WrapperNameCollector* collector);

    string                  name;
	string				    assetName;
    DateTimeArraySP         dates;
    DoubleArraySP           values;
    ObservationTypeArraySP  obsTypes;
    HolidayWrapper          holidays;
    string                  source;
};

typedef smartConstPtr<AssetHistory> AssetHistoryConstSP;
typedef smartPtr<AssetHistory>      AssetHistorySP;
typedef array<AssetHistorySP, AssetHistory> AssetHistoryArray;
typedef smartPtr<AssetHistoryArray> AssetHistoryArraySP;

// support for wrapper class
typedef MarketWrapper<AssetHistory> AssetHistoryWrapper;

DRLIB_END_NAMESPACE

#endif




