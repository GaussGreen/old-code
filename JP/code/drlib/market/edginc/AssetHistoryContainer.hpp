//----------------------------------------------------------------------------
//
//   Filename    : AssetHistoryContainer.hpp
//
//   Description : Holds all the asset histories for a given asset
//                 see AssetHistory QLib Design Proposal - Ian Stares Jan 2006
//
//   Author      : Ian Stares
//
//   Date        : January 19, 2006
//
//
//----------------------------------------------------------------------------

#ifndef ASSETHISTORYCONTAINER_HPP
#define ASSETHISTORYCONTAINER_HPP

#include "edginc/config.hpp"
#include "edginc/AssetHistory.hpp"
#include "edginc/ObservableHistoryContainer.hpp"

DRLIB_BEGIN_NAMESPACE

class MARKET_DLL AssetHistoryContainer : public CObject,
                              public virtual IObservableHistoryContainer {
public:
    static CClassConstSP const TYPE;

    // constructor which copies MarketObjects coming direct from cache
    AssetHistoryContainer(MarketObjectArraySP& histObjs);

    /** AssetHistoryContainer need to be immutable so clone just returns this*/
    virtual IObject* clone() const;

    /** returns the AssetHistory object corresponding to the source */
    AssetHistoryConstSP getAssetHistory(const string &source) const;
    
    // retrieve a single sample
    virtual double pastValue(const DateTime&            sampleDate,
                             const ObservationType*     obsType,
                             const string&              source,
                             const FixingType*          fixType,
                             const SamplingConvention*  sampleRule) const;

    // retrieve a single observation date
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                 const string&             source,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const;

    // is the given date a holiday for the relevant source
    virtual bool isHoliday(const DateTime& sampleDate,
                           const string&   source) const;
private:
    AssetHistoryContainer();
    AssetHistoryContainer(const AssetHistoryContainer &rhs);
    AssetHistoryContainer& operator=(const AssetHistoryContainer& rhs);

    const string sourceList() const;

    int findHistoryForSource(const string& source) const;

    static void load(CClassSP& clazz);
    static IObject* defaultAssetHistoryContainer();

    AssetHistoryArraySP       histories;
};

typedef smartPtr<AssetHistoryContainer> AssetHistoryContainerSP;
typedef smartConstPtr<AssetHistoryContainer> AssetHistoryContainerConstSP;

DRLIB_END_NAMESPACE

#endif




