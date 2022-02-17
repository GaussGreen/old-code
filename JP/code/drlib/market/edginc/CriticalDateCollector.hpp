//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CriticalDateCollector.hpp
//
//   Description : Critical date collector class
//
//   Author      : Andre Segger
//
//   Date        : 10 Jul 2001
//
//
//----------------------------------------------------------------------------

#ifndef CRITDATE_COLLECT_HPP
#define CRITDATE_COLLECT_HPP
#include "edginc/Collector.hpp"
#include "edginc/DividendList.hpp"
#include "edginc/Asset.hpp"

DRLIB_BEGIN_NAMESPACE

/** A class to collect important dates or events */
class MARKET_DLL CriticalDateCollector: public CObject, public virtual ICollector {
public:
    friend class CriticalDateCollHelper;
    static CClassConstSP const TYPE;

    CriticalDateCollector(const DateTime& startDate,
                          const DateTime& endDate,
                          bool            convertToDollars);

    void addDividend(const Dividend& newDiv);

    /** add specified dates to the collector. The type parameter is the type
        of the object supplying the dates. This could allow, for example,
        differentiation between equity and IR vols. Here benchmark dates means
        option expiry dates (eg for IR vols) */
    void addDates(const DateTimeArray& dates, CClassConstSP type);

    static DividendListSP collectAllDates(CAssetSP        asset,
                                          const DateTime& startDate,
                                          const DateTime& endDate,
                                          bool            convertToDollars);

    /** Collects vol 'critical' dates using accept method on
        supplied object. Only dates from objects of type clazz and
        that fall within startDate < date < endDate are included */
    static DateTimeArray collectVolDates(IObjectConstSP  object,
                                         CClassConstSP   clazz,
                                         const DateTime& startDate,
                                         const DateTime& endDate);


    // access methods
    DateTime getStartDate() const;
    DateTime getEndDate() const;
    bool requireDollarDivs() const;

private:
    static void load(CClassSP& clazz);
    CriticalDateCollector(const CriticalDateCollector &rhs);
    CriticalDateCollector& operator=(const CriticalDateCollector& rhs);

    DateTime            startDate; // $unregistered
    DateTime            endDate; // $unregistered
    bool                convertToDollars; // $unregistered
    DividendArray       dollarDivArray; // $unregistered
    DateTimeArray       volDates; // $unregistered
    CClassConstSP       requiredType; // $unregistered
};

typedef smartPtr<CriticalDateCollector>       CriticalDateCollectorSP;
typedef smartPtr<const CriticalDateCollector> CriticalDateCollectorConstSP;

DRLIB_END_NAMESPACE

#endif
