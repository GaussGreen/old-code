//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : IndexSpec.hpp
//
//   Description : specification for underlying/payoff index
//                 spot - as equity or fx
//                 fwdspot/futres - maturity date
//                 avgspot index - avg date list
//                 extreme index - max/min, sampling dates, intra-day or at close
//                 ir swap rate - tenor(5y), freq(semi), type (swap, cms etc), or reset/pay dates
//                 float rate - freq(3m), type (libor, libor in arrea etc), or reset/pay dates
//                 ir futures - freq(3m), type, maturity 
//                 credit index - tenor, 
//
//----------------------------------------------------------------------------

#ifndef INDEX_SPEC_HPP
#define INDEX_SPEC_HPP

#include "edginc/DateTime.hpp"
#include "edginc/Model.hpp"
#include "edginc/MarketFactor.hpp"
#include "edginc/AssetHistory.hpp"
#include "edginc/MarketObservable.hpp"
#include "edginc/Events.hpp"

DRLIB_BEGIN_NAMESPACE

/** Base interface class for IndexSpec */
class MARKET_DLL IIndexSpec : public virtual IProdCreator,
                              public virtual IMarketObservable,
                              public virtual FixingReqEvent::IEventHandler
{
public:
    /************************ methods ************************/
    /** returns index name */
    virtual string getName() const = 0;

    /** returns index factor */
    virtual IMarketFactorConstSP getFactor() const = 0;
    virtual IMarketFactorSP getFactor() = 0;

    /* IProdCreator:: */
    virtual double getValue(DateTime date, CashflowInfo &cfi ) const;
    
	/* IMarketObservable:: */
    virtual double pastValue(
	    const DateTime&             sampleDate,
	    const ObservationType*      obsType,   /* &ObservationExact() */
	    const ObservationSource*    source, /* IMarketObservable::getDefaultObsSource() */
	    const FixingType*           fixType,
	    const IObservationOverride* overrides,
	    const SamplingConvention*   sampleRule /* &UnadjustedConvention() */
        ) const;

    // IMarketObservable - retrieve a single observation date
    // Returns false if obs is to be omitted
    virtual bool observationDate(const DateTime&           sampleDate,
                                 const ObservationSource*  source,
                                 const SamplingConvention* sampleRule,
                                 DateTime*                 obsDate) const;
    
    // the IMarketObservable interface for retrieving past samples events
    virtual double addPastSampleEvent(const DateTime&             sampleDate,
                                    const ObservationType*      obsType,
                                    const ObservationSource*    source,
                                    const FixingType*           fixType,
                                    const IObservationOverride* overrides,
                                    const SamplingConvention*   sampleRule,
                                    PastSamplesCollector*        collector) const;

	/* IMarketObservable:: */
    virtual bool isHoliday(const DateTime& sampleDate, 
                           const ObservationSource* source) const;

protected:
    virtual AssetHistoryConstSP getAssetHistory(const string &source) const = 0;
};

/** Base class for IndexSpec */
class MARKET_DLL IndexSpec : public CObject, 
                  public virtual IIndexSpec
{
public:
    static CClassConstSP const TYPE;

    /************************ fields ************************/
    string name; // exported, name of this index spec

    DateTimeArray resetDates; // transient, populated via addResetDates()
    DateTime      today;  // transient, populated via setup()
    bool setupCalled;         // transient

    /************************ types ************************/
    struct MARKET_DLL EqualKey { // to be used with STL maps, sets,...
        bool operator()(const smartConstPtr<IndexSpec> a1, 
                        const smartConstPtr<IndexSpec> a2) const {
            return a1->isSame(a2.get());
        }
    };
    struct MARKET_DLL HashFcn { // to be used with STL maps, sets,...
        size_t operator()(const smartConstPtr<IndexSpec> a) const {
            return (size_t)a->hashCode();
        }
    };
    virtual int hashCode(void) const; ///< returns a hash value

    /************************ methods ************************/
    /* IProdCreator:: */
    virtual void addResetDates(const DateTimeArray &resetDates);

    /* IProdCreator:: */
    void setup(const IModel* model, const MarketData* market);

    /* FixingReqEvent::IEventHandler:: */
    void getEvents(
        const FixingReqEvent*, IModel* model,
        const DateTime& eDate, EventResults* events
    ) const;

    /* IIndexSpec:: */
    virtual string getName() const { return name; }
    virtual bool isSame(const IndexSpec* spec) const; ///< compare IndexSpec's

    IndexSpec(const CClassConstSP& type, const string &name) 
        : CObject(type), name(name), setupCalled(false) {}

protected:
    IndexSpec(const CClassConstSP& type=TYPE) : CObject(type), setupCalled(false) {}

private:
    static void load(CClassSP& clazz);
};

// smart pointer support
typedef smartPtr<IndexSpec> IndexSpecSP;
typedef smartConstPtr<IndexSpec> IndexSpecConstSP;
typedef array<IndexSpecSP, IndexSpec> IndexSpecArray;

DRLIB_END_NAMESPACE
#endif
