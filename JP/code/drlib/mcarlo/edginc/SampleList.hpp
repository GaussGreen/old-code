//   Filename    : SampleList.hpp
//
//   Description : SampleList representation
//
//   Author      : Stephen Hope
//
//   Date        : 27 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDGSAMPLELIST_HPP
#define EDGSAMPLELIST_HPP

#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/Asset.hpp"
#include "edginc/RefLevel.hpp"
#include "edginc/PastValues.hpp"
#include "edginc/MaturityPeriod.hpp"
#include "edginc/BadDayConvention.hpp"
#include "edginc/SamplingConvention.hpp"
#include "edginc/EventResults.hpp"

using namespace std;

DRLIB_BEGIN_NAMESPACE
class CVolProcessedBS;
class SimSeries;
class SampleList;
typedef smartConstPtr<SampleList> SampleListConstSP;
typedef smartPtr<SampleList> SampleListSP;
#ifndef QLIB_SAMPLELIST_CPP
EXTERN_TEMPLATE(class MCARLO_DLL_SP smartConstPtr<SampleList>);
EXTERN_TEMPLATE(class MCARLO_DLL_SP smartPtr<SampleList>);
EXTERN_TEMPLATE(IObjectSP MCARLO_DLL FieldGetSmartPtr<SampleListSP>(SampleListSP* t));
EXTERN_TEMPLATE(void MCARLO_DLL FieldSetSmartPtr<SampleListSP>(SampleListSP* t,
                                                    IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class MCARLO_DLL smartConstPtr<SampleList>);
INSTANTIATE_TEMPLATE(class MCARLO_DLL smartPtr<SampleList>);
INSTANTIATE_TEMPLATE(IObjectSP MCARLO_DLL FieldGetSmartPtr<SampleListSP>(SampleListSP* t));
INSTANTIATE_TEMPLATE(void MCARLO_DLL FieldSetSmartPtr<SampleListSP>(SampleListSP* t,
                                                         IObjectSP o));
#endif

/** Holds data corresponding to a set of samples for a single underlying.
    The samples can be weighted */
class MCARLO_DLL SampleList: public CObject,
                  virtual public IRefLevel,
                  virtual public IPastValues{
public:
    static CClassConstSP const TYPE;
    friend class SampleListHelper;

   /** Constructor */
    SampleList(const DateTimeArray&  dates,
               const DoubleArray&    values,
               const DoubleArray&    weights);

    /** Create a sample using a generated array of dates formed from first and
        last date and rules for generating the intervening dates. The
        resultant dates (exluding the specified first and last dates)
        are then business date adjusted. */
    static SampleListSP create(
        const DateTime&           valueDate, // to define historical samples
        const DateTime&           firstDate,
        const DateTime&           lastDate,
        const ExpiryConstSP&       matPeriod,
        bool                       stubAtEnd,
        const BadDayConventionConstSP& badDayConv,
        const HolidayConstSP&      hols,
        double                     histAvgToDate, // historic average to date
        double                     weightScaling); /* weights = weightScaling/
                                                      num samples*/

    /** Destructor */
    ~SampleList();

	virtual void GetMarket(const IModel* model, const CMarketDataSP market);
    /** SampleList validation */
    virtual void validatePop2Object();

    // PUBLIC METHODS

    /** How to get a IMCStateVarGen (because we implement IRefLevel) */
    virtual IStateVarGen* createStateVarGen(
        const IMultiFactors* multiFactors,
        const DateTime&      valueDate) const;

    /** returns the weighted sum of samples so far */
    double sumToDate(const DateTime& today)const;

    /** returns the average of the samples so far */
    double averageToDate(const DateTime& today)const;

    /** returns the expected life of a sample:
        expected life = 
        [ sum for i = 1 to N (weight(i) * time(i))]
        -----------------------------------------
        [sum for i = 1 to N (weight(i)]
        
        where N is number of future weights. */
    double expectedLife(const DateTime& today)const;

    /** returns the expected date of a sample (cf expectedLife) */
    DateTime expectedEndDate(const DateTime& today)const;

    /** returns the number of samples remaining in the future */
    int countFutureSamples(const DateTime& today)const;

    /** returns the weighting of samples that remain in the future */
    double futureWeight(const DateTime& today)const;

    /** return the first and last dates in a sample list */
    void getBoundingDates(DateTime& loDate,
                          DateTime& hiDate)const;

    /** Returns the first date in the sample list */
    const DateTime& getFirstDate() const;

    /** Returns the last date in the sample list */
    const DateTime& getLastDate() const;

    /** returns TRUE if all the weights sum to 1.0 
        (within a sensible tolerance) ? */
    bool weightsSumToOne()const;

    /** returns an estimate of the weighted sum of future samples */
    double futureSampleSum(const CAsset* asset,
                           const DateTime& today)const;

    /** returns the expected average */
    double expectedAverage(const CAsset* asset,
                           const DateTime& today)const;

    /** returns the average variance */
    double averageVariance(const CVolProcessedBS* procVol,
                           const DateTime& today,
                           bool  usePastWeights)const;

    /** returns the co-variance between two sample dates */
    double averageCovariance(const SampleList* other,
                             const CVolProcessedBS* procVol,
                             const DateTime& today)const;
    
    /** return average mode level.  Note Tested 
    double averageModeLevel(const CAsset* asset,
                            const CVolProcessedBS* procVol,
                            const DateTime& today)const;*/
    
    /** returns the value and weight corresponding to today,
	using the STAIRS interpolation on today in the dates array */
    bool getPreviousSample(const DateTime &today, 
                           DateTime&       previousDate, 
                           double&         value, 
                           double&         weight) const;
    
    /** populate future samples up to 'rollDate' with either the spot or the
        forward price of an asset
    */
    void roll(const Asset*    asset,
              const DateTime& today,
              const DateTime& rollDate,
              bool  useSpot = true);

	/** populate samples between date1 < date2, assuming that the spot level is
	    value1 and value2 on these dates, and using linear interpolation. */
	void rollLinear(const DateTime& date1, double value1,
					const DateTime& date2, double value2);

    /** Returns the dates in the sample */
    const DateTimeArray& getDates() const;

	const DoubleArray& getValues() const;

    /** Returns the weights in the sample. Need to review interaction
        of MC with this class */
    const DoubleArray& getWeights() const;

    /* methods below are implementation of IRefLevel and IPastValues
       interface */

    /** returns simulation dates (of all assets) which are
        strictly in the future */
    virtual DateTimeArray getFutureDates(const DateTime& valueDate) const;

    /** returns all the simulation dates (of all assets) */
    virtual const DateTimeArray& getAllDates() const;

    /** returns the date when a model should start simulating - will not
        return a date before today */
    virtual const DateTime& getSimStartDate(const DateTime& today) const;

    /** Returns all the simulation dates of the given asset */
    virtual const DateTimeArray& getDates(int iAsset) const;

    /** returns the number of future dates for given asset */
    virtual int numFutureDates(const DateTime& valueDate,
                               int             iAsset) const;
    
    /** returns the number of dates (past and future) for given asset */
    virtual int numDates(int iAsset) const;

    /** returns the number of assets ie 1 */
    virtual int getNumAssets() const;

    /** Returns true */
    virtual bool sameDatesPerAsset() const;

    /** Creates an object implementing the IMCPath interface */
    virtual IMCPath* createMCPath(const DateTime&    valueDate,
                                  const DoubleArray& spotsAtSimStart,
                                  const IPastValues* pastValues) const;

    /** This method is to be retired once we move to 'state variables' for the
        Monte Carlo. This is a hack to get around something that was designed in
        namely the fact that the ref level is responsible for choosing the
        sim start date */
    virtual void setSimStartDateToFirstRefDate();

    /** Returns the total number of past values */
    virtual int getNbPastValues(int iAsset) const;

    /** Returns the total number of past values 
        up to and including 'date' */
    virtual int getNbPastValues(const DateTime& valueDate,
                                int             iAsset) const;

    /** Returns an array holding the values on the requested dates. The
        length of the array is equal to the number of historic dates.
        An exception is thrown if there are no values for any of the dates
        or any historic values are 0. The supplied dates must be in
        increasing order */
    virtual DoubleArray getPastValues(const DateTimeArray& dates,
                                      int                  iAsset,
                                      const DateTime&      today) const;

    /** Rolls the PastValues through time populating any samples as they
        become historic */
    void roll(const Theta::Util&   thetaUtil,
              int                  iAsset,
              const Asset*         asset);

    /** rolls past values for all assets inside a MultiFactor */
    void roll(const Theta::Util&         thetaUtil,
              const IMultiMarketFactors* multiFactor);

    /** Records the dates for which this object needs points for 
        within a simulation */
    virtual void recordDates(SimSeries* simSeries) const;

    /** Generate an array of dates given a first and last date and rules
        for generating the intervening dates. The resultant dates
        (exluding the specified first and last dates) are then business
        date adjusted. Not sure where this method belongs but can go here
        for now */
    static DateTimeArraySP generateDates(
        const DateTime&                firstDate,
        const DateTime&                lastDate,
        const ExpiryConstSP&           matPeriod,
        bool                           stubAtEnd,
        const BadDayConventionConstSP& badDayConv,
        const HolidayConstSP&          hols);
        
    // Retrieves all past samples used by an instrument
    // For past dates the provided samples are used and if missing
    // samples are retrieved from the centralised asset history
    virtual void pastSamplesEvents(const DateTime& today,
                          const DateTimeArrayArray& requiredDates,
                          IMultiFactors* assets,
                          SamplingConvention* dateAdjust,
                          EventResults* events);

private:
    virtual void validate(const DateTime& today,
                          const DateTimeArrayArray& requiredDates,
                          IMultiFactors* assets,
                          SamplingConvention* dateAdjust);

    SampleList();
    SampleList& operator=(const SampleList& rhs);
    SampleList(const SampleList& rhs);
    class MCPath;
    friend class MCPath;
    
    DateTimeArray   dates;
    DoubleArray     values;
    DoubleArray     weights;
    bool            startSimAtFirstAvgDate; // to be removed
};



DRLIB_END_NAMESPACE

#endif
