//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : PastValues.hpp
//
//   Description : Records historic values for 1 or more assets
//
//   Author      : Mark A Robson
//
//   Date        : 19 Oct 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_PASTVALUES_HPP
#define EDR_PASTVALUES_HPP

#include "edginc/DateTime.hpp"
#include "edginc/CashFlow.hpp"
#include "edginc/Theta.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/ObservationMap.hpp"
#include "edginc/EventResults.hpp"

DRLIB_BEGIN_NAMESPACE
class IMultiMarketFactors;
class CDoubleMatrix;
class CAsset;
class SimSeries;

/** Defines interface which Provides views of simulation dates for
    multiple assets including past values */
class MARKET_DLL IPastValues: virtual public IObject{
public:
    static CClassConstSP const TYPE;
    
public:
    /** Returns the total number of past values */
    virtual int getNbPastValues(int iAsset) const = 0;

    /** Returns the total number of past values 
        up to and including 'date' */
    virtual int getNbPastValues(const DateTime& date,
                                int             iAsset) const = 0;

    /** Returns an array holding the values on the requested dates. The
        length of the array is equal to the number of historic dates.
        An exception is thrown if there are no values for any of the dates
        or any historic values are 0. The supplied dates must be in
        increasing order */
    virtual DoubleArray getPastValues(const DateTimeArray& dates,
                                      int                  iAsset,
                                      const DateTime&      today) const = 0;


    /** Returns the number of assets for which this object has data for */
    virtual int getNumAssets() const = 0;

    // validates and fills in past samples if possible
    // ensures all dates are there and all past dates have valid samples
    // For past dates the provided samples are used and if missing
    // samples are retrieved from the centralised asset history
    virtual void validate(const DateTime& today,
                          const DateTimeArrayArray& requiredDates,
                          IMultiFactors* assets,
                          SamplingConvention* dateAdjust) = 0;

    // Retrieves all past samples used by an instrument
    // For past dates the provided samples are used and if missing
    // samples are retrieved from the centralised asset history
    virtual void pastSamplesEvents(const DateTime& today,
                          const DateTimeArrayArray& requiredDates,
                          IMultiFactors* assets,
                          SamplingConvention* dateAdjust,
                          EventResults* events) = 0;

    /** Returns observation types if they exist */
    virtual ObservationTypeArray* getObsTypes(IMultiFactors* assets);

    /** Returns observation types if they exist */
    virtual ObservationSourceArray* getSources(IMultiFactors* assets);

    /** Indicates whether this PastValues object uses sampling rules to construct a history */
    virtual bool hasSamplingConvention();

    /** Rolls the PastValues through time populating any samples as they
        become historic */
    virtual void roll(const Theta::Util&   thetaUtil,
                      int                  iAsset,
                      const CAsset*        asset) = 0;

    /** rolls past values for all assets inside a MultiFactor */
    virtual void roll(const Theta::Util&         thetaUtil,
                      const IMultiMarketFactors* multiFactor) = 0;

    /** Instances of this class are thrown when a zero historic value 
        (or indeed no value at all) is tried to be used. Could argue
        whether this should derive from ModelException or not */
    class MARKET_DLL MissingSampleException: public ModelException{
    public:
        MissingSampleException(const string&    routine,
                               int              assetIndex,
                               const DateTime&  sampleDate);

        virtual ~MissingSampleException() throw ();

        /** returns the index of the asset which had the missing sample */
        int getAssetIndex() const;

        /** records the asset name whose sample was missing */
        void setAssetName(const string& name);

        /** creates a [deep] copy of the exception */
        virtual ModelException* clone() const;
        
        /** indicates whether this exception is derived from ModelException -
            used to drive whether this is stored as a 'cause' when a new
            exception is created */
        virtual bool isDerived() const;

        /** Returns the 'cause' exception if it's a MissingSampleException
            otherwise returns null. Equivalent to
            dynamic_cast<MissingSampleException*>e.getCause() */
        static MissingSampleException* getInstance(ModelException& e);

        MissingSampleException(
            const MissingSampleException& e);

    private:
        MissingSampleException(
            const ModelException& e,
            int                   assetIndex,
            const DateTime&       sampleDate);
        
        // fields
        int       assetIndex;
        DateTime  sampleDate;
    };

    class MARKET_DLL PositivePastValueCheck{
    public:
        static bool isValid(double value);
    }; 

    class MARKET_DLL NonNegativePastValueCheck{
    public:
        static bool isValid(double value);
    };

    /** Utility methods for creating various types of PastValues objects */
    class MARKET_DLL Util{
    public:
        /** Creates a PastValues where each asset has the same set of simulation
            dates. */
        static IPastValues* makeSimple(
            const DateTimeArray&  commonDates,
            const CDoubleMatrix&  pastValues);   // numCols = numAssets

        /** Creates a PastValues where each asset which has one historic
            past value  */
        static IPastValues* makeTrivial(
            const DateTime&      singleDate,
            const DoubleArray&   pastValues); // one entry per asset

        /** Creates a PastValues where each asset has the same set of sample
            dates using a NonNegativePastValue policy */
        static IPastValues* makeSimple(
            const DateTimeArray&       commonSimDates,
            const CDoubleMatrix&       pastValues,  // numCols = numAssets
            NonNegativePastValueCheck  notUsed);

        /** Creates a PastValues for a single asset which has one historic
            past value */
        static IPastValues* makeTrivial(
            const DateTime&      singleDate,
            const double         pastValue); // one entry per asset

        /** Creates a PastValues where each asset can have a different
            set of sample dates. */
        static IPastValues* makeGeneral(
            const CashFlowCluster&  samplesPerAsset);

    };

private:
    static void load(CClassSP& clazz);

#if 0
    // Simple class is not nested due to NT linking problems with nested templates
    // Implemented in PastValues.cpp
    template <class ValidPastValueCheck>
    class Simple;
#endif
    
    class General;
};
typedef smartConstPtr<IPastValues> IPastValuesConstSP;
typedef smartPtr<IPastValues> IPastValuesSP;
#ifndef QLIB_PASTVALUES_CPP
EXTERN_TEMPLATE(class MARKET_DLL_SP smartConstPtr<IPastValues>);
EXTERN_TEMPLATE(class MARKET_DLL_SP smartPtr<IPastValues>);
#else
INSTANTIATE_TEMPLATE(class MARKET_DLL smartConstPtr<IPastValues>);
INSTANTIATE_TEMPLATE(class MARKET_DLL smartPtr<IPastValues>);
#endif


DRLIB_END_NAMESPACE

#endif

