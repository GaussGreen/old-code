//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RefLevel.hpp
//
//   Description : Provides Different Ways of Handling "averaging in"
//
//   Author      : Mark A Robson
//
//   Date        : 26 Sep 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_REFLEVEL_HPP
#define EDR_REFLEVEL_HPP

#include "edginc/DateTime.hpp"
#include "edginc/DoubleMatrix.hpp"
#include "edginc/StateVariableCollector.hpp"
#include "edginc/StateVariableClient.hpp"
#include "edginc/HistoricalContext.hpp"
#include "edginc/Model.hpp"
#include "edginc/MarketData.hpp"

DRLIB_BEGIN_NAMESPACE
class IMultiFactors;
class IPastValues;

/** Provides Different Ways of Handling "averaging in" */
class MCARLO_DLL IRefLevel: virtual public IObject{
public:
    static CClassConstSP const TYPE;

	virtual void GetMarket(const IModel* model, const CMarketDataSP market) = 0;

    /** returns simulation dates (of all assets) which are
        strictly in the future */
    virtual DateTimeArray getFutureDates(const DateTime& valueDate) const = 0;

    /** returns all the simulation dates (of all assets) */
    virtual const DateTimeArray& getAllDates() const = 0;

    /** returns the date when a model should start simulating - will not
        return a date before today */
    virtual const DateTime& getSimStartDate(const DateTime& today) const = 0;

    /** returns the number of future dates for given asset */
    virtual int numFutureDates(const DateTime& valueDate,
                               int             iAsset) const = 0;

    /** returns the number of dates (past and future) for given asset */
    virtual int numDates(int iAsset) const = 0;

#if 0
    /** returns the number of assets */
    virtual int getNumAssets() const = 0;
#endif

    /** Returns true if each asset has the same set of simulation dates */
    virtual bool  sameDatesPerAsset() const = 0;

    /** Interface geared for interaction between MonteCarlo and RefLevel.
        A RefLevel can create an object implementing this interface */
    class MCARLO_DLL IMCPath{
    public:
        /** Main method. Returns refLevel given values for future path
            of specified asset (this is strictly after the simulation
            start date). The array must run from the next [future]
            date in the RefLevel to the last date */
        virtual double refLevel(int            iAsset,
                                const double*  futurePath) const = 0;

        /** Used for estimating refLevel when future path depends upon
            refLevel (ie vol interp level for LN) */
        virtual double refLevel(int                  iAsset,
                                const IMultiFactors* asset) const = 0;

        /** Returns an array of integers denoting which assets have a
            simulation date given the index into getAllDates(). The
            number of the assets (0, 1, 2, ..., numAssets-1) is the
            same as the order in which the series is supplied to the
            constructor */
        virtual const IntArray& assetsOnDate(int index) const = 0;

        //// utility method - returns the number of future points strictly
        //// after simulation start date
        virtual int numFutureDates(int iAsset) const = 0;

        virtual ~IMCPath(){}
    };
    typedef refCountPtr<IMCPath> IMCPathSP;

    /** Creates an object implementing the IMCPath interface. Note
        that the monte carlo will only supply values strictly after
        the simulation start */
    virtual IMCPath* createMCPath(const DateTime&    valueDate,
                                  const DoubleArray& spotsAtSimStart,
                                  const IPastValues* pastValues) const = 0;

    /** This method is to be retired once we move to 'state variables' for the
        Monte Carlo. This is a hack to get around something that was designed in
        namely the fact that the ref level is responsible for choosing the
        sim start date */
    virtual void setSimStartDateToFirstRefDate() = 0;

    /** Interface for the state variable for RefLevel. This is
        the type that products deal with in the payoff. The payoff obtains
        it by calling the getRefLevelSV() method below */
    class MCARLO_DLL IStateVar: public virtual IStateVariable,
                     public virtual IHistoricalContextGenerator {
    public:
        virtual ~IStateVar();

        /** Returns the reference level for given asset. For
            the PastPathGenerator then any future average values will
            be 'estimated' (current methodology is current spot). For
            a future path generator, then this must return 'estimated'
            values too (eg same values as the Past generator) until
            generatePath is invoked. At which point it will contain
            the correct value for the current path */
        virtual double refLevel(int iAsset) = 0;

        // stateless ref level
        virtual double refLevel2(int iAsset, IHistoricalContextSP history) = 0;
    };
    typedef smartPtr<IStateVar> IStateVarSP;

    /** A Generator of MC RefLevel Variables */
    class MCARLO_DLL IStateVarGen: public virtual IStateVariableGen,
    public virtual IStateVariableClient,
    public virtual VirtualDestructorBase
    {
    public:
        IStateVarGen(){}

        virtual ~IStateVarGen();

        /** Returns a RefLevel state variable which then provides
            access to the reference level. This is the method that
            products should call to get an IRefLevel::IStateVar. The
            previous IRefLevel::IStateVarSP (may be null) should be
            passed in. The return object may or may not be the same as
            oldStateVar.*/
        virtual IStateVarSP getRefLevelSV(IStateVarSP               oldStateVar,
                                          IStateVariableGen::IStateGen* pathGen) const = 0;

        /** Appends 'true' (ie non derived) state variable generators
            required to the supplied collector. Implementations typically call
            IStateVariableCollector::append */
        virtual void collectStateVars(
            IStateVariableCollectorSP svCollector) const = 0;

        /** Returns the dates used for computing the reference level
            of iAsset */
        virtual const DateTimeArray& getDates(int iAsset) const = 0;

        /** Returns the dates used for computing the reference level
            of iAsset */
        virtual const DateTimeArray& getAllDates() const = 0;
    };
    typedef smartPtr<IStateVarGen> IStateVarGenSP;

    /** How to get a IStateVariableGen from a RefLevel */
    virtual IStateVarGen* createStateVarGen(
        const IMultiFactors* multiFactors,
        const DateTime&      valueDate) const = 0;

    class MCARLO_DLL Util{
    public:
        /** Creates a RefLevel where each asset has the same set of simulation
            dates. */
        static IRefLevel* makeAverage(
            const DateTimeArray&  commonSimDates);

        /** Creates a RefLevel where each asset 'averages' in on a single
            day. SimStart date is today */
        static IRefLevel* makeTrivialAverage(const DateTime&   singleDate);

        /** Creates a Forward Starting reference level. SimStart is
            the fwdStartDate date */
        static IRefLevel* makeFwdStart(const DateTime&   fwdStartDate);

        /** creates a reference level which is always zero for all assets -
            essentially a dummy one to use if you want to handle
            the reference level yourself. It does not require any
            historic values */
        static IRefLevel* makeZero(const DateTime&   simStartDate);

        /** Creates a RefLevel where each asset has the same set of N simulation
            dates, and the level is the avg of the best/worst M out of the N. */
        static IRefLevel* makeAvgMofN(
            const DateTimeArray&  commonSimDates,
            bool                  isRefBest,
            int                   refNbKeep);
    };

private:
    static void load(CClassSP& clazz);

    class Average;
    class AvgMofN;
    class FwdStart;
    class Zero;
};
typedef smartConstPtr<IRefLevel> IRefLevelConstSP;
typedef smartPtr<IRefLevel> IRefLevelSP;
#ifndef QLIB_REFLEVEL_CPP
EXTERN_TEMPLATE(class MCARLO_DLL_SP smartConstPtr<IRefLevel>);
EXTERN_TEMPLATE(class MCARLO_DLL_SP smartPtr<IRefLevel>);
#else
INSTANTIATE_TEMPLATE(class MCARLO_DLL smartConstPtr<IRefLevel>);
INSTANTIATE_TEMPLATE(class MCARLO_DLL smartPtr<IRefLevel>);
#endif

DRLIB_END_NAMESPACE

#endif

