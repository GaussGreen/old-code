//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : RefLevel.cpp
//
//   Description : Handles idea of a reference level for a series of assets
//
//   Author      : Mark A Robson
//
//   Date        : 25 Sep 2001
//
//
//----------------------------------------------------------------------------
#include "edginc/config.hpp"
#define QLIB_REFLEVEL_CPP
#include "edginc/RefLevel.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/PastValues.hpp"
#include "edginc/Algorithm.hpp"
#include "edginc/Format.hpp"
#include "edginc/SVGenSpot.hpp"
#include "edginc/DateBuilder.hpp"
#include "edginc/ExplicitDates.hpp"

DRLIB_BEGIN_NAMESPACE
IRefLevel::IStateVar::~IStateVar(){}
IRefLevel::IStateVarGen::~IStateVarGen(){}

/** Class implementing IRefLevel where each asset has the same set of
    "averaging in" dates */
class IRefLevel::Average: public CObject,
                 virtual public IRefLevel{
public:
    class MCPath;
private:
    friend class  MCPath;
    DateTimeArray avgDates;
    bool          startSimAtFirstAvgDate; // to be removed
public:
    static CClassConstSP const TYPE;

    virtual void GetMarket(const IModel* model, const CMarketDataSP market) {
    }

    /** returns all the simulation dates (of all assets) */
    const DateTimeArray& getAllDates() const{
        return avgDates;
    }

    /** returns the date when a model should start simulating - will not
        return a date before today */
    virtual const DateTime& getSimStartDate(const DateTime& today) const{
        if (startSimAtFirstAvgDate){
            // start simulation on start date
            return (avgDates.front().isGreater(today)? avgDates.front(): today);
        }
        // otherwise start simulation immediately - not at first average date
        return today;
    }

    /** returns simulation dates (of all assets) which are strictly in
        the future */
    virtual DateTimeArray getFutureDates(const DateTime& valueDate) const{
        return valueDate.getFutureDates(avgDates);
    }

    /** Returns true if each asset has the same set of simulation dates */
    bool sameDatesPerAsset() const{
        return true;
    }

    /** returns the number of dates (past and future) for given asset */
    virtual int numDates(int iAsset) const{
        return avgDates.size();
    }

#if 0
    /** returns the number of assets */
    virtual int getNumAssets() const{
        return nbAssets;
    }
#endif

    /** returns the number of future dates for given asset */
    virtual int numFutureDates(const DateTime& valueDate,
                               int             iAsset) const{
        return valueDate.numFutureDates(avgDates);
    }

    class MCPath: public IRefLevel::IMCPath{
    private:
        const Average*  avgRefLevel;
        int             nbAssets;
        int             numDates; // for ease (same per asset)
        int             numFutureValues; // same per asset
        DoubleArray     sumSoFar;  // per asset
        IntArray        assetsIndexes;
    public:
        // constructor for when we start simulation at first avg date and
        // value date is before first avg date
        MCPath(const Average*     avgRefLevel,
               const DoubleArray& spotsAtSimStart):
            avgRefLevel(avgRefLevel),
            nbAssets(spotsAtSimStart.size()),
            numDates(avgRefLevel->numDates(0)),
            numFutureValues(numDates-1), // value at first date is known
            sumSoFar(nbAssets){
            // build up [trivial] assetsIndexes
            assetsIndexes = CIntArray(nbAssets);
            for (int i = 0; i < nbAssets; i++){
                assetsIndexes[i] = i;
            }

            // then calculate sumSoFar
            for (int a = 0; a < nbAssets; a++){
                sumSoFar[a] = spotsAtSimStart[a];
            }
        }

        // constructor for when we start simulation today
        MCPath(const Average*     avgRefLevel,
               const IPastValues* pastValues,
               const DateTime&    valueDate):
            avgRefLevel(avgRefLevel),
            nbAssets(pastValues->getNumAssets()),
            numDates(avgRefLevel->numDates(0)),
            numFutureValues(avgRefLevel->numFutureDates(valueDate, 0)),
            sumSoFar(nbAssets){
            // build up [trivial] assetsIndexes
            assetsIndexes = CIntArray(nbAssets);
            for (int i = 0; i < nbAssets; i++){
                assetsIndexes[i] = i;
            }

            // then calculate sumSoFar
            for (int a = 0; a < nbAssets; a++){
                DoubleArray past =
                    pastValues->getPastValues(avgRefLevel->avgDates,
                                              a,
                                              valueDate);
                for (int i = 0; i < past.size(); i++){
                    sumSoFar[a] += past[i];
                }
            }
        }
        /** Main method. Returns refLevel given values for future path of
            specified asset */
        virtual double refLevel(int            iAsset,
                                const double*  futurePath) const{
            double assetSum = sumSoFar[iAsset];
            for (int i = 0; i < numFutureValues; i++){
                assetSum += futurePath[i];
            }
            return (assetSum/numDates);
        }
        /** Used for estimating refLevel when future path depends upon
            refLevel (ie vol interp level for LN) */
        virtual double refLevel(int                  iAsset,
                                const IMultiFactors* mAsset) const{
            double assetSum = sumSoFar[iAsset];
            if (numFutureValues > 0){
                // extrapolate just using spot
                double spot = mAsset->assetGetSpot(iAsset);
                assetSum += numFutureValues * spot;
            }
            return (assetSum/numDates);
        }

        /** Returns an array of integers denoting which assets have a
            simulation date given the index into getAllDates(). The
            number of the assets (0, 1, 2, ..., numAssets-1) is the
            same as the order in which the series is supplied to the
            constructor */
        const CIntArray& assetsOnDate(int index) const{
            return assetsIndexes; // independent of index
        }

        //// utility method - returns the number of future points
        virtual int numFutureDates(int iAsset) const{
            return numFutureValues;
        }

    };

    IRefLevel::IMCPath* createMCPath(const DateTime&    valueDate,
                                     const DoubleArray& spotsAtSimStart,
                                     const IPastValues* pastValues) const{
        if (startSimAtFirstAvgDate && valueDate.isLess(avgDates.front())){
            return new MCPath(this, spotsAtSimStart);
        }
        // otherwise since our simStartDate starts immediately we don't need
        // to use spotsAtSimStart
        return new MCPath(this, pastValues, valueDate);
    }

    /** This method is to be retired once we move to 'state variables' for the
        Monte Carlo. This is a hack to get around something that was designed in
        namely the fact that the ref level is responsible for choosing the
        sim start date */
    virtual void setSimStartDateToFirstRefDate(){
        startSimAtFirstAvgDate = true;
    }

    /** Creates a IRefLevel where each asset has the same set of simulation
        dates. */
    Average(const DateTimeArray& avgDates): 
        CObject(TYPE), avgDates(avgDates), startSimAtFirstAvgDate(false){
        validatePop2Object();
    }

    virtual void validatePop2Object(){
        DateTime::ensureIncreasing(avgDates, "Average dates", true);
    }

    /**  State variable for this RefLevel. This is
         the type that products deal with in the payoff. */
    class StateVar: public virtual IStateVar{
        SVGenSpotSP             mcSpot;            //!< Spot path generator
        SVGenSpot::IStateVarSP  mcSpotSV;          //!< Spot path state variable
        vector<bool>         allDone;           //!< Whether reflevel has been computed
        vector<double>       sumSoFar;          //!< Sum so far for each asset
        vector<double>       spotsAtSimStart;   //!< Used to estimate a RefLevel
        int                  numDates;          //!< Number of reference dates
        bool                 isDoingPast;       //!< Whether we are doing past or future
    public:
        virtual ~StateVar(){}
        virtual bool doingPast() const{
            return isDoingPast;
        }
        //// update mc spot state variable (going from past to future)
        void update(IStateVariableGen::IStateGen* pathGen){
            mcSpotSV    = mcSpot->getSpotSV(pathGen);
            isDoingPast = mcSpotSV->doingPast();

            // Compute reflevels in case past won't be invoked
            for(unsigned int iAsset = 0; iAsset < allDone.size(); iAsset++) {
                refLevel(iAsset);
            }
        }

        class HistoricalContext : public IHistoricalContext {
            // this is used in the case where the product context can not be used
            // If the product has enough information about historical context, we
            // can declare something like
            // class product::HistoricalContext : ... public virtual statevar::HistoricalContext
            // that way, we can pass in product history context
            // in the tarn example, the statevar is storing its own copy of sumSoFar,
            // so i move it here
        public:
            vector<bool> allDone;
            vector<double> sumSoFar;
            vector<double> spotsAtSimStart;
            int numReferenceDates;
            int numDatesSoFar;

            HistoricalContext() {
                // determine whether any computation is necessary
                if (numDatesSoFar == numReferenceDates) {
                    for (unsigned i = 0; i < sumSoFar.size(); ++i)
                        allDone[i] = true;
                }
            }

            void deepCopyTo(IHistoricalContext* destination) const {
                HistoricalContext* dest = static_cast<HistoricalContext*>(destination);
                *dest = *this;
            }
        };

        virtual IHistoricalContextSP createHistoricalContext()
        {
            return IHistoricalContextSP(new HistoricalContext());
        }

        virtual IHistoricalContextSP getInitialHistoricalContext()
        {
            // xxx, revisit
            return IHistoricalContextSP(new HistoricalContext());
        }


        virtual double refLevel2(int iAsset, IHistoricalContextSP history)
        {
            // this function will be called per timepoint from the stateless payoff function
            // slotId represent current date
            // i guess an alternative is to pass in a DateTime object
            // note that this function can be called from a past or a future timepoint

            HistoricalContext* context =
                    static_cast<HistoricalContext*>(history.get());

            double sum;
            if (doingPast()) {
                sum = 0.0;
                context->allDone[iAsset] = false;
            } else {
                sum = context->sumSoFar[iAsset];
            }
            if (!context->allDone[iAsset]){
                // get curent spot
                double currentSpot = mcSpotSV->getSpotPrice(iAsset);
                context->sumSoFar[iAsset] += currentSpot;
                ++context->numDatesSoFar;

                // estimate the future levels
                if (context->numDatesSoFar < context->numReferenceDates){
                    sum += (context->numReferenceDates -
                            context->numDatesSoFar) * context->spotsAtSimStart[iAsset];
                } else {
                    // done for this path
                    context->allDone[iAsset] = true;
                }
            }
            return (sum/context->numReferenceDates);
        }

        /** Returns the reference level for given asset. When doing
            the the past any future average values will be 'estimated'
            (current methodology is current spot). For a future path
            generator, then this must return 'estimated' values too
            (eg same values as the Past generator) until generatePath
            is invoked. At which point it will contain the correct
            value for the current path */
        virtual double refLevel(int iAsset) {
            double sum;
            /** Switch off the allDone flag and start computation from
                scratch in case someone is calling the method multiple
                times in the past. */
            if(isDoingPast) {
                sum = 0.0;
                allDone[iAsset] = false;
            } else {
                sum = sumSoFar[iAsset];
            }
            if (!allDone[iAsset]){
                // then get SVPath
                const SVPath& path = mcSpotSV->path(iAsset);
                // evaluate sum for each asset
                int end = path.end();
                for (int iPath = path.begin(); iPath < end; iPath++){
                    sum += path[iPath];
                }
                if (isDoingPast){
                    sumSoFar[iAsset] = sum; // cache
                }
                // estimate the future levels
                if (numDates > end){
                    sum += (numDates - end) * spotsAtSimStart[iAsset];
                } else if (isDoingPast){
                    // all dates are in the past
                    allDone[iAsset] = true; // switch off for performance
                }
            }
            // then divide through
            return (sum/numDates);
        }

        /** constructor. */
        StateVar(const vector<double>& spotsAtSimStart,
                 SVGenSpotSP                  mcSpot,
                 IStateVariableGen::IStateGen* pathGen):
            mcSpot(mcSpot),
            allDone(spotsAtSimStart.size()),
            sumSoFar(spotsAtSimStart.size()),
            spotsAtSimStart(spotsAtSimStart),
            numDates(mcSpot->numDates()) {
                update(pathGen);
            }
    };

    /** A Generator of MC RefLevel Variables. Note that we need to tidy
        up these classes as there are 3 which are very similar. Suggests
        we create a base class for the SV and SVGen */
    class StateVarGen: public virtual IRefLevel::IStateVarGen{
        SVGenSpotSP       mcSpot; // this is our StateVar generator
        vector<double> spotsAtSimStart;
    public:
        virtual ~StateVarGen(){}

        StateVarGen(SVGenSpotSP mcSpot, const IMultiFactors* multiFactors,
                    int      numDates):
            mcSpot(mcSpot), spotsAtSimStart(multiFactors->NbAssets()){
            for (unsigned int iAsset = 0;
                 iAsset < spotsAtSimStart.size(); iAsset++){
                // simulation starts today
                spotsAtSimStart[iAsset] = multiFactors->assetGetSpot(iAsset);
            }
        }

        /** Create the corresponding State Variable for this State
            Variable Generator (NB Implies one state variable per
            generator). The previous IStateVariableSP (may be null)
            should be passed in.  The return object may or may not be
            the same as oldStateVar.  */
        virtual IStateVariableSP create(IStateVariableSP             oldStateVar,
                                     IStateVariableGen::IStateGen* pathGen) const{
            IStateVarSP myStateVar(&dynamic_cast<IStateVar&>(*oldStateVar));
            return getRefLevelSV(myStateVar, pathGen);
        }

        /** Returns a RefLevel state variable which then provides
            access to the reference level. This is the method that
            products should call to get an IRefLevel::IStateVar. The
            previous IRefLevel::IStateVarSP (may be null) should be
            passed in. The return object may or may not be the same as
            oldStateVar.*/
        virtual IStateVarSP getRefLevelSV(IStateVarSP               oldStateVar,
                                          IStateVariableGen::IStateGen* pathGen) const{
            if (oldStateVar.get()){
                // just update
                dynamic_cast<StateVar&>(*oldStateVar).update(pathGen);
                return oldStateVar;
            }
            return IStateVarSP(new StateVar(spotsAtSimStart, mcSpot, pathGen));
        }

        /** Appends 'true' (ie non derived) state variable generators
            required to the supplied collector. Implementations typically call
            IStateVariableCollector::append */
        virtual void collectStateVars(
            IStateVariableCollectorSP svCollector) const{
            svCollector->append(mcSpot.get());
        }

        /** Returns the dates used for computing the reference level
            of iAsset */
        const DateTimeArray& getDates(int iAsset) const {
            return mcSpot->getSimSeries()->getDates(iAsset);
        }

        /** Returns the dates used for computing the reference level
            of iAsset */
        const DateTimeArray& getAllDates() const {
            return mcSpot->getSimSeries()->getAllDates();
        }
    };

    /** How to get a IStateVariableGen from a RefLevel */
    virtual IStateVarGen* createStateVarGen(
        const IMultiFactors* multiFactors,
        const DateTime&      valueDate) const{
        SVGenSpotSP    mcSpot(new SVGenSpot(multiFactors->NbAssets(), avgDates));
        return new StateVarGen(mcSpot, multiFactors, avgDates.size());
    }

private:

    /* for reflection */
    Average(): CObject(TYPE), startSimAtFirstAvgDate(false){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Average, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRefLevel);
        EMPTY_SHELL_METHOD(defaultCon);
        FIELD(avgDates, "Averaging in dates");
        FIELD_NO_DESC(startSimAtFirstAvgDate);
        FIELD_MAKE_TRANSIENT(startSimAtFirstAvgDate);
    }

    static IObject* defaultCon(){
        return new Average();
    }
};

CClassConstSP const IRefLevel::Average::TYPE = CClass::registerClassLoadMethod(
    "IRefLevel::Average", typeid(IRefLevel::Average), load);


/** Creates a RefLevel where each asset has the same set of simulation
    dates. */
IRefLevel* IRefLevel::Util::makeAverage(
    const DateTimeArray&  commonSimDates){
    return new Average(commonSimDates);
}

/** Creates a RefLevel where each asset 'averages' in on a single
    day */
IRefLevel* IRefLevel::Util::makeTrivialAverage(
    const DateTime&      singleDate){
    DateTimeArray commonSimDates(1, singleDate);
    return makeAverage(commonSimDates);
}

/** Forward Starting RefLevel. Same date for all assets  */
class IRefLevel::FwdStart: public CObject, // avoid getting Average fields
                 virtual public IRefLevel{
private:
    DateTime      startDate;
    DateTimeArray allDates; // transient
public:
    static CClassConstSP const TYPE;

    virtual void GetMarket(const IModel* model, const CMarketDataSP market) {
    }

    /** returns all the simulation dates (of all assets) */
    const DateTimeArray& getAllDates() const{
        return allDates;
    }

    /** returns the date when a model should start simulating - will not
        return a date before today */
    virtual const DateTime& getSimStartDate(const DateTime& today) const{
        // start simulation on start date
        return (startDate.isGreater(today)? startDate: today);
    }

    /** returns simulation dates (of all assets) which are strictly in
        the future relative to supplied date */
    virtual DateTimeArray getFutureDates(const DateTime& date) const{
        return date.getFutureDates(allDates);
    }

    /** Returns true if each asset has the same set of simulation dates */
    bool sameDatesPerAsset() const{
        return true;
    }

    /** returns the number of dates (past and future) for given asset */
    virtual int numDates(int iAsset) const{
        return 1;
    }

#if 0
    /** returns the number of assets */
    virtual int getNumAssets() const{
        return nbAssets;
    }
#endif

    /** returns the number of future dates for given asset */
    virtual int numFutureDates(const DateTime& valueDate,
                               int             iAsset) const{
        return valueDate.numFutureDates(allDates);
    }


    class MCPath: public IRefLevel::IMCPath{
    private:
        DoubleArray     spotsAtStart;
        IntArray        assetsIndexes;
    public:
        // constructor
        MCPath(const DoubleArray& spotsAtStart,
               int                nbAssets): spotsAtStart(spotsAtStart){
            // build up [trivial] assetsIndexes
            assetsIndexes = CIntArray(nbAssets);
            for (int i = 0; i < nbAssets; i++){
                assetsIndexes[i] = i;
            }
        }
        /** Main method. Returns refLevel given values for future path of
            specified asset */
        virtual double refLevel(int            iAsset,
                                const double*  futurePath) const{
            return spotsAtStart[iAsset];
        }
        /** Used for estimating refLevel when future path depends upon
            refLevel (ie vol interp level for LN) */
        virtual double refLevel(int                  iAsset,
                                const IMultiFactors* mAsset) const{
            return spotsAtStart[iAsset];
        }

        /** Returns an array of integers denoting which assets have a
            simulation date given the index into getAllDates(). The
            number of the assets (0, 1, 2, ..., numAssets-1) is the
            same as the order in which the series is supplied to the
            constructor */
        const CIntArray& assetsOnDate(int index) const{
            return assetsIndexes; // independent of index
        }

        //// utility method - returns the number of future points after
        //// the sim start date
        virtual int numFutureDates(int iAsset) const{
            return 0;
        }

    };

    IRefLevel::IMCPath* createMCPath(const DateTime&    valueDate,
                                     const DoubleArray& spotsAtSimStart,
                                     const IPastValues* pastValues) const{
        DoubleArray spotsAtStart;
        if (valueDate.isGreaterOrEqual(startDate)){
            // spot at start is in the past, so use IPastValues
            int numAssets = pastValues->getNumAssets();
            spotsAtStart = DoubleArray(numAssets);
            for (int i = 0; i < numAssets; i++){
                spotsAtStart[i] = pastValues->getPastValues(allDates,
                                                            i,
                                                            valueDate)[0];
            }
        } else {
            spotsAtStart = spotsAtSimStart;
        }
        return new MCPath(spotsAtStart, pastValues->getNumAssets());
    }

    /** This method is to be retired once we move to 'state variables' for the
        Monte Carlo. This is a hack to get around something that was designed in
        namely the fact that the ref level is responsible for choosing the
        sim start date */
    virtual void setSimStartDateToFirstRefDate(){
        // nothing to do
    }

    /** Creates a IRefLevel where each asset has the same set of simulation
        dates. */
    FwdStart(const DateTime& startDate):
        CObject(TYPE), startDate(startDate){
        validatePop2Object();
    }

    virtual void validatePop2Object(){
        // build allDates
        allDates = DateTimeArray(1, startDate);
    }

    /**  State variable for this RefLevel. This is
         the type that products deal with in the payoff. */
    class StateVar: public virtual IStateVar{
        SVGenSpotSP             mcSpot;            //!< Spot path generator
        SVGenSpot::IStateVarSP  mcSpotSV;          //!< Spot path state variable
        vector<double>       estimatedLevels;   //!< The estimated level
        vector<double>       startLevels;       //!< The reference level
        vector<bool>         allDone;           //!< by asset
        bool                 isDoingPast;       //!< Whether we are doing past or future
    public:
        virtual ~StateVar(){}
        virtual bool doingPast() const{
            return mcSpotSV->doingPast();
        }
        //// update mc spot state variable (going from past to future)
        void update(IStateVariableGen::IStateGen* pathGen){
            mcSpotSV    = mcSpot->getSpotSV(pathGen);
            isDoingPast = mcSpotSV->doingPast();

            // Compute reflevels in case past won't be invoked
            for(unsigned int iAsset = 0; iAsset < allDone.size(); iAsset++) {
                refLevel(iAsset);
            }
        }


        virtual IHistoricalContextSP createHistoricalContext()
        {
            return IHistoricalContextSP();
        }


        virtual IHistoricalContextSP getInitialHistoricalContext()
        {
            // xxx, revisit
            return IHistoricalContextSP();
        }

        virtual double refLevel2(int iAsset, IHistoricalContextSP history)
        {
            return (0);
        }

        /** Returns the reference level for given asset. When doing
            the the past any future average values will be 'estimated'
            (current methodology is current spot). For a future path
            generator, then this must return 'estimated' values too
            (eg same values as the Past generator) until generatePath
            is invoked. At which point it will contain the correct
            value for the current path */
        virtual double refLevel(int iAsset) {
            /** Switch off the allDone flag and start computation from
                scratch in case someone is calling the method multiple
                times in the past. */
            if(isDoingPast) {
                allDone[iAsset] = false;
                startLevels[iAsset] = 0.0;
            }

            if (allDone[iAsset]){
                return startLevels[iAsset];
            }
            // else get SVPath
            const SVPath& path = mcSpotSV->path(iAsset);
            int pathEnd   = path.end();
            int pathBegin = path.begin();

            if((pathBegin == pathEnd) && isDoingPast) {
                // Reference level is in the future
                return estimatedLevels[iAsset];
            }

            // We can see tha path here
            if (isDoingPast) {
                // Path is in the past so save the value
                allDone[iAsset] = true;
                startLevels[iAsset] = path[0];
            }

            return (path[0]);
        }

        /** constructor. */
        StateVar(const vector<double>&     estimatedLevels,
                 SVGenSpotSP                  mcSpot,
                 IStateVariableGen::IStateGen* pathGen):
            mcSpot(mcSpot), estimatedLevels(estimatedLevels), startLevels(estimatedLevels.size()),
            allDone(startLevels.size()) {
                update(pathGen);
            }
    };

    /** A Generator of MC RefLevel Variables */
    class StateVarGen: public virtual IRefLevel::IStateVarGen{
        SVGenSpotSP              mcSpot; // this is our StateVar generator
        const IMultiFactors*  multiFactors;
        DateTime              valueDate;
        DateTime              startDate;
    public:
        virtual ~StateVarGen(){}

        StateVarGen(SVGenSpotSP mcSpot,
                    const IMultiFactors* multiFactors,
                    const DateTime&      valueDate,
                    const DateTime&      startDate):
        mcSpot(mcSpot), multiFactors(multiFactors), valueDate(valueDate),
        startDate(startDate) { }

        /** Create the corresponding State Variable for this State
            Variable Generator (NB Implies one state variable per
            generator). The previous IStateVariableSP (may be null)
            should be passed in.  The return object may or may not be
            the same as oldStateVar. */
        virtual IStateVariableSP create(IStateVariableSP             oldStateVar,
                                     IStateVariableGen::IStateGen* pathGen) const{
            IStateVarSP myStateVar(&dynamic_cast<IStateVar&>(*oldStateVar));
            return getRefLevelSV(myStateVar, pathGen);
        }

        /** Returns a RefLevel state variable which then provides
            access to the reference level. This is the method that
            products should call to get an IRefLevel::IStateVar. The
            previous IRefLevel::IStateVarSP (may be null) should be
            passed in. The return object may or may not be the same as
            oldStateVar. */
        virtual IStateVarSP getRefLevelSV(IStateVarSP               oldStateVar,
                                          IStateVariableGen::IStateGen* pathGen) const{
            if (oldStateVar.get()){
                // just update
                dynamic_cast<StateVar&>(*oldStateVar).update(pathGen);
                return oldStateVar;
            }

            int nbAssets = multiFactors->NbAssets();
            vector<double> estimatedLevels(nbAssets);
            for (int iAsset = 0; iAsset < nbAssets; iAsset++) {
                if (valueDate.isGreaterOrEqual(startDate)) {
                    // The state variable will extract it from the past
                    // and this field is not useful. We just put the current spot there
                    // spotLevels[iAsset] = mcSpot->getSpotSV(pathGen)->path(iAsset)[0];
                    estimatedLevels[iAsset] = multiFactors->assetGetSpot(iAsset);
                } else {
                    // Estimate reference level using the Fwd
                    estimatedLevels[iAsset] = multiFactors->assetFwdValue(iAsset, startDate);
                }
            }

            return IStateVarSP(new StateVar(estimatedLevels, mcSpot, pathGen));
        }
        /** Appends 'true' (ie non derived) state variable generators
            required to the supplied collector. */
        virtual void collectStateVars(
            IStateVariableCollectorSP svCollector) const{
            svCollector->append(mcSpot.get());
        }

        /** Returns the dates used for computing the reference level
            of iAsset */
        const DateTimeArray& getDates(int iAsset) const {
            return mcSpot->getSimSeries()->getDates(iAsset);
        }

        /** Returns the dates used for computing the reference level
            of iAsset */
        const DateTimeArray& getAllDates() const {
            return mcSpot->getSimSeries()->getAllDates();
        }
    };

    /** How to get a IStateVariableGen from a RefLevel */
    virtual IStateVarGen* createStateVarGen(
        const IMultiFactors* multiFactors,
        const DateTime&      valueDate) const{
        SVGenSpotSP    mcSpot(new SVGenSpot(multiFactors->NbAssets(), startDate));
        return new StateVarGen(mcSpot, multiFactors, valueDate, startDate);
    }

private:

    /* for reflection */
    FwdStart(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(FwdStart, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRefLevel);
        EMPTY_SHELL_METHOD(defaultCon);
        FIELD(startDate, "Start Date");
        FIELD(allDates, "Internal");
        FIELD_MAKE_TRANSIENT(allDates);
    }

    static IObject* defaultCon(){
        return new FwdStart();
    }
};

CClassConstSP const IRefLevel::FwdStart::TYPE =
CClass::registerClassLoadMethod(
    "IRefLevel::FwdStart", typeid(IRefLevel::FwdStart), load);

/** Creates a Forward Starting reference level. SimStart is
    the fwdStartDate date */
IRefLevel* IRefLevel::Util::makeFwdStart(
    const DateTime&   fwdStartDate){
    return new FwdStart(fwdStartDate);
}

/** Reference level that is always zero for all assets  */
class IRefLevel::Zero: public CObject,
                 virtual public IRefLevel{
private:
    DateTime      simStartDate;
    DateTimeArray allDates; // transient
public:
    static CClassConstSP const TYPE;

    virtual void GetMarket(const IModel* model, const CMarketDataSP market) {

    }

    /** returns all the simulation dates (of all assets) */
    const DateTimeArray& getAllDates() const{
        return allDates;
    }

    /** returns the date when a model should start simulating - will not
        return a date before today */
    virtual const DateTime& getSimStartDate(const DateTime& today) const{
        // start simulation on start date
        return (simStartDate.isGreater(today)? simStartDate: today);
    }

    /** returns simulation dates (of all assets) which are strictly in
        the future */
    virtual DateTimeArray getFutureDates(const DateTime& valueDate) const{
        return valueDate.getFutureDates(allDates);
    }

    /** Returns true if each asset has the same set of simulation dates */
    bool sameDatesPerAsset() const{
        return true;
    }

    /** returns the number of dates (past and future) for given asset */
    virtual int numDates(int iAsset) const{
        return 1;
    }

    /** returns the number of future dates for given asset */
    virtual int numFutureDates(const DateTime& valueDate,
                               int             iAsset) const{
        return valueDate.numFutureDates(allDates);
    }


    class MCPath: public IRefLevel::IMCPath{
    private:
        int             nbAssets;
        int             numFutureValues; // same per asset
        IntArray        assetsIndexes;
    public:
        // constructor
        MCPath(int nbAssets, int numFutureValues):
            nbAssets(nbAssets),
            numFutureValues(numFutureValues){
            // build up [trivial] assetsIndexes
            assetsIndexes = CIntArray(nbAssets);
            for (int i = 0; i < nbAssets; i++){
                assetsIndexes[i] = i;
            }
        }
        /** Main method. Returns refLevel given values for future path of
            specified asset */
        virtual double refLevel(int            iAsset,
                                const double*  futurePath) const{
            return 0.0; // by defn
        }

        /** Used for estimating refLevel when future path depends upon
            refLevel (ie vol interp level for LN) */
        virtual double refLevel(int                  iAsset,
                                const IMultiFactors* mAsset) const{
            return 0.0; // by defn
        }

        /** Returns an array of integers denoting which assets have a
            simulation date given the index into getAllDates(). The
            number of the assets (0, 1, 2, ..., numAssets-1) is the
            same as the order in which the series is supplied to the
            constructor */
        const CIntArray& assetsOnDate(int index) const{
            return assetsIndexes; // independent of index
        }

        //// utility method - returns the number of future points
        virtual int numFutureDates(int iAsset) const{
            return numFutureValues;
        }

    };

    IRefLevel::IMCPath* createMCPath(const DateTime&    valueDate,
                                     const DoubleArray& spotsAtSimStart,
                                     const IPastValues* pastValues) const{
        return new MCPath(pastValues->getNumAssets(),
                          valueDate.numFutureDates(allDates));
    }

    /** This method is to be retired once we move to 'state variables' for the
        Monte Carlo. This is a hack to get around something that was designed in
        namely the fact that the ref level is responsible for choosing the
        sim start date */
    virtual void setSimStartDateToFirstRefDate(){
        // nothing to do
    }

    /** Creates a IRefLevel where the level is always zero. */
    Zero(const DateTime& simStartDate):
        CObject(TYPE), simStartDate(simStartDate){
        validatePop2Object();
    }

    virtual void validatePop2Object(){
        // build transient array of all dates
        allDates = DateTimeArray(1, simStartDate);
    }

    /**  State variable for this RefLevel. This is
         the type that products deal with in the payoff. */
    class StateVar: public virtual IStateVar{
        bool isDoingPast;
    public:
        virtual ~StateVar(){}
        StateVar(bool isDoingPast): isDoingPast(isDoingPast){}

        virtual bool doingPast() const{
            return isDoingPast;
        }

        virtual IHistoricalContextSP createHistoricalContext()
        {
            return IHistoricalContextSP();
        }


        virtual IHistoricalContextSP getInitialHistoricalContext()
        {
            // xxx, revisit
            return IHistoricalContextSP();
        }

        virtual double refLevel2(int iAsset, IHistoricalContextSP history)
        {
            return (0);
        }


        /** Returns the reference level for given asset. For
            the PastPathGenerator then any future average values will
            be 'estimated' (current methodology is current spot). For
            a future path generator, then this must return 'estimated'
            values too (eg same values as the Past generator) until
            generatePath is invoked. At which point it will contain
            the correct value for the current path */
        virtual double refLevel(int iAsset) {
            return 0.0; // that's why it's a Zero RefLevel
        }
    };

    /** A Generator of MC RefLevel Variables */
    class StateVarGen: public virtual IRefLevel::IStateVarGen{
        DateTimeArray allDates;
    public:
        StateVarGen(const DateTime& valueDate):allDates(1, valueDate) {}

        virtual ~StateVarGen(){}

        /** Create the corresponding State Variable for this State
            Variable Generator (NB Implies one state variable per
            generator). The previous IStateVariableSP (may be null)
            should be passed in.  The return object may or may not be
            the same as oldStateVar.  */
        virtual IStateVariableSP create(IStateVariableSP             oldStateVar,
                                     IStateVariableGen::IStateGen* pathGen) const{
            IStateVarSP myStateVar(&dynamic_cast<IStateVar&>(*oldStateVar));
            return getRefLevelSV(myStateVar, pathGen);
        }

        /** Returns a RefLevel state variable which then provides
            access to the reference level. This is the method that
            products should call to get an IRefLevel::IStateVar. The
            previous IRefLevel::IStateVarSP (may be null) should be
            passed in. The return object may or may not be the same as
            oldStateVar.*/
        virtual IStateVarSP getRefLevelSV(IStateVarSP               oldStateVar,
                                          IStateVariableGen::IStateGen* pathGen) const{
            return IStateVarSP(new StateVar(pathGen->doingPast()));
        }

        /** Appends 'true' (ie non derived) state variable generators
            required to the supplied collector. */
        virtual void collectStateVars(
            IStateVariableCollectorSP svCollector) const{
            // no dependency
        }

        /** Returns the dates used for computing the reference level
            of iAsset */
        const DateTimeArray& getDates(int iAsset) const {
            return allDates;
        }

        /** Returns the dates used for computing the reference level
            of iAsset */
        const DateTimeArray& getAllDates() const {
            return allDates;
        }
    };

    /** How to get a IStateVariableGen from a RefLevel */
    virtual IStateVarGen* createStateVarGen(
        const IMultiFactors* multiFactors,
        const DateTime&      valueDate) const{
        return new StateVarGen(valueDate);
    }

private:

    /* for reflection */
    Zero(): CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(Zero, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRefLevel);
        EMPTY_SHELL_METHOD(defaultZero);
        FIELD(simStartDate, "Simulation Start Date");
        FIELD(allDates, "Internal");
        FIELD_MAKE_TRANSIENT(allDates);
    }

    static IObject* defaultZero(){
        return new Zero();
    }
};

CClassConstSP const IRefLevel::Zero::TYPE =
CClass::registerClassLoadMethod(
    "IRefLevel::Zero", typeid(IRefLevel::Zero), load);

/** creates a reference level which is always zero for all assets -
    essentially a dummy one to use if you want to handle
    the reference level yourself. It does not require any
    historic values */
IRefLevel* IRefLevel::Util::makeZero(const DateTime&   simStartDate){
    return new Zero(simStartDate);
}

/*************************************************************************/

/** Class implementing IRefLevel where each asset has the same set of
    N dates, and the level is the average of the best/worst M of the N */
class IRefLevel::AvgMofN: public CObject,
                 virtual public IRefLevel{
public:
    class MCPath;
    class StateVarGen;
    friend class StateVarGen;
private:
    friend class  MCPath;
    DateTimeArray avgDates;
    bool          isBest;
    int           nbKeep;
    bool          startSimAtFirstAvgDate; // to be removed
public:
    static CClassConstSP const TYPE;

    virtual void GetMarket(const IModel* model, const CMarketDataSP market) {

    }

    /** returns all the simulation dates (of all assets) */
    const DateTimeArray& getAllDates() const{
        return avgDates;
    }

    /** returns the date when a model should start simulating - will not
        return a date before today */
    virtual const DateTime& getSimStartDate(const DateTime& today) const{
        if (startSimAtFirstAvgDate){
            // start simulation on start date
            return (avgDates.front().isGreater(today)? avgDates.front(): today);
        }
        // otherwise start simulation immediately - not at first average date
        return today;
    }

    /** returns simulation dates (of all assets) which are strictly in
        the future */
    virtual DateTimeArray getFutureDates(const DateTime& valueDate) const{
        return valueDate.getFutureDates(avgDates);
    }

    /** Returns true if each asset has the same set of simulation dates */
    bool sameDatesPerAsset() const{
        return true;
    }

    /** returns the number of dates (past and future) for given asset */
    virtual int numDates(int iAsset) const{
        return avgDates.size();
    }

    /** returns the number of future dates for given asset */
    virtual int numFutureDates(const DateTime& valueDate,
                               int             iAsset) const{
        return valueDate.numFutureDates(avgDates);
    }

    class MCPath: public IRefLevel::IMCPath{
    private:
        const AvgMofN*      avgRefLevel;
        int                 nbAssets;
        int                 numDates; // for ease (same per asset)
        int                 numFutureValues; // same per asset
        DoubleMatrix        samplesSoFar; // [asset][sample date]
        DoubleArray         refIfPast;  // per asset. If all dates are past then we know the level
        mutable DoubleArray samplesForSort;  // [sample date]. Working area to avoid alloc in refLevel()
        IntArray            assetsIndexes;

        double calcAvgMofN(DoubleArray& samples) const {
            int    startIdx;
            int    endIdx;

            /* sorts best to worst so read appropriately */
            Algorithm::shellSort(samples);
            if (avgRefLevel->isBest)
            {
                startIdx = 0;
                endIdx = avgRefLevel->nbKeep;
            } else {
                startIdx = samples.size()-avgRefLevel->nbKeep;
                endIdx = samples.size();
            }
            double level = 0.0;
            for(int j=startIdx; j<endIdx; j++)
            {
                level += samples[j];
            }
            level /= avgRefLevel->nbKeep;
            return level;
        };

    public:
        // constructor for when we start simulation at first avg date and
        // value date is before first avg date
        MCPath(const AvgMofN*     avgRefLevel,
               const DoubleArray& spotsAtSimStart):
            avgRefLevel(avgRefLevel),
            nbAssets(spotsAtSimStart.size()),
            numDates(avgRefLevel->numDates(0)),
            numFutureValues(numDates - 1), // value at first date is known
            samplesSoFar(nbAssets, numDates),
            refIfPast(nbAssets),
            samplesForSort(numDates) {
            // build up [trivial] assetsIndexes
            assetsIndexes = CIntArray(nbAssets);
            for (int i = 0; i < nbAssets; i++){
                assetsIndexes[i] = i;
            }

            if (numFutureValues>0) {
                // can only record samples so far
                for (int a = 0; a < nbAssets; a++){
                    samplesSoFar[a][0] = spotsAtSimStart[a];
                }
            } else {
                // This will only happen if there is only one averaging date
                // Can calculate actual level
                for (int a = 0; a < nbAssets; a++){
                    refIfPast =spotsAtSimStart;
                }
            }
        }


        // constructor for when we start simulation today
        MCPath(const AvgMofN*     avgRefLevel,
               const IPastValues* pastValues,
               const DateTime&    valueDate):
            avgRefLevel(avgRefLevel),
            nbAssets(pastValues->getNumAssets()),
            numDates(avgRefLevel->numDates(0)),
            numFutureValues(avgRefLevel->numFutureDates(valueDate, 0)),
            samplesSoFar(nbAssets, numDates),
            refIfPast(nbAssets),
            samplesForSort(numDates) {
            // build up [trivial] assetsIndexes
            assetsIndexes = CIntArray(nbAssets);
            for (int i = 0; i < nbAssets; i++){
                assetsIndexes[i] = i;
            }

            if (numFutureValues>0) {
                // can only record samples so far
                for (int a = 0; a < nbAssets; a++){
                    DoubleArray past =
                        pastValues->getPastValues(avgRefLevel->avgDates,
                                                  a,
                                                  valueDate);
                    if (past.size()+numFutureValues != numDates) {
                        throw ModelException("IRefLevel::AvgMofN::MCPath",
                                             "Internal error : nb past (" + Format::toString(past.size()) +
                                             ") + numFutureValues (" + Format::toString(numFutureValues) +
                                             ") does not equal numDates (" + Format::toString(numDates) + ")!");
                    }
                    for (int i = 0; i < past.size(); i++){
                        samplesSoFar[a][i] = past[i];
                    }
                }
            } else {
                // can calculate actual level
                for (int a = 0; a < nbAssets; a++){
                    DoubleArray past =
                        pastValues->getPastValues(avgRefLevel->avgDates,
                                                  a,
                                                  valueDate);
                    refIfPast[a] = calcAvgMofN(past);
                }
            }
        }
        /** Main method. Returns refLevel given values for future path of
            specified asset */
        virtual double refLevel(int            iAsset,
                                const double*  futurePath) const{
            if (numFutureValues > 0){
                int    i;
                int    numPast = numDates - numFutureValues;
                for (i = 0; i < numPast; i++){
                    samplesForSort[i] = samplesSoFar[iAsset][i];
                }
                for (i = 0; i < numFutureValues; i++){
                    samplesForSort[numPast+i] = futurePath[i];
                }
                return calcAvgMofN(samplesForSort);
            }
            return refIfPast[iAsset];
        }
        /** Used for estimating refLevel when future path depends upon
            refLevel (ie vol interp level for LN) */
        virtual double refLevel(int                  iAsset,
                                const IMultiFactors* mAsset) const{
            if (numFutureValues > 0){
                // extrapolate just using spot
                int    i;
                int    numPast = numDates - numFutureValues;
                double spot = mAsset->assetGetSpot(iAsset);

                for (i = 0; i < numPast; i++){
                    samplesForSort[i] = samplesSoFar[iAsset][i];
                }
                for (i = 0; i < numFutureValues; i++){
                    samplesForSort[numPast+i] = spot;
                }
                return calcAvgMofN(samplesForSort);
            }
            return refIfPast[iAsset];
        }

        /** Returns an array of integers denoting which assets have a
            simulation date given the index into getAllDates(). The
            number of the assets (0, 1, 2, ..., numAssets-1) is the
            same as the order in which the series is supplied to the
            constructor */
        const CIntArray& assetsOnDate(int index) const{
            return assetsIndexes; // independent of index
        }

        //// utility method - returns the number of future points
        virtual int numFutureDates(int iAsset) const{
            return numFutureValues;
        }

    };

    IRefLevel::IMCPath* createMCPath(const DateTime&    valueDate,
                                     const DoubleArray& spotsAtSimStart,
                                     const IPastValues* pastValues) const{
        if (startSimAtFirstAvgDate && valueDate.isLess(avgDates.front())){
            return new MCPath(this, spotsAtSimStart);
        }
        // otherwise since our simStartDate starts immediately we don't need
        // to use spotsAtSimStart
        return new MCPath(this, pastValues, valueDate);
    }

    /** This method is to be retired once we move to 'state variables' for the
        Monte Carlo. This is a hack to get around something that was designed in
        namely the fact that the ref level is responsible for choosing the
        sim start date */
    virtual void setSimStartDateToFirstRefDate(){
        startSimAtFirstAvgDate = true;
    }

    /** Creates a IRefLevel where each asset has the same set of simulation
        dates. */
    AvgMofN(const DateTimeArray& avgDates,
            bool                 isBest,
            int                  nbKeep):
        CObject(TYPE), avgDates(avgDates),
        isBest(isBest), nbKeep(nbKeep), startSimAtFirstAvgDate(false) {
        validatePop2Object();
    }

    virtual void validatePop2Object(){
        DateTime::ensureIncreasing(avgDates, "average dates", true);
        if (nbKeep<1 || nbKeep>avgDates.size()) {
            throw ModelException("IRefLevel::AvgMofN::validatePop2Object",
                                 "nbKeep " + Format::toString(nbKeep) +
                                 " must be >0 and <= number of dates, " +
                                 Format::toString(avgDates.size()));
        }
    }

    /**  State variable for this RefLevel. This is
         the type that products deal with in the payoff. */
    class StateVar: public virtual IStateVar{
        SVGenSpotSP            mcSpot;          //!< Spot path generator
        SVGenSpot::IStateVarSP mcSpotSV;        //!< Spot path state variable
        DoubleArray         refIfPast;       /* per asset. If all dates are past
                                                then we know the level */
        bool                isBest;          // either best M of N or worst M of N
        int                 nbKeep;          // the value for M
        vector<bool>        allDone;         // by asset
        int                 numDates;
        vector<DoubleArray> samplesSoFar;    // [asset][sample date]
        mutable DoubleArray samplesForSort;  /* [sample date]. Working area to
                                                avoid alloc in refLevel() */
        vector<double>      spotsAtSimStart;
        bool                isDoingPast;     //!< Whether we are doing past or future

        double calcAvgMofN(DoubleArray& samples) const {
            int    startIdx;
            int    endIdx;

            /* sorts best to worst so read appropriately */
            Algorithm::shellSort(samples);
            if (isBest) {
                startIdx = 0;
                endIdx = nbKeep;
            } else {
                startIdx = samples.size() - nbKeep;
                endIdx = samples.size();
            }
            double level = 0.0;
            for(int j = startIdx; j < endIdx; j++) {
                level += samples[j];
            }
            level /= nbKeep;
            return level;
        }

    public:
        virtual ~StateVar(){}
        virtual bool doingPast() const{
            return isDoingPast;
        }
        //// update mc spot state variable (going from past to future)
        void update(IStateVariableGen::IStateGen* pathGen){
            mcSpotSV = mcSpot->getSpotSV(pathGen);
            isDoingPast = mcSpotSV->doingPast();

            // Compute reflevels in case past won't be invoked
            for(unsigned int iAsset = 0; iAsset < allDone.size(); iAsset++) {
                refLevel(iAsset);
            }
        }

        virtual IHistoricalContextSP createHistoricalContext()
        {
            return IHistoricalContextSP();
        }


        virtual IHistoricalContextSP getInitialHistoricalContext()
        {
            // xxx, revisit
            return IHistoricalContextSP();
        }

        virtual double refLevel2(int iAsset, IHistoricalContextSP history)
        {
            return (0);
        }


        /** Returns the reference level for given asset. When doing
            the the past any future average values will be 'estimated'
            (current methodology is current spot). For a future path
            generator, then this must return 'estimated' values too
            (eg same values as the Past generator) until generatePath
            is invoked. At which point it will contain the correct
            value for the current path */
        virtual double refLevel(int iAsset) {
            /** Switch off the allDone flag in case someone is
                computing the past multiple times. */
            if(isDoingPast) {
                allDone[iAsset] = false;
            };

            if (allDone[iAsset]) {
                // known
                return refIfPast[iAsset];
            }
            // else get SVPath
            const SVPath& path = mcSpotSV->path(iAsset);
            int begin = path.begin();
            int end   = path.end();

            // Copy over past values in preservation area
            int iPath;
            if(isDoingPast) {
                samplesSoFar[iAsset] = DoubleArray(0);
                for (iPath = path.begin(); iPath < end; iPath++) {
                    samplesSoFar[iAsset].push_back(path[iPath]);
                }
            }

            int iStep = 0;
            if(!isDoingPast) {
                // When doing future, get past values from preservation area
                for (; iStep < samplesSoFar[iAsset].size(); iStep++){
                    samplesForSort[iStep] = samplesSoFar[iAsset][iStep];
                }
            }

            // Add path from current state variable
            for (iPath = begin; iPath < end; iPath++, iStep++){
                samplesForSort[iStep] = path[iPath];
            }

            // See if the path is completed
            bool isCompleted = iStep == numDates;

            // estimate the remaining future levels
            for(; iStep < samplesForSort.size(); iStep++) {
                samplesForSort[iStep] = spotsAtSimStart[iAsset];
            }

            // Compute refLevel and cache it
            double refLevel = calcAvgMofN(samplesForSort);
            if (isDoingPast && isCompleted){
                allDone[iAsset] = true; // cache this number
                refIfPast[iAsset] = refLevel;
            }
            return refLevel;
        }

        /** constructor. numDates is total num of dates or zero
            if all in past */
        StateVar(bool                      isBest, //either best or worst of M of N
                 int                       nbKeep, // the value for M
                 const vector<double>&     spotsAtSimStart,
                 SVGenSpotSP                  mcSpot,
                 IStateVariableGen::IStateGen* pathGen):
        mcSpot(mcSpot),
        refIfPast(spotsAtSimStart.size(),0.0),
        isBest(isBest), nbKeep(nbKeep), allDone(spotsAtSimStart.size()),
        numDates(mcSpot->numDates()), samplesSoFar(spotsAtSimStart.size()),
        samplesForSort(numDates),
        spotsAtSimStart(spotsAtSimStart) {
            for(unsigned int iAsset = 0; iAsset < spotsAtSimStart.size(); iAsset++) {
                samplesSoFar[iAsset] = DoubleArray(0);
            }
            update(pathGen);
        }
    };

    /** A Generator of MC RefLevel Variables */
    class StateVarGen: public virtual IRefLevel::IStateVarGen{
        const AvgMofN*  avgRefLevel;
        SVGenSpotSP        mcSpot; // this is our StateVar generator
        vector<double>  spotsAtSimStart;
    public:
        virtual ~StateVarGen(){}

        StateVarGen(const AvgMofN* avgRefLevel,
                    SVGenSpotSP mcSpot, const IMultiFactors* multiFactors):
            avgRefLevel(avgRefLevel), mcSpot(mcSpot),
            spotsAtSimStart(multiFactors->NbAssets()){
            for (unsigned int iAsset = 0;
                 iAsset < spotsAtSimStart.size(); iAsset++){
                // simulation starts today
                spotsAtSimStart[iAsset] = multiFactors->assetGetSpot(iAsset);
            }
        }

        /** Create the corresponding State Variable for this State
            Variable Generator (NB Implies one state variable per
            generator). The previous IStateVariableSP (may be null)
            should be passed in.  The return object may or may not be
            the same as oldStateVar.  */
        virtual IStateVariableSP create(IStateVariableSP             oldStateVar,
                                     IStateVariableGen::IStateGen* pathGen) const{
            IStateVarSP myStateVar(&dynamic_cast<IStateVar&>(*oldStateVar));
            return getRefLevelSV(myStateVar, pathGen);
        }

        /** Returns a RefLevel state variable which then provides
            access to the reference level. This is the method that
            products should call to get an IRefLevel::IStateVar. The
            previous IRefLevel::IStateVarSP (may be null) should be
            passed in. The return object may or may not be the same as
            oldStateVar.*/
        virtual IStateVarSP getRefLevelSV(IStateVarSP               oldStateVar,
                                          IStateVariableGen::IStateGen* pathGen) const{
            if (oldStateVar.get()){
                // just update
                dynamic_cast<StateVar&>(*oldStateVar).update(pathGen);
                return oldStateVar;
            }
            // else create StateVar
            return IStateVarSP(new StateVar(avgRefLevel->isBest,
                                            avgRefLevel->nbKeep,
                                            spotsAtSimStart,
                                            mcSpot,
                                            pathGen));
        }

        /** Appends 'true' (ie non derived) state variable generators
            required to the supplied collector. */
        virtual void collectStateVars(
            IStateVariableCollectorSP svCollector) const{
            svCollector->append(mcSpot.get());
        }

        /** Returns the dates used for computing the reference level
            of iAsset */
        const DateTimeArray& getDates(int iAsset) const {
            return mcSpot->getSimSeries()->getDates(iAsset);
        }

        /** Returns the dates used for computing the reference level
            of iAsset */
        const DateTimeArray& getAllDates() const {
            return mcSpot->getSimSeries()->getAllDates();
        }
    };

    /** How to get a IStateVariableGen from a RefLevel */
    virtual IStateVarGen* createStateVarGen(
        const IMultiFactors* multiFactors,
        const DateTime&      valueDate) const{
        SVGenSpotSP    mcSpot(new SVGenSpot(multiFactors->NbAssets(), avgDates));
        return new StateVarGen(this, mcSpot, multiFactors);
    }

private:

    /* for reflection */
    AvgMofN(): CObject(TYPE), startSimAtFirstAvgDate(false){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(AvgMofN, clazz);
        SUPERCLASS(CObject);
        IMPLEMENTS(IRefLevel);
        EMPTY_SHELL_METHOD(defaultCon);
        FIELD(avgDates, "Averaging in dates");
        FIELD(nbKeep, "The M in M of N");
        FIELD(isBest, "Else worst");
        FIELD_NO_DESC(startSimAtFirstAvgDate);
        FIELD_MAKE_TRANSIENT(startSimAtFirstAvgDate);
    }

    static IObject* defaultCon(){
        return new AvgMofN();
    }
};

CClassConstSP const IRefLevel::AvgMofN::TYPE = CClass::registerClassLoadMethod(
    "IRefLevel::AvgMofN", typeid(IRefLevel::AvgMofN), load);


/** Creates a RefLevel where each asset has the same set of simulation
    dates. */
IRefLevel* IRefLevel::Util::makeAvgMofN(
    const DateTimeArray&  commonSimDates,
    bool                  isBest,
    int                   nbKeep) {
    return new AvgMofN(commonSimDates,
                       isBest, nbKeep);
}

/*************************************************************************/

void IRefLevel::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER_INTERFACE(IRefLevel, clazz);
    EXTENDS(IObject);
}


CClassConstSP const IRefLevel::TYPE = CClass::registerInterfaceLoadMethod(
    "IRefLevel", typeid(IRefLevel), load);



DRLIB_END_NAMESPACE
