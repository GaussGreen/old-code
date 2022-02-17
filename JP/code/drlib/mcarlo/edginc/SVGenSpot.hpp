//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SVGenSpot.hpp
//
//   Description : A Generator of MC Spot State Variables
//
//   Date        : March 2004
//
//
//----------------------------------------------------------------------------


#ifndef EDR_MCSPOT_HPP
#define EDR_MCSPOT_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/ElemStateVariableGen.hpp"
#include "edginc/SimSeries.hpp"
#include "edginc/SVPath.hpp"

DRLIB_BEGIN_NAMESPACE
class MCPath;
typedef smartPtr<MCPath> MCPathSP;
typedef vector<const MCPath*> MCPathArray;
typedef smartPtr<MCPathArray> MCPathArraySP;


/** A Generator of MC Spot State Variables */
class MCARLO_DLL MCPath : public virtual VirtualDestructorBase {
public:

    /** Interface for the state variable that MCPath produces. This is
        the type that products deal with in the payoff. The payoff obtains
        it by calling the getSpotSV() method below */
    class MCARLO_DLL IStateVar: public virtual IAdvanceableStateVariable{
    public:
        IStateVar();
        virtual ~IStateVar();

        /** Returns the number of assets */
        virtual int numAssets() const = 0;

        /** returns the path for the given asset. Use
            begin() and end() on SVPath to see which values are being
            simulated. Values before begin() are valid and will be
            historic values (Note the notion of multiple paths for one asset
            within a single simulated path has been dropped) */
        virtual const SVPath& path(int iAsset) const = 0;

        /** Returns an estimate of the largest value of the 'drift' we are
            likely to see between today and the last simulation
            date. Here, drift, means the ratio of the spot at the last
            simulation date and the current spot. (The drift depends upon
            market data (eg divs) and the randoms uses to simulate the
            path - we are trying to capture some number, x, so that we can
            say for 99% (say) of paths, the spot at the last date will be
            less than x * current spot. Used by quick greeks to minimise
            delta calculation */
        virtual double maxDriftEstimate(int iAsset) const = 0;

        // for stateless payoff
        virtual double getSpotPrice(int iAsset) const = 0;

        virtual void getSpotPrices(DoubleArray& prices) const = 0;
    };
    typedef smartPtr<IStateVar> IStateVarSP;
    typedef vector<IStateVarSP> IStateVarArray;

    virtual ~MCPath();

    /** Constructor - takes reference to supplied simSeries */
    MCPath(SimSeriesConstSP simSeries);

    /** Constructor - from numAssets and an array of dates (same dates per
        asset) */
    MCPath(int numAssets, const DateTimeArray& dates);

    /** Constructor - from numAssets and a single date */
    MCPath(int numAssets, const DateTime& date);

    /** Returns the sim series for this MCPath */
    SimSeriesConstSP getSimSeries() const;

    /** Returns the number of assets */
    int numAssets() const;
    /** Returns the number of dates (including all assets) */
    int numDates() const;
    /** Returns the number of dates for the specified asset */
    int numDates(int iAsset) const;

    /** Returns the union of all dates for a particular asset from all the MCPaths
        provided */
    template<class T> static
    DateTimeArraySP getAllDates(const T& mcPaths, int assetIdx) {
        vector<const DateTimeArray*> theDates(mcPaths.size());
        for (unsigned int i = 0; i < theDates.size(); i++){
            theDates[i] = &mcPaths[i]->simSeries->getDates(assetIdx);
        }
        return DateTimeArraySP(new DateTimeArray(DateTime::merge(theDates)));
    }

    /** Returns the union of all dates for all assets from all the MCPaths
        provided */
    template<class T> static
    DateTimeArraySP getAllDates(const T& mcPaths) {
        vector<const DateTimeArray*> theDates(mcPaths.size());
        for (unsigned int i = 0; i < theDates.size(); i++){
            theDates[i] = &mcPaths[i]->simSeries->getAllDates();
        }
        return DateTimeArraySP(new DateTimeArray(DateTime::merge(theDates)));
    }

    /** Creates a vector of IStateVars, with one IStateVar for each
        MCPath.  The allDatesPerAsset parameter must be the same sets
        of dates as returned by getAllDates. The beginIdx and endIdx
        are indexes into allDates and represent what the corresponding
        begin() and end() methods on each MCPath should return. The
        paths should be of size equal number of assets and contains
        the simulated path of each asset at 'allDates'.  The
        maxDriftEstimate contains a double for each asset which
        represents the value to be used for the maxDriftEstimate()
        method on MCPath::IStateVar */
    template<class T> static
    IStateVarArray createPaths(
        bool                         doingPast, // simulating the past?
        const T&                     mcPaths,
        const DateTimeArrayArray&    allDatesPerAsset,
        const vector<int>&           beginIdx,
        const vector<int>&           endIdx,
        const vector<const double*>& paths,
        const vector<double>&        maxDriftEstimates) {

        static const string method("MCPath::createPaths");
        try{
            // build empty vector
            vector<IStateVarSP> pathsByAsset(mcPaths.size());
            // then for each MCPath
            for (unsigned int i = 0; i < mcPaths.size(); i++){
                const MCPath& mcPath = *(mcPaths[i]);
                int numAssets = mcPath.numAssets();
                PathByAssetSP pathByAsset(new PathByAsset(numAssets,
                                                          maxDriftEstimates,
                                                          doingPast));
                pathsByAsset[i] = pathByAsset;
                // further optimisation remains to be done eg
                //bool sameDatesPerAsset = MCPath.simSeries->sameDatesPerAsset();
                // also optimise for case where the dates are the same in more
                // than one MCPath.
                // Then for each asset
                for (int j = 0; j < numAssets; j++){
                    const DateTimeArray& assetDates =
                        mcPath.simSeries->getDates(j); // get dates for this asset
                    if(!assetDates.size()) {
                        // No dates requested so create a dummy path with begin = end
                        pathByAsset->paths[j] =
                            SVPathSP(new SVPath(paths[j], 0, 0));
                    } else if (assetDates.size() == allDatesPerAsset[j].size()){
                        // easy!
                        pathByAsset->paths[j] =
                            SVPathSP(new SVPath(paths[j], beginIdx[j], endIdx[j]));
                    } else {
                        bool isTrivial; // this case will have been caught above
                        // This call creates a map which allows us to jump
                        // along the asset's dates
                        IntArray dateMap(
                            DateTime::createMapping(allDatesPerAsset[j],
                                                    assetDates,
                                                    isTrivial));
                        vector<int> offsets(assetDates.size()); /* what we'll pass
                                                                   to Path() */
                        const double* thePath = paths[j];

                        thePath += dateMap.front(); /* allows us to use trivial
                                                       mode on Path, perhaps */
                        // thisBegin will be used for what begin() will return
                        int thisBegin = assetDates.size();/* default value
                                                             (eg dates in past) */
                        bool setThisBegin = false; // not been set yet
                        int thisEnd = assetDates.size();//what end() will return
                        // loop over mapping so that we land on each asset date
                        int offsetIdx = 0; // where we are in offsets
                        for(int iStep = dateMap.front();
                            iStep < allDatesPerAsset[j].size();
                            offsetIdx++, iStep++, iStep += dateMap[iStep]){
                            offsets[offsetIdx] = iStep - dateMap.front();
                            if (iStep >= beginIdx[j] && !setThisBegin){
                                // we've crossed the start of the sim. Save idx
                                thisBegin = offsetIdx;
                                setThisBegin = true;
                            }
                            if (iStep >= endIdx[j] && thisEnd == assetDates.size()){
                                // we've crossed the end of the sim. Save idx
                                thisEnd = offsetIdx;
                            }
                        }
                        // typically this will be true most of the time
                        bool contiguous = offsets.back() == assetDates.size() - 1;
                        if (contiguous){
                            pathByAsset->paths[j] =
                                SVPathSP(new SVPath(thePath, thisBegin, thisEnd));
                        } else {
                            pathByAsset->paths[j] =
                                SVPathSP(new SVPath(thePath, offsets,
                                                thisBegin, thisEnd));
                        }
                    }
                }
            }
            return pathsByAsset;
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    /** Same as above for the case where allDates includes all dates from
        all assets */
    template<class T> static
    IStateVarArray createPaths(
        bool                         doingPast, // simulating the past?
        const T&                     mcPaths,
        DateTimeArraySP              allDates,
        int                          beginIdx,
        int                          endIdx,
        const vector<const double*>& paths,
        const vector<double>&        maxDriftEstimates) {


        int numAssets = paths.size();
        DateTimeArrayArray allDatesPerAsset(numAssets, *allDates);
        vector<int> beginIndices(numAssets, beginIdx);
        vector<int> endIndices(numAssets, endIdx);
        return createPaths<T>(doingPast, mcPaths, allDatesPerAsset, beginIndices,
                              endIndices, paths, maxDriftEstimates);
    }

    /** Essentially a wrapper for a vector of Paths but with a bit more to
        allow this to implement MCPath::IStateVar */
    class MCARLO_DLL PathByAsset: public virtual IStateVar{
        friend class MCPath;
        vector<SVPathSP> paths;
        vector<double> maxDrifts;
        bool           isDoingPast;
    public:
        ~PathByAsset() {}

        PathByAsset(const vector<SVPathSP>&  paths,
                    const vector<double>&  maxDrifts,
                    bool                   isDoingPast);

        /** Returns the number of assets */
        virtual int numAssets() const;

        /** Returns true if this state variable is being used to 'simulate'
            the past. This is a useful method for users of state variables -
            need to see how hard it is to implement */
        virtual bool doingPast() const;

        /** Returns the path for the specified asset. Note that the Path will
            go out of scope when this object is freed */
        const SVPath& path(int iAsset) const;

        SVPathSP getSV(int iAsset);

        /** Returns an estimate of the largest value of the 'drift' we are
            likely to see between today and the last simulation
            date. Here, drift, means the ratio of the spot at the last
            simulation date and the current spot. (The drift depends upon
            market data (eg divs) and the randoms uses to simulate the
            path - we are trying to capture some number, x, so that we can
            say for 99% (say) of paths, the spot at the last date will be
            less than x * current spot. Used by quick greeks to minimise
            delta calculation */
        virtual double maxDriftEstimate(int iAsset) const;

        virtual double getSpotPrice(int iAsset) const;

        virtual void getSpotPrices(DoubleArray& prices) const;

        virtual void advance();

        virtual void reset();

        virtual void prepare(bool mm)
        {
            for(size_t i=0; i<paths.size(); ++i)
            {
                if (paths[i].get())
                    paths[i]->prepare(mm);
            }
        }

    private:
        PathByAsset(int numAssets, const vector<double>& maxDrifts,
                    bool isDoingPast);
    };

    typedef smartPtr<PathByAsset> PathByAssetSP;

protected:
    ///// fields  ////////
    SimSeriesConstSP     simSeries;   // dates to simulate on

private:
};

#ifndef QLIB_MCSPOT_CPP
    EXTERN_TEMPLATE(class MCARLO_DLL_SP smartConstPtr<MCPath::IStateVar>);
    EXTERN_TEMPLATE(class MCARLO_DLL_SP smartPtr<MCPath::IStateVar>);
#else
    INSTANTIATE_TEMPLATE(class MCARLO_DLL smartConstPtr<MCPath::IStateVar>);
    INSTANTIATE_TEMPLATE(class MCARLO_DLL smartPtr<MCPath::IStateVar>);
#endif

class MCARLO_DLL SVGenSpot: virtual public MCPath,
              virtual public IElemStateVariableGen {
public:
    /** Constructor - takes reference to supplied simSeries */
    SVGenSpot(SimSeriesConstSP simSeries);

    /** Constructor - from numAssets and an array of dates (same dates per
        asset) */
    SVGenSpot(int numAssets, const DateTimeArray& dates);

    /** Constructor - from numAssets and a single date */
    SVGenSpot(int numAssets, const DateTime& date);

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.  The
        return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP             oldStateVar,
                                    IStateVariableGen::IStateGen* pathGen) const;

    /** Returns a MC Spot state variable which then provides access to
        the path etc. This is the method that products should call to
        get an MCPath::IStateVar. */
    IStateVarSP getSpotSV(IStateVariableGen::IStateGen* pathGen) const;

    void attachSVGen(IElemStateVariableGenVisitor*) const;
};

typedef smartPtr<SVGenSpot> SVGenSpotSP;
typedef vector<const SVGenSpot*> SVGenSpotArray;
typedef smartPtr<SVGenSpotArray> SVGenSpotArraySP;

class MCARLO_DLL MCQuadVar: virtual public MCPath,
                 virtual public IElemStateVariableGen {
public:
    /** Constructor - takes reference to supplied simSeries */
    MCQuadVar(SimSeriesConstSP simSeries);

    /** Constructor - from numAssets and an array of dates (same dates per
        asset) */
    MCQuadVar(int numAssets, const DateTimeArray& dates);

    /** Constructor - from numAssets and a single date */
    MCQuadVar(int numAssets, const DateTime& date);

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.  The
        return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP             oldStateVar,
                                    IStateVariableGen::IStateGen* pathGen) const;

    /** Returns a MC Spot state variable which then provides access to
        the path etc. This is the method that products should call to
        get an MCPath::IStateVar. */
    IStateVarSP getQuadVarSV(IStateVariableGen::IStateGen* pathGen) const;
    void attachSVGen(IElemStateVariableGenVisitor*) const;

};

typedef smartPtr<MCQuadVar> MCQuadVarSP;
typedef vector<const MCQuadVar*> MCQuadVarArray;
typedef smartPtr<MCQuadVarArray> MCQuadVarArraySP;

//////////////////////////////////////////////////////////////////////////

class MCARLO_DLL MCSqrtAnnualQuadVar: virtual public MCPath,
                                      virtual public IElemStateVariableGen {
public:
    /** Constructor - takes reference to supplied simSeries */
    MCSqrtAnnualQuadVar(SimSeriesConstSP simSeries);

    /** Constructor - from numAssets and an array of dates (same dates per
        asset) */
    MCSqrtAnnualQuadVar(int numAssets, const DateTimeArray& dates);

    /** Constructor - from numAssets and a single date */
    MCSqrtAnnualQuadVar(int numAssets, const DateTime& date);

    /** Create the corresponding State Variable for this State
        Variable Generator (from IStateVariableGen interface). The
        previous IStateVariableSP (may be null) should be passed in.  The
        return object may or may not be the same as oldStateVar. */
    virtual IStateVariableSP create(IStateVariableSP             oldStateVar,
                                    IStateVariableGen::IStateGen* pathGen) const;

    /** Returns a MC Spot state variable which then provides access to
        the path etc. This is the method that products should call to
        get an MCPath::IStateVar. */
    IStateVarSP getSqrtAnnualQuadVarSV(IStateVariableGen::IStateGen* pathGen) const;
    void attachSVGen(IElemStateVariableGenVisitor*) const;
};

typedef smartPtr<MCSqrtAnnualQuadVar> MCSqrtAnnualQuadVarSP;
typedef vector<const MCSqrtAnnualQuadVar*> MCSqrtAnnualQuadVarArray;
typedef smartPtr<MCSqrtAnnualQuadVarArray> MCSqrtAnnualQuadVarArraySP;

DRLIB_END_NAMESPACE

#endif
