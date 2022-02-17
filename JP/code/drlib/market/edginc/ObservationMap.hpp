//----------------------------------------------------------------------------
//
//   Filename    : ObservationMap.hpp
//
//   Description : Class to hold information about the actual observation dares
//                 i.e. when sample are actually observed and how these map to
//                 the sampling dates i.e. when sampling is attempted
//
//   Author      : Ian Stares   
//
//   Date        : July 3 2006
//

//
//----------------------------------------------------------------------------

#ifndef OBSMAP_HPP
#define OBSMAP_HPP

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/MultiFactors.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/SamplingConvention.hpp"
#include "edginc/ObservationSource.hpp"
#include "edginc/ObservationType.hpp"

DRLIB_BEGIN_NAMESPACE

/** for each asset this contains the details for each schedule
    note there may be several schedules being used per asset

    note by sample date we mean the date we attempt to retrieve a sample
    by observation date we mean the date we actually observe a sample
    Be careful though, if the ISDA rule is a move schedule rule the sample
    dates actually move as well with the holidays 
    See the ISDA spec for details */
class MARKET_DLL ObservationMap : public CObject {
public:
    static CClassConstSP const TYPE;

    // one set of sampling dates for all assets - the standard case
    void initialise(const DateTimeArray&          samplingDates,
                    const DateTimeArray&          refDates,
                    const IMultiFactors*          assets,
                    const ObservationSourceArray& sources,
                    const ObservationTypeArray&   obsTypes,
                    const SamplingConvention*     dateAdjust,
                    bool                          hasSamplingConvention,
                    const DateTime&               today);

    /** returns the modelling dates for today. This is a list containing:
        - for samples that have observed on/before today the sampling date
        - for the rest the observation dates at which samples will be known 
        - it does NOT contain the ref level dates*/
    const DateTimeCluster& getModellingDates() const;

    /** returns the modelling dates for today. This is a list containing:
        - for samples that have observed on/before today the sampling date
        - for the rest the observation dates at which samples will be known
        - it includes ref level dates */
    const DateTimeCluster& getAllModellingDates() const;

    /** returns the sample dates (NOT including RefLevel dates)
        note ISDA may cause some sampling dates to be omitted or 
        sampling dates to move with observations*/
    const DateTimeCluster& getSampleDates() const;

    /** what is the index of the modelling date corresponding to the
        given asset and sample date*/
    int getModellingIndex(int assetIdx, int sampleIdx) const;

    /** what is the index of the last sample date corresponding to the
        given asset and modelling date*/
    int getLastSampleIndex(int assetIdx, int simIdx) const;

    /** what is the index of the first sample date corresponding to the
        given asset and modelling date*/
    int getFirstSampleIndex(int assetIdx, int simIdx) const;

    ObservationMap();

private:
    static void load(CClassSP& clazz);
    static IObject* defaultObsMap();

    int                   numAssets;
    DateTimeClusterSP     observationDates;
    DateTimeClusterSP     sampleDates;
    DateTimeClusterSP     modelDates;
    DateTimeClusterSP     allDates;
    IntArrayArraySP       modelDateMap; // modelling date for each sample
    IntArrayArraySP       sampleMap; // last sample date for each modelling date 
};
typedef smartPtr<ObservationMap> ObservationMapSP;

DRLIB_END_NAMESPACE

#endif
