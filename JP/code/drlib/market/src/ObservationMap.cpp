//----------------------------------------------------------------------------
//
//   Filename    : Samples.cpp
//
//   Description : Class to hold all sampling information in one place 
//
//   Author      : Ian Stares   
//
//   Date        : July 3 2006
//

//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/ObservationMap.hpp"

DRLIB_BEGIN_NAMESPACE

ObservationMap::ObservationMap() : CObject(TYPE),
    numAssets(0),
    observationDates(0),
    sampleDates(0),
    modelDates(0),
    allDates(0),
    modelDateMap(0),
    sampleMap(0) {}

void ObservationMap::initialise(const DateTimeArray&          samplingDates,
                                const DateTimeArray&          refDates,
                                const IMultiFactors*          assets,
                                const ObservationSourceArray& sources,
                                const ObservationTypeArray&   obsTypes,
                                const SamplingConvention*     dateAdjust,
                                bool                          hasSamplingConvention,
                                const DateTime&               today) {
    static const string routine("ObservationMap::initialise");
    try {
        //set up the class members to be the right size
        numAssets = assets->NbAssets();
        modelDateMap = IntArrayArraySP(new IntArrayArray(numAssets));

        // only bother doing all the observation stuff if we're adjusting
        if (!dateAdjust->isUnadjusted()) {
            //set up the class members to be the right size
            allDates = DateTimeClusterSP(new DateTimeArrayArray(numAssets));
            observationDates = DateTimeClusterSP(new DateTimeArrayArray(numAssets));
            sampleDates = DateTimeClusterSP(new DateTimeArrayArray(numAssets));
            modelDates = DateTimeClusterSP(new DateTimeArrayArray(numAssets));
            sampleMap = IntArrayArraySP(new IntArrayArray(numAssets));
            DateTimeArray obsDates(numAssets);
            BoolArray sampled(numAssets);
            int assetIdx;
            // do ref dates - note we don't care about mappings and stuff
            // we just want to know which dates we'll need ref level samples for
            for (int dateIdx = 0; dateIdx < refDates.size(); ++dateIdx) {
                assets->observationDates(refDates[dateIdx], sources, dateAdjust,
                                        obsDates, sampled);
                for (assetIdx = 0; assetIdx < numAssets; ++assetIdx) {
                    if (sampled[assetIdx]) { // if date not omitted for this asset
                        if (dateAdjust->scheduleMoves() || obsDates[dateIdx] > today) {
                            ((*allDates)[assetIdx]).append(DateTimeSP(copy(&(obsDates[assetIdx]))));
                        } else {
                            ((*allDates)[assetIdx]).append(DateTimeSP(copy(&(refDates[dateIdx]))));
                        }
                    }
                }
            }
            // go through sampling dates one by one
            for (int dateIdx = 0; dateIdx < samplingDates.size(); ++dateIdx) {
                assets->observationDates(samplingDates[dateIdx], sources, dateAdjust,
                                        obsDates, sampled);
                for (assetIdx = 0; assetIdx < numAssets; ++assetIdx) {
                    if (sampled[assetIdx]) { // if date not omitted for this asset
                        if (dateAdjust->scheduleMoves()) {
                            ((*sampleDates)[assetIdx]).append(DateTimeSP(copy(&(obsDates[assetIdx]))));
                        } else {
                            ((*sampleDates)[assetIdx]).append(DateTimeSP(copy(&(samplingDates[dateIdx]))));
                        }
                        if (dateAdjust->scheduleMoves() || obsDates[assetIdx] > today) {
                            ((*modelDates)[assetIdx]).append(DateTimeSP(copy(&(obsDates[assetIdx]))));
                            ((*allDates)[assetIdx]).append(DateTimeSP(copy(&(obsDates[assetIdx]))));
                        } else {
                            ((*modelDates)[assetIdx]).append(DateTimeSP(copy(&(samplingDates[dateIdx]))));
                            ((*allDates)[assetIdx]).append(DateTimeSP(copy(&(samplingDates[dateIdx]))));
                        }
                        ((*observationDates)[assetIdx]).append(DateTimeSP(copy(&(obsDates[assetIdx]))));
                    }
                }
            }
            // build trimmed down modelDates arrays and the mapping
            // modelDateMap[j][i] gives index of sample/obs date for 
            // ith modelling date for jth asset
            for (assetIdx = 0; assetIdx < numAssets; ++assetIdx) {
                int numObsDates = (*observationDates)[assetIdx].size();
                (*modelDateMap)[assetIdx] = IntArraySP(new IntArray(numObsDates));
                (*sampleMap)[assetIdx] = IntArraySP(new IntArray(numObsDates));
                int modellingIdx = 0;
                int smpIdx = 0;
                (*(*modelDateMap)[assetIdx])[0] = 0;
                for (int i = 1; i < numObsDates; ++i) {
                    if ((*modelDates)[assetIdx][i-1] != (*modelDates)[assetIdx][i]) {
                        (*(*modelDateMap)[assetIdx])[i] = ++modellingIdx;
                        (*(*sampleMap)[assetIdx])[smpIdx++] = i-1;
                    } else {
                        (*(*modelDateMap)[assetIdx])[i] = modellingIdx;
                    }
                }
                // now clear dates and resize array
                DateTime::removeDuplicates((*modelDates)[assetIdx], false);
                int numModelDates = (*modelDates)[assetIdx].size();
                (*sampleMap)[assetIdx]->resize(numModelDates);
                // handle last date which might have been missed
                (*(*sampleMap)[assetIdx])[numModelDates-1] = numObsDates - 1;
           }

        } else {
            int assetIdx;
            DateTimeArray all = DateTime::merge(samplingDates, refDates);
            allDates = 
                DateTimeClusterSP(new DateTimeArrayArray(numAssets, all));
            observationDates = 
                DateTimeClusterSP(new DateTimeArrayArray(numAssets, samplingDates));
            sampleDates = observationDates;
            modelDates = observationDates;
            int numDates = (*observationDates)[0].size();
            for (assetIdx = 0; assetIdx < numAssets; ++assetIdx) {
                (*modelDateMap)[assetIdx] = IntArraySP(new IntArray(numDates));
                for (int i = 0; i < numDates; ++i) {
                    (*(*modelDateMap)[assetIdx])[i] = i;
                }
            }
            sampleMap = modelDateMap;
        }
#if 0
        // DISABLED
        // Setting the times in the dates to match the sampling convention has a few problems...
        // 1) The reflevel time is common to all assets, so if different assets have different
        //    timing rules, we have to be careful about setting up the values, i.e., what
        //    happens if the reflevel is EOD, the sampling convention is SOD and the monitoring
        //    starts on the same day as the reflevel?
        // 2) On the otherhand, parametric date builders are time neutral, and just populate the
        //    time with a default value.
        // Leave this until how this is all handled becomes clearer.
        if( hasSamplingConvention ) {
            // If we use a sampling convention, we must adjust all the times on our sampling
            // dates to fit the convention. This will ensure that we able to retrieve values
            // from the asset history.
            for (int assetIdx = 0; assetIdx < numAssets; ++assetIdx) {
                int time = DateTime::timeConvert(obsTypes[assetIdx]->indicativeTime());
                DateTime::setTimeOfDay((*observationDates)[assetIdx], time);
                DateTime::setTimeOfDay((*sampleDates)[assetIdx], time);
                DateTime::setTimeOfDay((*modelDates)[assetIdx], time);
                DateTime::setTimeOfDay((*allDates)[assetIdx], time);
            }
        }
#endif
 

    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** returns the releant dates for today. This is a list containing:
    - for samples that have observed on/before today the sampling date
    - for the rest the observation dates at which samples will be known*/
const DateTimeCluster& ObservationMap::getModellingDates() const {
    return *modelDates;
}

/** returns the releant dates for today. This is a list containing:
    - for samples that have observed on/before today the sampling date
    - for the rest the observation dates at which samples will be known
    - it includes ref level dates */
const DateTimeCluster& ObservationMap::getAllModellingDates() const {
    return *allDates;
}

/** returns the sample dates (NOT including RefLevel dates)
    note ISDA may cause some sampling dates to be omitted or 
    sampling dates to move with observations*/
const DateTimeCluster& ObservationMap::getSampleDates() const{
    return *sampleDates;
}

/** what is the index of the sim date corresponding to the
    given asset and sample date*/
int ObservationMap::getModellingIndex(int assetIdx, int sampleIdx) const {
    return (*(*modelDateMap)[assetIdx])[sampleIdx];
}

/** what is the index of the last sample date corresponding to the
    given asset and sim date*/
int ObservationMap::getLastSampleIndex(int assetIdx, int modellingIdx) const {
    return (*(*sampleMap)[assetIdx])[modellingIdx];
}

/** what is the index of the first sample date corresponding to the
    given asset and sim date*/
int ObservationMap::getFirstSampleIndex(int assetIdx, int modellingIdx) const {
    return modellingIdx == 0 ? 0 : 
                              (*(*sampleMap)[assetIdx])[modellingIdx-1] + 1;
}

void ObservationMap::load(CClassSP& clazz) {
   // clazz->setPublic(); // not visible publicly
    REGISTER(ObservationMap, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultObsMap);
    FIELD(numAssets, "The number of assets");
    FIELD_MAKE_TRANSIENT(numAssets);
    FIELD(observationDates, "Observation dates");
    FIELD_MAKE_TRANSIENT(observationDates);
    FIELD(sampleDates, "Sample dates");
    FIELD_MAKE_TRANSIENT(sampleDates);
    FIELD(modelDates, "Modelling dates");
    FIELD_MAKE_TRANSIENT(modelDates);
    FIELD(allDates, "Sample plus ref level dates");
    FIELD_MAKE_TRANSIENT(allDates);
    FIELD(modelDateMap, "Modelling date map");
    FIELD_MAKE_TRANSIENT(modelDateMap);
    FIELD(sampleMap, "Sample date map");
    FIELD_MAKE_TRANSIENT(sampleMap);
}

IObject* ObservationMap::defaultObsMap() {
    return new ObservationMap();
}

CClassConstSP const ObservationMap::TYPE = CClass::registerClassLoadMethod(
    "ObservationMap", typeid(ObservationMap), load);

DRLIB_END_NAMESPACE



