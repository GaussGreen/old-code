//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IndexedPerfList.cpp
//
//   Description : Class to sort doubles whilst providing access to which
//                 number ended up where
//
//   Author      : Mark A Robson
//
//   Date        : 29 November 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/IndexedPerfList.hpp"
#include "edginc/Maths.hpp"
#include <algorithm>

DRLIB_BEGIN_NAMESPACE

/** Creates IndexedPerfList holding numPerfs Cmpts. The assetIdx are
    initialised */
IndexedPerfList* IndexedPerfList::create(int numPerfs){
    return new IndexedPerfList(numPerfs);
}

/** Creates IndexedPerfList holding numPerfs Cmpts. The assetIdx are
    initialised */
IndexedPerfList::IndexedPerfList(int numPerfs): cmpts(numPerfs, (Cmpt*)0){
    try{
        for (int i = 0; i < numPerfs; i++){
            cmpts[i] = new Cmpt();
            cmpts[i]->assetIdx = i;
        }
    } catch (exception&){
        this->~IndexedPerfList();
        throw;
    }
}

/** Deletes components */
IndexedPerfList::~IndexedPerfList(){
    for (unsigned int i = 0; i < cmpts.size(); i++){
        delete cmpts[i];
    }
}

/** returns true if size == 0 */
bool IndexedPerfList::empty() const{
    return cmpts.empty();
}

/** Returns number of components */
int IndexedPerfList::size() const{
    return cmpts.size();
}

/** populate existing array of asset performances from a DoubleArray */
void IndexedPerfList::populate(const DoubleArray&  perfs){
    for (unsigned int i = 0; i < cmpts.size(); i++){
        cmpts[i]->assetIdx = i;
        cmpts[i]->perf = perfs[i];
    }
}

bool IndexedPerfList::sortMethod(const Cmpt* cmpt1, const Cmpt* cmpt2){
    return (cmpt1->perf > cmpt2->perf);
}

/** Sorts assets using perf field of components. Highest
    performance first */
void IndexedPerfList::sortByPerf(){
    sort(cmpts.begin(), cmpts.end(), sortMethod);
}

/** Calculates weighted sum using current order of IndexedPerf */
double IndexedPerfList::weightedSum(const DoubleArray&  weights) const{
    double sum = 0.0;
    for (unsigned int i = 0; i < cmpts.size(); i++){
        sum += weights[i]*cmpts[i]->perf;
    }
    return sum;
}

/** Store smallest absolute difference between each asset and its
    nearest neighbour (in terms of performance) whose performance is
    actually different (assets with equal performance are judged to be
    equally ranked) */
void IndexedPerfList::storeNearestPerfDiff(DoubleArray&  perfDiff) const{
    int size = cmpts.size();
    if (size < 2){
        // somewhat meaningless as there is no neighbour
        if (size > 0){
            perfDiff[0] = 0.0; // value shouldn't matter you would hope
        }
        return;
    }
    // deal with best and worst performers explicitly
    const Cmpt* best = cmpts.front();
    int i = 1; // start comparison with 2nd best asset
    do {
        perfDiff[best->assetIdx] = fabs(cmpts[i]->perf - best->perf);
        i++;
    } while (Maths::isZero(perfDiff[best->assetIdx]) && i < size);
    const Cmpt* worst = cmpts.back();
    i = size - 2; // start comparison with 2nd worst asset
    do {
        perfDiff[worst->assetIdx]= fabs(cmpts[i]->perf-worst->perf);
        i--;
    } while (Maths::isZero(perfDiff[worst->assetIdx]) && i >= 0);

    // then do remainder
    for (int iPerf = 1; iPerf < size-1; iPerf++){
        double thisPerf = cmpts[iPerf]->perf;
        double df1;
        // search amongst better performers
        int i = iPerf-1;
        do {
            df1 = fabs(cmpts[i]->perf - thisPerf);
            i--;
        } while (Maths::isZero(df1) && i >= 0);
        double df2;
        // search amongst worst performers
        i = iPerf+1;
        do {
            df2 = fabs(cmpts[i]->perf - thisPerf);
            i++;
        } while (Maths::isZero(df2) && i < size);
        // differences of zero don't count - the are equal in ranking
        int iAsset = cmpts[iPerf]->assetIdx;
        if (Maths::isZero(df1)){
            perfDiff[iAsset] = df2;
        } else if (Maths::isZero(df2)){
            perfDiff[iAsset] = df1;
        } else {
            perfDiff[iAsset] = Maths::min(df1, df2);
        }
    }      
}

DRLIB_END_NAMESPACE
