//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Hyb3CB.cpp
//
//   Description : temporary version of hyb3 base class work on the claim/bank and
//                 tree starting today without cluttering up the existing Hyb3
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Hyb3CB.hpp"
#include "edginc/IRVolPair.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/IndexSpecIR.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/Format.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/Results.hpp"
#include "edginc/MDFUtil.hpp"
#include "esl_log.h"

#include <stdio.h>

DRLIB_BEGIN_NAMESPACE


// export Hyb3 as an interface class with some member variables which means we
// don't have to supply default function implementations or export constructor
void Hyb3CB::load(CClassSP& clazz) {

    clazz->setPublic();
    REGISTER(Hyb3CB, clazz);
    SUPERCLASS(RateTree);

    FIELD(treeDataDebugFile, "Print tree debug data to named output file");
    FIELD_MAKE_OPTIONAL(treeDataDebugFile);
    FIELD(cetSmoothing,"");
    FIELD_MAKE_OPTIONAL(cetSmoothing);
    FIELD(zeroInterpStyle, "zero curve interpolation "
                 "style.  Defaults to LINEAR if not supplied");
    // ??? what about curveValueDate[2][3], currency etc.
    // ??? for tweaking what is copied
    FIELD_MAKE_OPTIONAL(zeroInterpStyle);
    FIELD(today, "today/reference date");
    FIELD_MAKE_TRANSIENT(today);
    FIELD(valueDate, "actual value date");
    FIELD_MAKE_TRANSIENT(valueDate);
}

// constructor - doesn't allocate any structures as derived models do that
Hyb3CB::Hyb3CB(const CClassConstSP &type) : RateTree(type),
treeBuilt(false), cetSmoothing(false), zeroInterpStyle(ZeroInterpStyle::LINEAR)
{
    // initialize the contents so that we have a reproductible 
    // behaviour if using uninitialized memory
    memset(&treeData, 0, sizeof(treeData));
    memset(&mktVolData, 0, sizeof(mktVolData));
    memset(&devData, 0, sizeof(devData));
    memset(&treeCurves, 0, sizeof(treeCurves));
}


// destructor - derived models should clear tree structures
Hyb3CB::~Hyb3CB() {
    if (treeBuilt) // cannot throw an exception from destructor
        fprintf(stderr,"%s: treeBuilt is still true\n",__FUNCTION__);
}


// interface reflection/class loading mechanism
CClassConstSP const Hyb3CB::TYPE = CClass::registerClassLoadMethod(
    "Hyb3CB", typeid(Hyb3CB), Hyb3CB::load);


bool Hyb3CBLoad(void) {
    return (Hyb3CB::TYPE != 0);
}


// register product discounting slice requirements
void Hyb3CB::registerZero(DateTime useDate, DateTime matDate, 
                          string curveName) {
    try {
        map<string, NamedClaimBank>::iterator cb;

        // any payments are dropped on or before the value date
        if (matDate <= valueDate)
            return;

        if (matDate < useDate) {
            throw ModelException("Cannot register zero with useDate "
            +useDate.toString()+" > matDate "+matDate.toString());
        }

        cb = claimBank.find(curveName);

        // if claimBank doesn't exist, create a new one and insert into map
        if (cb == claimBank.end()) {
            NamedClaimBank newClaimBank;
            newClaimBank.curveName = curveName;
            claimBank[curveName] = newClaimBank;
            cb = claimBank.find(curveName);
        }

        DateTime useDateL = max(today, useDate);
        cb->second.critZeroUseDates.push_back(useDateL.toIrDate());
        cb->second.critZeroMatDates.push_back(matDate.toIrDate());
        cb->second.critZeroLabel.push_back(string("ZERO/DISCOUNT"));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


// retrieve zero/discounting slice from the tree
void Hyb3CB::getZero(DateTime useDate, DateTime matDate, 
                     string curveName, TreeSliceSP &slice) {
    try {
        int sliceDim;
        map<string, NamedClaimBank>::iterator cb;
        DateTime currentDate = getCurrentDate();

        if (!dynamic_cast<TreeSliceRates*>(slice.get())) {
            slice = RateTree::createSimpleSlice(curveName);
            slice->name = "Zero "+useDate.toString()+" to "+matDate.toString();
        }
        TreeSliceRates& sliceL = dynamic_cast<TreeSliceRates&>(*slice);

        if (matDate <= valueDate)
            throw ModelException("Requested zero maturity date " + matDate.toString() +
                                 " must be greater than the model valueDate " +
                                 valueDate.toString());

        cb = claimBank.find(curveName);
        if (cb == claimBank.end())
            throw ModelException("Unable to find curveName " + curveName +
                                 "in the list of model registered claimBanks");

        if (useDate != currentDate)
            throw ModelException("requested zero observation date " + useDate.toString() +
                                 " must be the same as the model's current date " +
                                 currentDate.toString());

        CLAIM_BANK const* zeroBank = NULL;
        double* slicePtr = NULL;

        zeroBank = &cb->second.zeroBank;
        sliceDim = getCrvDim(curveName);
        sliceL.allocDim(sliceDim);  // reallocate slice dimension of required
        slicePtr = sliceL.getValuePtr();

        // ??? temporary fix as hyb3 zerobank function returns error if matDate == currentDate
        if (currentDate == matDate) {

            // set slice value as 1.0         
            if (Hyb3_SetSlice(slicePtr,
                              sliceDim,
                              1.0,
                              mTpIdx,
                              &treeData) != SUCCESS)
                throw ModelException("Hyb3_SetSlice failed: " + IrConverter::slog.pop());
        }
        else {
            double* zeroBkPt = NULL;

            // pointer to internal zero slice maturity - zeroBkPt is read only and is not freed
            zeroBkPt = Hyb3_ZbkReadZero((CLAIM_BANK*)zeroBank,
                                        matDate.toIrDate(),
                                        FALSE,  // ??? no interp, do we want to relax this
                                        currentDate.toIrDate(),
                                        mTpIdx,
                                        sliceDim,
                                        &treeData);

            if (zeroBkPt == NULL)
                throw ModelException("Hyb3_ZbkReadZero failed: " + IrConverter::slog.pop());

            // copy results of ReadZero function call to output slice
            if (Hyb3_CopySlice(slicePtr,
                               zeroBkPt,
                               sliceDim,
                               mTpIdx,
                               &treeData) != SUCCESS)
                throw ModelException("Hyb3_CopySlice failed: " + IrConverter::slog.pop());
        }

        slice->treeStep = range->treeStep;  // set internal slice timestep counter
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


// adjust parYield rate for any offsets where we are requesting/observing the
// rate no a date other than the reset date - uses deterministic ratio method
void Hyb3CB::dateAdjustIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, 
                               DateTime resetDate, TreeSlice& treeSlice) {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        string name = rateSpec.getName();
        char dcc = IrConverter::dccTo035A(*(rateSpec.dcc));
        char freq = IrConverter::toEslFreq(rateSpec.frequency.get());
        string curveName = rateSpec.getFactor()->getName();
        DateTime indexStartDate;
        double* sliceL = slice.getValuePtr();
        int currIdx = getCurrIdx(curveName);
        int curveIdx = getCrvIdx(curveName, currIdx);
        double ratio;

        if (currentDate == resetDate)
            return;  // no adjustment required if they're the same dates
        else if (currentDate < resetDate) {
            throw ModelException("Currently reset dates "+resetDate.toString()
                +" greater than the currentDate "+currentDate.toString()
                +" are not automatically adjusted - may be upgraded in the future");
        }

        if (rateSpec.fwdRateOffset.get())
            // ??? indexStartDate should be holiday adjusted
            indexStartDate = rateSpec.fwdRateOffset->toDate(resetDate);
        else
            indexStartDate = resetDate;  //  no offset

        // ??? re-investigate this later
        if (indexStartDate.daysDiff(resetDate) > 10) {
            throw ModelException("Reset date adjustment algorithm only supports small forward date "
                    "offsets (<=10 days) as the deterministic parYieldRatio function does "
                    "not currently account for forward starting rates.  Offset requested "
                    "in days = ", Format::toString(currentDate.daysDiff(resetDate)));
        }

        // if the reset date is before the currentDate when the rate is being asked for, 
        // approximate the rate by taking the ratio of the deterministic parYield rates 
        // between the actual reset date and the current date, then scale the stochastic 
        // payYield value by this ratio.  Essentially this is matching the first
        // moment.  We do not at this stage make volatility adjustments to the rate
        if (currentDate.daysDiff(resetDate) > 180) {
            throw ModelException("Current date ("+currentDate.toString()
                +") to resetDate ("+resetDate.toString()+") parYield approximate rate "
                    "adjustment is only supported if the date offset is < 180 days due to "
                    "issues around scaling the rate over too larger a time offset.  Requested "
                    "offset = "+Format::toString(currentDate.daysDiff(resetDate))
                    +" days for indexSpec "+rateSpec.name);
        }

        // returns ratio of (rate at resetDate) / (rate and current date)
        // ??? need a new function to take forward starting payYields rather than assuming 
        // the rate starts today - either spot offset or a significant forward starting offset
        if (ParYieldRatio(&ratio,
                        resetDate.toIrDate(), 
                        currentDate.toIrDate(),
                        0.0,  // no spread
                        &treeCurves[currIdx][curveIdx],
                        rateSpec.tenor->toMonths(),
                        dcc,
                        freq) != SUCCESS)
            throw ModelException(IrConverter::slog.pop());

        // scale stochastic rate by deterministic ratio of parYields
        if (Hyb3_MultiplyScalar(sliceL, 
                                slice.getDim(),
                                ratio, 
                                mTpIdx, 
                                &treeData) != SUCCESS) 
        {
            throw ModelException(IrConverter::slog.pop()
                +" Error scaling stochastic parYield by deterministic parYieldRatio ("
                +Format::toString(ratio)+") with currentDate = "+currentDate.toString()
                +", resetDate = "+resetDate.toString());
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb3CB::configureClaimBank(const CRIT_DATE* critDate, int nbCritDate) {

    try {
        int i;

        // add critical claim bank zero dates for each named bank
        map<string, NamedClaimBank>::iterator cb;

        // produce zero bank dates and optimize
        for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {

            string name = cb->first;
            try {
                int nbCrit = 0;
                long *critZeroMatDL = NULL;
                long *critZeroUseDL = NULL;
                // auto delete arrays
                IrConverter::AutoDeleteDRArray<long>
                    critZeroMatDLDelete(&critZeroMatDL, LONG, &nbCrit);
                IrConverter::AutoDeleteDRArray<long>
                    critZeroUseDLDelete(&critZeroUseDL, LONG, &nbCrit);

                // optimize non-critical dates
                int     nbMatDates = 0;
                long*   matDates   = NULL;
                int     nbUseDates = 0;
                long*   useDates   = NULL;
                int     j;

                // take copy of registered critical dates to use when printing out
                // the dates in the print function, along with the description labels
                cb->second.critZeroMatDatesUnsorted = cb->second.critZeroMatDates;
                cb->second.critZeroUseDatesUnsorted = cb->second.critZeroUseDates;

                // get pointer to required array to make debugging easier
                IRDate* optMatDates = &(cb->second.optZeroMatDates[0]);
                IRDate* optUseDates = &(cb->second.optZeroUseDates[0]);

                if (ZbkOptDates(cb->second.optZeroMatDates.size(),
                                optMatDates,
                                cb->second.optZeroUseDates.size(),
                                optUseDates,
                                nbCritDate,
                                const_cast<CRIT_DATE*>(critDate),
                                valueDate.toIrDate(),
                                &nbMatDates,
                                &matDates,
                                &nbUseDates,
                                &useDates) != SUCCESS)
                    throw ModelException(IrConverter::slog.pop());

                            
                // combine optimized dates with critical dates
                for (j=0; j<nbUseDates; ++j) {
                    cb->second.critZeroUseDates.push_back(useDates[j]);
                    cb->second.critZeroMatDates.push_back(matDates[j]);
                }

                // ??? for now, take a copy of the dateList to reuse cbkProcessDL function
                // which reallocs the arrays, so I have to use the C array rather than 
                // C++ vector
                nbCrit = cb->second.critZeroUseDates.size();
                if (cb->second.critZeroMatDates.size() != nbCrit) {
                    throw ModelException("Claim bank configuration error - number of registered critical "
                        "useDates "+Format::toString(nbCrit)
                        +" must equal number of registered critical matDates "
                        +Format::toString((int)cb->second.critZeroMatDates.size()));
                }

                critZeroMatDL = (long *) DR_Array (LONG, 0, nbCrit-1);
                critZeroUseDL = (long *) DR_Array (LONG, 0, nbCrit-1);

                if (critZeroMatDL == NULL || critZeroUseDL == NULL) {
                    throw ModelException("Unable to allocate memory for claim bank critical dates");
                }
                for (i = 0; i < nbCrit; i++) {
                    critZeroMatDL[i] = cb->second.critZeroMatDates[i];
                    critZeroUseDL[i] = cb->second.critZeroUseDates[i];
                }

                // sort and remove duplicates from the critical date list
                if (CbkProcessDL(&nbCrit,
                                &critZeroMatDL,
                                &nbCrit,
                                &critZeroUseDL) != SUCCESS)
                    throw ModelException(IrConverter::slog.pop());

                // copy result back to the claimBank vectors and resize
                cb->second.critZeroMatDates.resize(nbCrit);
                cb->second.critZeroUseDates.resize(nbCrit);

                for (i = 0; i < nbCrit; i++) {
                    cb->second.critZeroMatDates[i] = critZeroMatDL[i];
                    cb->second.critZeroUseDates[i] = critZeroUseDL[i];
                }

                // set the critical dates index for the last date
                cb->second.critDatesIdx = MAX(0, nbCrit - 1);
            }
            catch (exception& e) {
                throw ModelException(e, "In claim bank \""+name+"\"");
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


size_t Hyb3CB::getSliceSize(size_t nbDim) const {

    switch (nbDim) {
    case 1:
        return treeData.Width[0];
    case 2:
        return treeData.Width[0] *
               treeData.Width[1];
    case 3:
        return treeData.Width[0] *
               treeData.Width[1] *
               treeData.Width[2];
    default: 
        throw ModelException(__FUNCTION__, 
        "nbDim="+Format::toString(static_cast<int>(nbDim))+" must be 1,2 or 3");
    }
}


//------------------------------------------------------------------------------
// Get date on the timeline pointed by the index
//------------------------------------------------------------------------------
DateTime Hyb3CB::getDate(int dateIdx) const
{
    if (treeData.NbTP == -1)  // degenerate tree with no critical date (=0 price)
        return DateTime();

    if (dateIdx < 0 || dateIdx > treeData.NbTP) {
        throw ModelException(__FUNCTION__, 
        "Invalid time point index (" + Format::toString(dateIdx) 
        + "). Must be in the range 0..." + Format::toString(treeData.NbTP));
    }
    return DateTime::fromIrDate(treeData.TPDate[dateIdx]);
}

DateTimeArray Hyb3CB::getDates() const
{
    DateTimeArray dates( treeData.NbTP + 1 );
    for( int i = 0; i <= treeData.NbTP; ++i )
        dates[ i ] = DateTime::fromIrDate( treeData.TPDate[ i ] );
    return dates;
}

/** get last step (total num of steps) on time line */
int Hyb3CB::getLastStep() const
{
    return treeData.NbTP;
}

DateTime Hyb3CB::getCurveValueDate(string curveName) const {
    try {
        int currIdx = getCurrIdx(curveName);
        return currencyValueDate[currIdx];
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb3CB::update(int t, FDProduct::UpdateType type) {
    try {
        map<string, NamedClaimBank>::iterator cb;

        if (!treeBuilt)
            throw ModelException("Not able to update tree as the "
                                "tree has not been built");

        // set current range
        range->limits.bot1 = treeData.Bottom1[t];
        range->limits.top1 = treeData.Top1[t];
        range->limits.bot2 = treeData.Bottom2[t];
        range->limits.top2 = treeData.Top2[t];
        range->limits.bot3 = treeData.Bottom3[t];
        range->limits.top3 = treeData.Top3[t];
        range->treeStep = t;

        int T = treeData.NbTP;

        // Update tree to current time slice
        if (Hyb3_Lattice(&devData,
                        t,
                        T,
                        mktVolData,
                        &treeData) != SUCCESS)
            throw ModelException("Hyb3_Lattice Failed, " + IrConverter::slog.pop());

        mTpIdx = t;

        // update the registerd claim banks
        for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {

            int  addFlag = 0;
            long useDate = 0;
            long matDate = 0;
            int  idxActive = cb->second.critDatesIdx;
            string curveName = cb->second.curveName;
            long currentDate = getDate(mTpIdx).toIrDate();

            int zbDevMode;
            int currIdx = getCurrIdx(curveName);
            switch (treeData.TreeType) {
                case TTYPE_EQ1IR:
                    if (currIdx == 0)  // eqMode sets [0] as domestic index ??? change
                        zbDevMode = DISC_1D_NOCUPS;
                    else
                        throw ModelException("Internal model error - Tree is in EQ1IR mode "
                                            "and internal index must == 0 - value set as " +
                                            Format::toString(currIdx));
                    break;
                default:
                    throw ModelException("tree mode enum number " + 
                                        Format::toString(treeData.TreeType) + 
                                        " is not yet supported in qlib");
            }

            if (cb->second.critZeroMatDates.size() == 0)
                continue;

            if (idxActive >= 0 &&
                cb->second.critZeroMatDates[idxActive] >= currentDate) {

                addFlag = 1;
                useDate = cb->second.critZeroUseDates[idxActive];
                matDate = cb->second.critZeroMatDates[idxActive];

                cb->second.critDatesIdx--;
            }

            int curveIdx = getCrvIdx(curveName, currIdx);

            if (Hyb3_ZbkUpdate(&cb->second.zeroBank,
                               addFlag,
                               currentDate,
                               useDate,
                               mTpIdx,
                               treeData.NbTP,
                               curveIdx,
                               zbDevMode,
                               &devData,
                               &treeData) != SUCCESS)
                throw ModelException("Zero bank update failed for curve " + curveName + 
                                     " in DEV mode number " + Format::toString(zbDevMode) +
                                     ", " + IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


// ?? nothing to do in claim bank mode - remove when zerobankMode depricated
void Hyb3CB::updateFinalize(int t, FDProduct::UpdateType type) {

    if (!treeBuilt)
        throw ModelException(__FUNCTION__, "Not able to update tree as "
                             "the tree has not been built");
}


///////////////////////////////////
// roll forward tree from t+t to t
// maintaining internal zero banks
///////////////////////////////////
void Hyb3CB::updateStatePrice(int t, FDProduct::UpdateType type)
{
    int status = FAILURE;
    static char const* routine  = "Hyb3::update";
    TreeSliceRatesSP  tempSlice;

    int statePriceDim = -1; 
    //int T = mTreeData.NbTP;

    // set current range
    range->limits.bot1 = treeData.Bottom1[t];
    range->limits.top1 = treeData.Top1[t];
    range->limits.bot2 = treeData.Bottom2[t];
    range->limits.top2 = treeData.Top2[t];
    range->limits.bot3 = treeData.Bottom3[t];
    range->limits.top3 = treeData.Top3[t];
    range->treeStep = t;


    // allocate the slices on the first step
    if (type == FDProduct::FWD_0) {

        // Set the dimension of the state price (dependent on tree-type)
        switch (treeData.TreeType)
        {
        case TTYPE_FX2IR:
        case TTYPE_EQD2IR:
        case TTYPE_EQF2IR:
        case TTYPE_EQC2IR:
        case TTYPE_2IR2F1D:
            statePriceDim = 3;
            break;
        case TTYPE_EQ1IR:
        case TTYPE_2IR:
            statePriceDim = 2;
            break;
        default:
            throw ModelException("forward induction not implemented for "
                "model without EQ or FX dependency (i.e. CUPS2D etc)!");
            break;
        }

        statePriceSlice = DYNAMIC_POINTER_CAST<TreeSliceRates>(createSlice( ));
        statePriceSlice->name = "statePriceSlice";
        tempStatePriceSlice = DYNAMIC_POINTER_CAST<TreeSliceRates>(createSlice( ));
        tempStatePriceSlice->name = "tempStatePriceSlice";
        statePriceSlice->allocDim(statePriceDim);
        tempStatePriceSlice->allocDim(statePriceDim);
    }

    // swap the pointer around: we need to compute t+1 in order
    // to populate the FXRate correctly but need statePrice at t 
    // for express DEV calculation (might need to change in Lattice
    // in order to allow for forward induction as well - to do later)
    tempSlice           = tempStatePriceSlice;
    tempStatePriceSlice = statePriceSlice;
    statePriceSlice     = tempSlice;

    // check if we have the correct time
    statePriceSlice->testTreeStep() ; // should be this, but is protected currently...

    if (Hyb3_UpdateStatePrices( t,
        mktVolData,
        &treeData,
        &devData,                
        statePriceSlice->getValuePtr(),
        tempStatePriceSlice->getValuePtr() ) != SUCCESS )
    {
        DR_Error("Error updating the State Prices at time %d", t);
        goto RETURN;
    }

    // update the TimePoint Index
    mTpIdx = t;

    // update the next tree slices time step
    tempStatePriceSlice->treeStep = t+1;

    // StatePriceUpdate- function() : we need a state containing that info

    status = SUCCESS;

RETURN:

    if (status != SUCCESS)
        throw ModelException(routine, IrConverter::slog.pop());
}


/**  getStatePrices from the underlying model                               *
  *  supports the expressDEV (expressBACK) and the forward inductive scheme */
const TreeSlice & Hyb3CB::getStatePriceSlice(const int t) const
{
    // if we go forward in the tree, simply return the current state-price slice
    if (getInductionType() == FDModel::InductionType::FWD )
    {
        return *statePriceSlice; 
    }
    else if (getInductionType() == FDModel::InductionType::EXPRESSBACK )
    {
        map<int, TreeSliceSP>::const_iterator currentStatePriceSP = statePrices.find( t );
        if (currentStatePriceSP != statePrices.end() )
            return (*currentStatePriceSP->second);
    }
    else
        throw ModelException("state prices are only stored for forward and express backward"
                "induction schemes" );

    // return value to avoid warning of control path not returning values
    return *statePriceSlice;
};


void Hyb3CB::sliceExpand(TreeSlice& treeSlice, const string& curveName) const {
    try {
        int sliceDim = getCrvDim(curveName);
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
        slice.expand(sliceDim);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

/*
// Perform one backward induction step by moving slice from t+1 to t - no discounting
void Hyb3CB::sliceEv(TreeSlice& treeSlice) const {

    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
    // nothing to do at the back of the tree 
    if (mTpIdx == treeData.NbTP)
        return;   

    if (slice.getDim() != nbFactors)
        slice.allocDim(nbFactors);

    if (Fix3_Ev(slice.getValuePtr(), mTpIdx, treeData.NbTP, &devData, &treeData)
        != SUCCESS)
        throw ModelException(__FUNCTION__, IrConverter::IrConverter::slog.pop());
}
*/

void Hyb3CB::getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, DateTime resetDate, 
                        TreeSlice& treeSlice) {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        string name = rateSpec.getName();
        char dcc = IrConverter::dccTo035A(*rateSpec.dcc);
        char freq = IrConverter::toEslFreq(rateSpec.frequency.get());

        string curveName = rateSpec.getFactor()->getName();
        int curveDim = getCrvDim(curveName);
        if (slice.getDim() != curveDim)
            slice.allocDim(curveDim);

        DateTime indexStartDate;

        if (rateSpec.fwdRateOffset.get())
            // ??? indexStartDate should be holiday adjusted
            indexStartDate = rateSpec.fwdRateOffset->toDate(currentDate);
        else
            indexStartDate = currentDate;  //  no offset supplied

        if (indexStartDate < currentDate)
            throw ModelException("IR Index swap start date " + 
                                indexStartDate.toString() +
                                " can not be before the reset date " + currentDate.toString() +
                                ".  Index name = " + name);

        if (resetDate < today)
            throw ModelException("Tree can not calculate past rate index values - resetDate = "
                                + resetDate.toString() + " is before today = " + today.toString() +
                                ".  (NOTE: time of day may be the issue).  An internal  logic error " +
                                "exists as any past fixing requests should have been intercepted in the "
                                "index FDProduct code and this function should not have been called");

        if (currentDate < resetDate)
            throw ModelException("Currently future reset dates (" + resetDate.toString() + 
                                ") greater than the currentDate (" + currentDate.toString() + 
                                ") are not automatically adjusted to account for a forward observing " +
                                "offset.  This case requires future investigation");

        map<string, NamedClaimBank>::iterator cb = claimBank.find(curveName);
        if (cb == claimBank.end())
            throw ModelException("Unable to find named claim bank " + curveName +
                                " to calculate IR index payYield " + name);

        // annuity slice is not used, but required and populated for function call
        // so construct local slice
        TreeSliceRates annuitySlice(*range, "", -1);
        annuitySlice.allocDim(curveDim);

        double* sliceL = slice.getValuePtr();
        double* annuityL = annuitySlice.getValuePtr();

        // note: calculate the swap rate from the current date rather than the
        // reset date, and then adjust below if currentDate != resetDate
        if (Hyb3_ZbkParYield_t(sliceL,
                            annuityL,
                            curveDim,
                            &cb->second.zeroBank,
                            currentDate.toIrDate(),
                            indexStartDate.toIrDate(),
                            rateSpec.tenor->toMonths(),
                            dcc,
                            freq,
                            0.0,
                            mTpIdx,
                            &treeData) != SUCCESS)
            throw ModelException(IrConverter::slog.pop());

        slice.treeStep = range->treeStep;

        // adjust the stochastic rate if resetDate != currentDate - ie. estimating
        // the rate on a date that is not the actual reset date
        dateAdjustIRIndex(rateSpec, currentDate, resetDate, treeSlice);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP Hyb3CB::createMDF() const {
    MarketDataFetcherSP mdf(new RateTreeMDF(1));
    MDFUtil::setUseCurrencyBasis(*mdf, true);
    return mdf;
}


void Hyb3CB::populateIRParams(MKTVOL_DATA* mktvol, HYB3_TREE_DATA* tree, 
                             RateTree::IRModelParams* mp)
{
    mktvol->QLeft = 1.0 - mp->QLeft;
    mktvol->QRight = 1.0 - mp->QRight;
    mktvol->FwdShift = mp->FwdShift;
    mktvol->CetNbIter = mp->CetNbIter;

    for (int i = 0; i < 3; ++i)
    {
        mktvol->Alpha[i] = mp->FactorVol[i];
        mktvol->Beta[i] = mp->FactorMR[i];
        mktvol->Rho[i] = mp->FactorCorr[i];
    }

    // Check validity of input
    if (Hyb3_Param_Check(1,
                         mktvol,
                         tree) != SUCCESS) {        
        throw ModelException("Hyb3::PopulateIRParams", IrConverter::slog.pop());
    }
}


int Hyb3CB::getCrvIdx(const string& curveName, int currIdx) const {
    try {
        if(curveName.empty())
            throw ModelException("curveName not supplied");
        
        for (int i = 0; i < 3; i++) {
            if (curveName == this->curveName[currIdx][i])
                return i;
        }
        throw ModelException("Unable to find curveName " +
                            curveName + 
                            " in list of hyb3 registered curve names");
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

int Hyb3CB::getCurveIdx( const string & curveName ) const
{
    if (curveName.empty())
        return -1;

    int currIdx = getCurrIdx( curveName );
    int crvIdx = getCrvIdx( curveName, currIdx );
    return crvIdx * 10 + currIdx;
}


DRLIB_END_NAMESPACE
