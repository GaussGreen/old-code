//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : HybTree.hpp
//
//   Description : Hyb4
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/HybTree.hpp"
#include "edginc/TreeSlice.hpp"
#include "edginc/IRVolPair.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/IndexSpecIR.hpp"
#include "edginc/IndexSpecFX.hpp"
#include "edginc/IndexSpecEQ.hpp"
#include "edginc/SimpleEquity.hpp"
#include "edginc/FixedSettlement.hpp"
#include "edginc/RollingSettlement.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SRMFXVol.hpp"
#include "edginc/SRMEQVol.hpp"
#include "esl_log.h"
#include "esl_market.h"

#include "edginc/Format.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/Results.hpp"
#include "edginc/MDFUtil.hpp"

#include <stdio.h>

DRLIB_BEGIN_NAMESPACE

static int toDivType(int t) {
    switch (t) {
        case 0: return 'D'; // AMOUNT
        case 1: return 'Y'; // PERCENT
        case 2: return 'C'; // CONTINUOUS
        default: ASSERT(0); return -1;
    }
}

 // ??? reuse as much as existing for now - rewrite later   
 static void PopulateEQData(
    EQ_DATA    *eq_data,
    long       today,
    double     *eqSpotVolOverride,
    const SimpleEquity *eq,
    double     *EQCalibrationStrike,
    DividendListConstSP &dividendList)
{
    char *routine = "PopulateEQData";
    int i;
    
    if (eqSpotVolOverride==NULL)    /* 1.Disallow cut off      */
    {
        eq_data->EqCutOffFlag = FALSE;
        eq_data->EqCutOffLast = FALSE;
        eq_data->EqCutOffLevel   = EQSPOT_CUTOFF_R; /* Ratio to fwd vol level*/
        eq_data->EqBootStrapMode = EQ_NO_FAILURE_ALLOWED;
    }
    else                                         /* 3.Cut off at user level */
    {
        /* Cut off is allowed and there */
        /* is volatility cutoff data    */
        eq_data->EqCutOffFlag = TRUE;
        eq_data->EqCutOffLast = FALSE;
        eq_data->EqCutOffLevel = *eqSpotVolOverride;

        if (eq_data->EqCutOffLevel < TINY) /* The FAILURE return */
            throw ModelException(routine, "Unable to convert EQ cut off value!\n");

        eq_data->EqBootStrapMode = EQ_CONSTANT_SPOT_VOL;     
    }

    if (eq)
    {
        IrConverter::to_EQ_DATA(*eq_data, *eq);
        eq_data->ValueDate = today;

        // if no smile parameters in the environment, throw an error for now - sort out
        // properly when smile families are implemented
        if (eq_data->NbSmilePt == 0)
            throw ModelException("At least one equity smile date and a1,a2,a3 data set must "
                                 "be provided.  To turn smile off set a1=a2=a3 = 0.0");

        eq_data->Spot = eq->getSpot();

        // ??? Need to use a different vol representation, but not sure what
        const CVolBaseWrapper& eqVol = eq->getVol();
        const SRMEQVol *vol = dynamic_cast<const SRMEQVol*>(eqVol.get());
        if (!vol) 
            throw ModelException(routine, "Equity vol market structure must be of type SRM3EQVol");

        if (EQCalibrationStrike)
            throw ModelException(routine, "EQCalibrationStrike field not yet supported");
        /*
        LinearStrikeVolRequest volRequest(EQCalibrationStrike, 
                                            baseVolDates[0],
                                            baseVolDates.back(),
                                            false);

        DoubleArray vols(baseVolDates.size());
        CVolProcessedBS::calcVol(vol,
                                    &volRequest,
                                    eq,
                                    baseVolDates[0],
                                    baseVolDates,
                                    vols);
        */

        DateTimeArray compVolDates = vol->getCompVolDate();
        DoubleArray compVols = vol->getCompVol();

        eq_data->NbVol = compVolDates.size();
        for (i = 1; i <= eq_data->NbVol; i++)
        {
            eq_data->VolDate[i] = compVolDates[i-1].toIrDate();
            eq_data->Vol[i] = compVols[i-1];;
        }
        if (eq_data->NbVol > 0)
        {
            eq_data->VolDate[0] = today;
            eq_data->Vol[0] = eq_data->Vol[1];
        }

        DateTimeArray spotVolDates = vol->getSpotVolDate();
        DoubleArray spotVols = vol->getSpotVol();

        eq_data->NbInpSpotVol = spotVolDates.size();
        for (i = 0; i < eq_data->NbInpSpotVol; i++)
        {
            eq_data->InpSpotVolDate[i] = spotVolDates[i].toIrDate();
            eq_data->InpSpotVol[i] = spotVols[i];
        }

        // populate "static" equity data
        EquitySP equity = eq->getEquity();

        BorrowCurveSP borrow = equity->getBorrow();
        RhoBorrowPointwise rhoBorrow(0.0);  // ??? dummy sens
        ExpiryArrayConstSP borrowDates = equity->sensExpiries(&rhoBorrow);

        dividendList = equity->getDivList();
        const DividendArray& divs = dividendList->getArray();

        eq_data->NbBorrow = (*borrowDates).size();
        if (eq_data->NbBorrow > MAXNBDATE)
            throw ModelException(routine, "Nb of borrow curve points " +
                                Format::toString(eq_data->NbBorrow) + 
                                " exceeds maximum internal tree limit of " +
                                Format::toString(MAXNBDATE));

        // ??? have to check this properly
        for (i = 0; i < eq_data->NbBorrow; i++)
        {
            DateTime baseDate = equity->getValueDate();  // ??? no idea if this is correct date
            eq_data->BorrowDate[i] = (*borrowDates)[i]->toDate(baseDate).toIrDate();
            eq_data->Borrow[i] = borrow->interpAtDate(DateTime::fromIrDate(eq_data->BorrowDate[i]));
        }

        eq_data->NbFwd = divs.size();  // ??? need to check what dates in equity.sta are exactly
        if (eq_data->NbFwd > MAXNBEQDIV)
            throw ModelException(routine, "Number of dividends " +
                                Format::toString(eq_data->NbFwd) +
                                " exceeds maximum internal tree limit of " +
                                Format::toString(MAXNBEQDIV));
      
        // ??? first attempt at getting something to work
        for (i = 0; i < eq_data->NbFwd; i++)
        {
            eq_data->FwdDate[i] = divs[i].getExDate().toIrDate();
            eq_data->Fwd[i] = divs[i].getDivAmount();
            eq_data->FwdType[i] = toDivType(divs[i].getDivType()); 
        }

        const Settlement& settlement = equity->getSettlement();
        eq_data->SettleType = '0';
        {
            const FixedSettlement *fs = dynamic_cast<const FixedSettlement*>(&settlement);
            if (fs) {
                eq_data->SettleType = 'F';
                eq_data->NbSettle = fs->getFirstTradingDates().size();
                for (int i=0; i<eq_data->NbSettle; ++i) {
                    eq_data->LastTrading[i] = fs->getFirstTradingDates()[i].toIrDate();
                    eq_data->SettleDate[i] = fs->getSettlementDates()[i].toIrDate();
                }
            }
        }
        {
            const RollingSettlement *rs = dynamic_cast<const RollingSettlement*>(&settlement);
            if (rs) {
                eq_data->SettleType = 'R';
                eq_data->NbSettle = rs->getPeriod();
            }
        }
        if (eq_data->SettleType == '0') 
            throw ModelException("Unsupported settlement type "+settlement.getClass()->getName());
    }
    else
    {
        // if equity not used in product, set default values
        long dummyDate = Nxtday(today, 100);  // dummy first vol date

        eq_data->Spot = 0.0;
        eq_data->NbVol = 1;
        eq_data->ValueDate = today;
        eq_data->VolDate[0] = today;
        eq_data->VolDate[1] = dummyDate;
        eq_data->Vol[0] = 0.1;
        eq_data->Vol[1] = 0.1;
        eq_data->NbInpSpotVol = 0;
        eq_data->InpSpotVolDate[0] = today;
        eq_data->InpSpotVol[0] = 0.0;

        eq_data->NbBorrow = 1;
        eq_data->BorrowDate[0] = dummyDate;
        eq_data->Borrow[0] = 0.0;

        eq_data->NbFwd = 1;
        eq_data->FwdDate[0] = dummyDate;
        eq_data->Fwd[0] = 0.0;
        eq_data->FwdType[0] = 'D';

        eq_data->NbSettle = 0;
        eq_data->SettleType = 'R';  // ??? other option is fixed 'F'

        eq_data->NbSmilePt = 1;
        eq_data->SmileDate[0] = dummyDate;
        eq_data->a1[0] = 0.0;
        eq_data->a2[0] = 0.0;
        eq_data->a3[0] = 0.0;
    }

    // Check validity of input
    if (Hyb3_Eq_Check_Dyn_WType4(eq_data, today) != SUCCESS)
        throw ModelException(routine, IrConverter::slog.pop());

    if (Hyb3_Eq_Check_Sta_W(eq_data) != SUCCESS)
        throw ModelException(routine, IrConverter::slog.pop());
}

TreeSliceSP Hyb4temp::createSimpleSlice(
    const string & curveToDEV ) const
{
    TreeSliceRatesCompactSP newSlice( new TreeSliceRatesCompact(*range, curveToDEV, getCurveIdx(curveToDEV)));
    newSlice->treeStep = range->treeStep; //???? mTpIdx;
    return newSlice;
}

void Hyb4temp::load(CClassSP& clazz) {

    clazz->setPublic();
    REGISTER(Hyb4temp, clazz);
    SUPERCLASS(RateTree);

    FIELD(treeDataDebugFile, "Print tree debug data to named output file");
    FIELD_MAKE_OPTIONAL(treeDataDebugFile);
    FIELD(cetSmoothing,"");
    FIELD_MAKE_OPTIONAL(cetSmoothing);
    FIELD(zeroInterpStyle, "zero curve interpolation "
                 "style.  Defaults to LINEAR if not supplied");
    FIELD(nbFactorsDom,"Number of interest rate factors for domestic rate (default = 1)");
    FIELD_MAKE_OPTIONAL(nbFactorsDom);
    FIELD(nbFactorsFgn,"Number of interest rate factors for foreign rate (default = 1)");
    FIELD_MAKE_OPTIONAL(nbFactorsFgn);
    // ??? what about curveValueDate[2][3], currency etc.
    // ??? for tweaking what is copied
    FIELD_MAKE_OPTIONAL(zeroInterpStyle);
    FIELD(today, "today/reference date");
    FIELD_MAKE_TRANSIENT(today);
    FIELD(valueDate, "actual value date");
    FIELD_MAKE_TRANSIENT(valueDate);
}

// constructor - doesn't allocate any structures as derived models do that
Hyb4temp::Hyb4temp(const CClassConstSP &type) : RateTree(type),
treeBuilt(false), cetSmoothing(false), zeroInterpStyle(ZeroInterpStyle::LINEAR)
,nbFactorsDom(1),nbFactorsFgn(1)
{
    // initialize the contents so that we have a reproductible 
    // behaviour if using uninitialized memory
    memset(&treeData3, 0, sizeof(treeData3));
    memset(&treeData4, 0, sizeof(treeData4));
    memset(&mktVolData, 0, sizeof(mktVolData));
    memset(&devData, 0, sizeof(devData));
    memset(&treeCurves, 0, sizeof(treeCurves));
}


// destructor - derived models should clear tree structures
Hyb4temp::~Hyb4temp() {
    if (treeBuilt) // cannot throw an exception from destructor
        fprintf(stderr,"%s: treeBuilt is still true\n",__FUNCTION__);
}


// interface reflection/class loading mechanism
CClassConstSP const Hyb4temp::TYPE = CClass::registerClassLoadMethod(
    "Hyb4temp", typeid(Hyb4temp), Hyb4temp::load);


bool Hyb4tempLoad(void) {
    return (Hyb4temp::TYPE != 0);
}


// register product discounting slice requirements
void Hyb4temp::registerZero(DateTime useDate, DateTime matDate, 
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
void Hyb4temp::getZero(DateTime useDate, DateTime matDate, 
                     string curveName, TreeSliceSP &slice) {
    try {
        int sliceDim;
        map<string, NamedClaimBank>::iterator cb;
        DateTime currentDate = getCurrentDate();

        if (!dynamic_cast<TreeSliceRatesCompact*>(slice.get())) {
            slice = RateTree::createSimpleSlice(curveName);
            slice->name = "Zero "+useDate.toString()+" to "+matDate.toString();
        }
        TreeSliceRatesCompact& sliceL = dynamic_cast<TreeSliceRatesCompact&>(*slice);

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
            if (Hyb4_SetSlice(slicePtr,
                              sliceDim,
                              1.0,
                              mTpIdx,
                              &treeData4) != SUCCESS)
                throw ModelException("Hyb4_SetSlice failed: " + IrConverter::slog.pop());
        }
        else {
            double* zeroBkPt = NULL;

            // pointer to internal zero slice maturity - zeroBkPt is read only and is not freed
            zeroBkPt = Hyb4_ZbkReadZero((CLAIM_BANK*)zeroBank,
                                        matDate.toIrDate(),
                                        FALSE,  // ??? no interp, do we want to relax this
                                        currentDate.toIrDate(),
                                        mTpIdx,
                                        sliceDim,
                                        &treeData4);

            if (zeroBkPt == NULL)
                throw ModelException("Hyb4_ZbkReadZero failed: " + IrConverter::slog.pop());

            // copy results of ReadZero function call to output slice
            if (Hyb4_CopySlice(slicePtr,
                               zeroBkPt,
                               sliceDim,
                               mTpIdx,
                               &treeData4) != SUCCESS)
                throw ModelException("Hyb4_CopySlice failed: " + IrConverter::slog.pop());
        }

        slice->treeStep = range->treeStep;  // set internal slice timestep counter
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


// adjust parYield rate for any offsets where we are requesting/observing the
// rate no a date other than the reset date - uses deterministic ratio method
void Hyb4temp::dateAdjustIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, 
                               DateTime resetDate, TreeSlice& treeSlice) {
    try {
        TreeSliceRatesCompact& slice = dynamic_cast<TreeSliceRatesCompact&>(treeSlice);

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
        if (Hyb4_MultiplyScalar(sliceL, 
                                slice.getDim(),
                                ratio, 
                                mTpIdx, 
                                &treeData4) != SUCCESS) 
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

void Hyb4temp::configureClaimBank(const CRIT_DATE* critDate, int nbCritDate) {

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

size_t Hyb4temp::getSliceSize(size_t nbDim) const {

    if ((nbDim < 1) || (nbDim > 4))
        throw ModelException(__FUNCTION__, 
        "nbDim="+Format::toString(static_cast<int>(nbDim))+" must be 1, 2, 3, or 4");

    return treeData4.MaxIndex[nbDim-1] + 1;
}


//------------------------------------------------------------------------------
// Get date on the timeline pointed by the index
//------------------------------------------------------------------------------
DateTime Hyb4temp::getDate(int dateIdx) const
{
    if (treeData4.NbTP == -1)  // degenerate tree with no critical date (=0 price)
        return DateTime();

    if (dateIdx < 0 || dateIdx > treeData4.NbTP) {
        throw ModelException(__FUNCTION__, 
        "Invalid time point index (" + Format::toString(dateIdx) 
        + "). Must be in the range 0..." + Format::toString(treeData4.NbTP));
    }
    return DateTime::fromIrDate(treeData4.TPDate[dateIdx]);
}

DateTimeArray Hyb4temp::getDates() const
{
    DateTimeArray dates( treeData4.NbTP + 1 );
    for( int i = 0; i <= treeData4.NbTP; ++i )
        dates[ i ] = DateTime::fromIrDate( treeData4.TPDate[ i ] );
    return dates;
}

/** get last step (total num of steps) on time line */
int Hyb4temp::getLastStep() const
{
    return treeData4.NbTP;
}

DateTime Hyb4temp::getCurveValueDate(string curveName) const {
    try {
        int currIdx = getCurrIdx(curveName);
        return currencyValueDate[currIdx];
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb4temp::update(int t, FDProduct::UpdateType type) {
    try {
        map<string, NamedClaimBank>::iterator cb;

        if (!treeBuilt)
            throw ModelException("Not able to update tree as the "
                                "tree has not been built");

        range->limits.bot1 = treeData4.iMin[t];
        range->limits.top1 = treeData4.iMax[t];
        range->limits.bot2 = treeData4.jMin[t];
        range->limits.top2 = treeData4.jMax[t];
        range->limits.bot3 = treeData4.kMin[t];
        range->limits.top3 = treeData4.kMax[t];
        range->limits.bot4 = treeData4.LMin[t];
        range->limits.top4 = treeData4.LMax[t];

        range->offsets.offset1 = treeData4.NodeOffset0[t];
        range->offsets.offset2 = treeData4.NodeOffset1[t];
        range->offsets.offset3 = treeData4.NodeOffset2[t];
        range->offsets.offset4 = treeData4.NodeOffset3[t];

        range->treeStep = t;

        int T = treeData4.NbTP;

        // Update tree to current time slice
        if (Hyb4_Lattice(&devData,
                         &treeData4,
                         &latticeProg,
                         t,
                         T) != SUCCESS)
            throw ModelException("Hyb4_Lattice Failed, " + IrConverter::slog.pop());

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
            switch (treeData3.TreeType) {
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
                                        Format::toString(treeData3.TreeType) + 
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

            if (Hyb4_ZbkUpdate(&cb->second.zeroBank,
                               addFlag,
                               currentDate,
                               useDate,
                               mTpIdx,
                               treeData4.NbTP,
                               curveIdx,
                               zbDevMode,
                               &devData,
                               &treeData4) != SUCCESS)
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
void Hyb4temp::updateFinalize(int t, FDProduct::UpdateType type) {

    if (!treeBuilt)
        throw ModelException(__FUNCTION__, "Not able to update tree as "
                             "the tree has not been built");
}


///////////////////////////////////
// roll forward tree from t+t to t
// maintaining internal zero banks
///////////////////////////////////
void Hyb4temp::updateStatePrice(int t, FDProduct::UpdateType type)
{
}


/**  getStatePrices from the underlying model                               *
  *  supports the expressDEV (expressBACK) and the forward inductive scheme */
const TreeSlice & Hyb4temp::getStatePriceSlice(const int t) const
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


void Hyb4temp::sliceExpand(TreeSlice& treeSlice, const string& curveName) const {
    try {
        int sliceDim = getCrvDim(curveName);
        TreeSliceRatesCompact& slice = dynamic_cast<TreeSliceRatesCompact&>(treeSlice);
        slice.expand(sliceDim);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

/*
// Perform one backward induction step by moving slice from t+1 to t - no discounting
void Hyb4temp::sliceEv(TreeSlice& treeSlice) const {

    TreeSliceRatesCompact& slice = dynamic_cast<TreeSliceRatesCompact&>(treeSlice);
    // nothing to do at the back of the tree 
    if (mTpIdx == treeData4.NbTP)
        return;   

    if (slice.getDim() != nbFactors)
        slice.allocDim(nbFactors);

    if (Fix3_Ev(slice.getValuePtr(), mTpIdx, treeData4.NbTP, &devData, &treeData)
        != SUCCESS)
        throw ModelException(__FUNCTION__, IrConverter::IrConverter::slog.pop());
}
*/

void Hyb4temp::getIRIndex(const IndexSpecIR& rateSpec, DateTime currentDate, DateTime resetDate, 
                        TreeSlice& treeSlice) {
    try {
        TreeSliceRatesCompact& slice = dynamic_cast<TreeSliceRatesCompact&>(treeSlice);

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
        TreeSliceRatesCompact annuitySlice(*range, "", -1);
        annuitySlice.allocDim(curveDim);

        double* sliceL = slice.getValuePtr();
        double* annuityL = annuitySlice.getValuePtr();

        // note: calculate the swap rate from the current date rather than the
        // reset date, and then adjust below if currentDate != resetDate
        if (Hyb4_ZbkParYield_t(sliceL,
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
                            &treeData4) != SUCCESS)
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
MarketDataFetcherSP Hyb4temp::createMDF() const {
    MarketDataFetcherSP mdf(new RateTreeMDF(1));
    MDFUtil::setUseCurrencyBasis(*mdf, true);
    return mdf;
}

void Hyb4temp::init() {

    // initialise structures
    MktVol_Init(&mktVolData[FOREIGN]);
    MktVol_Init(&mktVolData[DOMESTIC]); 
    Hyb3_Tree_Init(&treeData3);
    Hyb4_Tree_Init(&treeData4);
    Hyb4_Dev_Init(&devData);
    Hyb4_Asset_IR_Init(&assetFor);
    Hyb4_Asset_IR_Init(&assetDom);
    Hyb4_Asset_FX_Init(&assetFx);
    Hyb4_Asset_EQ_Init(&assetEq);

    mTpIdx = -1;
//    mSmoothFlag = 0;  ???

    // initialiase the number of factors of the rates-modes
    mktVolData[FOREIGN].NbFactor = nbFactorsFgn;
    mktVolData[DOMESTIC].NbFactor = nbFactorsDom;

    treeBuilt = false;
}

void Hyb4temp::populateIRParams(MKTVOL_DATA* mktvol, HYB3_TREE_DATA* tree, 
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


int Hyb4temp::getCrvIdx(const string& curveName, int currIdx) const {
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

int Hyb4temp::getCurveIdx( const string & curveName ) const
{
    if (curveName.empty())
        return -1;

    int currIdx = getCurrIdx( curveName );
    int crvIdx = getCrvIdx( curveName, currIdx );
    return crvIdx * 10 + currIdx;
}

///////////////////////////////////////////////////////////

// helper functions to remove/integrate later on
static
void Hyb4CorrelOverwrites(
    FX_DATA            *fx_data,                    // (O) Fx data
    EQ_DATA            *eq_data,                    // (O) Eq data
    char               OverWriteString[6][MAXBUFF], // (I) Overwrite strings
    CORRELATION_DATA   *correlation_data)
{
    try {
        double temp[6];
                 
        temp[0] = correlation_data->CorrIR;
        temp[1] = correlation_data->CorrDomIRFX;
        temp[2] = correlation_data->CorrForIRFX;
        temp[3] = correlation_data->CorrDomIREQ;
        temp[4] = correlation_data->CorrForIREQ;
        temp[5] = correlation_data->CorrFXEQ;

        /*********************************************/
        /* Although the fx_data structure supports a */
        /* term-structure of correlations, we only   */
        /* use the first element in Build Tree       */
        /*********************************************/

        fx_data->Rho[0][0] = temp[0];
        fx_data->Rho[1][0] = temp[2];
        fx_data->Rho[2][0] = temp[1];
        fx_data->Rho[3][0] = temp[4];   
        fx_data->Rho[4][0] = temp[3];
        fx_data->Rho[5][0] = temp[5];   
                    
        if (strcmp (OverWriteString[0], "nil"))
        {                    
            fx_data->Rho[0][0] = atof (OverWriteString[0]);
        }
              
        if (strcmp (OverWriteString[1], "nil"))
        {
            fx_data->Rho[1][0] = atof (OverWriteString[1]);          
        }
              
        if (strcmp (OverWriteString[2], "nil"))
        {                    
            fx_data->Rho[2][0] = atof (OverWriteString[2]);                               
        }
              
        if (strcmp (OverWriteString[3], "nil"))
        {                    
            eq_data->Rho[3][0] = atof (OverWriteString[0]);
        }
              
        if (strcmp (OverWriteString[4], "nil"))
        {
            eq_data->Rho[4][0] = atof (OverWriteString[1]);          
        }
              
        if (strcmp (OverWriteString[5], "nil"))
        {                    
            eq_data->Rho[5][0] = atof (OverWriteString[2]);                               
        }

        if (Hyb3_Correl_Check_WType6 (fx_data, eq_data) != SUCCESS)
        {          
            throw ModelException(IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

static
void Hyb4FXSmileOverwrites(
    FX_DATA            *fx_data,                 
    char               OverWriteString[MAXBUFF], 
    FXVOLATILITY_DATA  *Volatility_data)
{
    try {
        long
            Year[2],
            Month[2],
            Day[2];
        int
            i;

        fx_data->Today = Volatility_data->ValueDate;
              
        fx_data->ValueDate = fx_data->Today;
        fx_data->SpotDays  = 0;

        fx_data->Spot = Volatility_data->FXSpotRate;
        // ignore BaseVolFreq entry

        fx_data->NbVol = Volatility_data->NbBaseVols;
        if (fx_data->NbVol > MAXNBDATE)
        {
            throw ModelException("Nb of vols exceeds max limit of "
                +Format::toString(MAXNBDATE)+" in (Fx_Input_W_WithSmile_DRI)");
        }
         
        Dsplit(fx_data->ValueDate, /* Split value date into month, day and year */
                &(Month[0]), 
                &(Day[0]), 
                &(Year[0]));

        for (i = 1; i <= fx_data->NbVol; i++)  /* Composite vols -> 1 offset*/
        {
            fx_data->VolDate[i] = Volatility_data->BaseVolDates[i-1];
            fx_data->FxVol[i] = Volatility_data->BaseVols[i-1]/100.0;

        }  /* for i */

        /* Repeat the first values for interpolation of pseudocomposite vols */
        /* in Get_TreeSpotVols                                               */

        if (fx_data->NbVol > 0)
        {
            fx_data->VolDate[0] = fx_data->ValueDate;
            fx_data->FxVol[0]   = fx_data->FxVol[1];
        }
         

        fx_data->NbInpSpotVol = Volatility_data->NbSpotVols;
        if (fx_data->NbInpSpotVol > MAXNBDATE)
        {
            throw ModelException("Nb of Inp Spot Vol exceeds max limit of "
                +Format::toString(MAXNBDATE)+" in ! (Fx_Input_W)");
        }

        if (fx_data->NbInpSpotVol > 0)
        {
        for (i = 0; i < fx_data->NbInpSpotVol; i++)   /* Spot vols -> 0 Offset*/
        {
                fx_data->InpSpotVolDate[i] = Volatility_data->SpotVolDates[i];
                fx_data->InpSpotVol[i] = Volatility_data->SpotVols[i]/100.0;
        }
        }


        /* deal with possible overlap between composite and spot vols */
        if ((fx_data->NbInpSpotVol > 0) && (fx_data->NbVol > 0))
        {
            long Idx;
            Idx = GetDLOffset(fx_data->NbVol,
                            &(fx_data->VolDate[1]),
                            fx_data->InpSpotVolDate[0],
                            CbkHIGHER);
            /*******************************************************/
            /* if all composite dates are < spotvoldate[0]         */
            /* then Idx = -999, in which case NbVol does not need  */
            /* to be changed.                                      */
            /*******************************************************/
            if (Idx >= 0L)
            {
                fx_data->NbVol = Idx;
            }        
        }
        
        if (strcmp (OverWriteString, "nil") == 0)    /* 1.Disallow cut off      */
        {

            fx_data->FxCutOffFlag = FALSE;
            fx_data->FxCutOffLast = FALSE;
            fx_data->FxCutOffLevel   = FXSPOT_CUTOFF_R; /* Ratio to fwd vol level*/
            fx_data->FxBootStrapMode = FX_NO_FAILURE_ALLOWED; 

        }
        else if (strcmp(OverWriteString,"last") == 0) /* 2.Cut off at last level */
        {
            /* Cut off is allowed and there */
            /* is volatility cutoff data    */
            fx_data->FxCutOffFlag = TRUE;
            fx_data->FxCutOffLast = TRUE;
            fx_data->FxCutOffLevel = FXSPOT_CUTOFF_R; /* Ratio to fwd vol level  */
            fx_data->FxBootStrapMode = FX_USE_LAST_LEVEL ;
        } 
        else                                         /* 3.Cut off at user level */
        {
            /* Cut off is allowed and there */
            /* is volatility cutoff data    */
            fx_data->FxCutOffFlag = TRUE;
            fx_data->FxCutOffLast = FALSE;
            fx_data->FxCutOffLevel = atof(OverWriteString);
            if (ABS(fx_data->FxCutOffLevel) < TINY) /* The FAILURE return */
            {
                throw ModelException("Unable to convert FX cut off value!");
            }
            fx_data->FxCutOffLevel /= 100.;
            fx_data->FxBootStrapMode = FX_CONSTANT_SPOT_VOL;
              
        }

        /* no equity smile yet */

        fx_data->NbSmilePt = 1;

        /* add arbitary number of days after the value date to ensure the smile date is
        > today as in some cases valueDate = Today. - chose to use 5 days but could have
        use any number */

        fx_data->SmileDate[0] = Nxtday(fx_data->ValueDate, 5);
        fx_data->a1[0] = 0.0;
        fx_data->a2[0] = 0.0;
        fx_data->a3[0] = 0.0;

        /* Check validity of input */ 
        if (Hyb3_Fx_Check_W (fx_data) != SUCCESS)
        {
            throw ModelException(IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

static
void Hyb4EQSmileOverwrites(
    EQ_DATA            *eq_data,                 
    char               OverWriteString[MAXBUFF], 
    EQVOLATILITY_DATA  *Volatility_data)
{
    try {
        long
            Year[2],
            Month[2],
            Day[2];
        int
            i;

        eq_data->ValueDate = Volatility_data->ValueDate;

        eq_data->Spot = Volatility_data->EQSpotRate;
        // ignore BaseVolFreq entry

        eq_data->NbVol = Volatility_data->NbBaseVols;
        if (eq_data->NbVol > MAXNBDATE)
        {
            throw ModelException("Nb of vols exceeds max limit of "
                +Format::toString(MAXNBDATE)+" in (Eq_Input_W_WithSmile_DRI)");
        }
         
        Dsplit(eq_data->ValueDate, /* Split value date into month, day and year */
                &(Month[0]), 
                &(Day[0]), 
                &(Year[0]));

        for (i = 1; i <= eq_data->NbVol; i++)  /* Composite vols -> 1 offset*/
        {
            eq_data->VolDate[i] = Volatility_data->BaseVolDates[i-1];
            eq_data->Vol[i] = Volatility_data->BaseVols[i-1]/100.0;

        }  /* for i */

        /* Repeat the first values for interpolation of pseudocomposite vols */
        /* in Get_TreeSpotVols                                               */

        if (eq_data->NbVol > 0)
        {
            eq_data->VolDate[0] = eq_data->ValueDate;
            eq_data->Vol[0]   = eq_data->Vol[1];
        }
         

        eq_data->NbInpSpotVol = Volatility_data->NbSpotVols;
        if (eq_data->NbInpSpotVol > MAXNBDATE)
        {
            throw ModelException("Nb of Inp Spot Vol exceeds max limit of "
                +Format::toString(MAXNBDATE)+" in ! (Fx_Input_W)");
        }

        if (eq_data->NbInpSpotVol > 0)
        {
        for (i = 0; i < eq_data->NbInpSpotVol; i++)   /* Spot vols -> 0 Offset*/
        {
                eq_data->InpSpotVolDate[i] = Volatility_data->SpotVolDates[i];
                eq_data->InpSpotVol[i] = Volatility_data->SpotVols[i]/100.0;
        }
        }


        /* deal with possible overlap between composite and spot vols */
        if ((eq_data->NbInpSpotVol > 0) && (eq_data->NbVol > 0))
        {
            long Idx;
            Idx = GetDLOffset(eq_data->NbVol,
                            &(eq_data->VolDate[1]),
                            eq_data->InpSpotVolDate[0],
                            CbkHIGHER);
            /*******************************************************/
            /* if all composite dates are < spotvoldate[0]         */
            /* then Idx = -999, in which case NbVol does not need  */
            /* to be changed.                                      */
            /*******************************************************/
            if (Idx >= 0L)
            {
                eq_data->NbVol = Idx;
            }        
        }
        
        if (strcmp (OverWriteString, "nil") == 0)    /* 1.Disallow cut off      */
        {

            eq_data->EqCutOffFlag = FALSE;
            eq_data->EqCutOffLast = FALSE;
            eq_data->EqCutOffLevel   = EQSPOT_CUTOFF_R; /* Ratio to fwd vol level*/
            eq_data->EqBootStrapMode = EQ_NO_FAILURE_ALLOWED; 

        }
        else if (strcmp(OverWriteString,"last") == 0) /* 2.Cut off at last level */
        {
            /* Cut off is allowed and there */
            /* is volatility cutoff data    */
            eq_data->EqCutOffFlag = TRUE;
            eq_data->EqCutOffLast = TRUE;
            eq_data->EqCutOffLevel = EQSPOT_CUTOFF_R; /* Ratio to fwd vol level  */
            eq_data->EqBootStrapMode = EQ_USE_LAST_LEVEL ;
        } 
        else                                         /* 3.Cut off at user level */
        {
            /* Cut off is allowed and there */
            /* is volatility cutoff data    */
            eq_data->EqCutOffFlag = TRUE;
            eq_data->EqCutOffLast = FALSE;
            eq_data->EqCutOffLevel = atof(OverWriteString);
            if (ABS(eq_data->EqCutOffLevel) < TINY) /* The FAILURE return */
            {
                throw ModelException("Unable to convert EQ cut off value!");
            }
            eq_data->EqCutOffLevel /= 100.;
            eq_data->EqBootStrapMode = EQ_CONSTANT_SPOT_VOL;
              
        }

        /* no equity smile yet */

        eq_data->NbSmilePt = 1;

        /* add arbitary number of days after the value date to ensure the smile date is
        > today as in some cases valueDate = Today. - chose to use 5 days but could have
        use any number */

        eq_data->SmileDate[0] = Nxtday(eq_data->ValueDate, 5);
        eq_data->a1[0] = 0.0;
        eq_data->a2[0] = 0.0;
        eq_data->a3[0] = 0.0;

        /* Check validity of input */ 
        if (Hyb3_Eq_Check_W (eq_data) != SUCCESS)
        {
            throw ModelException(IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

 /********************
 *
 * Hyb4 EQ Tree Mode
 *
**********************/


void Hyb4CB::load(CClassSP& clazz){

    clazz->setPublic();
    REGISTER(Hyb4CB, clazz);
    SUPERCLASS(Hyb4temp);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(ppy, "minimum number of tree nodes(points)/year");
    FIELD(maxStdDeviations, "number of standard deviations to trim the tree at");
    FIELD(forIRParams, "Collection of foreign interest rate model parameters");
    FIELD(domIRParams, "Collection of domestic interest rate model parameters");
    FIELD(FXSpotVolOverride,"Overrides market comp/spot vols for FX - sets constant spot vol for FX "
                            "diffusion (units = decimal NOT percentage)");
    FIELD_MAKE_OPTIONAL(FXSpotVolOverride);
    FIELD(EQSpotVolOverride,"Overrides market comp/spot vols for EQ - sets constant spot vol for EQ "
                            "diffusion (units = decimal NOT percentage)");
    FIELD_MAKE_OPTIONAL(EQSpotVolOverride);
    FIELD(legacyHyb4FXEQMode,"Override valueDate with today and do not extend yield curves back to today");
    FIELD_MAKE_OPTIONAL(legacyHyb4FXEQMode);
    // these fields are transient to allow tweaking to copy 
    // the model after ::getMarket() has been called 
    FIELD(corrIR, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrIR);
    FIELD(corrForIRFX, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrForIRFX);
    FIELD(corrDomIRFX, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrDomIRFX);
    FIELD(corrForIREQ, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrForIREQ);
    FIELD(corrDomIREQ, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrDomIREQ);
    FIELD(corrFXEQ, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrFXEQ);
}

CClassConstSP const Hyb4CB::TYPE = CClass::registerClassLoadMethod(
    "Hyb4CB", typeid(Hyb4CB), Hyb4CB::load);

// for class loading
bool Hyb4CBLoad(void) { return (Hyb4CB::TYPE != 0); }

Hyb4CB::~Hyb4CB() {
    clear();
}

void Hyb4CB::clear() {

    if (!treeBuilt)
        return;
    treeBuilt = false;

    map<string, NamedClaimBank>::iterator cb;
    for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {

        string curveName = cb->second.curveName;
        int curveDim;
        if (getCurrIdx(curveName) == FOREIGN)
            curveDim = 1;
        else
            curveDim = 2;

        Hyb4_CbkFree(&cb->second.zeroBank, &treeData4, curveDim);
    }

    Hyb4_Dev_Free(&devData, &treeData4);

    Hyb4_Asset_IR_Free(&assetFor, &treeData4);
    Hyb4_Asset_IR_Free(&assetDom, &treeData4);
    Hyb4_Asset_FX_Free(&assetFx , &treeData4);
    Hyb4_Asset_EQ_Free(&assetEq , &treeData4);
    
    Hyb4_Offset_Free(&treeData4);
    
    Hyb4_Tree_Limits_Unassign(&treeData4);
    
    Hyb4_Tree_Free(&treeData4);
    Hyb3_Tree_Free(&treeData3);

    for (int i=0; i<NBCRV; ++i)
        today2ValDateZeros[0][i] = -1.0;

    Hyb4CB::init();
}

IRGridPointAbsArraySP Hyb4CB::getSensitiveIRVolPoints(
        OutputNameConstSP  outputName,
        const CInstrument* inst) const {

    try {
        IRGridPointAbsArraySP volExposuresSP(new IRGridPointAbsArray);

        T_CURVE diffusionCurve;
        string volCalibIndex;
        IRVolRawSP volSP;
        CcyIRParamsSP ccyIRParams;

        // Note that checks below for matching outputName are just for efficiency
        // (e.g. to avoid unnecessary construction of diffusion T_CURVEs, vol
        // selection, and exposure processing.

        volSP = getIRVolRaw(forIRParams);
        if (outputName->equals(volSP->getBaseVol()->getName()) ||
            outputName->equals(volSP->getSwapVol()->getName()))
        {
            IrConverter::to_T_CURVE(
                diffusionCurve, forIRParams->curveToDiffuse.get(), false);

            volCalibIndex = forIRParams->volCalibIndex;
            ccyIRParams = forIRParams;
        } 
        else {
            volSP = getIRVolRaw(domIRParams);
            if (outputName->equals(volSP->getBaseVol()->getName()) ||
                outputName->equals(volSP->getSwapVol()->getName()))
            {
                IrConverter::to_T_CURVE(
                    diffusionCurve, domIRParams->curveToDiffuse.get(), false);

                volCalibIndex = domIRParams->volCalibIndex;
                ccyIRParams = domIRParams;
            }
            else 
                throw ModelException("outputName = \""+outputName->toString()
                + "\" does not match any vol object");
        }
        IRVolSelector volSelector(volSP, diffusionCurve, volCalibIndex,
                                  cetSmoothing, ccyIRParams);
        volSelector.getVolExposures(volExposuresSP, outputName);
        return volExposuresSP;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb4CB::getMarket(const MarketData*  market, IInstrumentCollectionSP instruments) {
}

void Hyb4CB::initModel() {

    try {
        int EqRateDim = 1;

        treeBuilt = true; 
        treeData3.FxMomentMatching = FALSE;

        CRIT_DATE* critDate=NULL;
        int nbCritDate = 0;
        IrConverter::AutoDeleteDRArray<CRIT_DATE>
            critDateDelete(&critDate, CRITDATE, &nbCritDate);

        // ??? note:  ZeroInterpTypeFlag is global variable defined in both Hyb3.lib
        // and fix3.lib - tidy up but have to check all Hyb3 products
        EslSetZeroInterpolation(zeroInterpStyle == ZeroInterpStyle::LINEAR ?
            ESL_INTERP_LINEAR : ESL_INTERP_FLATFWD);
        
        int i, j;

        map<string, NamedClaimBank>::iterator cb;

        // needed to setup timeline
        critDate = (CRIT_DATE *) DR_Array (CRITDATE, 0, 0);
        if (critDate == NULL)
            throw ModelException("Unable to allocate memory for critDate array");

        // today always critical date
        if (Add_To_DateList(&nbCritDate, 
                            &critDate,
                            today.toIrDate(),
                            0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
            throw ModelException("Add_To_DateList function failed : " + IrConverter::slog.pop());

        // valueDate always critical date
        if (Add_To_DateList(&nbCritDate, 
                            &critDate,
                            valueDate.toIrDate(),
                            0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
            throw ModelException("Add_To_DateList function failed : " + IrConverter::slog.pop());

        // add critical dates registered directly by products
        for (j = 0; j < critDates.size(); j++) {

            DateTime critDateL = critDates[j];
            if (critDateL > valueDate) {
                if (Add_To_DateList(&nbCritDate, 
                                    &critDate,
                                    critDateL.toIrDate(),
                                    0, 0, 0, 0, 0, 0, 0, 0, 0) != SUCCESS)
                    throw ModelException("Add_To_DateList function failed : " 
                                         + IrConverter::slog.pop());
            }
        }

        // store zero/discount factor between today and the value date to
        // in order to forward the price back to the value date
        for (i = 0; i < 2; i++) {
            for (j = 0; j < NBCRV; j++) {
                today2ValDateZeros[i][j] = ::pow(1.0 + treeCurves[i][j].Zero[0], 
                    -Daysact(today.toIrDate(), currencyValueDate[i].toIrDate())/365.0);

                // Adjust the zero dates from the value to today by scaling all of 
                // the zeroRates in the TCurves - domestic eq mode only uses index 0 of TCurve.
                if (!legacyHyb4FXEQMode)
                    RateTree::ExtendTreeCurve(&treeCurves[i][j], today);
            }
        }

        // configure and add claim banks
        configureClaimBank(critDate, nbCritDate);

        // add claim bank critical dates to list
        for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {

            // add combined critical zero bank dates
            int nbCrit = cb->second.critZeroUseDates.size();

            for (j=0; j<nbCrit; ++j) {
                long matDate = cb->second.critZeroMatDates[j];
                long useDate = cb->second.critZeroUseDates[j];
                if (Add_To_DateList(&nbCritDate,
                                    &critDate, 
                                    matDate,
                                    0, 0, 0, 0, 
                                    0, 0, 0, 0, 0) != SUCCESS)
                    throw ModelException("Add_To_DateList function failed : " 
                                         + IrConverter::slog.pop());

                if (Add_To_DateList(&nbCritDate, 
                                    &critDate,
                                    useDate,
                                    0,0, 0, 0, 
                                    0, 0, 0, 0, 0) != SUCCESS)
                    throw ModelException("Add_To_DateList function failed : " 
                                         + IrConverter::slog.pop());
            }
        }
    
        // remove duplicates, sort critical dates and store them to
        // help debugging by printing them out in the debug data file
        if (Sort_CritDate(nbCritDate, critDate) != SUCCESS)
            throw ModelException("Sort_CritDate failed: " + IrConverter::slog.pop());

        sortedCritDates.resize(1);
        sortedCritDates[0] = critDate[0].CritDate;
        for (i = 1; i < nbCritDate; i++) {
            if (sortedCritDates.back() == critDate[i].CritDate)
                continue;
            sortedCritDates.push_back(critDate[i].CritDate);
        }
        if (sortedCritDates.size()<=1) {
            range.reset(new TreeSliceRatesCompact::Range(0,0,0,0));
            treeData4.NbTP = -1;
            return;
        }
        // construct the timeline 
        if (Hyb3_Time_Line(today.toIrDate(), nbCritDate, critDate, 'I', &treeData3) 
            != SUCCESS)
            throw ModelException("Hyb3_Time_Line failed: " + IrConverter::slog.pop());

        // if printing out a debug file, turn on printing the default CET debug file also
        int printCET = FALSE;
        if (!treeDataDebugFile.empty())
            printCET = TRUE;

        string foreignIR  = forIRParams->curveToDiscount->getName();
        string domesticIR = domIRParams->curveToDiscount->getName();

        if (mEQName == domesticIR)
        {
            EqRateDim = 1;
            treeData3.TreeType = TTYPE_EQDFX2IR;
        }
        else if (mEQName == foreignIR)
        {
            EqRateDim = 0;
            treeData3.TreeType = TTYPE_EQFFX2IR;
        }
        else
        {
            throw ModelException("Currency of denomination of equity is missing: "
                                 + IrConverter::slog.pop());
        }

        if (Hyb3_Build_Tree(printCET, treeCurves, mktVolData, &fxData, &eqData, 
                            &treeData3) != SUCCESS)
            throw ModelException("Hyb3_Build_Tree failed: " + IrConverter::slog.pop());

        if (Hyb4_Tree_Alloc(&treeData4, &treeData3) != SUCCESS)
            throw ModelException("Hyb4_Tree_Alloc failed: " + IrConverter::slog.pop());

        treeData4.NbAsset   = 4;
        treeData4.NbAssetOn = 4;

        if (Hyb4_Tree_Limits_Assign(&treeData4, &treeData3) != SUCCESS)
            throw ModelException("Hyb4_Tree_Limits_Assign failed: " + IrConverter::slog.pop());

        if (Hyb4_Asset_IR_Alloc(&assetFor, &treeData4) != SUCCESS)
            throw ModelException("Hyb4_Asset_IR_Alloc failed: " + IrConverter::slog.pop());

        if (Hyb4_Asset_IR_Alloc(&assetDom, &treeData4) != SUCCESS)
            throw ModelException("Hyb4_Asset_IR_Alloc failed: " + IrConverter::slog.pop());

        if (Hyb4_Asset_FX_Alloc(&assetFx, &treeData4) != SUCCESS)
            throw ModelException("Hyb4_Asset_FX_Alloc failed: " + IrConverter::slog.pop());

        if (Hyb4_Asset_EQ_Alloc(&assetEq, &treeData4) != SUCCESS)
            throw ModelException("Hyb4_Asset_EQ_Alloc failed: " + IrConverter::slog.pop());

        if (Hyb4_CopyTreeFxEq(&treeData4, &treeData3, EqRateDim) != SUCCESS)
            throw ModelException("Hyb4_CopyTreeFxEq failed: " + IrConverter::slog.pop());

        if (Hyb4_Offset_Alloc(&treeData4) != SUCCESS)
            throw ModelException("Hyb4_Offset_Alloc failed: " + IrConverter::slog.pop());

        if (Hyb4_CopyAssetIR(&treeData4, &treeData3, &assetFor, &mktVolData[0], 0) != SUCCESS)
            throw ModelException("Hyb4_CopyAssetIR failed: " + IrConverter::slog.pop());

        assetFor.Fx = 2;

        if (Hyb4_CopyAssetIR(&treeData4, &treeData3, &assetDom, &mktVolData[1], 1) != SUCCESS)
            throw ModelException("Hyb4_CopyAssetIR failed: " + IrConverter::slog.pop());

        assetDom.Fx = -1;

        if (Hyb4_CopyAssetFX(&treeData4, &treeData3, &assetFx) != SUCCESS)
            throw ModelException("Hyb4_CopyAssetFX failed: " + IrConverter::slog.pop());

        assetFx.IrDom = 1;
        assetFx.IrFor = 0;

        if (Hyb4_CopyAssetEQ(&treeData4, &treeData3, &assetEq) != SUCCESS)
            throw ModelException("Hyb4_CopyAssetEQ failed: " + IrConverter::slog.pop());

        assetEq.IrDom = EqRateDim;

        treeData4.AssetType[0] = IRF;
        treeData4.asset[0] = (void *)&assetFor;

        treeData4.AssetType[1] = IRD;
        treeData4.asset[1] = (void *)&assetDom;

        treeData4.AssetType[2] = FX;
        treeData4.asset[2] = (void *)&assetFx;

        if (EqRateDim == 0)
        {
            treeData4.AssetType[3] = EQF;
        }
        else
        {
            treeData4.AssetType[3] = EQD;
        }
        treeData4.asset[3] = (void *)&assetEq;

        Hyb4_ConfigureLatticeCode(&treeData4, &latticeProg);

        range.reset(new TreeSliceRatesCompact::Range(
                    treeData4.MaxIndex[0],
                    treeData4.MaxIndex[1],
                    treeData4.MaxIndex[2],
                    treeData4.MaxIndex[3]));
        
        mTpIdx = treeData4.NbTP;

        // allocate claimBank internal memory
        for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {

            Hyb4_CbkInit(&cb->second.zeroBank);
            int nbCrit = cb->second.critZeroMatDates.size();
            string curveName = cb->second.curveName;

            if (nbCrit > 0) {
                int curveDim;
                // claim bank may be either 1D (foreign) 2D (domestic) or 
                // 3D (domestic equivalent of discounting foreign payment)
                // (ie FX*foreign)
                // if dimension already set (this case will be the foreign discounting
                // claimbank) use it, otherwise determine dimension and assign to cb
                if (cb->second.dimension == -1) { 
                    if (getCurrIdx(curveName) == FOREIGN)
                         curveDim = 1;
                    else
                         curveDim = 2;
                    cb->second.dimension = curveDim;
                }
                else
                    curveDim = cb->second.dimension;

                if (Hyb4_CbkAlloc(&cb->second.zeroBank, nbCrit, &treeData4, curveDim) != SUCCESS)
                    throw ModelException("Failed to allocate claim bank memory for curve " +
                                         cb->second.curveName);
            }
        }

        if (Hyb4_Dev_Alloc(&devData, &treeData4) != SUCCESS)
            throw ModelException("Hyb3_Dev_Alloc failed: " + IrConverter::slog.pop());

        print();
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb4CB::retrieveFactor() {
    try{
        int i,j;

        // ??? remove legacyHyb4FXEQMode field and export RateTree::CBLegacyPricingMode directly but
        // for now keep external interface the same and set the CBLegacy variable to the value of
        // exported flag in this model
        CBLegacyPricingMode = legacyHyb4FXEQMode;
    
        if (!discYC)
            throw ModelException("FDModel::discYC field is not populated - must contain the top product's "
                                 "discount curve. Internal error in FDModel/product code");

        currency[DOMESTIC] = discYC->getCcy();

        if (domIRParams->curveToDiffuse->getCcy() != currency[DOMESTIC])
            throw ModelException("Currency of supplied domestic curve to diffuse " + 
                                 domIRParams->curveToDiffuse->getCcy() +
                                 " must be the same as domestic currency defined by the product " +
                                 currency[DOMESTIC]);

        if (domIRParams->curveToDiscount->getCcy() != currency[DOMESTIC])
            throw ModelException("Currency of supplied domestic curve to diffuse " + 
                                 domIRParams->curveToDiffuse->getCcy() +
                                 " must be the same as domestic currency defined by the product " +
                                 currency[DOMESTIC]);

        currency[FOREIGN] = forIRParams->curveToDiffuse->getCcy();

        if (forIRParams->curveToDiscount->getCcy() != currency[FOREIGN])
            throw ModelException("Currency of supplied foreign curve to discount " + 
                                 forIRParams->curveToDiscount->getCcy() +
                                 " must be the same as foreign currency defined by curveToDiffuse" +
                                 currency[FOREIGN]);

        Hyb4CB::init();
        treeData3.Ppy = ppy;
        treeData3.NbSigmaMax = maxStdDeviations;

        for (i=0; i<2; ++i) {
            treeData3.CvDisc[i] = 1;
            treeData3.CvDiff[i] = 0;
            treeData3.CvIdx1[i] = 1;
            treeData3.CvIdx2[i] = 2;
        }

        // by convention set 0=diff curve, 1=discount curve, 2=other
        curveName[FOREIGN][0] = forIRParams->curveToDiffuse->getName();
        IrConverter::to_T_CURVE(treeCurves[FOREIGN][0], forIRParams->curveToDiffuse.get(), false);
        curveName[FOREIGN][1] = forIRParams->curveToDiscount->getName();
        IrConverter::to_T_CURVE(treeCurves[FOREIGN][1], forIRParams->curveToDiscount.get(), false);

        curveName[DOMESTIC][0] = domIRParams->curveToDiffuse->getName();
        IrConverter::to_T_CURVE(treeCurves[DOMESTIC][0], domIRParams->curveToDiffuse.get(), false);
        curveName[DOMESTIC][1] = domIRParams->curveToDiscount->getName();
        IrConverter::to_T_CURVE(treeCurves[DOMESTIC][1], domIRParams->curveToDiscount.get(), false);

        // for now only allow 1 FX to be registered
        if (fxFactors.size() != 1)
            throw ModelException("Hyb4 only supports one FX factor");
 
        // for now only allow 1 EQ to be registered
        if (eqFactors.size() != 1)
            throw ModelException("Hyb4 only supports one EQ factor");
 
        // check types and number of factors defined by instrument
        int nbYC=0;

        bool thirdCurve[2];
        thirdCurve[0] = false;
        thirdCurve[1] = false;
        for (i = 0; i < factors.size(); i++) {

            IMarketFactorConstSP &mf = factors[i];
            const YieldCurve* yc = dynamic_cast<const YieldCurve*>(mf.get());

            if (yc) {
                // for now must be either Foreign or Domestic
                if (yc->getCcy() != currency[FOREIGN] &&
                    yc->getCcy() != currency[DOMESTIC])
                    throw ModelException("Yield curve defined in instrument of currency " +
                                         yc->getCcy() +
                                         ". Hyb3 is configured as foreign = " + currency[DOMESTIC] +
                                         " and domestic = " + currency[FOREIGN]);

                for (j = 0; j < 2; j++) {
                    if (yc->getName() == curveName[j][0] ||
                        yc->getName() == curveName[j][1])
                        continue;
                    else {
                        if (yc->getCcy() == currency[j]) {
                            // define third currency or error if already defined
                            if (thirdCurve[j])
                                throw ModelException("Hyb3 engine supports maximum of three yield curves/currency "
                                                     "and supplied yield curve/factor " + yc->getName() +
                                                     "is a fourth for currency " + yc->getCcy());
                            curveName[j][2] = yc->getName();
                            IrConverter::to_T_CURVE(treeCurves[j][2], yc, false);
                            thirdCurve[j] = true;
                        }
                    }
                }
                nbYC++;
            }
            else {
                throw ModelException("Hyb4CB only supports yield curve market factors - not type " +
                                     mf->getClass()->getName());
            }
        }
        // if not supplied, set the third curve to equal the second/discount
        for (j = 0; j < 2; j++) {
            if (thirdCurve[j] == false) {
                curveName[j][2] = curveName[j][1];
                treeCurves[j][2] = treeCurves[j][1];
            }
        }

        // check curves all have the same valueDate for each currency
        for (i = 0; i < 2; i++) {
            long tmpValDate = treeCurves[i][0].ValueDate;
            if (legacyHyb4FXEQMode) {
                currencyValueDate[i] = getToday();
            }
            else {
                currencyValueDate[i] = DateTime::fromIrDate(tmpValDate);
            }
            for (j = 1; j < 3; j++) {
                if (treeCurves[i][j].ValueDate != tmpValDate)
                    throw ModelException("ValueDate " + 
                        DateTime::fromIrDate(treeCurves[0][i].ValueDate).toString() +
                        " defined in zero curve index [" + Format::toString(j) + 
                        "] is not the same as the valueDate defined in the other "
                        "zero curve(s) " + Format::toString((int)tmpValDate) +  
                        " for curve" + curveName[i][j]);
            }
        }

        // only FX and EQ and IR indexSpecs supported
        for (i = 0; i < indexSpecs.size(); i++) {
            string type = indexSpecs[i]->getClass()->getName();
            if (type != "IndexSpecIR" && 
                type != "IndexSpecFX" &&
                type != "IndexSpecEQ") {

                throw ModelException("Index spec type " + type + " is not supported by "
                                     "HybTree model - must be either FX or EQ or IR");
            }
        }

        const FXAsset* fx = dynamic_cast<FXAsset*>(fxFactors[0].get());
        if (!fx)
            throw ModelException("Internal error - FX market factor "
                                 "can not be cast to FXAsset");

        FXVOLATILITY_DATA fxVol;

        // populate foreign dimension with domesic values for
        populateTreeIRParams(FOREIGN);
        populateTreeIRParams(DOMESTIC); 

        strcpy(treeData3.Index[FOREIGN], forIRParams->volCalibIndex.c_str());
        strcpy(treeData3.Index[DOMESTIC], domIRParams->volCalibIndex.c_str());

        IrConverter::to_FXVOLATILITY_DATA(fxVol, fx);

        // populate equity vols and other equity data into tree

        const SimpleEquity* eq = dynamic_cast<SimpleEquity*>(eqFactors[0].get());

        double  value = 0;
        double *eqSpotVolOverride = NULL;
        
        if (EQSpotVolOverride.get())
        {
            value = EQSpotVolOverride->doubleValue() * 100;
            eqSpotVolOverride = &value;
        }

        PopulateEQData(&eqData,
                       today.toIrDate(),
                       eqSpotVolOverride,
                       eq,
                       NULL, // EQCalibrationStrike
                       dividendList);

        CORRELATION_DATA correlations;

        correlations.CorrIR      = corrIR;
        correlations.CorrForIRFX = corrForIRFX;
        correlations.CorrDomIRFX = corrDomIRFX;
        correlations.CorrForIREQ = corrForIREQ;
        correlations.CorrDomIREQ = corrDomIREQ;
        correlations.CorrFXEQ    = corrFXEQ;

        // ??? for now OwriteCorrel strings are set to nill, but may be used when model
        // supplies correlation overwrites

        char OwriteCorrel[6][MAXBUFF];
        strcpy(OwriteCorrel[0], "nil");
        strcpy(OwriteCorrel[1], "nil");
        strcpy(OwriteCorrel[2], "nil");
        strcpy(OwriteCorrel[3], "nil");
        strcpy(OwriteCorrel[4], "nil");
        strcpy(OwriteCorrel[5], "nil");

        Hyb4CorrelOverwrites(&fxData,
                             &eqData,
                             OwriteCorrel,
                             &correlations);
        
        char OwriteFxSpot[MAXBUFF];

        if (FXSpotVolOverride.get()) {  //  if supplied
            double value = FXSpotVolOverride->doubleValue();
            sprintf(OwriteFxSpot, "%lf", value * 100);
        }
        else
            strcpy(OwriteFxSpot, "nil");

        // ??? reuse this for now, but rewrite later on
        Hyb4FXSmileOverwrites(&fxData,
                              OwriteFxSpot,
                              &fxVol);
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}

// ??? check different modes
void Hyb4CB::sliceExpand(TreeSlice& treeSlice, const string& curveName) const {
    try {
        TreeSliceRatesCompact& slice = dynamic_cast<TreeSliceRatesCompact&>(treeSlice);

        int currIdx = getCurrIdx(curveName);

        // discounting a foreign product owned slice will always be in 3D mode
        if (currIdx == FOREIGN) {
            slice.expand(3);
            return;
        }
        if (slice.getDim() == 3)
            return;

        slice.expand(2);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb4CB::sliceDev(TreeSlice& treeSlice, int curveIdx) const {

    TreeSliceRatesCompact& slice = dynamic_cast<TreeSliceRatesCompact&>(treeSlice);
    try {
        string sliceName = slice.name;

        // nothing to do at the back of the tree 
        if (mTpIdx == treeData4.NbTP)
            return;

        // curveIdx == crvIdx * 10 + currIdx, see getCurveIdx
        int currIdx = curveIdx % 10;
        int crvIdx = curveIdx / 10;

        int devMode;
        int sliceDim = slice.getDim();

        if (currIdx == FOREIGN) {
            if (sliceDim == 1)
                throw ModelException("Hyb4 does not currently support discounting a "
                                     "foreign denominated curve - product must convert "
                                     "the foreign value to the domestic equivalent at the "
                                     "appropriate time and DEV in domestic only ");
            else if (sliceDim == 2)
                throw ModelException("Hyb4 does not currently support discounting a " +
                                    Format::toString(sliceDim) + " dimension Foreign currency "
                                    "denominated slice.  This usually arises when a foreign "
                                    "denominated component is fixing/reseting on a domestic "
                                    "index/rate (which Hyb4 considers a reverse cups).");
            else
                devMode = DISC_3D_CUPS;  // domestic eqivalent
        }
        else if (sliceDim == 4) {
            devMode = DISC_4D_CUPS;
        }
        else if (sliceDim == 3) {
            devMode = DISC_3D_CUPS;
        }
        else {
            devMode = DISC_2D_CUPS;
        }

        if (Hyb4_Dev(slice.getValuePtr(),
                     mTpIdx, 
                     treeData4.NbTP, 
                     crvIdx, 
                     devMode, 
                     const_cast<HYB4_DEV_DATA*>(&devData),
                     const_cast<HYB4_TREE_DATA*>(&treeData4)) != SUCCESS) {
            throw ModelException(IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, " slice name = " + slice.name);
    }
}


DateTime Hyb4CB::getCurveValueDate(string curveName) const {
    try {
        if (curveName.empty())
            throw ModelException(__FUNCTION__, "curveName not supplied");
        else {
            int currIdx = getCurrIdx(curveName);
            return currencyValueDate[currIdx];
        }
        return DateTime();
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb4CB::populateTreeIRParams(Currency ccy) {
    try {
        CcyIRParamsSP modelIRParams;
        IRVolRawSP volSP;

        if (ccy == FOREIGN) {   
            modelIRParams = forIRParams;
        }
        else {
            modelIRParams = domIRParams;
        }
        volSP = getIRVolRaw(modelIRParams);

        IRVolSelector volSelector(volSP, treeCurves[ccy][0],
                                  modelIRParams->volCalibIndex,
                                  cetSmoothing, modelIRParams);
        volSelector.getMktVolData(mktVolData[ccy]);

        // ??? temp "intermediate" structure to reuse logic from Hyb4 for now
        RateTree::IRModelParams mp;

        // retrieve model parameters from market cache
        IRExoticParamTable* pSmileTable = modelIRParams->smileTable.get();
        IRExoticParamTable* pModelTable = modelIRParams->modelTable.get();
        if (pSmileTable || pModelTable)
        {
            // Use new MarketTable object
            IrConverter::to_MODELPARAMETERS_DATA(mp, 
                                                modelIRParams->smileSet, 
                                                modelIRParams->modelSet,
                                                engineSet,
                                                *(engineTable.get()),
                                                *pSmileTable,
                                                *pModelTable);
        }
        else
        {
            // Use deprecated IRCalib object 
            IrConverter::to_MODELPARAMETERS_DATA(mp,
                modelIRParams->smileSet, 
                modelIRParams->modelSet,
                *(modelIRParams->irCalib.get()));
        }

        // populate tree structures from the model overwrite parameters or
        // underlying model parameters in the market data */
        populateIRParams(&mktVolData[ccy], &treeData3, &mp);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb4CB::update(int t, FDProduct::UpdateType type) {
    try {
        map<string, NamedClaimBank>::iterator cb;

        if (!treeBuilt)
            throw ModelException("Not able to update tree as the "
                                "tree has not been built");

        // set current range

        range->limits.bot1 = treeData4.iMin[t];
        range->limits.top1 = treeData4.iMax[t];
        range->limits.bot2 = treeData4.jMin[t];
        range->limits.top2 = treeData4.jMax[t];
        range->limits.bot3 = treeData4.kMin[t];
        range->limits.top3 = treeData4.kMax[t];
        range->limits.bot4 = treeData4.LMin[t];
        range->limits.top4 = treeData4.LMax[t];

        range->offsets.offset1 = treeData4.NodeOffset0[t];
        range->offsets.offset2 = treeData4.NodeOffset1[t];
        range->offsets.offset3 = treeData4.NodeOffset2[t];
        range->offsets.offset4 = treeData4.NodeOffset3[t];

        range->treeStep = t;

        int T = treeData4.NbTP;

        // Update tree to current time slice
        if (Hyb4_Lattice(&devData,
                         &treeData4,
                         &latticeProg,
                         t,
                         T) != SUCCESS)
            throw ModelException("Hyb4_Lattice Failed, " + IrConverter::slog.pop());

        mTpIdx = t;

        // update the registerd claim banks
        for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {

            int  addFlag = 0;
            long useDate = 0;
            long matDate = 0;
            int  idxActive = cb->second.critDatesIdx;
            string curveName = cb->second.curveName;  //  only used in error messages
            long currentDate = getDate(mTpIdx).toIrDate();

            int zbDevMode;
            // foreign currency discounter transformed to domestic equivalent
            if (cb->second.foreignDiscountedAsDomestic) {

                // ??? maybe remove these extra safety checks later on
                if ((treeData3.TreeType != TTYPE_FX2IR)    &&
                    (treeData3.TreeType != TTYPE_EQDFX2IR) &&
                    (treeData3.TreeType != TTYPE_EQFFX2IR))
                    throw ModelException("Internal consistency error - Hyb4 tree mode "
                        "must be TTYPE_FX2IR when discounting in a foreign currency");

                if (cb->second.dimension != 3)
                    throw ModelException("Internal consistency error - claim bank must be "
                        "three dimension slice when discounting in a foreign currency - "
                        "claimbank name " + cb->first + ", curveName " + curveName +
                        " dimension = " + Format::toString(cb->second.dimension));

                zbDevMode = DISC_3D_CUPS;

                if (cb->second.critZeroMatDates.size() == 0)
                    continue;

                if (idxActive >= 0 &&
                    cb->second.critZeroMatDates[idxActive] >= currentDate) {

                    addFlag = 1;
                    useDate = cb->second.critZeroUseDates[idxActive];
                    matDate = cb->second.critZeroMatDates[idxActive];

                    cb->second.critDatesIdx--;
                }

                // when discounting in domestic equivalent for a foreign curve, it is
                // assumed that we are always using the tree's domestic internal discount 
               // curve, which is always defined in Hyb4 as T_CURVE[DOM][1]
                int curveIdx = 1; 

                if (Hyb4_FXZbkUpdate(&cb->second.zeroBank,
                                     addFlag,
                                     currentDate,
                                     useDate,
                                     mTpIdx,
                                     treeData4.NbTP,
                                     curveIdx,
                                     zbDevMode,
                                     &devData,
                                     &treeData4) != SUCCESS)
                    throw ModelException("Zero bank update failed for curve " + curveName + 
                                        " in DEV mode number " + Format::toString(zbDevMode) +
                                        ", " + IrConverter::slog.pop());
            }
            else {

                int currIdx = getCurrIdx(curveName);
                if (currIdx == FOREIGN)  // foreign currency rate index 
                    zbDevMode = DISC_1D_NOCUPS;
                else  // domestic currency index or discounter
                    zbDevMode = DISC_2D_CUPS;

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

                if (Hyb4_ZbkUpdate(&cb->second.zeroBank,
                                   addFlag,
                                   currentDate,
                                   useDate,
                                   mTpIdx,
                                   treeData4.NbTP,
                                   curveIdx,
                                   zbDevMode,
                                   &devData,
                                   &treeData4) != SUCCESS)
                    throw ModelException("Zero bank update failed for curve " + curveName + 
                                        " in DEV mode number " + Format::toString(zbDevMode) +
                                        ", " + IrConverter::slog.pop());
            }
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb4CB::print()
{
    int i, j, k;

    map<string, NamedClaimBank>::iterator cb;

#if 0
    {
        FILE *stream = fopen("b.txt", "wb");

        fprintf(stream, "\nMKTVOL_DATA structure for Foreign IR\n\n");
        MktVol_PrintStructure(stream, &mktVolData[FOREIGN]);
        fprintf(stream, "\nMKTVOL_DATA structure for Domestic IR\n\n");
        MktVol_PrintStructure(stream, &mktVolData[DOMESTIC]);

        fprintf(stream, "\nTREE_DATA structure\n\n");
        PrintTree(stream, treeData3, 0);
        fclose(stream);

    }
#endif

    if (treeDataDebugFile.empty())
        return;

    FILE *stream = fopen(treeDataDebugFile.c_str(), "w");
    if (stream == NULL) {
        throw ModelException(__FUNCTION__, "Unable to open file " +
                             treeDataDebugFile + " for writing");
    }
    int count=0;

    fprintf(stream, "Model Settings:\n");
    fprintf(stream, "TodayDate: %s\n", today.toString().c_str());
    fprintf(stream, "valueDate: %s\n", valueDate.toString().c_str());
    fprintf(stream, "foreign currency value/spot date: %s\n", 
            currencyValueDate[FOREIGN].toString().c_str());
    fprintf(stream, "domestic currency value/spot date: %s\n", 
            currencyValueDate[DOMESTIC].toString().c_str());

    fprintf(stream, "\nFollowing critical dates registered with the engine:\n\n");
    for (i = 0; i < (int)sortedCritDates.size(); i++) {
        fprintf(stream, "%d   %ld\n", i, sortedCritDates[i]);
    }

    fprintf(stream, "\nFollowing claim banks registered with the engine:\n\n");

    for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {
        fprintf(stream, "%d, name: %s, curveName: %s, critDatesIdx %d\n",
                count,
                cb->first.c_str(),
                cb->second.curveName.c_str(),
                cb->second.critDatesIdx);

        fprintf(stream, "\n\tCritical Zero Dates used by tree\n");
        fprintf(stream, "\t[idx]  (zero Use/Obs Date)   (zero Maturity Date)\n");

        for (i = 0; i < (int)cb->second.critZeroMatDates.size(); i++) {
            fprintf(stream, "\t%d   %ld     %ld\n",
                    i,
                    cb->second.critZeroUseDates[i],
                    cb->second.critZeroMatDates[i]);
        }

        fprintf(stream, "\n\tOptional Zero Dates\n");
        fprintf(stream, "\t[idx]  (zero Use/Obs Date)   (zero Maturity Date)\n");

        for (i = 0; i < (int)cb->second.optZeroMatDates.size(); i++) {
            fprintf(stream, "\t%d   %ld     %ld\n",
                    i,
                    cb->second.optZeroUseDates[i],
                    cb->second.optZeroMatDates[i]);
        }
        fprintf(stream, "\n\tUnsorted Critical Zero Dates (as registered by products "
                        "before internal sorting and duplicate removal)\n");
        fprintf(stream, "\t[idx]  (zero Use/Obs Date)   (zero Maturity Date)    (label)\n");

        for (i = 0; i < (int)cb->second.critZeroMatDatesUnsorted.size(); i++) {
            fprintf(stream, "\t%d   %ld     %ld    %s\n",
                    i,
                    cb->second.critZeroUseDatesUnsorted[i],
                    cb->second.critZeroMatDatesUnsorted[i],
                    cb->second.critZeroLabel[i].c_str());
        }
    }
    
    fprintf(stream, "\nForeign Currency %s\n", currency[FOREIGN].c_str());
    fprintf(stream, "Domestic Currency %s\n\n", currency[DOMESTIC].c_str());
    fprintf(stream, "Nb  ForeignCurves  DomesticCurves\n");
    for (i = 0; i < 3; i++) {
        fprintf(stream, "%d  %s  %s\n", 
                i,
                curveName[0][i].c_str(),
                curveName[1][i].c_str());
    }

    fprintf (stream, "\n\n");

    // Slice sizes used for allocation
    fprintf(stream, "  W1   HW1    W2   HW2    W3   HW3\n");
    fprintf(stream, "%5d %5d %5d %5d %5d %5d\n\n\n",
                    treeData3.Width[0], treeData3.HalfWidth[0],
                    treeData3.Width[1], treeData3.HalfWidth[1],
                    treeData3.Width[2], treeData3.HalfWidth[2]);

    vector<string> currIdent(2);
    currIdent[0] = "FOREIGN";
    currIdent[1] = "DOMESTIC";

    for (i = 0; i < 2; i++) {
    
        fprintf (stream,"\n\n %s TIMELINE INFORMATION (IR)\n", currIdent[i].c_str());

        fprintf (stream,
                "Node      Date     Days  Max  Forward   "
                "Zero0   Discount0   Zero1   Discount1   "
                "Zero2   Discount2    IrMidNode    SpotVol \n");

        for (j = 0; j <= treeData3.NbTP; j++) {

            double Forward = 100. * treeData3.FwdRate[i][0][j]/treeData3.Length[j];

            int daysL = Daysact (treeData3.TPDate[0], treeData3.TPDate[j]);

            fprintf(stream,
                    "[%3d] \t%8ld  %3d  %3d  %7.4f  %7.4f   %8.6f  "
                    "%7.4f   %8.6f  %7.4f   %8.6f    %9.6f    %6.2f \n",
                    j,
                    treeData3.TPDate[j],
                    daysL,
                    treeData3.Top1[j],
                    Forward,
                    treeData3.ZeroRate[i][0][j] * 100.,
                    treeData3.ZeroCoupon[i][0][j],
                    treeData3.ZeroRate[i][1][j] * 100.,
                    treeData3.ZeroCoupon[i][1][j],
                    treeData3.ZeroRate[i][2][j] * 100.,
                    treeData3.ZeroCoupon[i][2][j],
                    exp(treeData3.IrZCenter[i][j]),
                    treeData3.SpotVol[i][j] * 100.);
            fflush(stream);
        }
    }

    // print out FX data
    fprintf (stream,"\nTIMELINE INFORMATION (FX)\n");
        
    fprintf (stream, "Node   Date      Forward    FxVol   "
                     " Spotvol   Rho's (Irf-Ird, Irf-FX, Ird-FX)\n");
                                        
    for (i = 0; i <= treeData3.NbTP; i++) {
        fprintf (stream, "[%3d]  %8ld  %12.6f    "
                        "%5.2f    %5.2f       %4.2f       %4.2f      %4.2f\n",
                        i, 
                        treeData3.TPDate[i],
                        treeData3.FwdFx[i],
                        treeData3.FxVol[i] * 100,
                        treeData3.SpotFxVol[i] * 100.,
                        treeData3.Rho[0][i],
                        treeData3.Rho[1][i],
                        treeData3.Rho[2][i]);
    }

    // Print out input zero curves
    for (i = 0; i < 2; i++) {
    
        fprintf (stream, "\n\n");
        fprintf (stream, "%s (%s)  Currency Curves (Index, COF, Risk)\n\n", 
                currIdent[i].c_str(), currency[i].c_str());

        for (k = 0; k < 3; k++) {

            fprintf (stream, "Maturity      Zero       Discount \n");

            for (j = 0; j < (treeCurves[i][k]).NbZero; j++) {

                double days = Daysact((treeCurves[i][k]).ValueDate, 
                                     (treeCurves[i][k]).ZeroDate[j]);
                // Discount factor up to the current date
                double discount = ::pow(1. + (treeCurves[i][k]).Zero[j], -days/365.);

                fprintf (stream,
                        "%ld   %9.6f     %8.6f \n",
                        (treeCurves[i][k]).ZeroDate[j],
                        (treeCurves[i][k]).Zero[j] * 100.,
                        discount);
            }
            fprintf (stream, "\n");
        }
    }

    fclose(stream);
}

// assuming call option
static double deltaToStrike(double delta, double vol, double fwdFX)
{
    double d1 = Normal_InvH(delta);
    double dummy = - d1 * vol + 0.5 * vol * vol;
    double impliedStrike = fwdFX * exp(dummy);

    return impliedStrike;
}

// ??? figure out how to turn this on/off but preferably not defining a full
// outputRequest for each one - maybe have a "debug info" output request
// that switches on all of this debug packet info
void Hyb4CB::recordOutput(Control* ctrl, Results* results) const {
    try {
        int i, j;
        
        if (!ctrl->requestsOutput(OutputRequest::DBG))
            return; // no debug info requested

        // record today and valueDates
        DateTimeSP todaySP(new DateTime(today));
        results->storeGreek(todaySP,
                            Results::DEBUG_PACKET,
                            OutputNameConstSP(new OutputName("TODAY")));

        DateTimeSP valueDateSP(new DateTime(valueDate));
        results->storeGreek(valueDateSP,
                            Results::DEBUG_PACKET,
                            OutputNameConstSP(new OutputName("VALUE_DATE")));

        DateTimeSP foreignValueDateSP(new DateTime(currencyValueDate[FOREIGN]));
        results->storeGreek(foreignValueDateSP,
                            Results::DEBUG_PACKET,
                            OutputNameConstSP(new OutputName("FOREIGN_IR_SPOT_DATE")));

        DateTimeSP domesticValueDateSP(new DateTime(currencyValueDate[DOMESTIC]));
        results->storeGreek(domesticValueDateSP,
                            Results::DEBUG_PACKET,
                            OutputNameConstSP(new OutputName("DOMESTIC_IR_SPOT_DATE")));

        // store tree critical dates
        DateTimeArraySP treeDates(new DateTimeArray(treeData3.NbTP+1));
        for (i = 0; i <= treeData3.NbTP; i++)
            (*treeDates)[i] = DateTime::fromIrDate(treeData3.TPDate[i]);

        results->storeGreek(treeDates,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("TREE_DATES")));

        // store number of days offsets from start of tree for each critical date
        IntArraySP daysOffsets(new IntArray(treeData3.NbTP+1));
        for (i = 0; i <= treeData3.NbTP; i++)
        {
            double days = Daysact(treeData4.TPDate[0], treeData3.TPDate[i]);
            (*daysOffsets)[i] = (int)days;
        }

        results->storeGreek(daysOffsets,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("TREE_DATE_OFFSETS")));

        // store FOREIGN ir rates
        CDoubleMatrixSP fgnZeros(new CDoubleMatrix(3, treeData3.NbTP+1));
        for (i = 0; i < 3; i++)
            for (j = 0; j <= treeData3.NbTP; j++)
                (*fgnZeros)[i][j] = treeData3.ZeroRate[0][i][j];

        results->storeGreek(fgnZeros,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FOREIGN_IR_ZEROS")));

        // store FOREIGN ir discount factors
        CDoubleMatrixSP fgnDiscFact(new CDoubleMatrix(3, treeData3.NbTP+1));
        for (i = 0; i < 3; i++)
            for (j = 0; j <= treeData3.NbTP; j++)
                (*fgnDiscFact)[i][j] = treeData3.ZeroCoupon[0][i][j];

        results->storeGreek(fgnDiscFact,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FOREIGN_IR_DISCOUNT_FACTORS")));

        // store DOMESTIC ir rates
        CDoubleMatrixSP domZeros(new CDoubleMatrix(3, treeData3.NbTP+1));
        for (i = 0; i < 3; i++)
            for (j = 0; j <= treeData3.NbTP; j++)
                (*domZeros)[i][j] = treeData3.ZeroRate[1][i][j];

        results->storeGreek(domZeros,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("DOMESTIC_IR_ZEROS")));

        // store DOMESTIC ir discount factors
        CDoubleMatrixSP domDiscFact(new CDoubleMatrix(3, treeData3.NbTP+1));
        for (i = 0; i < 3; i++)
            for (j = 0; j <= treeData3.NbTP; j++)
                (*domDiscFact)[i][j] = treeData3.ZeroCoupon[1][i][j];

        results->storeGreek(domDiscFact,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("DOMESTIC_IR_DISCOUNT_FACTORS")));

        // store foreign IR spotVols
        DoubleArraySP forIRSpotVol(new DoubleArray(treeData3.NbTP+1));
        for (i = 0; i <= treeData3.NbTP; i++)
            (*forIRSpotVol)[i] = treeData3.SpotVol[FOREIGN][i];

        results->storeGreek(forIRSpotVol,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FOREIGN_IR_SPOT_VOLS")));

        // store domestic IR spotVols
        DoubleArraySP domIRSpotVol(new DoubleArray(treeData3.NbTP+1));
        for (i = 0; i <= treeData3.NbTP; i++)
            (*domIRSpotVol)[i] = treeData3.SpotVol[DOMESTIC][i];

        results->storeGreek(domIRSpotVol,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("DOMESTIC_IR_SPOT_VOLS")));

        // store FX spot vols
        DoubleArraySP fxSpotVols(new DoubleArray(treeData3.NbTP+1));
        for (i = 0; i <= treeData3.NbTP; i++)
            (*fxSpotVols)[i] = treeData3.SpotFxVol[i];

        results->storeGreek(fxSpotVols,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FX_SPOT_VOLS")));

        // store FX composite vols
        DoubleArraySP fxCompVols(new DoubleArray(treeData3.NbTP+1));
        for (i = 0; i <= treeData3.NbTP; i++)
            (*fxCompVols)[i] = treeData3.FxVol[i];

        results->storeGreek(fxCompVols,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FX_COMP_VOLS")));

        // store forward FX 
        DoubleArraySP fxFwds(new DoubleArray(treeData3.NbTP+1));
        for (i = 0; i <= treeData3.NbTP; i++)
            (*fxFwds)[i] = treeData3.FwdFx[i];

        results->storeGreek(fxFwds,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FX_FWDS")));

        // Store local volatility function for various strikes
        CDoubleMatrixSP smileVols(new CDoubleMatrix(5, treeData3.NbTP+1));
        DoubleArraySP smileDeltas(new DoubleArray(5));
        (*smileDeltas)[0] = 0.1;
        (*smileDeltas)[1] = 0.25;
        (*smileDeltas)[2] = 0.5;
        (*smileDeltas)[3] = 0.75;
        (*smileDeltas)[4] = 0.9;

        for (i = 0; i <= treeData3.NbTP;i++)
        {
            double time = (*daysOffsets)[i] / 365.0;
            double vol = sqrt(time)* treeData3.FxVol[i]; // use comp Vol for estiamting delta
            double fxSpotVol = treeData3.SpotFxVol[i];   // spot Vol for local vol function

            for (j = 0; j < 5; j++)
            {
                double delta = (*smileDeltas)[j];
                double strikeFwdFx = deltaToStrike(delta, vol, treeData3.FwdFx[i]);
                double moneyness = strikeFwdFx/treeData3.FwdFx[i];
                double smile;

                if (Hyb3_Gfunc(&smile,
                            moneyness,
                            treeData3.A1C[i],
                            treeData3.A2C[i],
                            treeData3.A3C[i]) != SUCCESS)
                    throw ModelException("Hyb3_Gfunc failed: "+ IrConverter::slog.pop());

                smile *= fxSpotVol/moneyness;

                (*smileVols)[j][i] = smile;
            }
        }

        results->storeGreek(smileVols,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FX_SMILE_VOLS")));

        results->storeGreek(smileDeltas,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FX_SMILE_DELTAS")));
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void Hyb4CB::getFXIndex(TreeSlice& treeSlice) {
    try {
        if ((treeData3.TreeType != TTYPE_EQDFX2IR) &&
            (treeData3.TreeType != TTYPE_EQFFX2IR))
            throw ModelException("Hyb4 tree type must be in 4d mode to "
                "calculate an FXIndex - tree mode is currently set to " + 
                Format::toString(treeData3.TreeType));

        TreeSliceRatesCompact& slice = dynamic_cast<TreeSliceRatesCompact&>(treeSlice);

        if (slice.getDim() != 3)
            slice.allocDim(3);

        // for now simply copy the FX Spot slice into the provided slice
        if (Hyb4_CopySlice(slice.getValuePtr(),
                           devData.AssetValue2[0],
                           3,       // FX always 2D slice in this mode
                           mTpIdx,  // current step
                           &treeData4) != SUCCESS) {
            throw ModelException("Hyb4_CopySlice failed", IrConverter::slog.pop());
        }
        slice.treeStep = range->treeStep;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb4CB::getFXIndex(const IndexSpecFX& fxSpec, DateTime currentDate, DateTime resetDate,
                          TreeSlice& treeSlice) {
    try {
        // ??? We assume that the FX is always defined by the 
        // foreign and domestic discount curves - double check this here
        // to be doubly sure - remove this check later

/*
        ??? fxSpec.baseCurve/riskCurve are blank strings - not populated.  Don't know
            why they are not populated by the indexSpec - check if this should be the case
            or not
       if (fxSpec.baseCurve != curveName[DOMESTIC][1])
            throw ModelException("Model domestic (discount) curveName " + 
                                 curveName[DOMESTIC][1] + " must be the same curve "
                                 "as the base curve defined in the indexSpecFX " +
                                 fxSpec.baseCurve);
        if (fxSpec.riskCurve != curveName[FOREIGN][1])
            throw ModelException("Model foreign (discount) curveName " + 
                                 curveName[FOREIGN][1] + " must be the same curve "
                                 "as the risk curve defined in the indexSpecFX " +
                                 fxSpec.riskCurve);
*/
        TreeSliceRatesCompact& slice = dynamic_cast<TreeSliceRatesCompact&>(treeSlice);

        if (slice.getDim() != 3)
            slice.allocDim(3);

        // for now simply copy the FX Spot slice into the provided slice
        if (Hyb4_CopySlice(slice.getValuePtr(),
                           devData.AssetValue2[0],
                           3,       // FX always 3D slice in this mode
                           mTpIdx,  // current step
                           &treeData4) != SUCCESS) {
            throw ModelException("Hyb4_CopySlice failed", IrConverter::slog.pop());
        }

        // adjust for any date offsets
        dateAdjustFXIndex(fxSpec, currentDate, resetDate, treeSlice);
        slice.treeStep = range->treeStep;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb4CB::getEQIndex(TreeSlice& treeSlice) {
    try {
        if ((treeData3.TreeType != TTYPE_EQDFX2IR) &&
            (treeData3.TreeType != TTYPE_EQFFX2IR))
            throw ModelException("Hyb4 tree type must be in 4d mode to "
                "calculate an FXIndex - tree mode is currently set to " + 
                Format::toString(treeData3.TreeType));

        TreeSliceRatesCompact& slice = dynamic_cast<TreeSliceRatesCompact&>(treeSlice);

        if (slice.getDim() != 4)
            slice.allocDim(4);

        // for now simply copy the FX Spot slice into the provided slice
        if (Hyb4_CopySlice(slice.getValuePtr(),
                           devData.AssetValue3[0],
                           4,       // FX always 2D slice in this mode
                           mTpIdx,  // current step
                           &treeData4) != SUCCESS) {
            throw ModelException("Hyb4_CopySlice failed", IrConverter::slog.pop());
        }
        slice.treeStep = range->treeStep;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb4CB::getEQIndex(const IndexSpecEQ& eqSpec, DateTime currentDate, DateTime resetDate,
                          TreeSlice& treeSlice) {
    try {
        // ??? We assume that the FX is always defined by the 
        // foreign and domestic discount curves - double check this here
        // to be doubly sure - remove this check later

/*
        ??? fxSpec.baseCurve/riskCurve are blank strings - not populated.  Don't know
            why they are not populated by the indexSpec - check if this should be the case
            or not
       if (fxSpec.baseCurve != curveName[DOMESTIC][1])
            throw ModelException("Model domestic (discount) curveName " + 
                                 curveName[DOMESTIC][1] + " must be the same curve "
                                 "as the base curve defined in the indexSpecFX " +
                                 fxSpec.baseCurve);
        if (fxSpec.riskCurve != curveName[FOREIGN][1])
            throw ModelException("Model foreign (discount) curveName " + 
                                 curveName[FOREIGN][1] + " must be the same curve "
                                 "as the risk curve defined in the indexSpecFX " +
                                 fxSpec.riskCurve);
*/
        TreeSliceRatesCompact& slice = dynamic_cast<TreeSliceRatesCompact&>(treeSlice);

        if (slice.getDim() != 4)
            slice.allocDim(4);

        if (Hyb4_CopySlice(slice.getValuePtr(),
                           devData.AssetValue3[0],
                           4,
                           mTpIdx,
                           &treeData4) != SUCCESS) {
            throw ModelException("Hyb4_CopySlice failed", IrConverter::slog.pop());
        }

        // adjust for any date offsets
        dateAdjustEQIndex(eqSpec, currentDate, resetDate, treeSlice);
        slice.treeStep = range->treeStep;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb4CB::dateAdjustEQIndex(const IndexSpecEQ& fxSpec, DateTime currentDate, 
                               DateTime resetDate, TreeSlice& treeSlice) {}

// use deterministic FX ratio trick to adjust the stochastic FX slice
void Hyb4CB::dateAdjustFXIndex(const IndexSpecFX& fxSpec, DateTime currentDate, 
                                 DateTime resetDate, TreeSlice& treeSlice) {
    try {
        TreeSliceRatesCompact& slice = dynamic_cast<TreeSliceRatesCompact&>(treeSlice);

        if (slice.getDim() != 3)
            slice.allocDim(3);

        if (currentDate == resetDate)
           return;  // no adjustment required if they're the same dates
        else if (currentDate < resetDate)
            throw ModelException("Currently reset dates " + resetDate.toString() +
                                 " greater than the currentDate " + currentDate.toString() +
                                 " are not automatically adjusted");

        // if the reset date is before the currentDate when the rate is being asked for, 
        // approximate the rate by taking the ratio of the deterministic fwd FX ratios
        // between the actual reset date and the current date, then scale the stochastic 
        // FX value by this ratio.  Essentially this is matching the first
        // moment.  We do not at this stage try to make volatility adjustments
        int nbDaysDiff = currentDate.daysDiff(resetDate);
        if (nbDaysDiff > 180)
            throw ModelException("Current date " + currentDate.toString() + " to resetDate " +
                  resetDate.toString() + "FX approximate rate adjustment is only supported if "
                  "the date offset is < 180 days due to issues around scaling the rate over "
                  "too larger a time offset.  Requested offset of " + Format::toString(nbDaysDiff) +
                  "days for indexSpec " + fxSpec.name + " is too large ");

        double ratio;
        if (Hyb3_FX_Fwd_Ratio(&ratio, 
                              &treeCurves[DOMESTIC][1],
                              &treeCurves[FOREIGN][1],
                              resetDate.toIrDate(),
                              currentDate.toIrDate()) != SUCCESS)
            throw ModelException("Hyb3_FX_Fwd_Ratio failed: " + IrConverter::slog.pop());

        if (Hyb4_MultiplyScalar(slice.getValuePtr(),
                                3,
                                ratio,
                                mTpIdx,
                                &treeData4) != SUCCESS)
            throw ModelException("Hyb4_MultiplyScalar failed: " + IrConverter::slog.pop());
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


// Hyb4 eq mode always uses index [0] for domestic yield curve
int Hyb4CB::getCurrIdx(const string& curveName) const {
    
    if(curveName.empty())
        throw ModelException(__FUNCTION__, "curveName is not supplied");
    
    // loop through all curves, find a name match and double check that there
    // are no other duplicate curve names in the other currency
    // ??? no need to do this every time - so ensure setup is robust
    int i, j, index=-1;
    bool match=false;
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 3; j++) {
            if (curveName == this->curveName[i][j]) {
                if (match == true)
                    throw ModelException(__FUNCTION__, "Curve name " + curveName +
                                         "appears as both a domestic and foreign currency "
                                         "in Hyb4 tree registered curve names");                
                index = i;
                match=true;
                break;  // if match found check the other currency as well
            }
        }
    }

    if (index == -1)
        throw ModelException(__FUNCTION__, "Unable to find curveName " +
                             curveName + 
                             " in list of Hyb4 registered yield curve names");
    
    return index;
}


// interest rate curves are always one/first dimension
int Hyb4CB::getCrvDim(const string& curveName) const {
    try {
        int currIdx = getCurrIdx(curveName);
        if (currIdx == FOREIGN)
            return 1;
        else
            return 2;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// forward the price back to the valueDate as the convention is to report values to
// the valueDate
double Hyb4CB::getPrice0(const TreeSlice& price) const {
    try {
        double vdPrice;
        const TreeSliceLayer * layer = dynamic_cast< const TreeSliceLayer * >( &price );
        if (layer)
            vdPrice = layer->getPrice0();
        else
            vdPrice = price.getCentre();

        // ??? old method:
        // now forward to valueDate using precalculate discount factors between today/valueDate
        // tree is setup so curveToDiffuse = [1], which is what we are reporting the price in
        // so we use that curve's discountFactor to forward back to valueDate
        // vdPrice /= today2ValDateZeros[DOMESTIC][1];

        // ??? new method
        // as today=valueDate in the current implementation, no need to forward value
        // from today to valueDate
        return vdPrice;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// register product discounting slice requirements - for now the rule with
// pricing on Hyb4FX is that we don't want to discount any foreign denominated
// payments in foreign but instead convert to domestic equivalent by FX at payment date
// and discount in domestic from there.  To achieve this within the engine, we
// use a 3D domestic discounting slice, and essentially discount 1.0*FX(matDate) 
// in the claim bank
void Hyb4CB::registerZero(DateTime useDate, DateTime matDate, 
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

        string curveNameL = curveName;
        int currIdx = getCurrIdx(curveNameL);
        bool foreignAsDomestic = false;

        
        if (currIdx == FOREIGN) {
            // T_CURVE index 1 is always defined as the tree discount
            // curve for each currency.  FX is defined as DOM_disc/FOR_disc so
            // to transform foreign discount curve to domestic, the foreign curve
            // must be the tree discount curve, and we will then discount the result
            // in the domestic tree discount curve only.
            if (curveName != Hyb4CB::curveName[FOREIGN][1])
                throw ModelException("engine only supports discounting foreign currency "
                    "curve if foreign curve name " + curveName + " is the foreign discount "
                    "curve defined in the tree " + Hyb4CB::curveName[FOREIGN][1]);
            curveNameL = foreignZeroName(curveName);  // change curve name to internal name
            foreignAsDomestic = true;
        }

        cb = claimBank.find(curveNameL);

        // if claimBank doesn't exist, create a new one and insert into map
        if (cb == claimBank.end()) {
            NamedClaimBank newClaimBank;
            // store original curve name inside claim bank structure - for foreign discount
            // curve this will be different from the claim bank's key name in the map but
            // allows us to keep the original curve name for error messages etc.
            newClaimBank.curveName = curveName;  
            claimBank[curveNameL] = newClaimBank;
            cb = claimBank.find(curveNameL);
        }

        cb->second.foreignDiscountedAsDomestic = foreignAsDomestic;

        DateTime useDateL = max(today, useDate);
        cb->second.critZeroUseDates.push_back(useDateL.toIrDate());
        cb->second.critZeroMatDates.push_back(matDate.toIrDate());

        // foreign discounters are always 3D domstic curve equivalent
        if (foreignAsDomestic) {
            cb->second.critZeroLabel.push_back(string("FOREIGN DENOMINATED ZERO/DISCOUNT"));
            cb->second.dimension = 3;
        }
        else {
            cb->second.critZeroLabel.push_back(string("ZERO/DISCOUNT"));
            cb->second.dimension = 2;
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// retrieve zero/discounting slice from the tree
void Hyb4CB::getZero(DateTime useDate, DateTime matDate, 
                       string curveName, TreeSliceSP &slice) {
    try {
        int bankSliceDim;
        map<string, NamedClaimBank>::iterator cb;
        DateTime currentDate = getCurrentDate();
        string curveNameL = curveName;

        if (!dynamic_cast<TreeSliceRatesCompact*>(slice.get())) {
            slice = RateTree::createSimpleSlice(curveName);
            slice->name = "Zero "+useDate.toString()+" to "+matDate.toString();
        }
        TreeSliceRatesCompact& sliceL = dynamic_cast<TreeSliceRatesCompact&>(*slice);

        if (matDate <= valueDate)
            throw ModelException("Requested zero maturity date " + matDate.toString() +
                                 " must be greater than the model valueDate " +
                                 valueDate.toString());

        int currIdx = getCurrIdx(curveName);
        if (currIdx == FOREIGN) {
            curveNameL = foreignZeroName(curveName);
        }

        cb = claimBank.find(curveNameL);
        if (cb == claimBank.end())
            throw ModelException("Unable to find curveName " + curveName + " (which is "
                                 "represented internally as domestic discounting curveName " +
                                 curveNameL + ") in the list of model registered claimBanks");

        if (useDate != currentDate)
            throw ModelException("requested zero observation date " + useDate.toString() +
                                 " for curveName " + curveName + "must be the same as the model's "
                                 "current date " + currentDate.toString());

        CLAIM_BANK const* zeroBank = NULL;
        double* slicePtr = NULL;

        zeroBank = &cb->second.zeroBank;
        bankSliceDim = cb->second.dimension;
        bool foreignAsDomestic = cb->second.foreignDiscountedAsDomestic;

        // double check
        if (currIdx == FOREIGN && foreignAsDomestic == true) {
            if (bankSliceDim != 3)
                throw ModelException("Internal consistency error - claimBank curve " + curveNameL +
                                    " dimension must equal 3 when representing a foreign curve "
                                    " discounting slice - slice dimension currently set to " +
                                    Format::toString(bankSliceDim));
        }
        else {  // domestic must be 2D
            if (bankSliceDim != 2)
                throw ModelException("Internal consistency error - claimBank curve " + curveNameL +
                                    " dimension must equal 2 when representing a domestic curve "
                                    " discounting slice - slice dimension currently set to " +
                                    Format::toString(bankSliceDim));
        }

        if (currentDate == matDate) {
            if (currIdx == FOREIGN) {
                // whether native or domestic equivalent mode, value = FX at today
                sliceL.allocDim(3);
                slicePtr = sliceL.getValuePtr();
                if (Hyb4_CopySlice(slicePtr,
                                   devData.AssetValue2[0],
                                   3,
                                   mTpIdx,
                                   &treeData4) != SUCCESS)
                    throw ModelException("Hyb4_CopySlice failed: " + IrConverter::slog.pop());
            }
            else {  // domestic
                sliceL.allocDim(2);
                slicePtr = sliceL.getValuePtr();
                if (Hyb4_SetSlice(slicePtr,
                                  2,
                                  1.0,  // set slice value as 1.0
                                  mTpIdx,
                                  &treeData4) != SUCCESS)
                    throw ModelException("Hyb4_SetSlice failed: " + IrConverter::slog.pop());
            }
        }
        else {
            double* zeroBkPt = NULL;

            // pointer to internal zero slice maturity - zeroBkPt is read only and should not
            // be freed here
            zeroBkPt = Hyb4_ZbkReadZero((CLAIM_BANK*)zeroBank,
                                        matDate.toIrDate(),
                                        FALSE,  // ??? no interp, do we want to relax this
                                        currentDate.toIrDate(),
                                        mTpIdx,
                                        bankSliceDim,
                                        &treeData4);

            if (zeroBkPt == NULL)
                throw ModelException("Hyb3_ZbkReadZero failed: " + IrConverter::slog.pop());

            sliceL.allocDim(bankSliceDim);
            slicePtr = sliceL.getValuePtr();
            // deep copy results of ReadZero function call to output slice
            if (Hyb4_CopySlice(slicePtr,
                                zeroBkPt,
                                bankSliceDim,
                                mTpIdx,
                                &treeData4) != SUCCESS)
                throw ModelException("Hyb4_CopySlice failed: " + IrConverter::slog.pop());
        }

        slice->treeStep = range->treeStep;  // set internal slice timestep counter
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

string Hyb4CB::foreignZeroName(string curveName) {
    string name = curveName + "-FOREIGN_AS_DOMESTIC_ZERO";
    return name;
}

DRLIB_END_NAMESPACE
