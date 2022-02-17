//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Hyb3CBFX.hpp
//
//   Description : hyb3 domestic FX mode FD model class
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Hyb3CBFX.hpp"
#include "edginc/TreeSlice.hpp"
#include "edginc/IRVolPair.hpp"
#include "edginc/VolRequestRaw.hpp"
#include "edginc/IrConverter.hpp"
#include "edginc/IndexSpecIR.hpp"
#include "edginc/IndexSpecFX.hpp"
#include "edginc/Maths.hpp"
#include "edginc/SRMFXVol.hpp"
#include "esl_log.h"
#include "esl_market.h"


DRLIB_BEGIN_NAMESPACE

// helper functions to remove/integrate later on
static
void HYB3CorrelOverwrites(
    FX_DATA            *fx_data,                    // (O) Fx data
    char               OverWriteString[3][MAXBUFF], // (I) Overwrite strings
    CORRELATION_DATA   *correlation_data)
{
    try {
        double temp[3];
                 
        temp[0] = correlation_data->CorrIR;
        temp[1] = correlation_data->CorrDomIRFX;
        temp[2] = correlation_data->CorrForIRFX;


        /*********************************************/
        /* Although the fx_data structure supports a */
        /* term-structure of correlations, we only   */
        /* use the first element in Build Tree       */
        /*********************************************/
        fx_data->Rho[0][0] = temp[0];   /* Flat correlation curves */
        fx_data->Rho[1][0] = temp[2];
        fx_data->Rho[2][0] = temp[1];   /* Reverse here: Rho[1]=temp[2] and Rho[2]=temp[1] */
                    
        if (strcmp (OverWriteString[0], "nil"))                      /* Correlation overwrite */
        {                    
            fx_data->Rho[0][0] = atof (OverWriteString[0]);/* same comment as above*/
        }
              
        if (strcmp (OverWriteString[1], "nil"))
        {
            fx_data->Rho[1][0] = atof (OverWriteString[1]);          
        }
              
        if (strcmp (OverWriteString[2], "nil"))
        {                    
            fx_data->Rho[2][0] = atof (OverWriteString[2]);                               
        }
              

        if (Hyb3_Correl_Check_WType3 (fx_data) != SUCCESS)                  /* Check validity of input */
        {          
            throw ModelException(IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


static
void HYB3FXSmileOverwrites(
    FX_DATA            *fx_data,                 // (O) Fx data
    char               OverWriteString[MAXBUFF], // (I) FX Owrite str
    FXVOLATILITY_DATA  *Volatility_data,         // (I) FXVols file name
    FXSMILE_DATA       *Smile_data,              // (I) FX smile file name
    char               NbFXSmileOWS[MAXBUFF],    // (I)
    char               FXSmileParamOWS[MAXNBDATE][MAXBUFF])
{
    try {
        long
            Year[2],
            Month[2],
            Day[2];
        int
            i,
            readerror;     /* Reading error status  */

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
         
        // Read FX Smile data
        if (strstr(NbFXSmileOWS, "nil") != NULL)
        {

            // if smile data structure is NULL, turn off the FX smile by setting all
            // the smile parameters to 0.0
            if (Smile_data == NULL)
            {
                fx_data->NbSmilePt = 1;

                /* add arbitary number of days after the value date to ensure the smile date is
                > today as in some cases valueDate = Today. - chose to use 5 days but could have
                use any number */
                fx_data->SmileDate[0] = Nxtday(fx_data->ValueDate, 5);
                fx_data->a1[0] = 0.0;
                fx_data->a2[0] = 0.0;
                fx_data->a3[0] = 0.0;
            }
            else
            {
                fx_data->NbSmilePt = Smile_data->NbSmileDates;
                if ((fx_data->NbSmilePt < 0) ||
                    (fx_data->NbSmilePt > MAXNBDATE))
                {
                    throw ModelException("Fx_Input_W: Nb FX smile param lines out of range!");
                }

                for (i=0; i < fx_data->NbSmilePt; i++)
                {
                    fx_data->SmileDate[i] = Smile_data->SmileDates[i];
                    fx_data->a1[i] = Smile_data->A1[i];
                    fx_data->a2[i] = Smile_data->A2[i];
                    fx_data->a3[i] = Smile_data->A3[i];
                }  
            }
        }
        else
        {
            /* read the overwrite strings */
            readerror = sscanf (NbFXSmileOWS, 
                                "%d \n", 
                                &(fx_data->NbSmilePt));

            if (readerror != 1)
            {        
                throw ModelException("Fx_Input_W: Cannot read Nb FX Smile Param OWS");
            }

            if ((fx_data->NbSmilePt < 0) ||
                (fx_data->NbSmilePt > MAXNBDATE))
            {
                throw ModelException("Fx_Input_W: Nb FX smile param lines OWS out of range!");
            }

            for (i = 0; i < fx_data->NbSmilePt; i++)
            {
                readerror = sscanf (FXSmileParamOWS[i], 
                                    "%ld\t%lf\t%lf\t%lf \n",
                                    &(fx_data->SmileDate[i]),
                                    &(fx_data->a1[i]),
                                    &(fx_data->a2[i]),
                                    &(fx_data->a3[i]));

                if (readerror != 4)
                {          
                    throw ModelException("Fx_Input_W: Cannot read FX smile params OWS");
                }
            }
        }


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

 /********************
 *
 * Hyb3 FX Tree Mode
 *
**********************/


void Hyb3CBFX::load(CClassSP& clazz){

    clazz->setPublic();
    REGISTER(Hyb3CBFX, clazz);
    SUPERCLASS(Hyb3CB);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(ppy, "minimum number of tree nodes(points)/year");
    FIELD(maxStdDeviations, "number of standard deviations to trim the tree at");
    FIELD(forIRParams, "Collection of foreign interest rate model parameters");
    FIELD(domIRParams, "Collection of domestic interest rate model parameters");
    FIELD(FXSpotVolOverride,"Overrides market comp/spot vols for FX - sets constant spot vol for FX "
                            "diffusion (units = decimal NOT percentage)");
    FIELD_MAKE_OPTIONAL(FXSpotVolOverride);
    FIELD(legacyHyb3FXMode,"Override valueDate with today and do not extend yield curves back to today");
    FIELD_MAKE_OPTIONAL(legacyHyb3FXMode);
    
    // these fields are transient to allow tweaking to copy 
    // the model after ::getMarket() has been called 
    FIELD(fxIndexSpecSupplied, "")
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(fxIndexSpecSupplied);  // ??? 
    FIELD(corrIR, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrIR);
    FIELD(corrForIRFX, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrForIRFX);
    FIELD(corrDomIRFX, "");
    FIELD_MAKE_TRANSIENT_BUT_TWEAKABLE(corrDomIRFX);
}

CClassConstSP const Hyb3CBFX::TYPE = CClass::registerClassLoadMethod(
    "Hyb3CBFX", typeid(Hyb3CBFX), Hyb3CBFX::load);

// for class loading
bool Hyb3CBFXLoad(void) { return (Hyb3CBFX::TYPE != 0); }

Hyb3CBFX::~Hyb3CBFX() {
    clear();
}

// reset the underlying structures to NULL etc.
void Hyb3CBFX::reset() {

    // initialise structures
    MktVol_Init(&mktVolData[FOREIGN]);
    MktVol_Init(&mktVolData[DOMESTIC]); 

    Hyb3_Tree_Init(&treeData);
    Hyb3_Dev_Init(&devData);

    mTpIdx = -1;
    //    mSmoothFlag = 0;  ???

    // initialiase the number of factors of the rates-modes
    mktVolData[FOREIGN].NbFactor = forIRParams->nbFactors;
    mktVolData[DOMESTIC].NbFactor = domIRParams->nbFactors;

    treeBuilt = false;
}

// ??? AK: why can't we use the sliceDim of the ClaimBanks to clear the map?
void Hyb3CBFX::clear() {

    if (!treeBuilt)
        return;
    treeBuilt = false;

    map<string, NamedClaimBank>::iterator cb;
    for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {

        string curveName = cb->second.curveName;
        int curveDim;
        if (getCurrIdx(curveName) == FOREIGN)
            curveDim = mktVolData[FOREIGN].NbFactor;
        else
            curveDim = mktVolData[FOREIGN].NbFactor+mktVolData[DOMESTIC].NbFactor;

        Hyb3_CbkFree(&cb->second.zeroBank, &treeData, curveDim);
    }

    Hyb3_Dev_Free(&devData, &treeData);
    Hyb3_Tree_Free(&treeData);

    for (int i=0; i<NBCRV; ++i)
        today2ValDateZeros[0][i] = -1.0;

    // make sure that the data is reset to initial state
    reset();
}

IRGridPointAbsArraySP Hyb3CBFX::getSensitiveIRVolPoints(
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

void Hyb3CBFX::getMarket(const MarketData*  market, IInstrumentCollectionSP instruments) {
    try {
        today = market->GetReferenceDate();

        // ??? need final decision by business what to do here
        // set valueDate = today for now, but today may be SOD or EOD but we
        // want valueDate to be EOD no matter what time today is - this means
        // payments are always dropped anytime today
        valueDate = DateTime(today.getDate(), DateTime::END_OF_DAY_TIME);

        forIRParams->getData(this, market);
        domIRParams->getData(this, market);

        // in FX mode, must provide curve to discount
        if (forIRParams->curveToDiscount.isEmpty())
            throw ModelException("Must supply foreign curve to discount for Hyb3CBFX model");
        if (domIRParams->curveToDiscount.isEmpty())
            throw ModelException("Must supply domestic curve to discount for Hyb3CBFX model");

        // top level instrument must supply discount factor curve in order to define the
        // pricing/PnL currency of the trade.  This is retrieved through a complicated tangle
        // of code that eventually populates a marketDataFetcher by KComponent implementing the
        // IInstrumentCollection interface function:
        // virtual string discountYieldCurveName() const = 0
        // which returns the name of the discount curve defined in the instrument
        string discYCName = getDomesticYCName();
        if (discYCName.empty())
            throw ModelException("Instrument must supply a discount Curve");

        string modelDiscYCName = domIRParams->curveToDiscount->getName();
        if (discYCName != modelDiscYCName)
            throw ModelException("Name of curveToDiscount field set in Hyb3FX model (" +
                modelDiscYCName + ") must be same as defined on the instrument (" + 
                discYCName + ")");

        discYC = domIRParams->curveToDiscount.getSP();
        
        // Get Engine
        if (!engineTable.isEmpty())
            engineTable.getData(this, market);

        // recurse all the potential components in the instrument and look for the type
        // IMarketFactor::TYPE, which represent marketData type fields (yield curves, 
        // index Specs etc) and store them for retrieval later on
        class RetrieveFactors : public ObjectIteration::IAction {
            FDModel * model;
            const MarketData * market;
            IMarketFactorArray & factors;

        public:
            RetrieveFactors( FDModel * model, const MarketData * market, IMarketFactorArray & factors ) :
                model( model ), market( market ), factors( factors ) { factors.clear(); }

            virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {

                IMarketFactor* factor = dynamic_cast<IMarketFactor*>(state.getObject().get());
                string name = factor->getName();
                string type = factor->getClass()->getName();

                if (FXAsset::TYPE->isInstance(factor)) 
                    return true;

                int i;
                for (i = 0; i < factors.size(); ++i ) {
                    if(name == factors[i]->getName() && 
                        type == factors[i]->getClass()->getName())
                        break;
                }
                if (i >= factors.size() && model->acceptFactor(factor)) {
                    factors.push_back( IMarketFactorSP::attachToRef( factor ) );
                    return false;
                }
                else
                    return true;
            }
        };
        RetrieveFactors retrieveFactors( this, market, factors );
        ObjectIteration iteration( IMarketFactor::TYPE );
        iteration.recurse( retrieveFactors, instruments );

        fxIndexSpecSupplied = false;

        // recurse all the potential components in the instrument and look for the type
        // IndexSpec::TYPE, and store them for later retrieval.  Store the IndexSpecFX
        // types in a separate array to simplify logic of determining which mode the tree
        // should be configured for
        class RetrieveIndexSpecs : public ObjectIteration::IAction {

            IndexSpecArray & indexSpecs;
            IMarketFactorArray & fxFactors;
            bool & fxIndexSpecSupplied;

        public:
            RetrieveIndexSpecs(IndexSpecArray & indexSpecs, IMarketFactorArray & fxFactors, 
                bool & fxIndexSpecSupplied) :
                indexSpecs(indexSpecs), fxFactors(fxFactors), fxIndexSpecSupplied(fxIndexSpecSupplied)
                {indexSpecs.clear(); fxFactors.clear();}

            virtual bool invoke(ObjectIteration::State& state, IObjectConstSP obj) {

                int i;
                IndexSpec * indexSpec = dynamic_cast< IndexSpec * >( state.getObject().get() );

                string name = indexSpec->getName();
                string type = indexSpec->getClass()->getName();

                for (i = 0; i < indexSpecs.size(); ++i ) {
                    if( name == indexSpecs[i]->getName() && 
                        type == indexSpecs[i]->getClass()->getName())
                        break;
                }
                if (i >= indexSpecs.size()) {

                    indexSpecs.push_back( IndexSpecSP::attachToRef(indexSpec));
                    // if FX, get the FX market factor at the same time
                    IndexSpecFX * fxIndexSpec = dynamic_cast<IndexSpecFX*>(indexSpec);
                    if (fxIndexSpec && !fxIndexSpec->sameCurve) 
                    {
                        fxIndexSpecSupplied = true;

                        IMarketFactor *factor = fxIndexSpec->getFactor().get();
                        if (factor && !factor->getName().empty()) {
                            // test for duplicates
                            int j;
                            for (j = 0; j < fxFactors.size(); ++j) {
                                if( name == fxFactors[ j ]->getName())
                                break;
                            }
                            if( j >= fxFactors.size() )
                                fxFactors.push_back(fxIndexSpec->getFactor());
                        }
                    }
                    return false;
                }
                else
                    return true;
            }
        };
        RetrieveIndexSpecs retrieveIS (indexSpecs, fxFactors, fxIndexSpecSupplied);
        ObjectIteration retrieveISIter(IndexSpec::TYPE);
        retrieveISIter.recurse(retrieveIS, instruments);

        string domCurrency = discYC->getCcy();

        // FX between foreign and domestic discount curves must be supplied
        // ??? ignore for now any others, but check later if supplied that FX spot is the same
        string fxName = market->getFXName(forIRParams->curveToDiscount->getName(),
                                        domIRParams->curveToDiscount->getName());
        mFXName = fxName;
        MarketObjectSP fxObj = market->GetData(fxName, FXAsset::TYPE);

        FXAsset *fx = dynamic_cast<FXAsset*>(fxObj.get());
        if (!fx)
            throw ModelException("Internal error - Unable to cast an FXAsset::TYPE MarketFactor "
                                "to a FXAsset* pointer, which should never happen");

        fx->getMarket(this, market);

        IMarketFactor* fxFact = dynamic_cast<IMarketFactor*>(fxObj.get());
        string name = fxFact->getName();
        string type = fxFact->getClass()->getName();
        int j;
        for (j = 0; j < fxFactors.size(); ++j) {
            if( name == fxFactors[ j ]->getName())
            break;
        }
        if (j >= fxFactors.size())
            fxFactors.push_back( IMarketFactorSP::attachToRef( fxFact ) );

        // retrieve correlations, into structure of forIR-domIR, forIR-FX, domIR-FX
        double nan = sqrt(-1.0);
        corrIR = nan;
        corrForIRFX = nan;
        corrDomIRFX = nan;

        // forIR/domIR
        // check all yield curves of different currencies and model curves - must be same correlation
        // value or at least one provided
        string factor1, factor2;

        factor1 = forIRParams->curveToDiscount->getName();
        factor2 = domIRParams->curveToDiscount->getName();
        if (market->hasCorrelationData(factor1, factor2)) {

            string corrName = market->getCorrelationName(factor1, factor2);

            CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                        forIRParams->curveToDiscount->getClass(),
                                                        domIRParams->curveToDiscount->getClass(),
                                                        Correlation::TYPE,
                                                        market);
            corrIR = CorrelationSP::dynamicCast(corrBase)->getCorrelation(); 
        }

        factor1 = forIRParams->curveToDiffuse->getName();
        factor2 = domIRParams->curveToDiffuse->getName();
        if (market->hasCorrelationData(factor1, factor2)) {

            string corrName = market->getCorrelationName(factor1, factor2);

            CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                        forIRParams->curveToDiscount->getClass(),
                                                        domIRParams->curveToDiscount->getClass(),
                                                        Correlation::TYPE,
                                                        market);
            double result = CorrelationSP::dynamicCast(corrBase)->getCorrelation();
            if (!Maths::finite(corrIR))
                corrIR = result;
            else {
                if (!Maths::equals(corrIR, result))
                    throw ModelException("Correlation defined for curveToDiffuse factors " +
                                        factor1 + " and " + factor2 + "(value = " + 
                                        Format::toString(result) +
                                        ") does not equal previously calculated correlation for " +
                                        "model discount curves");
            }
        }
        
        // ??? more to do here to check all other correlations of yield curve to yield curve
        if (!Maths::finite(corrIR))
            throw ModelException("??? tmp error - correlation between forIR and domIR not found");


        // foreignIR/FX correlation
        // ??? just check curve to diffuse and curve to discount FXs for now, but have to check all
        // of them later

        string foreignIR = forIRParams->curveToDiscount->getName();
        string domesticIR = domIRParams->curveToDiscount->getName();

        if (market->hasFXData(foreignIR, domesticIR)) {
            string fxName = market->getFXName(foreignIR, domesticIR);

            if (market->hasCorrelationData(foreignIR, fxName)) {
                string corrName = market->getCorrelationName(foreignIR, fxName);

                CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                            forIRParams->curveToDiscount->getClass(),
                                                            FXAsset::TYPE,
                                                            Correlation::TYPE,
                                                            market);
                corrForIRFX = CorrelationSP::dynamicCast(corrBase)->getCorrelation();
            }
        }

        if (!Maths::finite(corrForIRFX))
            throw ModelException("Correlation not found in market for foreignIR ("
                                 + foreignIR + "), FX (" + fxName + ") on discount curve");

        // domesticIR/FX correlation
        // ??? just check curve to diffuse and curve to discount FXs for now, but have to check all
        // of them later

        foreignIR = forIRParams->curveToDiscount->getName();
        domesticIR = domIRParams->curveToDiscount->getName();

        if (market->hasFXData(foreignIR, domesticIR)) {

            string fxName = market->getFXName(foreignIR, domesticIR);
            if (market->hasCorrelationData(domesticIR, fxName)) {

                string corrName = market->getCorrelationName(domesticIR, fxName);

                CorrelationBaseSP corrBase = getCorrelation(corrName,
                                                            domIRParams->curveToDiscount->getClass(),
                                                            FXAsset::TYPE,
                                                            Correlation::TYPE,
                                                            market);
                corrDomIRFX = CorrelationSP::dynamicCast(corrBase)->getCorrelation();
            }
        }

        if (!Maths::finite(corrForIRFX))
            throw ModelException("Correlation not found in market for domesticIR ("+
                                 domesticIR+"), FX ("+fxName+") on discount curve");

    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb3CBFX::initModel() {

    try {
        treeBuilt = true; 
        treeData.FxMomentMatching = TRUE;

        CRIT_DATE* critDate=NULL;
        int nbCritDate = 0;
        IrConverter::AutoDeleteDRArray<CRIT_DATE>
            critDateDelete(&critDate, CRITDATE, &nbCritDate);

        // ??? note:  ZeroInterpTypeFlag is global variable defined in both hyb3.lib
        // and fix3.lib - tidy up but have to check all hyb3 products
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
                if (!legacyHyb3FXMode)
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
            range.reset(new TreeSliceRates::Range(0,0,0,0,0,0));
            treeData.NbTP = -1;
            return;
        }
        // construct the timeline 
        if (Hyb3_Time_Line(today.toIrDate(), nbCritDate, critDate, 'I', &treeData) 
            != SUCCESS)
            throw ModelException("Hyb3_Time_Line failed: " + IrConverter::slog.pop());

        // if printing out a debug file, turn on printing the default CET debug file also
        int printCET = FALSE;
        if (!treeDataDebugFile.empty())
            printCET = TRUE;
        if (Hyb3_Build_Tree(printCET, treeCurves, mktVolData, &fxData, &eqData, 
                            &treeData) != SUCCESS)
            throw ModelException("Hyb3_Build_Tree failed: " + IrConverter::slog.pop());

        // tree may be in 2D CUPS or 3D mode
        switch(treeData.TreeType) {

            case TTYPE_2IR: {
                range.reset(new TreeSliceRates::Range(
                            -treeData.HalfWidth[0],
                            treeData.HalfWidth[0],
                            -treeData.HalfWidth[1],
                            treeData.HalfWidth[1]));
                break;
            }
            case TTYPE_2IR2F1D: // Hyb2+1 mode
            case TTYPE_FX2IR: {
                range.reset(new TreeSliceRates::Range(
                            -treeData.HalfWidth[0],
                            treeData.HalfWidth[0],
                            -treeData.HalfWidth[1],
                            treeData.HalfWidth[1],
                            -treeData.HalfWidth[2],
                            treeData.HalfWidth[2]));
                break;
            }
            default:
                throw ModelException("Internal error - incorrect/unknown TreeType flag "
                                     "set in HYB3_TREE_DATA structure. Flag value = " +
                                     Format::toString(treeData.TreeType));
        }
        
        mTpIdx = treeData.NbTP;

        // allocate claimBank internal memory
        for (cb = claimBank.begin(); cb != claimBank.end(); cb++) {

            Hyb3_CbkInit(&cb->second.zeroBank);
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
                         curveDim = mktVolData[FOREIGN].NbFactor;
                    else
                         curveDim = mktVolData[FOREIGN].NbFactor + mktVolData[DOMESTIC].NbFactor;
                    cb->second.dimension = curveDim;
                }
                else
                    curveDim = cb->second.dimension;

                if (Hyb3_CbkAlloc(&cb->second.zeroBank, nbCrit, &treeData, curveDim) != SUCCESS)
                    throw ModelException("Failed to allocate claim bank memory for curve " +
                                         cb->second.curveName);
            }
        }

        if (Hyb3_Dev_Alloc(&devData, &treeData) != SUCCESS)
            throw ModelException("Hyb3_Dev_Alooc failed: " + IrConverter::slog.pop());

        print();
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


void Hyb3CBFX::retrieveFactor() {
    try{
        int i,j;

        // ??? remove legacyHyb3FXMode field and export RateTree::CBLegacyPricingMode directly but
        // for now keep external interface the same and set the CBLegacy variable to the value of
        // exported flag in this model
        CBLegacyPricingMode = legacyHyb3FXMode;
    
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

        // reset the data structures to initial values
        reset();

        treeData.Ppy = ppy;
        treeData.NbSigmaMax = maxStdDeviations;

        for (i=0; i<2; ++i) {
            treeData.CvDisc[i] = 1;
            treeData.CvDiff[i] = 0;
            treeData.CvIdx1[i] = 1;
            treeData.CvIdx2[i] = 2;
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
            throw ModelException("Hyb3 only supports one FX factor");
 
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
                throw ModelException("Hyb3BCFX only supports yield curve market factors - not type " +
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
            if (legacyHyb3FXMode) {
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

        // only FX and IR indexSpecs supported
        for (i = 0; i < indexSpecs.size(); i++) {
            string type = indexSpecs[i]->getClass()->getName();
            if (type != "IndexSpecIR" && 
                type != "IndexSpecFX") {

                throw ModelException("Index spec type " + type + " is not supported by "
                                     "Hyb3CBFX model - must be either FX or IR");
            }
        }

        const FXAsset* fx = dynamic_cast<FXAsset*>(fxFactors[0].get());
        if (!fx)
            throw ModelException("Internal error - FX market factor "
                                 "can not be cast to FXAsset");

        /*
        // ??? doesn't work as riskCcy has not got market - still just yield curve wrappers
        // dom ycs in fx assset must be domestic discount ycs
        string fxDomYCName = fx->getBaseCcy()->getName();
        if (fxDomYCName != domIRParams->curveToDiscount->getName())
            throw ModelException("Domestic yield curve " + fxDomYCName + " defined for FX Asset " + 
                                 fx->getName() +
                                 " by convention must be the model's domestic curveToDiscount yield curve " +
                                 domIRParams->curveToDiscount->getName());

        string fxForYCName = fx->getRiskCcy()->getName();
        if (fxForYCName != forIRParams->curveToDiscount->getName())
            throw ModelException("Foreign yield curve " + fxForYCName + " defined for FX Asset " + 
                                 fx->getName() +
                                 " by convention must be the model's foreign curveToDiscount yield curve " +
                                 forIRParams->curveToDiscount->getName());
        */
        bool fxMode = false;

        if (fxIndexSpecSupplied == true) {
            fxMode = true;
        }
        else  {
            vector<bool> factorInd(factors.size());
            for (i = 0; i < factors.size(); i++)
                factorInd[i] = false;

            for (i = 0; i < indexSpecs.size(); i++) {
                IndexSpecSP is = indexSpecs[i];
                const IndexSpecIR* ir = dynamic_cast<const IndexSpecIR*>(is.get());

                if (ir) {
                    int j;
                    // find factor associated with indexSpec and mark with indicator
                    for (j = 0; j < factors.size(); j++) {
                       string factorName = factors[j]->getName();
                       string irFactorName = ir->getFactor()->getName();
                       if (factorName == irFactorName) {
                            factorInd[j] = true;  // match found - factor derived from indexSpec
                            break;
                       }
                    }
                }
            }
            // now for all IR factors not marked as true (ie., derived from index spec), 
            // check currency and if any is foreign, must be in 3D mode
            for (i = 0; i < factors.size(); i++) {
                IMarketFactorConstSP &mf = factors[i];

                const YieldCurve* yc = dynamic_cast<const YieldCurve*>(mf.get());
                if (yc) {
                    if (factorInd[i] == false) {
                        string ccy=yc->getCcy();  // help debugging
                        if (ccy != discYC->getCcy()) {
                            fxMode = true;
                            break;
                        }
                    }
                }
            }
        }
        
        if (fxMode)
            treeData.TreeType = TTYPE_FX2IR;
        else
        {
            if ( mktVolData[FOREIGN].NbFactor == 1 ) 
                treeData.TreeType = TTYPE_2IR;
            else if (  mktVolData[FOREIGN].NbFactor == 2 )
                treeData.TreeType = TTYPE_2IR2F1D;
            else
                throw ModelException( "Internal error - more than two factors per IR asset not supported in Hyb3");
        }

        FXVOLATILITY_DATA fxVol;
        CORRELATION_DATA correlations;

        // only support a maximum of three distinct ycs, so check if more were supplied
        // ??? not sure if factors has repeats or these are sorted out before hand
        /*if (factors.size() > 3)
            throw ModelException("A maximum of three distinct " + discYC->getCcy() +
                                 " currency curves are supported by hyb3 - " +
                                 Format::toString(factors.size()) + " supplied"); */

        // populate foreign dimension with domesic values for
        populateTreeIRParams(FOREIGN);
        populateTreeIRParams(DOMESTIC); 

        strcpy(treeData.Index[FOREIGN], forIRParams->volCalibIndex.c_str());
        strcpy(treeData.Index[DOMESTIC], domIRParams->volCalibIndex.c_str());

        IrConverter::to_FXVOLATILITY_DATA(fxVol, fx);

        correlations.CorrIR = corrIR;
        correlations.CorrForIRFX = corrForIRFX;
        correlations.CorrDomIRFX = corrDomIRFX;


        // AK: do we need to overwrite corr-information here ???
        // currently hard coded information, needs to done properly
        mktVolData[FOREIGN].CorrSwapSt    = Nxtmth( today.toIrDate(), 120L, 1L);
        mktVolData[FOREIGN].CorrSwapMat   = Nxtmth(mktVolData[FOREIGN].CorrSwapSt,  120L, 1L);
        mktVolData[FOREIGN].CorrSwapDCC   = '3';
        mktVolData[FOREIGN].CorrSwapFreq  = 'A';
        mktVolData[DOMESTIC].CorrSwapSt    = Nxtmth( today.toIrDate(), 120L, 1L);
        mktVolData[DOMESTIC].CorrSwapMat   = Nxtmth(mktVolData[DOMESTIC].CorrSwapSt,  120L, 1L);
        mktVolData[DOMESTIC].CorrSwapDCC   = '3';
        mktVolData[DOMESTIC].CorrSwapFreq  = 'A';


        // ??? for now OwriteCorrel strings are set to nill, but may be used when model
        // supplies correlation overwrites
        char OwriteCorrel[3][MAXBUFF];
        strcpy(OwriteCorrel[0], "nil");
        strcpy(OwriteCorrel[1], "nil");
        strcpy(OwriteCorrel[2], "nil");

        HYB3CorrelOverwrites(&fxData,
                             OwriteCorrel,
                             &correlations);
        
        FXSMILE_DATA dummyFXSmile;
        IrConverter::to_FXSMILE_DATA(dummyFXSmile, fx);

        if (dummyFXSmile.NbSmileDates == 0)
            throw ModelException("At least one FX smile date and a1,a2,a3 data set must "
                                 "be provided.  To turn smile off set a1=a2=a3 = 0.0. "
                                 "Check that smileExpiry, smileA1... in "
                                 "FXAsset::SRMFX::Vol are not empty.");

        char OwriteFxSpot[MAXBUFF];
        char NbFXSmileOWS[MAXBUFF];
        char FXSmileParamOWS[MAXNBDATE][MAXBUFF];

        if (FXSpotVolOverride.get()) {  //  if supplied
            double value = FXSpotVolOverride->doubleValue();
            sprintf(OwriteFxSpot, "%lf", value * 100);
        }
        else
            strcpy(OwriteFxSpot, "nil");

        strcpy(NbFXSmileOWS, "nil");  // read smile params from dummyFXSmile structure

        // ??? reuse this for now, but rewrite later on
        HYB3FXSmileOverwrites(&fxData,
                              OwriteFxSpot,
                              &fxVol,
                              &dummyFXSmile,
                              NbFXSmileOWS,
                              FXSmileParamOWS);
    }
    catch (exception& e){
        throw ModelException(e, __FUNCTION__);
    }
}

// ??? check different modes
void Hyb3CBFX::sliceExpand(TreeSlice& treeSlice, const string& curveName) const {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

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

void Hyb3CBFX::sliceDev(TreeSlice& treeSlice, int curveIdx) const {

    TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);
    try {
        string sliceName = slice.name;

        // nothing to do at the back of the tree 
        if (mTpIdx == treeData.NbTP)
            return;

        // curveIdx == crvIdx * 10 + currIdx, see getCurveIdx
        int currIdx = curveIdx % 10;
        int crvIdx = curveIdx / 10;

        int devMode;
        int sliceDim = slice.getDim();

        if (currIdx == FOREIGN) {
            if (sliceDim == 1)
                throw ModelException("Hyb3 does not currently support discounting a "
                                     "foreign denominated curve - product must convert "
                                     "the foreign value to the domestic equivalent at the "
                                     "appropriate time and DEV in domestic only ");
            else if (sliceDim == 2) // ??? what about Hyb2+1
                throw ModelException("Hyb3 does not currently support discounting a " +
                                    Format::toString(sliceDim) + " dimension Foreign currency "
                                    "denominated slice.  This usually arises when a foreign "
                                    "denominated component is fixing/reseting on a domestic "
                                    "index/rate (which Hyb3 considers a reverse cups).");
            else
                devMode = DISC_3D_CUPS;  // domestic eqivalent
        }
        else if (sliceDim == 3) {
            if (treeData.TreeType == TTYPE_2IR2F1D)
                devMode = DISC_3D_2IR2F1D_CUPS;
            else
                devMode = DISC_3D_CUPS;
        }
        else {
            devMode = DISC_2D_CUPS;
        }

        if (Hyb3_Dev(slice.getValuePtr(),
                     mTpIdx, 
                     treeData.NbTP, 
                     crvIdx, 
                     devMode, 
                     const_cast<HYB3_DEV_DATA*>(&devData),
                     const_cast<HYB3_TREE_DATA*>(&treeData)) != SUCCESS) {
            throw ModelException(IrConverter::slog.pop());
        }
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__, " slice name = " + slice.name);
    }
}


DateTime Hyb3CBFX::getCurveValueDate(string curveName) const {
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

void Hyb3CBFX::populateTreeIRParams(Currency ccy) {
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

        // ??? temp "intermediate" structure to reuse logic from hyb3 for now
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
        populateIRParams(&mktVolData[ccy], &treeData, &mp);
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb3CBFX::update(int t, FDProduct::UpdateType type) {
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
            string curveName = cb->second.curveName;  //  only used in error messages
            long currentDate = getDate(mTpIdx).toIrDate();

            int zbDevMode;
            // foreign currency discounter transformed to domestic equivalent
            if (cb->second.foreignDiscountedAsDomestic) {

                // ??? maybe remove these extra safety checks later on
                if (treeData.TreeType != TTYPE_FX2IR)
                    throw ModelException("Internal consistency error - hyb3 tree mode "
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
               // curve, which is always defined in hyb3 as T_CURVE[DOM][1]
                int curveIdx = 1; 

                if (Hyb3_FXZbkUpdate(&cb->second.zeroBank,
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
            else {

                int currIdx = getCurrIdx(curveName);
                if (currIdx == FOREIGN)  // foreign currency rate index 
                {
                    if (mktVolData[FOREIGN].NbFactor == 1)
                        zbDevMode = DISC_1D_NOCUPS;
                    else if (mktVolData[FOREIGN].NbFactor == 2)
                        zbDevMode = DISC_2D_1IR2F_NOCUPS;
                    else
                        throw ModelException("foreign market data has more than 2 factors "
                                             "which is not currently supported in Hyb3CB");
                }
                else  // domestic currency index or discounter
                {
                    if (mktVolData[FOREIGN].NbFactor == 1)
                        zbDevMode = DISC_2D_CUPS;
                    else if (mktVolData[FOREIGN].NbFactor == 2)
                        zbDevMode = DISC_3D_CUPS;
                    else
                        throw ModelException("foreign market data has more than 2 factors "
                        "which is not currently supported in Hyb3CB");
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
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb3CBFX::print()
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
        PrintTree(stream, treeData, 0);
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
                    treeData.Width[0], treeData.HalfWidth[0],
                    treeData.Width[1], treeData.HalfWidth[1],
                    treeData.Width[2], treeData.HalfWidth[2]);

    vector<string> currIdent(2);
    currIdent[0] = "FOREIGN";
    currIdent[1] = "DOMESTIC";

    for (i = 0; i < 2; i++) {
    
        fprintf (stream,"\n\n %s TIMELINE INFORMATION (IR)\n", currIdent[i].c_str());

        fprintf (stream,
                "Node      Date     Days  Max  Forward   "
                "Zero0   Discount0   Zero1   Discount1   "
                "Zero2   Discount2    IrMidNode    SpotVol \n");

        for (j = 0; j <= treeData.NbTP; j++) {

            double Forward = 100. * treeData.FwdRate[i][0][j]/treeData.Length[j];

            int daysL = Daysact (treeData.TPDate[0], treeData.TPDate[j]);

            fprintf(stream,
                    "[%3d] \t%8ld  %3d  %3d  %7.4f  %7.4f   %8.6f  "
                    "%7.4f   %8.6f  %7.4f   %8.6f    %9.6f    %6.2f \n",
                    j,
                    treeData.TPDate[j],
                    daysL,
                    treeData.Top1[j],
                    Forward,
                    treeData.ZeroRate[i][0][j] * 100.,
                    treeData.ZeroCoupon[i][0][j],
                    treeData.ZeroRate[i][1][j] * 100.,
                    treeData.ZeroCoupon[i][1][j],
                    treeData.ZeroRate[i][2][j] * 100.,
                    treeData.ZeroCoupon[i][2][j],
                    exp(treeData.IrZCenter[i][j]),
                    treeData.SpotVol[i][j] * 100.);
            fflush(stream);
        }
    }

    // print out FX data
    fprintf (stream,"\nTIMELINE INFORMATION (FX)\n");
        
    fprintf (stream, "Node   Date      Forward    FxVol   "
                     " Spotvol   Rho's (Irf-Ird, Irf-FX, Ird-FX)\n");
                                        
    for (i = 0; i <= treeData.NbTP; i++) {
        fprintf (stream, "[%3d]  %8ld  %12.6f    "
                        "%5.2f    %5.2f       %4.2f       %4.2f      %4.2f\n",
                        i, 
                        treeData.TPDate[i],
                        treeData.FwdFx[i],
                        treeData.FxVol[i] * 100,
                        treeData.SpotFxVol[i] * 100.,
                        treeData.Rho[0][i],
                        treeData.Rho[1][i],
                        treeData.Rho[2][i]);
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
void Hyb3CBFX::recordOutput(Control* ctrl, Results* results) const {
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
        DateTimeArraySP treeDates(new DateTimeArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*treeDates)[i] = DateTime::fromIrDate(treeData.TPDate[i]);

        results->storeGreek(treeDates,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("TREE_DATES")));

        // store number of days offsets from start of tree for each critical date
        IntArraySP daysOffsets(new IntArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
        {
            double days = Daysact(treeData.TPDate[0], treeData.TPDate[i]);
            (*daysOffsets)[i] = (int)days;
        }

        results->storeGreek(daysOffsets,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("TREE_DATE_OFFSETS")));

        // store FOREIGN ir rates
        CDoubleMatrixSP fgnZeros(new CDoubleMatrix(3, treeData.NbTP+1));
        for (i = 0; i < 3; i++)
            for (j = 0; j <= treeData.NbTP; j++)
                (*fgnZeros)[i][j] = treeData.ZeroRate[0][i][j];

        results->storeGreek(fgnZeros,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FOREIGN_IR_ZEROS")));

        // store FOREIGN ir discount factors
        CDoubleMatrixSP fgnDiscFact(new CDoubleMatrix(3, treeData.NbTP+1));
        for (i = 0; i < 3; i++)
            for (j = 0; j <= treeData.NbTP; j++)
                (*fgnDiscFact)[i][j] = treeData.ZeroCoupon[0][i][j];

        results->storeGreek(fgnDiscFact,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FOREIGN_IR_DISCOUNT_FACTORS")));

        // store DOMESTIC ir rates
        CDoubleMatrixSP domZeros(new CDoubleMatrix(3, treeData.NbTP+1));
        for (i = 0; i < 3; i++)
            for (j = 0; j <= treeData.NbTP; j++)
                (*domZeros)[i][j] = treeData.ZeroRate[1][i][j];

        results->storeGreek(domZeros,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("DOMESTIC_IR_ZEROS")));

        // store DOMESTIC ir discount factors
        CDoubleMatrixSP domDiscFact(new CDoubleMatrix(3, treeData.NbTP+1));
        for (i = 0; i < 3; i++)
            for (j = 0; j <= treeData.NbTP; j++)
                (*domDiscFact)[i][j] = treeData.ZeroCoupon[1][i][j];

        results->storeGreek(domDiscFact,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("DOMESTIC_IR_DISCOUNT_FACTORS")));

        // store foreign IR spotVols
        DoubleArraySP forIRSpotVol(new DoubleArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*forIRSpotVol)[i] = treeData.SpotVol[FOREIGN][i];

        results->storeGreek(forIRSpotVol,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FOREIGN_IR_SPOT_VOLS")));

        // store domestic IR spotVols
        DoubleArraySP domIRSpotVol(new DoubleArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*domIRSpotVol)[i] = treeData.SpotVol[DOMESTIC][i];

        results->storeGreek(domIRSpotVol,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("DOMESTIC_IR_SPOT_VOLS")));

        // store FX spot vols
        DoubleArraySP fxSpotVols(new DoubleArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*fxSpotVols)[i] = treeData.SpotFxVol[i];

        results->storeGreek(fxSpotVols,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FX_SPOT_VOLS")));

        // store FX composite vols
        DoubleArraySP fxCompVols(new DoubleArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*fxCompVols)[i] = treeData.FxVol[i];

        results->storeGreek(fxCompVols,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FX_COMP_VOLS")));

        // store forward FX 
        DoubleArraySP fxFwds(new DoubleArray(treeData.NbTP+1));
        for (i = 0; i <= treeData.NbTP; i++)
            (*fxFwds)[i] = treeData.FwdFx[i];

        results->storeGreek(fxFwds,
                            Results::DEBUG_PACKET, 
                            OutputNameConstSP(new OutputName("FX_FWDS")));

        // Store local volatility function for various strikes
        CDoubleMatrixSP smileVols(new CDoubleMatrix(5, treeData.NbTP+1));
        DoubleArraySP smileDeltas(new DoubleArray(5));
        (*smileDeltas)[0] = 0.1;
        (*smileDeltas)[1] = 0.25;
        (*smileDeltas)[2] = 0.5;
        (*smileDeltas)[3] = 0.75;
        (*smileDeltas)[4] = 0.9;

        for (i = 0; i <= treeData.NbTP;i++)
        {
            double time = (*daysOffsets)[i] / 365.0;
            double vol = sqrt(time)* treeData.FxVol[i]; // use comp Vol for estiamting delta
            double fxSpotVol = treeData.SpotFxVol[i];   // spot Vol for local vol function

            for (j = 0; j < 5; j++)
            {
                double delta = (*smileDeltas)[j];
                double strikeFwdFx = deltaToStrike(delta, vol, treeData.FwdFx[i]);
                double moneyness = strikeFwdFx/treeData.FwdFx[i];
                double smile;

                if (Hyb3_Gfunc(&smile,
                            moneyness,
                            treeData.A1C[i],
                            treeData.A2C[i],
                            treeData.A3C[i]) != SUCCESS)
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


void Hyb3CBFX::getFXIndex(TreeSlice& treeSlice) {
    try {
        if (treeData.TreeType != TTYPE_FX2IR)
            throw ModelException("Hyb3 tree type must be TTYPE_FX2IR to "
                "calculate an FXIndex - tree mode is currently set to " + 
                Format::toString(treeData.TreeType));

        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        if (slice.getDim() != 3)
            slice.allocDim(3);

        // for now simply copy the FX Spot slice into the provided slice
        if (Hyb3_CopySlice(slice.getValuePtr(),
                           devData.FxSpot,
                           3,       // FX always 2D slice in this mode
                           mTpIdx,  // current step
                           &treeData) != SUCCESS) {
            throw ModelException("Hyb3_CopySlice failed", IrConverter::slog.pop());
        }
        slice.treeStep = range->treeStep;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

void Hyb3CBFX::getFXIndex(const IndexSpecFX& fxSpec, DateTime currentDate, DateTime resetDate,
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
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

        if (slice.getDim() != 3)
            slice.allocDim(3);

        // for now simply copy the FX Spot slice into the provided slice
        if (Hyb3_CopySlice(slice.getValuePtr(),
                           devData.FxSpot,
                           3,       // FX always 2D slice in this mode
                           mTpIdx,  // current step
                           &treeData) != SUCCESS) {
            throw ModelException("Hyb3_CopySlice failed", IrConverter::slog.pop());
        }

        // adjust for any date offsets
        dateAdjustFXIndex(fxSpec, currentDate, resetDate, treeSlice);
        slice.treeStep = range->treeStep;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// use deterministic FX ratio trick to adjust the stochastic FX slice
void Hyb3CBFX::dateAdjustFXIndex(const IndexSpecFX& fxSpec, DateTime currentDate, 
                                 DateTime resetDate, TreeSlice& treeSlice) {
    try {
        TreeSliceRates& slice = dynamic_cast<TreeSliceRates&>(treeSlice);

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

        if (Hyb3_MultiplyScalar(slice.getValuePtr(),
                                3,
                                ratio,
                                mTpIdx,
                                &treeData) != SUCCESS)
            throw ModelException("Hyb3_MultiplyScalar failed: " + IrConverter::slog.pop());
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}


// hyb3 eq mode always uses index [0] for domestic yield curve
int Hyb3CBFX::getCurrIdx(const string& curveName) const {
    
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
                                         "in hyb3 tree registered curve names");                
                index = i;
                match=true;
                break;  // if match found check the other currency as well
            }
        }
    }

    if (index == -1)
        throw ModelException(__FUNCTION__, "Unable to find curveName " +
                             curveName + 
                             " in list of hyb3 registered yield curve names");
    
    return index;
}


// interest rate curves are always one/first dimension
// ??? unless it's a Hyb2+1 tree
int Hyb3CBFX::getCrvDim(const string& curveName) const {
    try {
        int currIdx = getCurrIdx(curveName);
        if (currIdx == FOREIGN)
            return mktVolData[FOREIGN].NbFactor;
        else
            return mktVolData[FOREIGN].NbFactor+mktVolData[DOMESTIC].NbFactor;
    }
    catch (exception& e) {
        throw ModelException(e, __FUNCTION__);
    }
}

// forward the price back to the valueDate as the convention is to report values to
// the valueDate
double Hyb3CBFX::getPrice0(const TreeSlice& price) const {
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
// pricing on Hyb3FX is that we don't want to discount any foreign denominated
// payments in foreign but instead convert to domestic equivalent by FX at payment date
// and discount in domestic from there.  To achieve this within the engine, we
// use a 3D domestic discounting slice, and essentially discount 1.0*FX(matDate) 
// in the claim bank
void Hyb3CBFX::registerZero(DateTime useDate, DateTime matDate, 
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
            if (curveName != Hyb3CB::curveName[FOREIGN][1])
                throw ModelException("engine only supports discounting foreign currency "
                    "curve if foreign curve name " + curveName + " is the foreign discount "
                    "curve defined in the tree " + Hyb3CB::curveName[FOREIGN][1]);
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
void Hyb3CBFX::getZero(DateTime useDate, DateTime matDate, 
                       string curveName, TreeSliceSP &slice) {
    try {
        int bankSliceDim;
        map<string, NamedClaimBank>::iterator cb;
        DateTime currentDate = getCurrentDate();
        string curveNameL = curveName;

        if (!dynamic_cast<TreeSliceRates*>(slice.get())) {
            slice = RateTree::createSimpleSlice(curveName);
            slice->name = "Zero "+useDate.toString()+" to "+matDate.toString();
        }
        TreeSliceRates& sliceL = dynamic_cast<TreeSliceRates&>(*slice);

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
                if (Hyb3_CopySlice(slicePtr,
                                   devData.FxSpot,
                                   3,
                                   mTpIdx,
                                   &treeData) != SUCCESS)
                    throw ModelException("Hyb3_CopySlice failed: " + IrConverter::slog.pop());
            }
            else {  // domestic
                int domSliceDim = mktVolData[FOREIGN].NbFactor + mktVolData[DOMESTIC].NbFactor;
                sliceL.allocDim(domSliceDim); 
                slicePtr = sliceL.getValuePtr();
                if (Hyb3_SetSlice(slicePtr,
                                  domSliceDim,
                                  1.0,  // set slice value as 1.0
                                  mTpIdx,
                                  &treeData) != SUCCESS)
                    throw ModelException("Hyb3_SetSlice failed: " + IrConverter::slog.pop());
            }
        }
        else {
            double* zeroBkPt = NULL;

            // pointer to internal zero slice maturity - zeroBkPt is read only and should not
            // be freed here
            zeroBkPt = Hyb3_ZbkReadZero((CLAIM_BANK*)zeroBank,
                                        matDate.toIrDate(),
                                        FALSE,  // ??? no interp, do we want to relax this
                                        currentDate.toIrDate(),
                                        mTpIdx,
                                        bankSliceDim,
                                        &treeData);

            if (zeroBkPt == NULL)
                throw ModelException("Hyb3_ZbkReadZero failed: " + IrConverter::slog.pop());

            sliceL.allocDim(bankSliceDim);
            slicePtr = sliceL.getValuePtr();
            // deep copy results of ReadZero function call to output slice
            if (Hyb3_CopySlice(slicePtr,
                                zeroBkPt,
                                bankSliceDim,
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

string Hyb3CBFX::foreignZeroName(string curveName) {
    string name = curveName + "-FOREIGN_AS_DOMESTIC_ZERO";
    return name;
}

DRLIB_END_NAMESPACE
