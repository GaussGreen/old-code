//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD1DRetLN.hpp
//
//   Description : one factor finite difference engine for 1 LN
//                 asset factors   dX_i = drift_i*dt + sigma_i * dW_i with i=1,
//
//   Author      : Xiaolan Zhang
//
//   Date        : Feb 14, 2006
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/FD1DRetSolver.hpp"
#include "edginc/FD1DRetLN.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/UtilFuncs.hpp"

#include "edginc/AssetUtil.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/ProtEquity.hpp"


DRLIB_BEGIN_NAMESPACE

////////////////////////// FD1DRetLN /////////////////////////
void FD1DRetLN::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FD1DRetLN, clazz);
    SUPERCLASS(FD1DRet);
    EMPTY_SHELL_METHOD(defaultFD1DRetLN);
}

//----------------------------------------------------------------
// helpers
CClassConstSP const FD1DRetLN::TYPE = CClass::registerClassLoadMethod(
    "FD1DRetLN", typeid(FD1DRetLN), load);

//----------------------------------------------------------------

FD1DRetLN::FD1DRetLN() : FD1DRet(TYPE), 
                       volType(VolSurface::TYPE->getName()){

    whichChgVarDim1 = LOG_X;  //1: S; 2: log(S); 5: log(fwd)
}

//----------------------------------------------------------------

FD1DRetLN::~FD1DRetLN(){}

//----------------------------------------------------------------

void FD1DRetLN::validatePop2Object()
{
    static const string method = "FD1DRetLN::validatePop2Object";
    try {
        FD1DRet::validatePop2Object();

        if (isFwdInduction)
        {
            throw ModelException(method, "Fwd Induction hasn't implemented yet! ");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

MarketDataFetcherSP FD1DRetLN::createMDF() const {
    return MarketDataFetcherSP(new MDFAssetVol(volType));
}

IModel::WantsRiskMapping FD1DRetLN::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** retrieving market data from MDF */    
void FD1DRetLN::retrieveFactor(){
    static const string method = "FD1DRetLN::retrieveFactor";
    try{
        FD1DRet::retrieveFactor();

        //temporary
        underlying = CAssetConstSP::dynamicCast( factors[0] );
        if ((StruckEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ||
            (ProtEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ){
                throw ModelException("FD1DRetLN::retrieveFactor", "Protected and struck haven't been implemented yet!");
        }

        CAssetConstSP plainAsset = underlying;

        volRequests.resize(factors.size());
        drifts.resize(factors.size()); 
        vols.resize(factors.size()); 
        stepForward.resize(factors.size());
        
        //init vol
        const IFDProductLN* prodLN = 
            dynamic_cast<const IFDProductLN*>(prod.get());
        if (!prodLN){
            throw ModelException(method, "Product does not implement "
                                 "IFDProductLN interface");
        }

        int iAsset;
        for (iAsset=0; iAsset < factors.size(); iAsset++) {
            // for each asset have an array of vol requests
            volRequests[iAsset] = 
                prodLN->getVolInterp(iAsset);
        }
        // interpolate the vol using our LN request
        //should retrieve vol from Underlier, not from mAsset since Underlier maybe changed!
        volLN1 = CVolProcessedBSSP(underlying->getProcessedVol(volRequests[0].get()));

        // get time metric, use 1st asset for now
        timeMetric = volLN1->GetTimeMetric();
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

void FD1DRetLN::generateDriftAndVols(const DateTimeArray& futurePathDates) {
    static const string routine("FD1DRetLN::"
                                "generateDriftAndVols");
    try{
        CDoubleArray fwdVar(timeLine->NumOfStep+1); // should be -1?
        
        for(int iAsset=0;iAsset<factors.size();iAsset++) {
            // interpolate the vol using our LN request
            // XXX factor vs asset
            CVolProcessedSP vol(dynamic_cast<const CAsset *>(
                factors[iAsset].get())->getProcessedVol(volRequests[iAsset].get())); 
            CVolProcessedBSSP volBS(CVolProcessedBSSP::dynamicCast(vol));
            volBS->CalcVar(futurePathDates, volBS->forward, fwdVar);

            // drifts[iStep] is drift FROM iStep point TO iStep+1
            // vols[iStep] is vol FROM iStep point TO iStep+1  (to double check)
            double totalVar = 0.0; // used for max drift calc
            for (int iStep = 0; iStep < timeLine->NumOfStep; iStep ++)
            {
                //drift is the mu, and PDECoef will do the corresponding chg of variable
                double drift = stepForward[iAsset][iStep+1]/
                    stepForward[iAsset][iStep];

                //drifts = mu in SDE
                drifts[iAsset][iStep] = log(drift);

                totalVar += fwdVar[iStep];
                //theDrift *= Maths::max(1.0, drift);

                double vol = sqrt(fwdVar[iStep] );
                vol = vol / sqrt(timeLine->TradeTime[iStep+1]-timeLine->TradeTime[iStep]);
                vols[iAsset][iStep] = vol;                
            }
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}

CVolProcessedBSConstSP FD1DRetLN::getProcessedVol()
{
    return volLN1; //to review, for now, just use the vol in mAsset[0], it's to getTimeMetric
}

//----------------------------------------------------------------

void FD1DRetLN::initModel()
{
    //call the base;
    FD1DRet::initModel();
}

void FD1DRetLN::finaliseModel(CControl*    control){

    //timeLine should be already be set up
    int iAsset;
    for (iAsset=0; iAsset < factors.size(); iAsset++) {
        stepForward[iAsset] = DoubleArray(timeLine->NumOfStep+1, 0.0);
        dynamic_cast<const CAsset *>(
            factors[iAsset].get())->fwdValue(timeLine->StepDates, stepForward[iAsset]);
        drifts[iAsset]   = DoubleArray(timeLine->NumOfStep+1);
        vols[iAsset]     = DoubleArray(timeLine->NumOfStep+1);
    }

    // precompute drift and vols could possibly have a class
    // "process" which did this - for now leave as local functions
    generateDriftAndVols(timeLine->StepDates);
       
    FD1DRet::finaliseModel(control);    
}

//----------------------------------------------------------------
/** update payoffIndex          */

void FD1DRetLN::payoffIndexUpdate (int& step, FDProduct::UpdateType type){
    payoffIndex->update(step, type);
}

//----------------------------------------------------------------

/** calculate the stuff which specifies the grid */
/**calculate FD grid and its underlying*/
/**alpha : truncation1D*/
void FD1DRetLN::setFdBounds(double& volForBound1, double alpha1, double& outLowB1, double& outUpB1){

    //get the vol from the underlying to set the FDBound
    volForBound1 = vols[0][timeLine->NumOfStep-1];  //to check if it's NumOfStep or NumOfStep-1
          
    double longest_T = timeLine->TradeTime[timeLine->NumOfStep];
    double temp1 = alpha1* volForBound1 * sqrt(longest_T);
    double undSpot1 = underlying->fwdValue(timeLine->StepDates[0]);

    switch (whichChgVarDim1) {
        case X:{
            //S
            outLowB1 = undSpot1 * exp(-temp1);
            outUpB1 = undSpot1* exp(temp1);
        }        
        break;
        case LOG_X: {
            //log(S)
            outLowB1 = log(undSpot1) - temp1;
            outUpB1 = log(undSpot1) + temp1;
        }
        break;
        case LOG_FWDX: {
            //log(fwd)
            outLowB1 = log(undSpot1) - temp1;
            outUpB1 = log(undSpot1) + temp1;
        }
        break;
    }
}

    /**----------------------------------------------------------------
        calculate vol for current step for a set of spot levels
        returns number of vol calculated - one (flat for all node) or num 
    ----------------------------------------------------------------*/
//to review and test, for dblbarrier
int FD1DRetLN::GetStepVol(int step, vector<double>& vol, const double* s, int start, int end){

    vol.resize(end-start+1, vols[0][step]);  //vols[iAsset][step]

    return (end-start+1);
}

//----------------------------------------------------------------
//  setup the coefficients of the PDE needed to be solved
//  Assume the most general PDE is 
//  U_t + a*U_x + c*U_xx + f*U + g = 0;
//  this ft is called at each time step, so they're time depd.
//----------------------------------------------------------------

void FD1DRetLN::pdeCoeff(int step, double** coeff, 
                        int bot1, int top1){

    int ix;

    if (step == timeLine->NumOfStep){
        throw ModelException("FD1DRetLN::pdeCoeff", 
            "pdeCoeff shouldn't be called at last step!.");
    }
            
    double mu1 = drifts[0][step];
    double dt_trading =timeLine->TradeYrFrac[step+1];;
    double sig1 = vols[0][step];   //mAsset[0]
    double sig12 = 0.5 * sig1 * sig1 * dt_trading;
    double df = discYC->pv(timeLine->StepDates[step], timeLine->StepDates[step+1]);
    double r = -log(df);

    // fpr accessing und level, need review !!!
    payoffIndex->update(step, FDProduct::BWD);
    double* und1 = payoffIndex->getValue( step ).getValues();

    //fill all coeff of PDE
    for (ix = bot1; ix <= top1; ix++){
            coeff[1][ix] = sig12 ;    //c*dt, coeff of U_xx
            coeff[2][ix] = -r;    //f*dt, coeff of U
            coeff[3][ix]= 0.0;  //g*dt, 

            switch (whichChgVarDim1) {
                case X:{//S
                    coeff[0][ix] = mu1 * und1[ix];  //a *dt, coeff of U_x
                    coeff[1][ix] *= und1[ix] * und1[ix];   //c*dt, coeff of U_xx
                }        
                break;
                case LOG_X:{//log(S)
                    coeff[0][ix] = mu1 - sig12 ;  //a *dt, coeff of U_x
                }
                break;
                case LOG_FWDX:{//log(fwd)
                    coeff[0][ix] = - sig12 ;  //a *dt, coeff of U_x
                }
                break;
                default:{
                    throw ModelException("FD1DRetLN::pdeCoeff", "This change of variable is not available");
            }
        }
    }
}

//----------------------------------------------------------------

FDProductSP FD1DRetLN::makeProduct(const IProdCreatorSP & creator)
{
    const IndexSpecEQ * spec = dynamic_cast< const IndexSpecEQ * >( creator.get() );
    if( spec )
    {
        if( spec->getFactor()->getName() == factors[0]->getName() )
        {
            return payoffIndex = FDModel::makeProduct( IProdCreatorSP( new Spot( spec->getFactor()->getName() ) ) );
        }
        else
        {
            throw ModelException( "FD1DRetLN::makeProduct",
                "IndexSpec source is '" + spec->getFactor()->getName() + "', "
                "while '" + factors[0]->getName() + "' "
                "' is expected" );
          }
    }

    return FDModel::makeProduct( creator );
}

//----------------------------------------------------------------

bool FD1DRetLNLoad(){
    return (FD1DRetLN::TYPE && true);
}

DRLIB_END_NAMESPACE
