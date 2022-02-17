#include "edginc/config.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/FD1DSolver.hpp"
#include "edginc/FD1DLN.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/VolSurface.hpp"
#include "edginc/UtilFuncs.hpp"

#include "edginc/AssetUtil.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/ProtEquity.hpp"


DRLIB_BEGIN_NAMESPACE

////////////////////////// FD1DLN /////////////////////////
void FD1DLN::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FD1DLN, clazz);
    SUPERCLASS(FD1D);
    EMPTY_SHELL_METHOD(defaultFD1DLN);
}

//----------------------------------------------------------------
// helpers
CClassConstSP const FD1DLN::TYPE = CClass::registerClassLoadMethod(
    "FD1DLN", typeid(FD1DLN), load);

//----------------------------------------------------------------

FD1DLN::FD1DLN( const CClassConstSP & type ) :
    FD1D( type ),
    volType(VolSurface::TYPE->getName())
{
    changeOfVar = LOG_X;  //1: S; 2: log(S); 5: log(fwd)
}

//----------------------------------------------------------------

FD1DLN::~FD1DLN(){}

//----------------------------------------------------------------

void FD1DLN::validatePop2Object()
{
    static const string method = "FD1DLN::validatePop2Object";
    try
    {
        FD1D::validatePop2Object();

        if (isFwdInduction)
            throw ModelException(method, "Fwd Induction hasn't implemented yet! ");
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

MarketDataFetcherSP FD1DLN::createMDF() const {
    return MarketDataFetcherSP(new MDFAssetVol(volType));
}

IModel::WantsRiskMapping FD1DLN::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** retrieving market data from MDF */    
void FD1DLN::retrieveFactor(){
    static const string method = "FD1DLN::retrieveFactor";
    try{
        FD1D::retrieveFactor();

        if ((StruckEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ||
            (ProtEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ){
                throw ModelException("FD1DLN::retrieveFactor", "Protected and struck haven't been implemented yet!");
        }

        volRequests.resize(factors.size());
        drifts.resize(factors.size()); 
        vols.resize(factors.size()); 
        
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

void FD1DLN::generateDriftAndVols(const DateTimeArray& futurePathDates)
{
    static const string routine("FD1DLN::generateDriftAndVols");
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
                double drift = stepForward[iStep+1]/
                    stepForward[iStep];

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

CVolProcessedBSConstSP FD1DLN::getProcessedVol()
{
    return volLN1; //to review, for now, just use the vol in mAsset[0], it's to getTimeMetric
}

//----------------------------------------------------------------

void FD1DLN::initModel()
{
    //call the base;
    FD1D::initModel();
}

void FD1DLN::finaliseModel(CControl*    control)
{
    //timeLine should be already be set up
    int iAsset;
    for (iAsset=0; iAsset < factors.size(); iAsset++) {
        dynamic_cast<const CAsset *>(
            factors[iAsset].get())->fwdValue(timeLine->StepDates, stepForward);
        drifts[iAsset]   = DoubleArray(timeLine->NumOfStep+1);
        vols[iAsset]     = DoubleArray(timeLine->NumOfStep+1);
    }

    // precompute drift and vols could possibly have a class
    // "process" which did this - for now leave as local functions
    generateDriftAndVols(timeLine->StepDates);
       
    // call parent last
    FD1D::finaliseModel( control );
}

//----------------------------------------------------------------

/** calculate the stuff which specifies the grid */
/**calculate FD grid and its underlying*/
/**alpha : truncation1D*/
void FD1DLN::setFdBounds(
    double alpha1,
    double & outLowB1, double & outUpB1 ) const
{
    //get the vol from the underlying to set the FDBound
    double volForBos = vols[0][timeLine->NumOfStep-1];  //to check if it's NumOfStep or NumOfStep-1
          
    double longest_T = timeLine->TradeTime[timeLine->NumOfStep];
    double temp1 = alpha1* volForBos * sqrt(longest_T);
    double undSpot1 = underlying->fwdValue(timeLine->StepDates[0]);

    switch (changeOfVar) {
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
int FD1DLN::GetStepVol(int step, vector<double>& vol, const double* s, int start, int end) const
{
    vol.resize(end-start+1, vols[0][step]);  //vols[iAsset][step]

    return (end-start+1);
}

//----------------------------------------------------------------
//  setup the coefficients of the PDE needed to be solved
//  Assume the most general PDE is 
//  U_t + a*U_x + c*U_xx + f*U + g = 0;
//  this ft is called at each time step, so they're time depd.
//----------------------------------------------------------------

void FD1DLN::pdeCoeff(
    int step,
    TreeSliceEQ & a,
    TreeSliceEQ & c,
    TreeSliceEQ & f,
    TreeSliceEQ & g ) const
{
    int ix;

    if (step == timeLine->NumOfStep){
        throw ModelException("FD1DLN::pdeCoeff", 
            "pdeCoeff shouldn't be called at last step!.");
    }
            
    double mu1 = drifts[0][step];
    double dt_trading =timeLine->TradeYrFrac[step+1];;
    double sig1 = vols[0][step];   //mAsset[0]
    double sig12 = 0.5 * sig1 * sig1 * dt_trading;
    double df = discYC->pv(timeLine->StepDates[step], timeLine->StepDates[step+1]);
    double r = -log(df);

    int bot1, top1;
    xValue->getCalcRange( bot1, top1 );
    double * s = xValue;

    //fill all coeff of PDE
    for (ix = bot1; ix <= top1; ix++){
            c[ix] = sig12 ;    //c*dt, coeff of U_xx
            f[ix] = -r;    //f*dt, coeff of U
            g[ix]= 0.0;  //g*dt, 

            switch (changeOfVar) {
                case X:{//S
                    a[ix] = mu1 * s[ix];  //a *dt, coeff of U_x
                    c[ix] *= s[ix] * s[ix];   //c*dt, coeff of U_xx
                }        
                break;
                case LOG_X:{//log(S)
                    a[ix] = mu1 - sig12 ;  //a *dt, coeff of U_x
                }
                break;
                case LOG_FWDX:{//log(fwd)
                    a[ix] = - sig12 ;  //a *dt, coeff of U_x
                }
                break;
                default:{
                    throw ModelException("FD1DLN::pdeCoeff", "This change of variable is not available");
            }
        }
    }
}

//----------------------------------------------------------------

DRLIB_END_NAMESPACE
