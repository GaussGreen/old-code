//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DLN.hpp
//
//   Description : two factor finite difference engine for two LN
//                 asset factors   dX_i = drift_i*dt + sigma_i * dW_i with i=1,2
//
//   Author      : Xiaolan Zhang
//
//   Date        : May 27, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FD2DLN.hpp"
#include "edginc/LatticeProdEDR.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/IInstrumentCollection.hpp"
#include "edginc/AssetUtil.hpp"
#include "edginc/StruckEquity.hpp"
#include "edginc/ProtEquity.hpp"
#include "edginc/VolSurface.hpp"

DRLIB_BEGIN_NAMESPACE

//const double DEFAULT_LBOUND2 = 0.000025; 

//const double DEFAULT_P_UND2 = 0.01;


////////////////////////// FD2DLN /////////////////////////
void FD2DLN::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FD2DLN, clazz);
    SUPERCLASS(FD2D);
    EMPTY_SHELL_METHOD(defaultFD2DLN);
}

//----------------------------------------------------------------
// helpers
CClassConstSP const FD2DLN::TYPE = CClass::registerClassLoadMethod(
    "FD2DLN", typeid(FD2DLN), load);

//----------------------------------------------------------------

FD2DLN::FD2DLN() : FD2D(TYPE), 
                       volType(VolSurface::TYPE->getName()){

    whichChgVarDim1 = LOG_X;  //1: S; 2: log(S); 5: log(fwd)
    whichChgVarDim2 = whichChgVarDim1;  //1: S; 2: log(S); 5: log(fwd)
}

//----------------------------------------------------------------

FD2DLN::~FD2DLN(){}

//----------------------------------------------------------------

void FD2DLN::validatePop2Object()
{
    static const string method = "FD2DLN::validatePop2Object";
    try {
        FD2D::validatePop2Object();

        if (whichChgVarDim2 != whichChgVarDim1){
            throw ModelException(method, 
                                 "FD2DLN, only support whichChgVarDim1 = whichChgVarDim2");
        }
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


MarketDataFetcherSP FD2DLN::createMDF() const {
    return MarketDataFetcherSP(new MDFAssetVol(volType));
}

IModel::WantsRiskMapping FD2DLN::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

/** Invoked after instrument has got its market data. */
void FD2DLN::getMarket(const MarketData* market, IInstrumentCollectionSP instruments)
{
    FDModel::getMarket(market, instruments);

    if( factors.size() != 2 )
        throw ModelException( "FD2DLN::getMarket", "only 2 factors suppported" );
}

/** retrieving market data from MDF */    
void FD2DLN::retrieveFactor(){
    static const string method = "FD2DLN::retrieveFactor";
    try{
        FD2D::retrieveFactor();

        underlying = CAssetConstSP::dynamicCast( factors[0] );
        underlying2 = CAssetConstSP::dynamicCast( factors[1] );
        if( ! underlying || ! underlying2 )
            throw ModelException( "FD2DLN::retrieveFactor", "only asset underliers suppported" );

        if ((StruckEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ||
            (ProtEquity::TYPE->isInstance(underlying) && !(AssetUtil::isBasket(underlying))) ){
                throw ModelException("FD2DLN::retrieveFactor", "Protected and struck haven't been implemented yet!");
        }

        if ((StruckEquity::TYPE->isInstance(underlying2) && !(AssetUtil::isBasket(underlying2))) ||
            (ProtEquity::TYPE->isInstance(underlying2) && !(AssetUtil::isBasket(underlying2))) ){
                throw ModelException("FD2DLN::retrieveFactor", "Protected and struck haven't been implemented yet!");
        }

        volRequests.resize(factors.size());
        drifts.resize(factors.size()); 
        vols.resize(factors.size()); 
        //volLN.resize(factors.size());
        //maxDrifts(prod->getNumAssets(), 1.0),
        //fwds(prod->getfactors.size()());
        stepForward.resize(factors.size());
        
        //init vol
        const IFDProductLN* prodLN = 
            dynamic_cast<const IFDProductLN*>(prod.get());
        if (!prodLN){
            throw ModelException(method, "Product does not implement "
                                 "IFDProductLN interface");
        }

        //get correlation
        rho = (*correls)[0][1];
        
        int iAsset;
        for (iAsset=0; iAsset < factors.size(); iAsset++) {
            // for each asset have an array of vol requests
            volRequests[iAsset] = 
                prodLN->getVolInterp(iAsset);
        }
        // interpolate the vol using our LN request
        //should retrieve vol from Underlier, not from mAsset since Underlier maybe changed!
        volLN1 = CVolProcessedBSSP(underlying->getProcessedVol(volRequests[0].get()));
        volLN2 = CVolProcessedBSSP(underlying2->getProcessedVol(volRequests[1].get()));

        // get time metric, use 1st asset for now
        timeMetric = volLN1->GetTimeMetric();
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

void FD2DLN::generateDriftAndVols(const DateTimeArray& futurePathDates) {
    static const string routine("FD2DLN::"
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
                //drift *= exp(-0.5 * fwdVar[iStep]);
                //drifts[iAsset][iStep] = drift;

                //drifts = mu in SDE
                drifts[iAsset][iStep] = log(drift);

                totalVar += fwdVar[iStep];
                //theDrift *= Maths::max(1.0, drift);

                double vol = sqrt(fwdVar[iStep] );
                vol = vol / sqrt(timeLine->TradeTime[iStep+1]-timeLine->TradeTime[iStep]);
                vols[iAsset][iStep] = vol;                
            }
            // here we're simulating what happens in generatePath
            // want theDrift to represent [reasonable] worst
            // case drift. A value of 2 for MAX_DRIFT_RAND_NUMBER
            // means that we catch >96% of paths
//            theDrift *= exp(sqrt(totalVar)*MAX_DRIFT_RAND_NUMBER);
//            if (theDrift > maxDrifts[iAsset]){
//                maxDrifts[iAsset] = theDrift;
//            }
        }
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}


CVolProcessedBSConstSP FD2DLN::getProcessedVol()
{
    return volLN1; //to review, for now, just use the vol in mAsset[0], it's to getTimeMetric
}

//----------------------------------------------------------------

/** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
double FD2DLN::getPrice0( const TreeSlice & price ) const
{
    const string method = "FD2DLN::getPrice0";

    // to do : interpolate for centre point, otherwise same grid tweak won't work !
    // to do : integrate the price if isRandomInitialVol == true

    try {

        int  botS1, topS1, botS2, topS2;
        const TreeSlice & s1 = payoffIndex->getValue( 0 );
        s1.getCalcRange( botS1, topS1 );
        const TreeSlice & s2 = payoffIndex2->getValue( 0 );
        s2.getCalcRange( botS2, topS2 );

        double undSpot1 = underlying->getSpot();
        double undSpot2 = underlying2->getSpot();

        return interpF2(undSpot1, 
                        undSpot2, 
                        s1.getValues(),
                        s2.getValues(),
                        price.getValues(),
                        topS1 - botS1 + 1,
                        topS2 - botS2 + 1);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

void FD2DLN::initModel()
{
    //call the base;
    FD2D::initModel();
}

void FD2DLN::finaliseModel(CControl*    control){

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
       
    FD2D::finaliseModel(control);    
}

//----------------------------------------------------------------
/** update payoffIndex          */

void FD2DLN::payoffIndexUpdate (int& step, FDProduct::UpdateType type){
    payoffIndex->update(step, type);
    payoffIndex2->update(step, type);
}

//----------------------------------------------------------------

/** calculate the stuff which specifies the grid */
/**calculate FD grid and its underlying*/
/**alpha : truncation1D*/
void FD2DLN::setFdBounds(double& volForBound1, double alpha1, double& outLowB1, double& outUpB1,
                         double& volForBound2, double alpha2, double& outLowB2, double& outUpB2){

    //get two vol from 2 underlying to set the FDBound
    volForBound1 = vols[0][timeLine->NumOfStep-1];  //to check if it's NumOfStep or NumOfStep-1
    volForBound2 = vols[1][timeLine->NumOfStep-1];  
          
    double longest_T = timeLine->TradeTime[timeLine->NumOfStep];

    double temp1 = alpha1* volForBound1 * sqrt(longest_T);
    double temp2 = alpha2* volForBound2 * sqrt(longest_T);

    double undSpot1 = underlying->getSpot();
    double undSpot2 = underlying2->getSpot();


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

    switch (whichChgVarDim2) {
        case X:{
            //S
            outLowB2 = undSpot2 * exp(-temp2);
            outUpB2 = undSpot2* exp(temp2);
        }        
        break;
        case LOG_X: {
            //log(S)
            outLowB2 = log(undSpot2) - temp2;
            outUpB2 = log(undSpot2) + temp2;
        }
        break;
        case LOG_FWDX: {
            //log(fwd)
            outLowB2 = log(undSpot2) - temp2;
            outUpB2 = log(undSpot2) + temp2;
        }
        break;
    }

}

//----------------------------------------------------------------
//  setup the coefficients of the PIDE needed to be solved
//  Assume the most general PDE is 
//  U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = 0;
//  this ft is called at each time step, so they're time depd.
//----------------------------------------------------------------

void FD2DLN::pdeCoeff(int step, double*** coeff, 
                        int bot1, int top1, int bot2, int top2 ){

    int ix;
    int iy;

    if (step == timeLine->NumOfStep){
        throw ModelException("FD2DLN::pdeCoeff", 
            "pdeCoeff shouldn't be called at last step!.");
    }
            
    double mu1 = drifts[0][step];
    double mu2 = drifts[1][step];

    double dt_trading =timeLine->TradeYrFrac[step+1];;

    double sig1 = vols[0][step];   //mAsset[0]
    double sig2 = vols[1][step];  //mAsset[1]
    double sig12 = sig1 * sig1 * dt_trading;
    double sig22 = sig2 * sig2* dt_trading;

    double df = discYC->pv(timeLine->StepDates[step], timeLine->StepDates[step+1]);
    double r = -log(df);

    // fpr accessing und level, need review !!!
    payoffIndex->update(step, FDProduct::BWD);
    double* und1 = payoffIndex->getValue( step ).getValues();

    payoffIndex2->update(step, FDProduct::BWD);
    double* und2 = payoffIndex2->getValue( step ).getValues();
    
    //fill all coeff of PDE
    for (ix = bot1; ix <= top1; ix++){
        for (iy = bot2; iy <= top2; iy++){

            coeff[2][ix][iy] = 0.5 * sig12 ;  //c
            coeff[3][ix][iy] = 0.5 * sig22 ;  //d
            coeff[4][ix][iy] = rho * sig1 * sig2 * dt_trading;  //e
            coeff[5][ix][iy] = -r;  //f

            switch (whichChgVarDim1) {
                case X:{//S
                    coeff[0][ix][iy] = mu1 * und1[ix];  //a
                    coeff[1][ix][iy] = mu2 * und2[ix];  //b
                    coeff[2][ix][iy] *= und1[ix] * und1[ix];  //c
                    coeff[3][ix][iy] *= und2[ix] * und2[ix];  //d
                    coeff[4][ix][iy] *= und1[ix] * und2[ix];  //e
                }        
                break;
                case LOG_X:{//log(S)
                    coeff[0][ix][iy] = mu1 - 0.5 * sig12 ;  //a
                    coeff[1][ix][iy] = mu2 - 0.5 * sig22 ;  //b
                }
                break;
                case LOG_FWDX:{//log(fwd)
                    coeff[0][ix][iy] = - 0.5 * sig12 ;  //a
                    coeff[1][ix][iy] = - 0.5 * sig22 ;  //b
                }
                break;
                default:{
                    throw ModelException("FD2DLN::pdeCoeff", "This change of variable is not available");
                }
            }
        }
    }
}

//----------------------------------------------------------------

FDProductSP FD2DLN::makeProduct(const IProdCreatorSP & creator)
{
    const IndexSpecEQ * spec = dynamic_cast< const IndexSpecEQ * >( creator.get() );
    if( spec )
    {
        if( spec->getFactor()->getName() == factors[0]->getName() )
        {
            return payoffIndex = FDModel::makeProduct( IProdCreatorSP( new Spot( spec->getFactor()->getName() ) ) );
        }
        else if( spec->getFactor()->getName() == factors[1]->getName() )
        {
            return payoffIndex2 = FDModel::makeProduct( IProdCreatorSP( new Spot2( spec->getFactor()->getName() ) ) );
        }
        else
        {
            throw ModelException( "FD2DLN::makeProduct",
                "IndexSpec source is '" + spec->getFactor()->getName() + "', "
                "while '" + factors[0]->getName() + "' "
                "or '" + factors[1]->getName() + "' are expected" );
        }
    }

    return FDModel::makeProduct( creator );
}

//----------------------------------------------------------------

bool FD2DLNLoad(){
    return (FD2DLN::TYPE && true);
}

DRLIB_END_NAMESPACE
