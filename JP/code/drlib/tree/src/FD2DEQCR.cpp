//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FD2DEQCR.hpp
//
//   Description : two factor finite difference engine for credit
//                 asset factor dX = drift*dt + sigma*dW
//                 credit factor dp = k(p_0 - p)dt + vol*sqrt(p)*dZ
//                 Z and W have a correlation rho.
//
//   Author      : Bruno Melka
//
//   Date        : 10 Jun 2005
//
//----------------------------------------------------------------------------

//#include <math.h>
#include "edginc/config.hpp"
#include "edginc/FD2DEQCR.hpp"
#include "edginc/mathlib.hpp"
#include "edginc/EquityBase.hpp"
#include "edginc/ATMVolRequest.hpp"
#include "edginc/MRSpotVolRequest.hpp"
#include "edginc/MRSpotVolProcessed.hpp"
#include "edginc/MDFAssetVol.hpp"
#include "edginc/VolSVCJ.hpp"

DRLIB_BEGIN_NAMESPACE

const double DEFAULT_LBOUND2 = 0.000025; 

const double DEFAULT_P_UND2 = 0.01;


////////////////////////// FD2DEQCR /////////////////////////
void FD2DEQCR::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(FD2DEQCR, clazz);
    SUPERCLASS(FD2D);
    EMPTY_SHELL_METHOD(defaultFD2DEQCR);
}

//----------------------------------------------------------------
// helpers
CClassConstSP const FD2DEQCR::TYPE = CClass::registerClassLoadMethod(
    "FD2DEQCR", typeid(FD2DEQCR), load);

//----------------------------------------------------------------

FD2DEQCR::FD2DEQCR() : FD2D(TYPE){

    isVariableGrid = false;  

    whichChgVarDim1 = LOG_X;  
    whichChgVarDim2 = X;

}

//----------------------------------------------------------------

FD2DEQCR::~FD2DEQCR(){}

//----------------------------------------------------------------


/** Create a MarketDataFetcher which will be used for retrieving market data etc */
MarketDataFetcherSP FD2DEQCR::createMDF() const {
    return MarketDataFetcherSP(new MDFAssetVol(VolSVCJ::TYPE->getName(),
                                               VolSurface::TYPE->getName()));
}

//----------------------------------------------------------------

IModel::WantsRiskMapping FD2DEQCR::wantsRiskMapping() const {
    return riskMappingIrrelevant;
}

//----------------------------------------------------------------
/** retrieving market data from MDF */
void FD2DEQCR::retrieveFactor(){
    static const string method = "FD2DEQCR::retrieveFactor";
    try{
        FD2D::retrieveFactor();

        underlying = CAssetConstSP::dynamicCast( factors[0] );
        if( ! underlying )
            throw ModelException( "FD2DEQCR::retrieveFactor", "only asset underlier suppported" );

        //credit spread

        // copy the vol that's in the multiasset
        //vol = VolSVCJSP::dynamicCast(VolRequestRaw::copyVolBase(*mAsset, 0));
        // get time metric
        //timeMetric = vol->GetTimeMetric();
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

MarketObjectSP FD2DEQCR::GetMarket(const MarketData*    market,
                                    const string&        name,
                                    const CClassConstSP& type) const {
    
    static const string method("FD2DEQCR::GetMarket");
    try {

        if (CAsset::TYPE->isAssignableFrom(type)){
            MarketObjectSP asset = market->GetData(name, type);
            // need to remove const but just go round for now
            const_cast<FD2DEQCR*>(this) ->underlying = CAssetConstSP::dynamicCast(asset);
            return asset;
        }
        else if (YieldCurve::TYPE->isAssignableFrom(type)){
            MarketObjectSP curve = market->GetData(name, type);
            // need to remove const but just go round for now
            const_cast<FD2DEQCR*>(this) ->discYC = YieldCurveConstSP::dynamicCast(curve);
            return curve;
        }
        else if (ICDSParSpreads::TYPE->isAssignableFrom(type)){
            MarketObjectSP spreads = market->GetData(name, type);
            // need to remove const but just go round for now
            const_cast<FD2DEQCR*>(this) ->credit = ICDSParSpreadsConstSP::dynamicCast(spreads);
            return spreads;
        }
        else {
            // then just invoke parent with overridden type
            return market->GetData(name, type);
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

/** interpolate/integrate for prices at t=0 with asset=asset(t=0) etc. */
double FD2DEQCR::getPrice0( const TreeSlice & price ) const
{
    const string method = "FD2DEQCR::getPrice0";
    try {
        const TreeSlice & s = payoffIndex->getValue( 0 );
        int botS, topS;
        s.getCalcRange(botS, topS);

        vector<double> c(dim2);
        computeUndLevel2(0, dim2-1, 0, &*c.begin());

        double undSpot1 = underlying->getSpot();
        double undSpot2 = init_CIR;

        return interpF2(undSpot1,
                        undSpot2,
                        s.getValues(),
                        &*c.begin(),
                        price.getValues(),
                        topS - botS + 1,
                        dim2);

    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

//----------------------------------------------------------------

/** get vol etc */
void FD2DEQCR::initModel(){

    // retrieve vol factor
    const EquityBase* eq = dynamic_cast<const EquityBase*>(underlying.get());
    if (eq !=0){
        ATMVolRequestSP volRequest(new ATMVolRequest());
        volLN = CVolProcessedBSSP::dynamicCast(CVolProcessedSP(eq->getProcessedVol(volRequest.get())));
        // get time metric
        timeMetric = volLN->GetTimeMetric();
    }
    else {
        throw ModelException("FD2DEQCR::initModel", "only equity base is currently supported!");
    }
    // call base method to set up time line
    FD2D::initModel();
}

//----------------------------------------------------------------

void FD2DEQCR::finaliseModel(CControl*    control){

    if (whichChgVarDim2 != X) {
        throw ModelException( "FD2DEQCR::finaliseModel", 
                             "whichChgVarDim2 should be X (working variable is X)!");
    }

    // compute vols
    vols.resize(timeLine->StepDates.size());
    volLN->CalcVol(timeLine->StepDates, CVolProcessedBS::forward, vols);

    MRSpotVolRequest betaRequest;
    CVolProcessed* volcds = credit->getProcessedVol(&betaRequest);
    smartPtr<MRSpotVolProcessed> volCdsData(&dynamic_cast<MRSpotVolProcessed&>(*volcds));
    DoubleArray volTemp(timeLine->StepDates.size());
    volCdsData->spotVol(timeLine->StepDates[0], timeLine->StepDates, volTemp);

    vol_CIR = volTemp[0];
    mr_CIR = volCdsData->meanReversion();

    CashFlowArraySP defRatesDates;
    defRatesDates = credit->defaultRates()->getCleanSpreadCurve();

    init_CIR = 365.0 * log(credit->defaultRates()->calcDefaultPV(timeLine->StepDates[0], timeLine->StepDates[1]))
                / (timeLine->StepDates[0].getDate() - timeLine->StepDates[1].getDate());

    createTermStructureCIR(defRatesDates);

    FD2D::finaliseModel(control);

    // compute fwds
    //underlying->fwdValue(timeLine->StepDates, stepForward1);

    int iAsset = 0;
    stepForward[iAsset] = DoubleArray(timeLine->NumOfStep+1, 0.0);
    dynamic_cast<const CAsset *>(
        factors[iAsset].get())->fwdValue(timeLine->StepDates, stepForward[iAsset]);
}

//----------------------------------------------------------------

void FD2DEQCR::createTermStructureCIR(CashFlowArraySP defRatesDates) {

    lr_CIR.resize(timeLine->StepDates.size());
    DoubleArray lr_CIR_on_curve(defRatesDates->size());
    double gamma = sqrt(mr_CIR * mr_CIR + 2.0 * vol_CIR * vol_CIR);
    int prevDate = timeLine->StepDates[0].getDate();
    
    for (int i = 0; i < defRatesDates->size(); i++) {
        int newDate = (*defRatesDates)[i].date.getDate();
        double spread = (*defRatesDates)[i].amount;
        double theta = (newDate - prevDate) / 365.0;
        double chi = exp(gamma * theta) - 1;
        double denom = chi * (gamma + mr_CIR) + 2.0 * gamma;
        chi *= 2.0 / denom;
        double phi = 2.0 * gamma * exp((gamma + mr_CIR) * theta * 0.5);
        phi /= denom;
        lr_CIR_on_curve[i] = (init_CIR * chi - spread * theta) * vol_CIR * vol_CIR;
        lr_CIR_on_curve[i] /= (2.0 * mr_CIR * log(phi));
        prevDate = newDate;
    }

    int k = 0;
    for (int j = 0; j < timeLine->StepDates.size(); j++) {
        if (timeLine->StepDates[j].isGreater((*defRatesDates)[k].date)) {
            k++;
        }
        lr_CIR[j] = lr_CIR_on_curve[k];
    }
}


//----------------------------------------------------------------


/** calculate the stuff which specifies the grid */
//calculate FD grid and its underlying
//alpha : truncation1D

void FD2DEQCR::setFdBounds(double& volForBound1, double alpha1, double& outLowB1, double& outUpB1,
                double& volForBound2, double alpha2, double& outLowB2, double& outUpB2){

    volForBound1 = vols[0];            // vol for spot    
    volForBound2 = vol_CIR;            // vol for credit spread        
    
    //set up default boundaries for FD, can be modified by product at init()    
    double c_LR = lr_CIR[0];        // credit long term
    double k = mr_CIR;                // credit mean reversion
    double sigmaCredit = vol_CIR;    // vol credit spread
    
    //maybe need to combine model and products
    //hard code for now, to change!!!!!!!!!!!!!!!!
    double longest_T = timeLine->TradeTime[timeLine->NumOfStep];

    double temp1 = alpha1* volForBound1 * sqrt(longest_T);
    double temp2 = alpha2* volForBound1 * sqrt(longest_T);

    double undSpot1 = underlying->getSpot();
    double undSpot2 = init_CIR;

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

    if (dim2 == 1) {
        //lBound2 = uBound2 = undSpot2 ;
        outUpB2=outLowB2=undSpot2;
    }
    else {
        if (longest_T < 1/k) {
            temp2 = alpha2 * sigmaCredit * sqrt(undSpot2) * sqrt(longest_T);
            outLowB2 = max(undSpot2 - temp2, DEFAULT_LBOUND2);
            outUpB2 = undSpot2 + temp2;        
        }
        else { //long term
            //hard code for now
            double guessVsup = 0.2; 
            double p_bar = 1 - DEFAULT_P_UND2;
            double deltaVup = 0.01;

            double alpha = 2.0 * k * c_LR / (sigmaCredit * sigmaCredit);
            double beta = alpha / c_LR;

            //double lBound;
            double uBound;
            double cdf;

            outLowB2 = DEFAULT_LBOUND2; 
          
            //obtain the upper boundary
            cdf = 0.0;
            while (cdf < p_bar ){
                guessVsup += deltaVup;                
                uBound = guessVsup * beta;
                cdf = imsl_d_gamma_cdf(uBound, alpha);                
            }
            outUpB2 = guessVsup;        

            //due to same issue as mentioned above,
            //we include this part to make sure that V0 is in the Grid
            double middle = (outLowB2 + outUpB2) / 2.0;

            if ((undSpot2 < outLowB2) || (undSpot2 > outUpB2)){

                if(undSpot2 < middle){
                    double dis = middle - undSpot2;

                    //only make sure that we have at least some points around V0    
                    outLowB2 = Maths::max(outLowB2 - 2.0 * dis, 0.0);
                }else{
                    double dis = undSpot2 - middle ;
                    
                    outUpB2 = Maths::max(outUpB2 + 2.0 * dis, 0.0);
                }
            }
        }       
    }
}

//----------------------------------------------------------------
//  setup the coefficients of the PIDE needed to be solved
//  Assume the most general PDE is 
//  U_t + a*U_x + b*U_y + c*U_xx + d*U_yy + e*U_xy + f*U + Jump(x,y) = 0;
//  this ft is called at each time step, so they're time depd.
//----------------------------------------------------------------

void FD2DEQCR::pdeCoeff(int step, double*** coeff, 
                        int bot1, int top1, int bot2, int top2 ){

    double sigmaSpot = vols[step];        // vol of spot 
    double rho   = correl;                // correlation            
    double c_LR = lr_CIR[step];            // long term default intensity
    double k = mr_CIR;                    // credit mean reversion
    double sigmaCredit = vol_CIR;        // vol of credit
   
    int ix;
    int iy;

    if (step == timeLine->NumOfStep){
        throw ModelException("FD2DEQCR::pdeCoeff", 
            "pdeCoeff shouldn't be called at last step !");
    }
            
    double mu = 0;
    double dt_trading = timeLine->TradeYrFrac[step+1];
    double df = discYC->pv(timeLine->StepDates[step], timeLine->StepDates[step+1]);
    double r = -log(df);

    if (!Maths::isZero(stepForward[0][step])){
        mu = log(stepForward[0][step+1]/stepForward[0][step]);    //[iAsset][step]
    }else{
        throw ModelException("FD2DEQCR::pdeCoeff", 
                             "division by zero.");
    }



//    if(stepForward1[step] ==0 ){
//        // hard coded for now !!!
//        throw ModelException("FD2DEQCR::pdeCoeff", 
//        "division by zero");
//    }else{
//        mu = log(stepForward1[step+1]/stepForward1[step]);    
//    }

    // fpr accessing underlying level, need review !!!
    payoffIndex->update(step, FDProduct::BWD);
    double* und = payoffIndex->getValue( step ).getValues();

    //fill all coeff of PDE
    for (ix = bot1; ix <= top1; ix++){
        for (iy = bot2; iy <= top2; iy++){

            coeff[1][ix][iy] = k * (c_LR - gridLevel2[iy]) * dt_trading;  //b

            coeff[2][ix][iy] = 0.5 * sigmaSpot * sigmaSpot * dt_trading;  //c

            coeff[3][ix][iy] = 0.5 * gridLevel2[iy] * sigmaCredit * sigmaCredit * dt_trading;  //d

            coeff[4][ix][iy] = rho * sigmaSpot * sigmaCredit * sqrt(gridLevel2[iy]) * dt_trading;  //e

            coeff[5][ix][iy] = -(r + dt_trading * gridLevel2[iy]);  //f

            switch (whichChgVarDim1) {
                case X:{//S
                    coeff[0][ix][iy] = mu * und[ix];  //a
                    coeff[2][ix][iy] *= und[ix] * und[ix];  //c
                    coeff[4][ix][iy] *= und[ix];  //e
                }        
                break;
                case LOG_X:{//log(S)
                    coeff[0][ix][iy] = mu - 0.5 * sigmaSpot * sigmaSpot * dt_trading;  //a
                }
                break;
                case LOG_FWDX:{//log(fwd)
                    coeff[0][ix][iy] = - 0.5 * sigmaSpot * sigmaSpot * dt_trading;  //a
                }
                break;
                default:{
                    throw ModelException("NewFD1FLV::pdeCoeff", "This change of variable is not available");
                }
            }
        }
    }
}

//----------------------------------------------------------------


bool FD2DEQCRLoad(){
    return (FD2DEQCR::TYPE && true);
}

DRLIB_END_NAMESPACE
