//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SRMInitialGuess.cpp
//
//   Description : Initial Guess for SRM
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Function.hpp"
#include "edginc/Asset.hpp"
#include "edginc/MarketDataFetcher.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/VolProcessedBS.hpp"
#include "edginc/LinearStrikeTSVolRequest.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/VolProcessedDVF.hpp"
#include "edginc/Optimizer.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Maths.hpp"
#include "edginc/VanillaGrid.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/Spline.hpp"
#include "edginc/IRVol.hpp"
#include "edginc/MarketDataFetcherSRM.hpp"
#include "edginc/SwapMaturityVolRequest.hpp"
#include "edginc/SRMInitialGuess.hpp"
#include "edginc/SRMUtil.hpp"
#include "edginc/SRMConstants.hpp"
#include "edginc/SRMRatesUtil.hpp"
#include "edginc/SRMRatesHJMUtil.hpp"
#include "edginc/SRMEquityUtil.hpp"

DRLIB_BEGIN_NAMESPACE

// Class for objective function
class SRMIV_ObjFunc: public MFunctionND{
public:
    virtual void operator()(const CDoubleArray&  x,
                            CDoubleArray&        f) const 
    {
        static const string method = "SRMIV_ObjFunc::operator()";

        if ((x[0]!=x[0])
            ||(x[1]!=x[1])
            ||(x[2]!=x[2])){

            throw ModelException(method, "optimization failed");
        }

        int nbStrike = normStrike.size();
        double s = 0;
        for (int iStrike = 0; iStrike < nbStrike; iStrike++){        
            double A = x[0];
            double B = x[1];
            double Cx = x[2]*normStrike[iStrike];
            double v = 1;
            double v0;
            if (iMat==0){
                v0 = 0;
            } else {
                v0 = atmImpVar[iMat-1]/atmImpVar[iMat];
            }
            double g;
            if (Cx==0){
                g = 1;
            } else {   
                g = 1 
                    + A*A*(1 - (tanh(Cx*v)-tanh(Cx*v0))/(Cx*(v-v0)))
                    + B*B*(1 - 4*(atan(exp(Cx*v))-atan(exp(Cx*v0)))/(Cx*(v-v0))+(tanh(Cx*v)-tanh(Cx*v0))/(Cx*(v-v0)))
                    + 2*A*(log(cosh(Cx*v)/cosh(Cx*v0))/(Cx*(v-v0)))
                    + 2*B*(1-2*(atan(exp(Cx*v))-atan(exp(Cx*v0)))/(Cx*(v-v0)))
                    + 2*A*B*(log(cosh(Cx*v)/cosh(Cx*v0))/(Cx*(v-v0))+(1/cosh(Cx*v)-1/cosh(Cx*v0))/(Cx*(v-v0)));
            }
            double atmVarInc;
            if (iMat==0){
                atmVarInc = atmImpVar[0];                          
            } else {
                atmVarInc = (atmImpVar[iMat]-atmImpVar[iMat-1]);
            }
            s += Maths::square(impVarInc[iStrike]/atmVarInc-g);
        }
        // "insure" constraints
        double l1 = 1+x[1]-fabs(x[0]);
        double l2 = 1+2*x[1]-Maths::square(x[0]);

        double objf =  s/nbStrike + exp(-20.0*l1) + exp(-20.0*l2); 
        if (objf+1==objf){
            f[0] = 1.0e+30;
        } else {
            f[0] = objf;
        }
    }

    SRMIV_ObjFunc(const RangeArray&  defRanges,
                  int                nbVars,
                  int                nbFuncs,
                  DoubleArray impVarInc,
                  DoubleArray atmImpVar,
                  DoubleArray normStrike,
                  int iMat):
        MFunctionND(nbVars, nbFuncs, defRanges),
        impVarInc(impVarInc),
        atmImpVar(atmImpVar),
        normStrike(normStrike),
        iMat(iMat){}
    
private:

    DoubleArray impVarInc;
    DoubleArray atmImpVar;
    DoubleArray normStrike;
    int iMat;
};

class SRMIV_ObjFuncLoc: public MFunctionND{
public:
    virtual void operator()(const CDoubleArray&  x,
                            CDoubleArray&        f) const 
    {
        static const string method = "SRMIV_ObjFuncLoc::operator()";

        if ((x[0]!=x[0])
            ||(x[1]!=x[1])
            ||(x[2]!=x[2])){

            throw ModelException(method, "optimization failed");
        }

        int nbStrike = normStrike.size();
        double s = 0;
        for (int iStrike = 0;iStrike<nbStrike;iStrike++){
            double Cx = x[2]*normStrike[iStrike];
            double g = 1 + x[0]*tanh(Cx) + x[1]*(1-1/cosh(Cx));
            s += Maths::square(locVolN[imat][iStrike]-g);
        }
        // "insure" constraints
        double l1 = 1+x[1]-fabs(x[0]);
        double l2 = 1+2*x[1]-Maths::square(x[0]);

        double objf =  s/nbStrike + exp(-20.0*l1) + exp(-20.0*l2); 
        if (objf+1==objf){
            f[0] = 1.0e+30;
        } else {
            f[0] = objf;
        }
    }

    SRMIV_ObjFuncLoc(const RangeArray&  defRanges,
                     int                nbVars,
                     int                nbFuncs,
                     DoubleMatrix locVolN,
                     DoubleArray normStrike,
                     int imat):
        MFunctionND(nbVars, nbFuncs, defRanges),
        locVolN(locVolN),
        normStrike(normStrike),
        imat(imat){}

private:
    
    CDoubleMatrix locVolN;
    DoubleArray normStrike;
    int imat;
};

// This function calculates the square difference between input function and 
// hyperbolic local volatility model
SRMInitialGuess::SRMInitialGuess():CObject(TYPE){}

// constructor used by the addin
SRMInitialGuess::SRMInitialGuess(
    MarketDataSP    market,
    CAssetWrapper   asset,
    ExpiryArraySP   expiries,
    const DoubleArray&              strike,
    const string&                   volType,
    const string&                   methodType,
    const string&                   volMap):
    CObject(TYPE),
    market(market),
    asset(asset),
    expiries(expiries),
    strike(strike),
    volType(volType),
    methodType(methodType),
    volMap(volMap){

    MarketDataFetcherSP mdf(new MarketDataFetcherLN(volType));
    NonPricingModel model(mdf);
    this->asset.getData(&model, market);
        
    computeGuess();
}

// constructor used by the calibrator
SRMInitialGuess::SRMInitialGuess(
    MarketDataSP market,
    Calibrator::ObjFuncSP objFunc,
    const Calibrator::InstanceIDArray& ids,
    const string& methodType,
    const string& volType):
    CObject(TYPE),
    market(market),
    methodType(methodType),
    volType(volType){
    static const string method = "SRMInitialGuess::SRMInitialGuess";
    try {
        // check we calibrate all the variables
        if (ids.size()!=4){
            throw ModelException("SRM initial guess needs all the parameters to be calibrated");
        }       
                
        if (!CString::equalsIgnoreCase(methodType,COMPMET)
            &&!CString::equalsIgnoreCase(methodType,LOCMET)
            &&!CString::equalsIgnoreCase(methodType,IMPMET)
            &&!CString::equalsIgnoreCase(methodType,DUPIREMET)
            &&!CString::equalsIgnoreCase(methodType,GENMET)){
            throw ModelException("SRM initial guess method is not recognized");
        }       
                
        // chek volType
        if (volType.empty()){
            this->volType = VOLPREF;
        }

        // get name
        CClassConstSP clazz = CClass::forName("Calibrator::InstanceID");
        CFieldConstSP fieldName = clazz->getDeclaredField("name");
        volName = fieldName->getString(ids[0]);
        
        // spotvol or compvol
        CFieldConstSP fieldFieldName = clazz->getDeclaredField("fieldName");
        string name;
        int iId;
        volMap = "";
        for (iId=0; iId<ids.size();iId++){
            name = fieldFieldName->getString(ids[iId]);
            if (CString::equalsIgnoreCase(name,COMPVOL)){
                volMap = COMPVOL;
            }
            if (CString::equalsIgnoreCase(name,SPOTVOL)){
                volMap = SPOTVOL;
            }
        }
        if (CString::equalsIgnoreCase(volMap,"")){
            throw ModelException();
        }
        
        // get expiries
        Calibrator::InstanceIDSP copyId = Calibrator::InstanceID::dbArray("SRMEQ::Vol", volName,volMap);
        expiries = Calibrator::InstanceID::getExpiries(copyId,objFunc);
        

        // get asset and strike
        VanillaGrid::LeastSquareSimpleSP objFuncLS;
        try{
            objFuncLS = VanillaGrid::LeastSquareSimpleSP::dynamicCast(objFunc);
        } catch (exception&) {
            Calibrator::ObjFuncLeastSquareComboSP objFuncLSCombo = Calibrator::ObjFuncLeastSquareComboSP::dynamicCast(objFunc);
            objFuncLS = VanillaGrid::LeastSquareSimpleSP::dynamicCast(objFuncLSCombo->getObjFuncArray()[0]);
        }
                
        string assetName = objFuncLS->getAsset()->getName();
        asset = CAssetWrapper(assetName);
        MarketDataFetcherSP mdf(new MarketDataFetcherLN(this->volType));
        NonPricingModel model(mdf);
        asset.getData(&model, market);

        CDoubleMatrix instStrikesUsed;
        instStrikesUsed = *(objFuncLS->getWeightMatrix()->getInstStrikesUsed());
        int nbStrike = instStrikesUsed.numRows();
        strike = DoubleArray(nbStrike);
        int iStrike;
        for (iStrike=0; iStrike<nbStrike; iStrike++){
            strike[iStrike] = instStrikesUsed[0][iStrike];
        }
        
        // Add validate strike method 
        validateStrike(objFuncLS);

        // compute initial guess
        computeGuess();
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void SRMInitialGuess::validateStrike(VanillaGrid::LeastSquareSimpleSP objFuncLS){
    static const string method = "SRMInitialGuess::validateStrike";

    CDoubleMatrix instStrikesUsed;
    instStrikesUsed = *(objFuncLS->getWeightMatrix()->getInstStrikesUsed());
    int nbStrike = instStrikesUsed.numRows();
    CDoubleMatrix weights;
    weights = *(objFuncLS->getWeightMatrix()->getWeights());

    DateTimeArray maturities = *(objFuncLS->getWeightMatrix()->getDateTimeMaturities());
        
    int nbMat = maturities.size();
    DoubleArray fwdPrice(nbMat);
        
    const double MIN_RATIO = 0.01;
    vector<double> maxWeight(nbMat); // i.e., initialised to zero

    int k = 0,
        i = 0,
        imat = 0;

    for( i = 0 ; i < nbMat ; i++)
    {
        for( k = 0 ; k < nbStrike ; k++) 
        { 
            if( weights[i][k] > maxWeight[i])
            {
                maxWeight[i] = weights[i][k]; 
            }
        }
    }

    int minMat = nbMat;

    for(imat = 0 ; imat < nbMat ; imat++)
    {
        if( maxWeight[imat] >= 0.001)
        {
            minMat = imat;
            break;
        }
    }
        
    if (minMat == nbMat)
    {
        throw ModelException (method , "Invalid weight matrix : some weights must be greater than 0.1%" );
    }

    int t = minMat;
    fwdPrice[t] = asset->fwdValue(maturities[t]);
    double fwd = fwdPrice[t];

    int minIdx = 0, 
        maxIdx = 0;
                
    for(k = 0 ; k < nbStrike ; k++)
    {
        if( weights[t][k] > MIN_RATIO * maxWeight[t])
        {
            minIdx = k;
            break;
        }
    }
        
    for(k = minIdx ; k < nbStrike ; k++)
    {
        if( weights[t][k] > MIN_RATIO * maxWeight[t])
        {
            maxIdx = k;
        }
    }
                
    if( instStrikesUsed[t][minIdx] >= fwd )
    {
        throw ModelException (method , " All strikes are above the forward " );
    }
        
    if( instStrikesUsed[t][maxIdx] <= fwd )
    {
        throw ModelException (method, " All strikes are below the forward " );
    }

}
void SRMInitialGuess::tweakOneParameter(double shift, DoubleArraySP guessParam, 
                                        DoubleArray& shiftGuess, DoubleArraySP& tweakGuess){
    try{
        int nbMat = shiftGuess.size();

        for(int i=0;i<nbMat;++i){
            shiftGuess[i] = shift * (*guessParam)[i];
            (*tweakGuess)[i] = (*guessParam)[i] + shiftGuess[i];
        }
    }
    catch(exception& e){
        throw ModelException(e, "SRMInitialGuess::tweakOneParameter");
    }
    
}

void SRMInitialGuess::getAdjustedGuess(CDoubleMatrix sensitivities, CDoubleMatrix valGuess, DoubleArraySP& adjustedGuess){
    
    try{
        int nbMat = sensitivities.numCols();
        int nbStrikes = sensitivities.numRows();

        for(int t=0;t<nbMat;++t){
            double A = 0.;
            for(int k=0;k<nbStrikes;++k){
                A += Maths::square(sensitivities[t][k]);
            }
            double B = 0.;
            for(int k=0;k<nbStrikes;++k){
                B += sensitivities[t][k] * valGuess[t][k];
            }
            if(A > 0.){
                (*adjustedGuess)[t] -= B/A;
            }
        }
    }
    catch(exception& e){
        throw ModelException(e, "SRMInitialGuess::getAdjustedGuess");
    }
                               
}

void SRMInitialGuess::getObjFuncArray(CDoubleMatrix valGuess,DoubleArraySP& objFuncVals, int& count1, int& count2,int& count3){

    count1 = 0;
    count2 = 0;
    count3 = 0;
    int nbMat = valGuess.numCols();
    int nbStrikes = valGuess.numRows();
   
    for(int t=0;t<nbMat;++t){
        for(int k=0;k<nbStrikes;++k){
            (*objFuncVals)[t] += Maths::square(valGuess[t][k]);
        }
      
    }

    for(int t=0;t<nbMat;++t){
        if(((*objFuncVals)[t] <= 4.5E-06) && ((*objFuncVals)[t] > 4.0E-06)){ ++count1;}
        if(((*objFuncVals)[t] <= 5.0E-06) && ((*objFuncVals)[t] > 4.5E-06)){ ++count2;}
        if(((*objFuncVals)[t] <= 5.5E-06) && ((*objFuncVals)[t] > 5.0E-06)){ ++count3;}
    }
}

void SRMInitialGuess::nbrOfExpiriesToImprove(CDoubleMatrix valGuess,double tol,int& count,DoubleArraySP& myRes){

    int nbMat = valGuess.numCols();
    int nbStrikes = valGuess.numRows();
    myRes = DoubleArraySP(new DoubleArray(nbMat));
    count = 0;

    for(int iMat=0;iMat<nbMat;++iMat){
        for(int k=0;k<nbStrikes;++k){
            (*myRes)[iMat] += Maths::square(valGuess[iMat][k]);
        }
        if( (*myRes)[iMat] > tol){++count;}
    }
}

int SRMInitialGuess::getCenterIndex(VanillaGrid::LeastSquareSimpleSP objFuncLS){
    static const string method = "SRMInitialGuess::getCenterIndex";
    
    try{
        CDoubleMatrix instStrikesUsed;
        instStrikesUsed = *(objFuncLS->getWeightMatrix()->getInstStrikesUsed());
        int nbStrike = instStrikesUsed.numRows();
        CDoubleMatrix weights;
        weights = *(objFuncLS->getWeightMatrix()->getWeights());

        DateTimeArray maturities = *(objFuncLS->getWeightMatrix()->getDateTimeMaturities());
        int nbMat = maturities.size();

        /* get index for the first maturity with non-zero weights */
        vector<double> maxWeight(nbMat); // i.e., initialised to zero

        int k = 0,
            i = 0,
            imat = 0;

        for( i = 0 ; i < nbMat ; i++)
        {
            for( k = 0 ; k < nbStrike ; k++) 
            { 
                if( weights[i][k] > maxWeight[i])
                {
                    maxWeight[i] = weights[i][k]; 
                }
            }
        }

        int minMat = nbMat;

        for(imat = 0 ; imat < nbMat ; imat++)
        {
            if( maxWeight[imat] >= 0.001)
            {
                minMat = imat;
                break;
            }
        }
        /*****/

        double fwdPrice = asset->fwdValue(maturities[minMat]);
        int index = 0;
        for(k = 0; k < nbStrike; ++k){
            if(instStrikesUsed[minMat][k] < fwdPrice){
                index += 1;
            }
        }

        if(index == 0){
            throw ModelException(method, "empty weightMatrix");
        }

        return index;
    }
    catch(exception& e){
        throw ModelException(e, "SRMInitialGuess::getCenterIndex");
    }
}


void SRMInitialGuess::saveGuess(Calibrator::ObjFuncSP objFunc){
    static const string method = "SRMInitialGuess::saveGuess";
    try {
        // format guess
        int nbVar = 4;
        int iMat;
        int nbMat = expiries->size();
        DoubleArray x(nbVar*nbMat);
        for (iMat=0; iMat<nbMat;iMat++){
            x[iMat] = (*atmVol)[iMat];
            x[iMat+nbMat] = (*smileA1)[iMat];
            x[iMat+2*nbMat] = (*smileA2)[iMat];
            x[iMat+3*nbMat] = (*smileA3)[iMat];
        }

        // create id array
        Calibrator::InstanceIDArray idsGuess(nbVar);
        idsGuess[0] = Calibrator::InstanceID::dbArray("SRMEQ::Vol", volName,volMap);
        idsGuess[1] = Calibrator::InstanceID::dbArray("SRMEQ::Vol", volName,"smileA1");
        idsGuess[2] = Calibrator::InstanceID::dbArray("SRMEQ::Vol", volName,"smileA2");
        idsGuess[3] = Calibrator::InstanceID::dbArray("SRMEQ::Vol", volName,"smileA3");
        for (int iId=0; iId<nbVar;iId++){
            Calibrator::InstanceID::initialise(idsGuess[iId], objFunc);
        }

        // save guess
        IObjectSP adjGroup(objFunc->getAdjustableGroup());
        Calibrator::InstanceID::applyAdjustment(idsGuess, adjGroup, x);
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}


void SRMInitialGuess::computeGuess(){
    static const string method = "SRMInitialGuess::computeGuess";
    try{                                                            
        CDoubleMatrixSP result;
        if (CString::equalsIgnoreCase(methodType,LOCMET)) {
            result = smileLocVol();
        } else if (CString::equalsIgnoreCase(methodType,IMPMET)){
            result = smileImpVol();
        } else if (CString::equalsIgnoreCase(methodType,GENMET)) {
            result = smileGen();
        } else if (CString::equalsIgnoreCase(methodType,DUPIREMET)) {
            result = smileDupireRates();
        } else if (CString::equalsIgnoreCase(methodType,COMPMET)){
            // composite method
            try{
                result = smileDupireRates();
                validateGuess(result);
            }catch(exception&){
                try{
                    result = smileLocVol();
                    validateGuess(result);
                }catch(exception&){
                    try{
                        result = smileImpVol();
                        validateGuess(result);
                    }catch(exception&){
                        result = smileGen();
                    }
                }
            }
        }
        else {
            throw ModelException("methodType is not recognized");
        }
        mapGuess(result);
    } 
    catch(exception& e){
        throw ModelException(e, method);
    }
}


CDoubleMatrixSP SRMInitialGuess::smileImpVol(){
    static const string method = "SRMInitialGuess::smileImpVol";
    try {
        DateTime refDate = market->GetReferenceDate();
        DateTime lastDate = expiries->back()->toDate(refDate);
                
        int nbMaturities = expiries->size();
        int nbStrike = strike.size();
                
        CDoubleMatrixSP result(new CDoubleMatrix(nbMaturities,4));
        DoubleArray atmImpVol(nbMaturities);
        DoubleArray atmImpVar(nbMaturities);
        DoubleArray normStrike(nbStrike);
        DoubleMatrix impVar(nbMaturities, nbStrike);
        DoubleArray bootVar(nbStrike);
        DoubleArray impVarInc(nbStrike);
                
        /*get atm implied volatilities*/
        int iStrike, iMat;
        LinearStrikeTSVolRequest req(asset->fwdValue(refDate), refDate, lastDate, false);
        CVolProcessedBSSP atmvolproc(asset->getProcessedVol(&req));
        for (iMat=0;iMat<nbMaturities;iMat++){
            DateTime mat = (*expiries)[iMat]->toDate(refDate);   
            atmImpVol[iMat]= atmvolproc->CalcVol(refDate, mat);
            atmImpVar[iMat]= refDate.yearFrac(mat)* ::pow(atmImpVol[iMat],2);
            // some validation
            if (!Maths::isPositive(atmImpVar[iMat])){
                throw ModelException("variance is negative");                           
            }
            if (iMat>0){
                if (!Maths::isPositive(atmImpVar[iMat]-atmImpVar[iMat-1])){
                    throw ModelException("variance is not increasing");     
                }
            }
        }
                
        /*get implied volatilities*/
        for (iStrike=0;iStrike<nbStrike;iStrike++){
            LinearStrikeTSVolRequest req(strike[iStrike], refDate, lastDate, false);
            CVolProcessedBSSP volproc(asset->getProcessedVol(&req));
            for (iMat=0;iMat<nbMaturities;iMat++){
                DateTime mat = (*expiries)[iMat]->toDate(refDate);
                impVar[iMat][iStrike]= refDate.yearFrac(mat)* ::pow(volproc->CalcVol(refDate, mat),2);
            }
        }

        /*normalize strikes*/
        double spotPrice = asset->fwdValue(refDate);
        for (iStrike=0;iStrike<nbStrike;iStrike++){
            normStrike[iStrike] = log(strike[iStrike]/spotPrice);
        }
                
        /*set the optimizer*/
        DoubleArray x(3);
        DoubleArray guess(3);
        guess[0] = -1.0;
        guess[1] = 1.0;
        guess[2] = 1.0;
        RangeArray ranges(3);
        ranges[0] = RangeSP(new Range (OpenBoundary(-50.0), OpenBoundary(50.0)));
        ranges[1] = RangeSP(new Range (OpenBoundary(0.01), OpenBoundary(100.0)));
        ranges[2] = RangeSP(new Range (OpenBoundary(0.01), OpenBoundary(100.0)));

        /*launch bootstrap and optimization*/
        int iMatBoot;
        for (iMat=0;iMat<nbMaturities;iMat++){
            /*bootstrap*/
            for (iStrike=0; iStrike<nbStrike; iStrike++){
                if (iMat==0){
                    bootVar[iStrike] = 0;
                } else {   
                    bootVar[iStrike] = 0;
                    for (iMatBoot = 0; iMatBoot<iMat; iMatBoot++){
                        double v = atmImpVar[iMatBoot]/atmImpVar[iMat];
                        double v0;
                        if (iMatBoot==0){
                            v0 = 0;
                        } else {
                            v0 = atmImpVar[iMatBoot-1]/atmImpVar[iMatBoot];
                        }
                        double A = (*result)[iMatBoot][1];
                        double B = (*result)[iMatBoot][2];
                        double Cx = (*result)[iMatBoot][3]*normStrike[iStrike];
                        double g;
                        if (Cx==0){
                            g = 1;
                        } else {   
                            g = 1 
                                + A*A*(1 - (tanh(Cx*v)-tanh(Cx*v0))/(Cx*(v-v0)))
                                + B*B*(1 - 4*(atan(exp(Cx*v))-atan(exp(Cx*v0)))/(Cx*(v-v0))+(tanh(Cx*v)-tanh(Cx*v0))/(Cx*(v-v0)))
                                + 2*A*(log(cosh(Cx*v)/cosh(Cx*v0))/(Cx*(v-v0)))
                                + 2*B*(1-2*(atan(exp(Cx*v))-atan(exp(Cx*v0)))/(Cx*(v-v0)))
                                + 2*A*B*(log(cosh(Cx*v)/cosh(Cx*v0))/(Cx*(v-v0))+(1/cosh(Cx*v)-1/cosh(Cx*v0))/(Cx*(v-v0)));
                        }
                        if (iMatBoot==0) {
                            bootVar[iStrike] += atmImpVar[iMatBoot]*g;
                        } else {
                            bootVar[iStrike] += (atmImpVar[iMatBoot]-atmImpVar[iMatBoot-1])*g;
                        }
                    }
                }
                impVarInc[iStrike] = impVar[iMat][iStrike] - bootVar[iStrike];
            }
                        
            /*optimization*/
            SRMIV_ObjFunc srmobjFunc(ranges,3,1,impVarInc,atmImpVar,normStrike,iMat);
            QuasiNewton qn;
            qn.minimize(srmobjFunc,guess,x);
            (*result)[iMat][1] = x[0];
            (*result)[iMat][2] = x[1]; 
            (*result)[iMat][3] = x[2];
        }
                
        if (CString::equalsIgnoreCase(volMap,COMPVOL)){
            for (iMat=0; iMat<nbMaturities; iMat++){
                (*result)[iMat][0] = atmImpVol[iMat];
            }
        } else {
            (*result)[0][0] = atmImpVol[0];
            for (iMat=1; iMat<nbMaturities; iMat++){
                (*result)[iMat][0] = sqrt((atmImpVar[iMat]-atmImpVar[iMat-1])/
                                          ((*expiries)[iMat-1]->toDate(refDate)).yearFrac((*expiries)[iMat]->toDate(refDate)));
            }
        }
        return result;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

CDoubleMatrixSP SRMInitialGuess::smileLocVol() {
    static const string method = "SRMInitialGuess::smileLocVol";
    try {
        DateTime refDate = market->GetReferenceDate();
        DateTime lastDate = expiries->back()->toDate(refDate);
                
        int nbMaturities = expiries->size();
        int nbStrike = strike.size();
                
        CDoubleMatrixSP result(new CDoubleMatrix(nbMaturities,4));
        DoubleArray atmlocVol(nbMaturities);
        DoubleMatrix locVol(nbMaturities, nbStrike);
        DoubleMatrix locVolN(nbMaturities, nbStrike);
        DoubleArray normStrike(nbStrike);

        LocVolRequest lvolreq(refDate,false,false,false,true,0.005,0.001,0.01); 
        IVolProcessedSP proc(asset->getProcessedVol(&lvolreq));
        CVolProcessedDVFSP lvolproc(CVolProcessedDVFSP::dynamicCast(proc));
                
        /*get local volatilities*/
        int iStrike, iMat;
        double spotPrice = asset->fwdValue(refDate);
        for (iMat=0;iMat<nbMaturities;iMat++){
            DateTime mat = (*expiries)[iMat]->toDate(refDate);
            atmlocVol[iMat]= lvolproc->CalcLocVol( mat, spotPrice);
            if (!Maths::isPositive(atmlocVol[iMat])){
                throw ModelException("Local vol is negative");                          
            }
            for (iStrike=0;iStrike<nbStrike;iStrike++){
                locVol[iMat][iStrike]= lvolproc->CalcLocVol( mat, strike[iStrike]);
            }
        }

        /*normalize local volatility*/
        for (iMat=0;iMat<nbMaturities;iMat++){
            for (iStrike=0;iStrike<nbStrike;iStrike++){
                locVolN[iMat][iStrike] = locVol[iMat][iStrike]/atmlocVol[iMat];
            }
        }
                
        /*normalize strikes*/
        for (iStrike=0;iStrike<nbStrike;iStrike++){
            normStrike[iStrike] = log(strike[iStrike]/spotPrice);
        }

        /*set the optimizer*/
        DoubleArray x(3);
        DoubleArray guess(3);
        guess[0] = -1.0;
        guess[1] = 1.0;
        guess[2] = 1.0;
        RangeArray ranges(3);
        ranges[0] = RangeSP(new Range (OpenBoundary(-50.0), OpenBoundary(50.0)));
        ranges[1] = RangeSP(new Range (OpenBoundary(0.01), OpenBoundary(100.0)));
        ranges[2] = RangeSP(new Range (OpenBoundary(0.01), OpenBoundary(100.0)));

        /*optimization*/
        for (iMat=0;iMat<nbMaturities;iMat++){
            SRMIV_ObjFuncLoc srmobjFuncLoc(ranges,3,1,locVolN,normStrike,iMat);
            QuasiNewton qn;
            qn.minimize(srmobjFuncLoc,guess,x);
            (*result)[iMat][1] = x[0];
            (*result)[iMat][2] = x[1]; 
            (*result)[iMat][3] = x[2]; 
        }
                        
        if (CString::equalsIgnoreCase(volMap,COMPVOL)){
            LinearStrikeTSVolRequest req(spotPrice, refDate, lastDate, false);
            CVolProcessedBSSP volproc(asset->getProcessedVol(&req));
            for (iMat=0;iMat<nbMaturities;iMat++){
                DateTime mat = (*expiries)[iMat]->toDate(refDate);
                (*result)[iMat][0] = volproc->CalcVol(refDate, mat);
            }
        } else {
            for (iMat=0; iMat<nbMaturities; iMat++){
                (*result)[iMat][0] = atmlocVol[iMat];
            }
        }
        return result;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

CDoubleMatrixSP SRMInitialGuess::smileGen(){
    static const string method = "SRMInitialGuess::smileGen";
    try{
        DateTime refDate = market->GetReferenceDate();
        DateTime lastDate = expiries->back()->toDate(refDate);
        int nbMaturities = expiries->size();
        CDoubleMatrixSP result(new CDoubleMatrix(nbMaturities,4));
        double C0 = 0.5;
        double C1 = 10;
        double C2 = 2;
        int iMat;
        for (iMat=0; iMat<nbMaturities; iMat++){
            (*result)[iMat][1] = -0.5;
            (*result)[iMat][2] = 1.0;
            (*result)[iMat][3] = C0+C1*exp(-C2*refDate.yearFrac((*expiries)[iMat]->toDate(refDate)));
        }

        DoubleArray atmImpVol(nbMaturities);
        DoubleArray atmImpVar(nbMaturities);
        LinearStrikeTSVolRequest req(asset->fwdValue(refDate), refDate, lastDate, false);
        CVolProcessedBSSP atmvolproc(asset->getProcessedVol(&req));
        for (iMat=0; iMat<nbMaturities; iMat++){
            DateTime mat = (*expiries)[iMat]->toDate(refDate);   
            atmImpVol[iMat]= atmvolproc->CalcVol(refDate, mat);
            atmImpVar[iMat]= refDate.yearFrac(mat)* ::pow(atmImpVol[iMat],2);
        }

        if (CString::equalsIgnoreCase(volMap,COMPVOL)){
            for (iMat=0; iMat<nbMaturities; iMat++){
                (*result)[iMat][0] = atmImpVol[iMat];
            }
        } else {
            (*result)[0][0] = atmImpVol[0];
            for (iMat=1; iMat<nbMaturities; iMat++){
                (*result)[iMat][0] = sqrt((atmImpVar[iMat]-atmImpVar[iMat-1])/
                                          ((*expiries)[iMat-1]->toDate(refDate)).yearFrac((*expiries)[iMat]->toDate(refDate)));
            }
        }
        return result;
    }
    catch(exception& e){
        throw ModelException(e, method);
    }
}

void SRMInitialGuess::validateGuess(CDoubleMatrixSP result)
{
    static const string method = "SRMInitialGuess::validateGuess";
                
    try
    {
        int nbMaturities = expiries->size();
        int iMat, iRes;
                
        // Numbers ?
        for (iMat=0; iMat<nbMaturities; iMat++)
        {
            for (iRes = 0; iRes < 4; iRes++)
            {
                if (!Maths::finite((*result)[iMat][iRes]))
                {
                    throw ModelException();
                }

                // somewhat arbitrary bound
                if (fabs((*result)[iMat][iRes]) > 1000.)
                {
                    throw ModelException();
                }
            }                                       
        }

        // Constraints ?
        for (iMat=0; iMat<nbMaturities; iMat++)
        {
            double B = (*result)[iMat][2];

            if (!Maths::isPositive(B))
            {
                throw ModelException();
            }
        }
    }
    catch(exception& e)
    {
        throw ModelException(e, method);
    }
}

void SRMInitialGuess::mapGuess(CDoubleMatrixSP result)
{
    int nbMat = expiries->size();
    atmVol = DoubleArraySP(new DoubleArray(nbMat));
    smileA1 = DoubleArraySP(new DoubleArray(nbMat));
    smileA2 = DoubleArraySP(new DoubleArray(nbMat));
    smileA3 = DoubleArraySP(new DoubleArray(nbMat));

    // map guess
    for (int iMat=0; iMat<nbMat; iMat++)
    {
        (*atmVol)[iMat] = (*result)[iMat][0];
        (*smileA1)[iMat] = (*result)[iMat][1]*(*result)[iMat][3];
        (*smileA2)[iMat] = (*result)[iMat][2]*(*result)[iMat][3]*(*result)[iMat][3]; 
        (*smileA3)[iMat] = fabs((*result)[iMat][1]) + fabs((*result)[iMat][2]); 
    }
                
    // collar
    DateTime refDate = market->GetReferenceDate();
    for (int iMat=0; iMat<nbMat; iMat++)
    {
        double t = refDate.yearFrac( ((*expiries)[iMat])->toDate(refDate) );
        double maxCurv = Maths::max(2., 1000./(1.+t*t));  

        (*atmVol)[iMat]  = Maths::collar((*atmVol)[iMat],  0.9, 0.02);
        (*smileA1)[iMat] = Maths::collar((*smileA1)[iMat], 14.5, -14.5);
        (*smileA2)[iMat] = Maths::collar((*smileA2)[iMat], maxCurv, 0.01);
        (*smileA3)[iMat] = Maths::collar((*smileA3)[iMat], 5.0, 0.01);
    }
}

void SRMInitialGuess::validate(CDoubleMatrixSP result){
    static const string method = "SRMInitialGuess::validate";

    try{
        int nbMat = expiries->size();

        for(int iMat = 0; iMat < nbMat; ++iMat){
            (*result)[iMat][0] = Maths::collar((*result)[iMat][0], 14.5, -14.5);
            (*result)[iMat][1] = Maths::collar((*result)[iMat][1], 9000., 0.001);
            (*result)[iMat][2] = Maths::collar((*result)[iMat][2], 5., 0.01);
            (*result)[iMat][3] = Maths::collar((*result)[iMat][3], 0.9, 0.02);
        }

        //adjust parameters
        for(int iMat = 0; iMat < nbMat; ++iMat){
            (*smileA1)[iMat] = (*result)[iMat][0];
            (*smileA2)[iMat] = (*result)[iMat][1];
            (*smileA3)[iMat] = (*result)[iMat][2];
            (*atmVol)[iMat]  = (*result)[iMat][3];
        }
    }catch(exception&){
        throw ModelException(method, "SRMInitialGuess::validate failed !");
    }
        
}

 
ExpiryArraySP SRMInitialGuess::getExpiries(){
    return expiries;
}

int SRMInitialGuess::getNbExpiries(){
    return expiries->size();
}

string SRMInitialGuess::getVolMap(){
    return volMap;
}
DoubleArraySP SRMInitialGuess::getAtmVol(){
    return atmVol;
}

void SRMInitialGuess::adjustAtmVol(DoubleArraySP guessAtmVal){
    atmVol = guessAtmVal;
}

void SRMInitialGuess::adjustSmileA1(DoubleArraySP guessSmileA1){
    smileA1 = guessSmileA1;
}

void SRMInitialGuess::adjustSmileA2(DoubleArraySP guessSmileA2){
    smileA2 = guessSmileA2;
}


DoubleArraySP SRMInitialGuess::getSmileA1(){
    return smileA1;
}

DoubleArraySP SRMInitialGuess::getSmileA2(){
    return smileA2;
}

DoubleArraySP SRMInitialGuess::getSmileA3(){
    return smileA3;
}

int SRMInitialGuess::getNbStrikes(){
    return strike.size();
}

CClassConstSP const SRMInitialGuess::TYPE =
CClass::registerClassLoadMethod("SRMInitialGuess", typeid(SRMInitialGuess), load);

// Invoked when Class is 'loaded'
void SRMInitialGuess::load(CClassSP& clazz){
    REGISTER(SRMInitialGuess, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultSRMInitialGuess);
    FIELD(asset, "");
    FIELD(market, "");
    FIELD(expiries, "");
    FIELD(strike, "");
    FIELD(volType, "");
    FIELD_MAKE_OPTIONAL(volType);
    FIELD(methodType, "");
    FIELD_MAKE_OPTIONAL(methodType);
    FIELD(volMap, "");        
    FIELD_MAKE_OPTIONAL(volMap);
}

IObject* SRMInitialGuess::defaultSRMInitialGuess(){
    return new SRMInitialGuess();
}


const string SRMInitialGuess::SPOTVOL = "spotVol";
const string SRMInitialGuess::COMPVOL = "compVol";
const string SRMInitialGuess::COMPMET = "Composite";
const string SRMInitialGuess::IMPMET = "Implied Vol";
const string SRMInitialGuess::LOCMET = "Local Vol";
const string SRMInitialGuess::DUPIREMET = "RatesDupire";
const string SRMInitialGuess::GENMET = "Generic";
const string SRMInitialGuess::VOLPREF = "VolPreferred";


class SRMInitialGuessAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    MarketDataSP market;
    CAssetWrapper asset;
    ExpiryArraySP expiries;
    DoubleArray strike;
    string volType;
    string methodType;
    string volMap;

private:
    CDoubleMatrixSP computeGuessAddin(){
        static const string method = "SRMInitialGuessAddin::computeGuessAddin";
        try{
                                
            SRMInitialGuess guess(market,
                                  asset,
                                  expiries,
                                  strike,
                                  volType,
                                  methodType,
                                  volMap);                      

            int nbMat = guess.getNbExpiries();
            CDoubleMatrixSP result(new CDoubleMatrix(nbMat,4));
            DoubleArraySP atmVol = guess.getAtmVol();
            DoubleArraySP smileA1 = guess.getSmileA1();
            DoubleArraySP smileA2 = guess.getSmileA2();
            DoubleArraySP smileA3 = guess.getSmileA3();

            for (int iMat=0; iMat<nbMat;iMat++){
                (*result)[iMat][0]= (*atmVol)[iMat];
                (*result)[iMat][1]= (*smileA1)[iMat];
                (*result)[iMat][2]= (*smileA2)[iMat];
                (*result)[iMat][3]= (*smileA3)[iMat];
            }
            return result;
        } catch(exception& e){
            throw ModelException(e, method);
        }
    }


    static IObjectSP callcomputeGuessAddin(SRMInitialGuessAddin* params) {
        return params->computeGuessAddin();
    }

    // for reflection
    SRMInitialGuessAddin():
        CObject(TYPE), 
        volType(SRMInitialGuess::VOLPREF),
        methodType(SRMInitialGuess::COMPMET),
        volMap(SRMInitialGuess::SPOTVOL) {}

    // Invoked when Class is 'loaded'
    static void load(CClassSP& clazz)
    {
        REGISTER(SRMInitialGuessAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultSRMInitialGuessAddin);
        FIELD(asset, "");
        FIELD(market, "");
        FIELD(expiries, "");
        FIELD(strike, "");
        FIELD(volType, "");
        FIELD_MAKE_OPTIONAL(volType);
        FIELD(methodType, "");
        FIELD_MAKE_OPTIONAL(methodType);
        FIELD(volMap, "");        
        FIELD_MAKE_OPTIONAL(volMap);

        Addin::registerClassObjectMethod(
            "SRMInitialGuess",
            Addin::MARKET,
            "compute initial values for SRM calibration",
            TYPE,
            false,
            Addin::expandMulti,
            (Addin::ObjMethod*)callcomputeGuessAddin);
    }

    static IObject* defaultSRMInitialGuessAddin(){
        return new SRMInitialGuessAddin();
    }
};

CClassConstSP const SRMInitialGuessAddin::TYPE =
CClass::registerClassLoadMethod("SRMInitialGuessAddin", typeid(SRMInitialGuessAddin), load);


bool SRMInitialGuessLoad() {
    return (SRMInitialGuessAddin::TYPE != 0);
}

DRLIB_END_NAMESPACE
