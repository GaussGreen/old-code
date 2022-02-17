//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : MertonLVCalib.cpp
//
//   Description : MertonLVCalib
//
//   Author      : Francois Lu
//
//   Date        : 2 Feb 2004
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/VolMertonLV.hpp"
#include "edginc/Vanilla.hpp"
#include "edginc/LinearStrikeVolRequest.hpp"
#include "edginc/Black.hpp"
#include "edginc/NonPricingModel.hpp"

DRLIB_BEGIN_NAMESPACE

#if 0
// debug
#define DEBUG_FILE_NAME 
#ifdef DEBUG_FILE_NAME
FILE * fileMertonLV=fopen("C:\\temp\\MertonLV.txt","w");
#endif
#endif

/////////////////////

class MertonLVCalib;
typedef smartPtr<MertonLVCalib> MertonLVCalibSP;

class MertonLVCalib : public CInstrument, public CClosedFormLN::IIntoProduct
{
public:
    static CClassConstSP const TYPE;
    static void load(CClassSP& clazz);
    static IObject* defaultMertonLVCalib(){
        return new MertonLVCalib();
    }
    // override base implementation if required
    virtual void GetMarket(const IModel*, const CMarketDataSP);   
    virtual void Validate();
	
	void setParams(const CDoubleArray&);

	//	compute the forwards
	void setForwards();


	DateTime getValueDate() const {
		return MarketSurface->getBaseDate();
	}

    /** Returns the name of the instrument's discount currency. */
    virtual string discountYieldCurveName() const;

    /** Implementation of ClosedForm::IntoProduct interface */
    virtual CClosedFormLN::IProduct* createProduct(CClosedFormLN* model) const;

	friend void initialisationCallPriceTarget(	//compute the prices of the calls in the market
		vector<vector<double> >&,
		const MertonLVCalib*);
	friend void costFunction(					//penalty function described as in the isml description 
		const CDoubleArray&,					//for non_least_square optimization method
		CDoubleArray&, 
		const MertonLVCalib*, 
		const vector< vector<double> >&);
	friend void calcImpliedMertonVolSurface(	//compute the implied Merton Vol surface
		const vector<vector<double> >&,
		const MertonLVCalib*, vector<vector<double> >&, 
		const CDoubleArray&);
	friend void recalibration(					//to finish the calibration by an other method than levenberg_Maquard
		const vector<vector<double> >& mktPrices, 
		const MertonLVCalib* instTemp, 
		DoubleArray& calibratedParams);

private:
    DoubleMatrix    WeightMatrix;   //Weight Matrix for Calibration
    vector<double> forwardValues; // $unregistered
	// for product to access instrument data
    friend class MertonLVCalibClosedFormProd;  
    MertonLVCalib();
    CAssetWrapper           asset;
    YieldCurveWrapper       discount;
    InstrumentSettlementSP  instSettle;

	bool wantRecalibrateMertonLVCalib;

	// transient field
	VolSurfaceSP		MarketSurface;
	VolMertonLVSP		MertonObject;
};

// helpers
void MertonLVCalib::load(CClassSP& clazz){
    clazz->setPublic();
    REGISTER(MertonLVCalib, clazz);
    SUPERCLASS(CInstrument);
    IMPLEMENTS(CClosedFormLN::IIntoProduct);
    EMPTY_SHELL_METHOD(defaultMertonLVCalib);
    FIELD(asset, "Underlying of option");
    FIELD(discount, "Discount curve");
    FIELD(instSettle, "instrument settlement");

	// calibration params
    FIELD(WeightMatrix,  "Weight Matrix for Calibration.");

    FIELD(MarketSurface,     "market implied vol surface");
	FIELD_MAKE_TRANSIENT(MarketSurface);
    FIELD(MertonObject,     "Merton vol parameters");
	FIELD_MAKE_TRANSIENT(MertonObject);
   
	FIELD(wantRecalibrateMertonLVCalib,  "want Recalibrate MertonLVCalib?");
}

// constructor
MertonLVCalib::MertonLVCalib(): CInstrument(TYPE) {}

void MertonLVCalib::Validate()
{
}


void MertonLVCalib::GetMarket(const IModel*       model, 
								const CMarketDataSP market)
{
    static const string method = "MertonLVCalib::GetMarket";
    try {
        CAsset::getAssetMarketData(model, market.get(), "V", discount, asset);
        discount.getData(model, market);
        instSettle->getMarket(model, market.get());
        // asset and discount curve must not be null 
        if (!asset) {
            throw ModelException(method, "Asset is NULL.");
        }
        if (!discount) {
            throw ModelException(method,  "Discount curve is NULL.");
        }
        // since the assets have changed, we have to validate the instrument
        validatePop2Object();

		// get vol surface specifically for later use
		NonPricingModel nonPricingModel;
		MarketSurface = VolSurfaceSP::dynamicCast(market->GetData(&nonPricingModel, 
										asset->getName(),VolSurface::TYPE));
        if (!MarketSurface) {
            throw ModelException(method, "Asset does not have vol surface.");
        }

		// get Merton+jump vol
		MertonObject = VolMertonLVSP::dynamicCast(market->GetData(&nonPricingModel, 
										asset->getName(), VolMertonLV::TYPE));
        if (!MarketSurface) {
            throw ModelException(method, "Asset does not have vol surface.");
        }

		// validate surface and weight sizes
		const DoubleMatrix& volMatrix = MarketSurface->getVolMatrix();
		if (volMatrix.numRows() != WeightMatrix.numRows())
		{
            throw ModelException(method, "weight matrix has "
										+ Format::toString(WeightMatrix.numRows())
										+ " rows but vol surface has "
										+ Format::toString(volMatrix.numRows())
										+ " rows.");
		}
		if (volMatrix.numCols() != WeightMatrix.numCols())
		{
            throw ModelException(method, "weight matrix has "
										+ Format::toString(WeightMatrix.numCols())
										+" columns but vol surface has "
										+ Format::toString(volMatrix.numCols())
										+" columns.");
		}
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

CClassConstSP const MertonLVCalib::TYPE = CClass::registerClassLoadMethod(
    "MertonLVCalib", typeid(MertonLVCalib), load);

class MertonLVCalibClosedFormProd: public CClosedFormLN::IProduct{
private:
    const MertonLVCalib*  inst; // a reference

public:
    MertonLVCalibClosedFormProd(const MertonLVCalib* instr): inst(instr){}

    void price(CClosedFormLN*   model,
               Control*        control, 
               CResults*       results)const;
	void calibrate() const;


};

/*gives the BS price of a call or put with a modification of the vol and of the forward*/
static double priceBSModified(const DateTime& valueDate,
								const DateTime& matDates,
								bool isCall,
								double strike,
								const InstrumentSettlement* instSettle,
								const Asset* asset,
								const YieldCurve* discount,
								double foward,
								double fwrdModification,
								double varianceModification)
{

    double absStrike = strike;
    double fwdPrice = foward*fwrdModification;
    double variance = varianceModification;
    double discFactor = instSettle->pv(valueDate,
                                       matDates, 
                                       discount, 
                                       asset);
    double strikePremium = Black::price(isCall, fwdPrice, absStrike, discFactor, variance);
	return strikePremium;
}

/*compute the price of a call or put option that follow a Merton process.*/
static double priceMerton(const DateTime& valueDate,
						  const DateTime& matDates,
						  double dt,
						  bool isCall,
						  double strike,
						  const InstrumentSettlement* instSettle,
						  const Asset* asset,
						  const YieldCurve* discount,
						  double forward,
						  double ATMVol,
						  double JumpRate,
						  double JumpMean,
						  double JumpWidth)
{
	const double tolerance=0.001;///if the fraction of price added is less than tolerance time the already estimated price, we stop. 
	const int maxIteration=50;///Maximum of iterations. 
	double MertonPrice=0;
	double increment;
	double factor=1;
	double varianceModification = JumpWidth*JumpWidth;
	double fwrdModification = JumpRate*dt*(1-exp(JumpMean+0.5*varianceModification));
	int i=0;
	do{
		increment=factor*priceBSModified(valueDate,
										 matDates,
										 isCall,
										 strike,
										 instSettle,
										 asset,
										 discount,
										 forward,
										 exp(fwrdModification+i*(JumpMean+0.5*varianceModification)),
										 ATMVol*ATMVol*dt+i*varianceModification);
		MertonPrice+=increment;
		i++;
		factor*=JumpRate*dt/i;
	} while(increment>tolerance*MertonPrice && (i<maxIteration));

	return exp(-JumpRate*dt)*MertonPrice;
}

/*given the parameters of the jumps, this method compute the volatility of the diffusion in a Merton process 
to mach the price of the call or put of price "price"*/
static double impliedMertonVol(const DateTime& valueDate,
								const DateTime& matDates,
								double dt,
								bool isCall,
								double strike,
								const InstrumentSettlement* instSettle,
								const Asset* asset,
								const YieldCurve* discount,
								const CDoubleArray& calibratedParams,
								double price,
								double forward)
{
	bool overprice=false;
	double result=0.1;
	double precision=0.00005;
	double precision2=0.05;
	double Call_Merton;
	int Maxiteration=100;
	int iteration=0;
	Call_Merton=priceMerton(valueDate,
								matDates,
								dt,
								isCall,
								strike,
								instSettle,
								asset,
								discount,
								forward,
								precision,
								calibratedParams[1],
								calibratedParams[2],
								calibratedParams[3]);
	if(Call_Merton>price) return 0.0;
	else{
		while (precision2>precision && (iteration<Maxiteration)){
			iteration++;
			Call_Merton=priceMerton(valueDate,
									matDates,
									dt,
									isCall,
									strike,
									instSettle,
									asset,
									discount,
									forward,
									result,
									calibratedParams[1],
									calibratedParams[2],
									calibratedParams[3]);
			if (Call_Merton<price) {
				result+=precision2;
				if (overprice) precision2/=2.0;
			}
			else {
				overprice=true;
				result-=precision2;
				precision2/=2.0;
			}
		}
		return result;
	}
}

/*fills mktPrices with the prices of the call options given the implied volatility surface*/
void initialisationCallPriceTarget(vector<vector<double> >& mktPrices,const MertonLVCalib*  instTemp)
{
	DateTime valueDate = instTemp->MarketSurface->getBaseDate();
	const DateTimeArray& matDates = instTemp->MarketSurface->getDates();
	const DoubleArray& surfaceStrikes = instTemp->MarketSurface->getStrikes();
	bool isCall = true;

	double variance;
	int i, j;
	for (i=0; i<(int)mktPrices.size(); i++){
		for (j=0; j<(int)mktPrices[i].size(); j++){
			LinearStrikeVolRequestSP volRequest(new LinearStrikeVolRequest(
                                    surfaceStrikes[j], 
                                    valueDate, 
                                    matDates[i],
                                    false));
			CVolProcessedSP vol(instTemp->MarketSurface->getProcessedVol(volRequest.get(),instTemp->asset.get()));
			CVolProcessedBSSP volBS = CVolProcessedBSSP::dynamicCast(vol);
			variance = volBS->CalcVar(valueDate, matDates[i]);
			mktPrices[i][j]=priceBSModified(valueDate,
											matDates[i],
											isCall,
											surfaceStrikes[j],
											instTemp->instSettle.get(),
											instTemp->asset.getSP().get(),
											instTemp->discount.getSP().get(),
											instTemp->forwardValues[i],
											1.0,
											variance);
		}
	}
}

/*fills impliedMertonVolSurface with the implied Merton vols (given the parameters of the jumps)
 given the matrices of the call prices*/
void calcImpliedMertonVolSurface(const vector<vector<double> >& mktPrices, 
								 const MertonLVCalib* instTemp, 
								 vector<vector<double> >& impliedMertonVolSurface, 
                                 const CDoubleArray& calibratedParams){
    DateTime valueDate = instTemp->MarketSurface->getBaseDate();
    const DateTimeArray& matDates = instTemp->MarketSurface->getDates();
    const DoubleArray& surfaceStrikes = instTemp->MarketSurface->getStrikes();
    bool isCall = true;

    int i, j;
	
    // trading yr fraction
    vector<double> dt(matDates.size());
    TimeMetricConstSP timeMetric = instTemp->MarketSurface->getTimeMetric();

    ///debug///
//	fprintf(fileMertonLV, "TEST, ");
//	for (i=0; i<(int)mktPrices[1].size(); i++) fprintf(fileMertonLV, "%f, ",surfaceStrikes[i]);
//	fprintf(fileMertonLV, "\n");

    for (i=0; i<(int)dt.size() ; i++)
    {
        dt[i]= timeMetric->yearFrac(valueDate, matDates[i]);
    }

    for (i=0; i<(int)mktPrices.size(); i++){
//		fprintf(fileMertonLV, "%f, ",dt[i]);
        for (j=0; j<(int)mktPrices[i].size(); j++){
            impliedMertonVolSurface[i][j]=impliedMertonVol(valueDate,
                                                           matDates[i],
                                                           dt[i],
                                                           isCall,
                                                           surfaceStrikes[j],
                                                           instTemp->instSettle.get(),
                                                           instTemp->asset.getSP().get(),
                                                           instTemp->discount.getSP().get(),
                                                           calibratedParams,
                                                           mktPrices[i][j],
                                                           instTemp->forwardValues[i]);
//			fprintf(fileMertonLV, "%f, ",mktPrices[i][j]);
        }
//		fprintf(fileMertonLV, "\n");
    }
//	fclose(fileMertonLV);
}

/*compute the penalty function in the format required by isml for non_least_square optimization method*/
void costFunction(const CDoubleArray& params, CDoubleArray& sumError, const MertonLVCalib*  instTemp, const vector< vector<double> >&mktPrices)
{
	int i, j;

	DateTime valueDate = instTemp->MarketSurface->getBaseDate();
	const DateTimeArray& matDates = instTemp->MarketSurface->getDates();
	const DoubleArray& surfaceStrikes = instTemp->MarketSurface->getStrikes();
	bool isCall = true;

	// trading yr fraction
	vector<double> dt(matDates.size());
	TimeMetricConstSP timeMetric = instTemp->MarketSurface->getTimeMetric();
	for (i=0; i<(int)dt.size() ; i++)
	{
        dt[i]= timeMetric->yearFrac(valueDate, matDates[i]);
	}
	double priceTemp;
	int index=0;

	for (i=0; i<(int)mktPrices.size(); i++) // num of maturities
	{
		for (j=0; j<(int)mktPrices[i].size(); j++) // num of strikes
		{				
			if(instTemp->WeightMatrix[j][i]!=0){
				priceTemp = priceMerton(valueDate,
										 matDates[i],
										 dt[i],
										 isCall,
										 surfaceStrikes[j],
										 instTemp->instSettle.get(),
										 instTemp->asset.getSP().get(),
										 instTemp->discount.getSP().get(),
										 instTemp->forwardValues[i],
										 params[0],
										 params[1],
										 params[2],
										 params[3]);
				sumError[index]=sqrt(instTemp->WeightMatrix[j][i])*(priceTemp-mktPrices[i][j]);
				index++;
			}
		}
	}
}

class MertonLVOptimizer: public MFunctionND{
public:
    virtual void operator()(const CDoubleArray&  x,
                           CDoubleArray& f) const;
	CDoubleArray guessOptimizer;
	const MertonLVCalib* instOptimizer; 
	vector<vector<double> >mktPricesOptimizer;
	MertonLVOptimizer(const MertonLVCalib* , 
						vector<vector<double> >, 
						int NumberVars, 
						int MumberFunctions);
};

MertonLVOptimizer::MertonLVOptimizer(const MertonLVCalib* instTemp,
									 vector<vector<double> >mktPricesTemp,
									 int NumberVars,
									 int MumberFunctions):MFunctionND(NumberVars, MumberFunctions)
{
	instOptimizer=instTemp;
	mktPricesOptimizer=mktPricesTemp;
}

void MertonLVOptimizer::operator()(const CDoubleArray&  x,
										CDoubleArray&  f) const
{	
	costFunction(x, f, instOptimizer, mktPricesOptimizer);
}
 
void MertonLVCalib::setParams(const CDoubleArray& calibratedParams){
    MertonObject->ATMVol=calibratedParams[0];
    MertonObject->JumpRate=calibratedParams[1];
    MertonObject->JumpMean=calibratedParams[2];
    MertonObject->JumpWidth=fabs(calibratedParams[3]);

}

void MertonLVCalib::setForwards(){
	int numberMat, i;
	const DateTimeArray& matDates = MarketSurface->getDates();
	numberMat=matDates.size();
	forwardValues.resize(numberMat);
	for (i=0; i<numberMat;i++){
		forwardValues[i]=asset->fwdValue(matDates[i]);
	}
}


/** Returns the name of the instrument's discount currency */
string MertonLVCalib::discountYieldCurveName() const {
    return discount.getName();
}


/////////////*****			core of the calibration			*********///////////////
void MertonLVCalibClosedFormProd::price(CClosedFormLN*  model,
                                   Control*        control, 
                                   CResults*       results)const
{
    static const string method = "MertonLVCalibClosedFormProd::price";

    try {
		int i,j;
		int numberFunction=0;
		for (i=0; i<inst->WeightMatrix.numRows(); i++){
			for(j=0; j<inst->WeightMatrix.numCols(); j++){
				if (inst->WeightMatrix[j][i]) numberFunction++;
			}
		}
		vector<vector<double> >mktPrices(inst->WeightMatrix.numRows());
		vector<vector<double> >impliedMertonVolSurface(inst->WeightMatrix.numRows());
		for (i=0; i<(int)mktPrices.size(); i++)
		{
			mktPrices[i].resize(inst->WeightMatrix.numCols());
			impliedMertonVolSurface[i].resize(inst->WeightMatrix.numCols());
		}
		const_cast<MertonLVCalib*>(inst)->setForwards();
		initialisationCallPriceTarget(mktPrices, inst);
		CDoubleArray guess(4);        // initial guess
		CDoubleArray calibratedParams(4); // results
		guess[0]=inst->MertonObject->ATMVol;
		guess[1]=inst->MertonObject->JumpRate;
		guess[2]=inst->MertonObject->JumpMean;
		guess[3]=inst->MertonObject->JumpWidth;
		
		if(inst->wantRecalibrateMertonLVCalib){
			MertonLVOptimizer optimizer(inst,mktPrices,4,numberFunction);
			LevenbergMarquardt().minimize(optimizer, guess, calibratedParams);

			const_cast<MertonLVCalib*>(inst)->setParams(calibratedParams);
			calcImpliedMertonVolSurface(mktPrices, inst, impliedMertonVolSurface, calibratedParams);
		}
		else calcImpliedMertonVolSurface(mktPrices, inst, impliedMertonVolSurface, guess);

		DoubleArrayArraySP test(new DoubleArrayArray(inst->WeightMatrix.numRows()));
		for (i=0; i< inst->WeightMatrix.numRows();i++)
		{
			(*test)[i].resize(inst->WeightMatrix.numCols());
			for(j=0; j<inst->WeightMatrix.numCols();j++)
				(*test)[i][j] = impliedMertonVolSurface[i][j];
		}
		IObjectSP result = test;
        OutputNameConstSP outputName(new OutputName("ImpliedMertonVolSurface"));
		results->storeGreek(result, Results::DEBUG_PACKET, outputName);
		
		OutputNameConstSP valueATMVol(new OutputName("ATMVol"));
		results->storeGreek(IObjectSP(CDouble::create(inst->MertonObject->ATMVol)), Results::DEBUG_PACKET, valueATMVol);
		
		OutputNameConstSP valueJumpRate(new OutputName("JumpRate"));
		results->storeGreek(IObjectSP(CDouble::create(inst->MertonObject->JumpRate)), Results::DEBUG_PACKET, valueJumpRate);
		
		OutputNameConstSP valueJumpMean(new OutputName("JumpMean"));
		results->storeGreek(IObjectSP(CDouble::create(inst->MertonObject->JumpMean)), Results::DEBUG_PACKET, valueJumpMean);
		
		OutputNameConstSP valueJumpWidth(new OutputName("JumpWidth"));
		results->storeGreek(IObjectSP(CDouble::create(inst->MertonObject->JumpWidth)), Results::DEBUG_PACKET, valueJumpWidth);

    }
    catch (exception& e) {
        throw ModelException(&e, method);
    }
}
	
/** Implementation of ClosedForm::IntoProduct interface */
CClosedFormLN::IProduct* MertonLVCalib::createProduct(
    CClosedFormLN* model) const{
    return new MertonLVCalibClosedFormProd(this);
}

bool MertonLVCalibLoad()
{
	return (true && MertonLVCalib::TYPE);
}

DRLIB_END_NAMESPACE
