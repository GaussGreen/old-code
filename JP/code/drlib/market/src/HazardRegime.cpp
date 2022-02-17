#include "edginc/config.hpp"
#define QLIB_HAZARDREGIME_CPP
#include "edginc/HazardRegime.hpp"
#include "edginc/Addin.hpp"
#include "edginc/CreditSpreadCurve.hpp"
#include "edginc/CDSParSpreads.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/Results.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/MultiRegimeFactor.hpp"

DRLIB_BEGIN_NAMESPACE

void HazardRegime::validatePop2Object(){
	const static string method = "HazardRegime::validatePop2Object";

	try{
		validate();

        // assume the last state is default state
        for (int i=0; i<states.size(); i++){
            if (!Maths::isZero(generator[states.size()-1][i]-states[i])){
                throw ModelException(method, "the last column of generator should be equal to states.");
            }
            if (!Maths::isZero(generator[i][states.size()-1])){
                throw ModelException(method, "the last row of generator should all be zero.");
            }
        }

		if (recoveryRate < 0. || recoveryRate > 1.){
			throw ModelException(method, "recoveryRate is not between 0% and 100%.");
		}

		int nRows = generator.numRows();
		int nCols = generator.numCols();
		if (nRows != nCols){
			throw ModelException(method, "generator must be a square matrix.");
		}
        int size = (nRows-1)*(nRows-2);
		generatorEntryForCalib.resize(size);
        rowIdx.resize(size);
        colIdx.resize(size);
		int iEntry = 0;
		for (int irow = 0; irow < nRows-1; ++irow){
			for (int icol = 0; icol < nCols-1; ++icol){
				if (irow != icol){
					generatorEntryForCalib[iEntry] = generator[icol][irow];
					rowIdx[iEntry] = irow;
					colIdx[iEntry] = icol;
					iEntry++;
				}
			}
		}

        size = nRows - 2;
        jumpEntryForCalib.resize(size);
        rowIdxJump.resize(size);
        colIdxJump.resize(size);
		for (int i = 0; i < nRows-2; ++i){
			jumpEntryForCalib[i] = jumpsOnLogSpot[i][i+1];
			rowIdxJump[i] = i+1;
			colIdxJump[i] = i;
		} 
	} catch (exception& e) {
		throw ModelException(e, method);
	}
}

double	HazardRegime::calcRiskyFactor(double tradYear) const {
	DoubleMatrix	tmpMatrix(nbStates-1, nbStates-1);
    for (int iRow = 0; iRow < nbStates-1; iRow++){
        for (int iCol = 0; iCol < nbStates-1; iCol++){
            tmpMatrix[iCol][iRow] = generator[iCol][iRow];
        }
        tmpMatrix[iRow][iRow] += generator[nbStates-1][iRow];
        tmpMatrix[iRow][iRow] -= (1. - recoveryRate) * states[iRow];
    }
    DoubleMatrix tmpMatrix2 = tmpMatrix.exp(tradYear, -1.0, 7);
    double			factor = 0.;
	for (int iCol = 0; iCol < nbStates-1; ++iCol){
		factor += tmpMatrix2[iCol][initStateIdx];
	}
	return factor;
}

double HazardRegime::calcCreditSpread(double tradYear) const {
	double factor = calcRiskyFactor(tradYear);
	double spread = - log(factor)/tradYear;
    return spread;
}

DoubleArray HazardRegime::calcCreditSpread(DoubleArray tradYears) const {
	DoubleArray		creditSpread(tradYears.size());
	int i = 0;
	for (; i < tradYears.size(); ++i){
		creditSpread[i] = -log(calcRiskyFactor(tradYears[i]))/tradYears[i];
	}
	return creditSpread;
}

// Linetsky's dependency power
double HazardRegime::calcEquityCreditDependency(int stateIdx, int whichMethod) const{
    const static string method = "HazardRegime::calcEquityCreditDependency";
        
    int i;
    double sum1 = 0.;
    double sum2 = 0;

    switch (whichMethod){
        // defined as p = d(lnS)/d(lnCS) = (dS/S)/(dCS/CS). 
        // = (E[dS]/S)/(E[dCS]/CS)
        case 0:
            for (i = 0; i < nbStates-1; i++){
                sum1 += generator[i][stateIdx] * (states[i] - states[stateIdx]);
            }
            for (i=0; i < nbStates-1; i++){
                sum2 += generator[i][stateIdx] * (exp(jumpsOnLogSpot[i][stateIdx])-1.0);
            }
            return sum2 / (sum1 / states[stateIdx]);
        case 1:
            // = E[dlnS]/E[dlnCS]
            for (i = 0; i < nbStates-1; i++){
                sum1 += generator[i][stateIdx] * (log(states[i]) - log(states[stateIdx]));
            }
            for (i=0; i < nbStates-1; i++){
                sum2 += generator[i][stateIdx] * jumpsOnLogSpot[i][stateIdx];
            }
            return sum2 / sum1;
        case 2:
            // = E[dlnS] / (E[dCS]/CS)
            for (i = 0; i < nbStates-1; i++){
                sum1 += generator[i][stateIdx] * (states[i] - states[stateIdx]);
            }
            for (i=0; i < nbStates-1; i++){
                sum2 += generator[i][stateIdx] * jumpsOnLogSpot[i][stateIdx];
            }
            return sum2 / (sum1 / states[stateIdx]);
        default:
            throw ModelException(method, "whichMethod must be 0, 1 or 2");
            break;
    }
}

IObject* HazardRegime::defaultCtor(){
    return new HazardRegime();
}

/** Invoked when Class is 'loaded' */
void HazardRegime::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(HazardRegime, clazz);
    SUPERCLASS(RegimeFactor);
    IMPLEMENTS(Calibrator::IAdjustable);
    EMPTY_SHELL_METHOD(defaultCtor);

	FIELD(recoveryRate, "recovery rate");
    FIELD(jumpAsymmetry, "jump size can be different when credit jumps from state i to j or from j to i");
    FIELD_MAKE_OPTIONAL(jumpAsymmetry);

	FIELD(generatorEntryForCalib, "independent entries of generator written in array form");
	FIELD_MAKE_TRANSIENT(generatorEntryForCalib);
    FIELD(rowIdx, "independent entries of generator written in array form");
	FIELD_MAKE_TRANSIENT(rowIdx);
    FIELD(colIdx, "independent entries of generator written in array form");
	FIELD_MAKE_TRANSIENT(colIdx);

    FIELD(jumpEntryForCalib, "independent entries of jumpOnLogSpot written in array form");
	FIELD_MAKE_TRANSIENT(jumpEntryForCalib);
    FIELD(rowIdxJump, "independent entries of jumpOnLogSpot written in array form");
	FIELD_MAKE_TRANSIENT(rowIdxJump);
    FIELD(colIdxJump, "independent entries of jumpOnLogSpot written in array form");
	FIELD_MAKE_TRANSIENT(colIdxJump);

	// fields to be calibrated
	Calibrator::IAdjustable::registerField(
        clazz, "generatorEntryForCalib",
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));

    Calibrator::IAdjustable::registerField(
        clazz, "jumpEntryForCalib",
        new Range(ClosedBoundary(0.0), Infinity(Infinity::Plus)));
}

HazardRegime::HazardRegime(): 
RegimeFactor(TYPE){}

/** Called after (calibrator) adjustments have been made to fields */
void HazardRegime::update(){
    for (int iRow = 0; iRow < generator.numRows(); iRow++){
        generator[generator.numCols()-1][iRow] = states[iRow];
    }
	for (int i=0; i < generatorEntryForCalib.size(); ++i){
		generator[colIdx[i]][rowIdx[i]] = generatorEntryForCalib[i];
	}
    for (int iRow = 0; iRow < generator.numRows(); iRow++){
        generator[iRow][iRow] = 0.;
        double sum = 0.;
        for (int iCol = 0; iCol < generator.numCols(); iCol++){
            sum += generator[iCol][iRow];
        }
        generator[iRow][iRow] = -sum;
    }

    // assume jump additive 
    for (int i=0; i < jumpEntryForCalib.size(); ++i){
		jumpsOnLogSpot[colIdxJump[i]][rowIdxJump[i]] = jumpEntryForCalib[i];
        jumpsOnLogSpot[rowIdxJump[i]][colIdxJump[i]] = jumpAsymmetry * jumpEntryForCalib[i];
	}
    for (int iRow = 0; iRow < jumpsOnLogSpot.numRows()-1; iRow++){
        for (int iCol = iRow-2; iCol >= 0; iCol--){
            jumpsOnLogSpot[iCol][iRow] = jumpsOnLogSpot[iCol+1][iRow] + jumpsOnLogSpot[iCol][iCol+1];
        }
        for (int iCol = iRow+2; iCol < jumpsOnLogSpot.numRows()-1; iCol++){
            jumpsOnLogSpot[iCol][iRow] = jumpsOnLogSpot[iCol-1][iRow] + jumpsOnLogSpot[iCol][iCol-1];
        }
    }
}

/** Called after adjustments have been made to fields (eg calibrator) */
void HazardRegime::fieldsUpdated(const CFieldArray& fields){
    update();
}

DoubleArray	HazardRegime::getGeneratorEntryForCalib() const {
	return generatorEntryForCalib;
}

DoubleArray	HazardRegime::getJumpEntryForCalib() const {
	return jumpEntryForCalib;
}

string HazardRegime::getName() const {
	return name;
}

double HazardRegime::getRecoveryRate() const {
    return recoveryRate;
}

double	HazardRegime::getGenerator(int fromState, int toState) const {
	return generator[toState][fromState];
}

double	HazardRegime::getJumpsOnLogSpot(int fromState, int toState) const {
	return jumpsOnLogSpot[toState][fromState];
}

double	HazardRegime::getState(int iState) const {
	return states[iState];
}

void    HazardRegime::setJumpAsymmetry(double value){
    jumpAsymmetry = value;
}

CClassConstSP const HazardRegime::TYPE =
CClass::registerClassLoadMethod("HazardRegime", typeid(HazardRegime), load);

DEFINE_TEMPLATE_TYPE(HazardRegimeWrapper);

/* external symbol to allow class to be forced to be linked in */
bool HazardRegimeLinkIn(){
    return (HazardRegime::TYPE != 0);  
}

class HazardRegimeAddin: public CObject{
    static CClassConstSP const TYPE;

    HazardRegimeSP		regime;
	DoubleArray			tradYears;

	static IObjectSP calc(HazardRegimeAddin* params) {
        return params->calcCreditSpread();
    }

    IObjectSP calcCreditSpread(){
		DoubleArraySP output(new DoubleArray(tradYears.size()));
        *output = regime->calcCreditSpread(tradYears);
        return IObjectSP(output);
    }

    HazardRegimeAddin(): CObject(TYPE) {}

    static void load(CClassSP& clazz){
        REGISTER(HazardRegimeAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultHazardRegimeAddin);

        FIELD(regime, "");
		FIELD(tradYears, "");

        Addin::registerClassObjectMethod("CALC_CREDIT_SPREADS",
                                          Addin::MARKET,
                                          "returns credit spread curve",
                                          TYPE,
                                          false,
                                          Addin::expandMulti,
                                          (Addin::ObjMethod*)calc);
    }

    static IObject* defaultHazardRegimeAddin(){
        return new HazardRegimeAddin();
    }
};

CClassConstSP const HazardRegimeAddin::TYPE = CClass::registerClassLoadMethod(
	"HazardRegimeAddin", typeid(HazardRegimeAddin), HazardRegimeAddin::load);

class HazardRegimeLeastSquareFit: public Calibrator::ObjFuncLeastSquare{
public:
	static CClassConstSP const TYPE;

	virtual void validatePop2Object(){
		static const string method = "HazardRegimeLeastSquareFit::validatePop2Object";
		try {
			if (tradYears.size()!=creditSpreadTgt.size()){
				throw ModelException(method, "tradYears and creditSpreadTgt must have the same size.");
			}
		} catch (exception& e){
			throw ModelException(e, method);
		}
	}

	void getMarket(MarketData* market){
	}

	IObjectSP getAdjustableGroup(){
		return IObjectSP::attachToRef(this);
	}

	int getNbFuncs() const{
		return tradYears.size();
	}

	void calcValue(CDoubleArray& funcvals) const {
		static const string method = "RegimeFactorHazardRegimeLeastSquareFit::calcValue";
		try{
			int iFunc = 0;
			for (int iTradYear = 0; iTradYear < tradYears.size(); ++iTradYear){
				funcvals[iFunc++] = (hazard->calcCreditSpread(tradYears[iTradYear]) - creditSpreadTgt[iTradYear])*10000.;
			}
		}
		catch(exception& e){
			throw ModelException(e, method);
		}
	}

	// for reflection
	HazardRegimeLeastSquareFit():
	Calibrator::ObjFuncLeastSquare(TYPE){}

	HazardRegimeLeastSquareFit(MarketDataSP		    market,
                                HazardRegimeWrapper hazard,
								DoubleArray		    tradYears,
								DoubleArray		    creditSpreadTgt):
	Calibrator::ObjFuncLeastSquare(TYPE),
	market(market),
    hazard(hazard),
	tradYears(tradYears),
	creditSpreadTgt(creditSpreadTgt){
		validatePop2Object();
	}

private:
	static void load(CClassSP& clazz){
		REGISTER(HazardRegimeLeastSquareFit, clazz);
		SUPERCLASS(Calibrator::ObjFuncLeastSquare);
		EMPTY_SHELL_METHOD(defaultHazardRegimeLeastSquareFit);

		FIELD(market, "");
        FIELD(hazard, "");
		FIELD(tradYears, "");
		FIELD(creditSpreadTgt, "");
	}

	static IObject* defaultHazardRegimeLeastSquareFit(){
		return new HazardRegimeLeastSquareFit();
	}

	//registered fields
	MarketDataSP            market;
    HazardRegimeWrapper		hazard;    
	DoubleArray				tradYears;
	DoubleArray				creditSpreadTgt;
};
typedef smartPtr<HazardRegimeLeastSquareFit> HazardRegimeLeastSquareFitSP;

class HazardRegimeDependencyFit: public Calibrator::ObjFuncLeastSquare{
public:
	static CClassConstSP const TYPE;

	virtual void validatePop2Object(){
		static const string method = "HazardRegimeDependencyFit::validatePop2Object";
		try {
			if (targets.size()!=hazard->getNbStates()-1){
				throw ModelException(method, "size of targets must be equal to nbStates - 1 of hazard.");
			}
		} catch (exception& e){
			throw ModelException(e, method);
		}
	}

	void getMarket(MarketData* market){
	}

	IObjectSP getAdjustableGroup(){
		return IObjectSP::attachToRef(this);
	}

	int getNbFuncs() const{
		return targets.size();
	}

	void calcValue(CDoubleArray& funcvals) const {
		static const string method = "HazardRegimeDependencyFit::calcValue";
		try{
			int iFunc = 0;
			for (int i = 0; i < targets.size(); ++i){
				funcvals[iFunc++] = (hazard->calcEquityCreditDependency(i,whichMethod) - targets[i])*10000.;
			}
		}
		catch(exception& e){
			throw ModelException(e, method);
		}
	}

	// for reflection
	HazardRegimeDependencyFit():
	Calibrator::ObjFuncLeastSquare(TYPE){}

	HazardRegimeDependencyFit(MarketDataSP		    market,
                                HazardRegimeWrapper hazard,
								DoubleArray		    targets,
                                int                 whichMethod):
	Calibrator::ObjFuncLeastSquare(TYPE),
	market(market),
    hazard(hazard),
	targets(targets),
    whichMethod(whichMethod){
		validatePop2Object();
	}

private:
	static void load(CClassSP& clazz){
		REGISTER(HazardRegimeDependencyFit, clazz);
		SUPERCLASS(Calibrator::ObjFuncLeastSquare);
		EMPTY_SHELL_METHOD(defaultHazardRegimeDependencyFit);

		FIELD(market, "");
        FIELD(hazard, "");
		FIELD(targets, "");
        FIELD(whichMethod, "");
	}

	static IObject* defaultHazardRegimeDependencyFit(){
		return new HazardRegimeDependencyFit();
	}

	//registered fields
	MarketDataSP            market;
    HazardRegimeWrapper		hazard;    
	DoubleArray				targets;
    int                     whichMethod;
};
typedef smartPtr<HazardRegimeDependencyFit> HazardRegimeDependencyFitSP;

class HazardRegimeCalibAddin: public CObject{
public:
    static CClassConstSP const TYPE;

    MarketDataSP            market;
    HazardRegimeWrapper     hazard;
    YieldCurveWrapper       discount;
    CreditCurveWrapper		creditCurve;        
	string					creditType;
    OptimizerNDSP           optimizer;
    DoubleArray             dependencyPowers;        // size equals to numbers of non-defaulting states
    double                  jumpAsymmetry;
    int                     whichMethod;

    // transcient fields
    ExpiryArray     		expiries;
	DoubleArray				tradYears;
	DoubleArray				creditSpreads;          //target, continuously compounding
    DoubleArray             creditSpreadResults;    //results, continuously compounding

    virtual void validatePop2Object(){
        static const string method = "HazardRegimeCalibAddin::validatePop2Object";
        try {
			if (!market.get()){
				throw ModelException(method, "market is not supplied.");
			}

			// target ...
            NonPricingModel model;
		    creditCurve.getData(&model, market);
            discount.getData(&model, market);

            MultiRegimeFactorWrapper    multi(hazard.getName());
            multi.getData(&model, market);
            hazard = HazardRegimeWrapper(multi.get()->getHazardRegime());

            DateTime	valueDate;
		    market->GetReferenceDate(valueDate);

            // create a risky curve
            IYieldCurveSP riskyCurve = creditCurve.get()->makeRiskyCurve(*discount.get());
            if (!riskyCurve.get()){
                throw ModelException("ConvBond::getRiskyCurve", "failed to create risky curve");
            }

            expiries = discount->getExpiries().obj();
            int iDate, numDates = expiries.size();
			DateTimeArray   dates(numDates);
            tradYears.resize(numDates);
            creditSpreads.resize(numDates);
            creditSpreadResults.resize(numDates);

			for (iDate = 0; iDate < numDates; iDate++){
                dates[iDate] = expiries[iDate]->toDate(valueDate);
                tradYears[iDate] = (dates[iDate].getDate()-valueDate.getDate())/365.;
                creditSpreads[iDate] = riskyCurve->pv(valueDate, dates[iDate])/discount->pv(valueDate, dates[iDate]);
                creditSpreads[iDate] = -log(creditSpreads[iDate])/tradYears[iDate];
			}           
            if (hazard->getNbStates()-1 != dependencyPowers.size()){
                throw ModelException(method, "size of dependencyPowers should be equal to numbers of non-default states");
            }
            hazard->setJumpAsymmetry(jumpAsymmetry);
        } catch (exception& e){
            throw ModelException(e, method);
        }
    }

    CResultsSP calibHazardRegime(){
        static const string method("HazardRegimeCalibAddin::calibHazardRegime");

        // calibrator
        Calibrator          calibrator(optimizer);

        // objFunc
        HazardRegimeLeastSquareFit      objFunc(market,
                                                hazard,
												tradYears,
											creditSpreads);

        // if there's a market, give objFunc a chance to get it
        // otherwise, trust (!) that objFunc already has all the market
        // data it needs
        if (market.get()){
            objFunc.getMarket(market.get());
        }

        // further validation irrespective of getMarket having been
        // called or not
        objFunc.validate();

        // fit the initial credit 
        int initStateIdx = hazard.get()->getInitStateIdx();
        int nbStates = hazard.get()->getNbStates();
        if ((initStateIdx == 0)||(initStateIdx == nbStates-2)){
            throw ModelException(method, "initstate can not be either the lowest or the highest non-default state");
        }

        DoubleArray states = hazard.get()->getStates();
        double initState = creditSpreads[0]/(1.-hazard.get()->getRecoveryRate());
        double ratio = initState / states[initStateIdx];
        for (int i=0; i<states.size()-1; i++){
            states[i] *= ratio;
        }
        hazard.get()->setStates(states);

        // ids
		DoubleArray		generatorEntryForCalib = hazard->getGeneratorEntryForCalib();
        string          name = hazard->getName();

        Calibrator::InstanceIDArray ids;
        ids.reserve(generatorEntryForCalib.size());
		for (int i=0; i < generatorEntryForCalib.size(); i++){
            double initGuess = 1.e-6;
            if (Maths::isPositive(generatorEntryForCalib[i])){
                initGuess = generatorEntryForCalib[i];
            }
   
            ids.push_back(Calibrator::InstanceIDSP(
				new Calibrator::InstanceIDDb("HazardRegime",
                                            name,
                                            "generatorEntryForCalib",
                                            true,
                                            initGuess,
                                            i,
                                            false)));
		}

        // run the calibrator
        CResultsSP  results = calibrator.run(objFunc, ids);

        // merge other results ...
        CResultsSP res(new Results());

        res->storeScalarGreek(initState, 
                            "Calibrator",
                            OutputNameSP(new OutputName("HazardRegime","initState")));

        res->storeGreek(IObjectSP(states.clone()),
                        "Calibrator",
                        OutputNameSP(new OutputName("HazardRegime","states")));

        for (int iTradYear = 0; iTradYear < tradYears.size(); iTradYear++){
            creditSpreadResults[iTradYear] = hazard->calcCreditSpread(tradYears[iTradYear]);
        }
        res->storeGreek(IObjectSP(tradYears.clone()),
                        "Calibrator",
                        OutputNameSP(new OutputName("creditSpreads", "tradYears")));
        res->storeGreek(IObjectSP(creditSpreads.clone()),
                        "Calibrator",
                        OutputNameSP(new OutputName("creditSpreads", "target")));
        res->storeGreek(IObjectSP(creditSpreadResults.clone()),
                        "Calibrator",
                        OutputNameSP(new OutputName("creditSpreads", "fitting")));

        // calib to equity credit dependency 
        CResultsSP  results2 = calibEquityCreditDependency();

        res->storeGreek(IObjectSP::constCast(results2->retrieveGreek("Calibrator",OutputNameSP(new OutputName("calibratedVars")))),
                        "Calibrator",
                        OutputNameSP(new OutputName("calibratedVars2")));

        res->storeGreek(IObjectSP::constCast(results2->retrieveGreek("Calibrator",OutputNameSP(new OutputName("dependencyPowers")))),
                        "Calibrator",
                        OutputNameSP(new OutputName("dependencyPowers")));

        results->merge("Calibrator",res.get());

        return results;
    }

    CResultsSP calibEquityCreditDependency(){
        static const string method("HazardRegimeCalibAddin::calibEquityCreditDependency");

        // calibrator
        Calibrator      calibrator(optimizer);

        // objFunc
        HazardRegimeDependencyFit      objFunc(market,
                                                hazard,
                                                dependencyPowers,
                                                whichMethod);

        // if there's a market, give objFunc a chance to get it
        // otherwise, trust (!) that objFunc already has all the market
        // data it needs
        if (market.get()){
            objFunc.getMarket(market.get());
        }

        // further validation irrespective of getMarket having been
        // called or not
        objFunc.validate();

        // ids
		DoubleArray		jumpEntryForCalib = hazard->getJumpEntryForCalib();
        string          name = hazard->getName();

        Calibrator::InstanceIDArray ids;
        ids.reserve(jumpEntryForCalib.size());
		for (int i=0; i < jumpEntryForCalib.size(); i++){
            ids.push_back(Calibrator::InstanceIDSP(
				new Calibrator::InstanceIDDb("HazardRegime",
                                            name,
                                            "jumpEntryForCalib",
                                            true,
                                            jumpEntryForCalib[i],
                                            i,
                                            false)));
		}

        // run the calibrator
        CResultsSP  results = calibrator.run(objFunc, ids);

        DoubleArray     powers(hazard->getNbStates()-1);
        for (int i = 0; i < powers.size(); i++){
            powers[i] = hazard->calcEquityCreditDependency(i, whichMethod);
        }
        results->storeGreek(IObjectSP(powers.clone()),
                        "Calibrator",
                        OutputNameSP(new OutputName("dependencyPowers")));
        return results;
    }

    static IObjectSP calib(HazardRegimeCalibAddin* params) {
		return params->calibHazardRegime();
    }

    static IObjectSP calibDependency(HazardRegimeCalibAddin* params) {
		return params->calibEquityCreditDependency();
    }

    HazardRegimeCalibAddin():
    CObject(TYPE) {}

private:

    static void load(CClassSP& clazz){
        REGISTER(HazardRegimeCalibAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultHazardRegimeCalibAddin);

        FIELD(market, "");
        FIELD(hazard, "");
        FIELD(discount, "");
        FIELD(creditCurve, "");
		FIELD(creditType, "");
        FIELD(optimizer, "");

        FIELD(dependencyPowers, "");
        FIELD_MAKE_OPTIONAL(dependencyPowers);
        FIELD(jumpAsymmetry, "");
        FIELD_MAKE_OPTIONAL(jumpAsymmetry);
        FIELD(whichMethod, "");
        FIELD_MAKE_OPTIONAL(whichMethod);

        // transients
        FIELD(expiries, "");
        FIELD_MAKE_TRANSIENT(expiries);
        FIELD(tradYears, "");
        FIELD_MAKE_TRANSIENT(tradYears);
        FIELD(creditSpreads, "");
        FIELD_MAKE_TRANSIENT(creditSpreads);
        FIELD(creditSpreadResults, "");
        FIELD_MAKE_TRANSIENT(creditSpreadResults);

        Addin::registerClassObjectMethod("CALIB_HAZARD_REGIME",
                                          Addin::MARKET,
                                          "Calibrate HazardRegime to credit spread curve",
                                          TYPE,
                                          true,
                                          Addin::returnHandle,
                                          (Addin::ObjMethod*)calib);

        Addin::registerClassObjectMethod("CALIB_EQUITY_CREDIT_DEPENDENCY",
                                          Addin::MARKET,
                                          "Calibrate equity credit dependency (Linetsky)",
                                          TYPE,
                                          true,
                                          Addin::returnHandle,
                                          (Addin::ObjMethod*)calibDependency);

    }

    static IObject* defaultHazardRegimeCalibAddin(){
        return new HazardRegimeCalibAddin();
    }
};

CClassConstSP const HazardRegimeLeastSquareFit::TYPE =
CClass::registerClassLoadMethod("HazardRegimeLeastSquareFit", typeid(HazardRegimeLeastSquareFit), load);

CClassConstSP const HazardRegimeDependencyFit::TYPE =
CClass::registerClassLoadMethod("HazardRegimeDependencyFit", typeid(HazardRegimeDependencyFit), load);

CClassConstSP const HazardRegimeCalibAddin::TYPE =
CClass::registerClassLoadMethod("HazardRegimeCalibAddin", typeid(HazardRegimeCalibAddin), load);

DRLIB_END_NAMESPACE
