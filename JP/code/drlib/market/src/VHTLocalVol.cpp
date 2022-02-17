//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : VHTLocalVol.cpp
//
//   Description : Parameterized Local Volatility
//
//   Author      : Manos Venardos
//
//   Date        : 2 November 2006
//
//   $Id: $
//----------------------------------------------------------------------------

#include "edginc/config.hpp"

#include "edginc/VolBaseParam.hpp"
#include "edginc/VolProcessedDVF.hpp"
#include "edginc/VolatilityDVF.hpp"
#include "edginc/Calibrator.hpp"
#include "edginc/VolRequestTime.hpp"
#include "edginc/LocVolRequest.hpp"
#include "edginc/ILocalVolGrid.hpp"
#include "edginc/MarketDataFetcherLN.hpp"
#include "edginc/NonPricingModel.hpp"
#include "edginc/Addin.hpp"

DRLIB_BEGIN_NAMESPACE

FORWARD_DECLARE(VHTLocalVolCalculator);

FORWARD_DECLARE_REF_COUNT(VHTLocalVolGrid);


/** Parameterized Vol */
class VHTLocalVol: public VolBaseParam,
                   virtual public IVolProcessed,
                   virtual public IVolatilityBS,
                   virtual public IVolatilityDVF,
                   public virtual Calibrator::IAdjustable {
private:

    /** Our VolParam class. This class exists primarily to keep the
        ComputeImpVol and spotVolSurfaceFromStrikes methods out of the
        CVolBaseParamSurface interface */
    class VHTLocalVolParam: public CVolParam{
    public:
        static CClassConstSP const TYPE;
        /** Default constructor */
        VHTLocalVolParam();

        /** Destructor */
        ~VHTLocalVolParam();

        /** Lattice based method for Implied Vol */
        virtual void ComputeImpVol(const CVolBase*          vol,
                                   const CLatticeDouble&    strikes,
                                   const DateTimeArray&     maturities,
                                   CLatticeDouble&          impV) const;

        /** Creates a "Matrix Vol Surface" on provided Strikes  */
        virtual VolSurface* spotVolSurfaceFromStrikes(const CVolBase* vol, const CDoubleArray&   strikes) const;

    private:
        static void load(CClassSP& clazz);
        
        static IObject* defaultCtor();
    };

    friend class VHTLocalVolParam;
    friend class VolProcessedVHTLocalVol;
    friend class LVGenerator;

public:

    static CClassConstSP const TYPE;

    /** Destructor */
    ~VHTLocalVol();

    /** Basic validation */
    void validatePop2Object();

    /** Method that builds a CVolParam. */
    virtual CVolParam* createVolParam() const;

    /** Needed for IAdjustable interface. Returns market data name for vol */
    virtual string getName() const;

    /** Entry point for volatility interpolation */
    IVolProcessed* getProcessedVol(const CVolRequest* volRequest,
                                   const CAsset*      asset) const;

    void getMarket(const IModel* model, const MarketData* market);

private:
    /** create a volatility calculator with a caching of forward values
        not implemented here */
    VHTLocalVolGrid* CreateVolCalculator(
        const DateTimeArray& maturity, bool isIntraDayInterp = true);

    /** Lattice method for vol interpolation */
    void computeImpVol(const CLatticeDouble&      strikes,
                       const DateTimeArray&       maturities,
                       CLatticeDouble&            impV) const;

    /** Computes implied vol for a given strike, maturity */
    double calcVol(double strike, double tradYear) const;

    /** Creates a "Matrix Vol Surface" on provided Strikes */
    VolSurface* spotVolSurfaceFromStrikes(const CDoubleArray& strikes) const;

    // MANDATORY FIELDS
    DateTimeArray                           benchmarkDates;         //!< Benchmark dates
    DoubleArray                             atmVols;                //!< ATMVols
    DoubleArray                             vhtParamA;              //!< Parameter A
    DoubleArray                             vhtParamB;              //!< Parameter B
    DoubleArray                             vhtParamC;              //!< Parameter C
    DoubleArray                             vhtParamD;              //!< Parameter D
    
    // TRANSIENT FIELDS
    DateTime                                baseDate;               //!< Base date
    VHTLocalVolCalculatorSP                 calc;                   //!< Calculator object


    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz);

    VHTLocalVol();

    static IObject* defaultCtor();

protected:
    /** Build the parameterised vol and cache any values */
    void buildCache();

    /** Called after adjustments have been made to fields (eg calibrator) */
    virtual void fieldsUpdated(const CFieldArray& fields);

    /** Called after (calibrator) adjustments have been made to fields */
    void update();
};

typedef smartPtr<VHTLocalVol> VHTLocalVolSP;


//////////////////////////////////////////////////////////////////////////


/** Interpolates vol at a given strike maturity */
class VHTLocalVolCalculator: public CObject {
public:
    static CClassConstSP const TYPE;

    /** Full constructor */
    VHTLocalVolCalculator(const DoubleArray& tradYears,
                          const DoubleArray& atmFwdVars,
                          const DoubleArray& vhtParamA,
                          const DoubleArray& vhtParamB,
                          const DoubleArray& vhtParamC,
                          const DoubleArray& vhtParamD):
    CObject(TYPE), tradYears(tradYears), atmFwdVars(atmFwdVars), 
    vhtParamA(vhtParamA), vhtParamB(vhtParamB), vhtParamC(vhtParamC), vhtParamD(vhtParamD) {
        static const string method("VHTLocalVolInterpolant::VHTLocalVolInterpolant");
        try {
            validatePop2Object();
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Validation */
    virtual void validatePop2Object() {
        static const string method("VHTLocalVolInterpolant::validatePop2Object");
        try {
            // Validate dimensions
            ASSERT(tradYears.size() == atmFwdVars.size());
            ASSERT(tradYears.size() == vhtParamA.size());
            ASSERT(tradYears.size() == vhtParamB.size());
            ASSERT(tradYears.size() == vhtParamC.size());
            ASSERT(tradYears.size() == vhtParamD.size());
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

    /** Volatility interpolation method */
    double calcLocalVolSquared(double x, double tradYear) const {
        static const string method("VHTLocalVolCalculator:calcLocalVolSquared");
        try{
            if (Maths::isZero(tradYear)){
                return atmFwdVars[0];
            }
            
            // Determine tradYearsIdx such that: tradYears[tradYearsIdx-1] <= tradYear < tradYears[tradYearsIdx]
			int tradYearsIdx = 0;
			while (tradYearsIdx < tradYears.size() && tradYears[tradYearsIdx] <= tradYear){
                tradYearsIdx++;
            }
            if(tradYearsIdx == tradYears.size()) {
                --tradYearsIdx;
            }
            
            // Get relevant parameters
            double atmFwdVar = atmFwdVars[tradYearsIdx];
            double A = vhtParamA[tradYearsIdx];
            double B = vhtParamB[tradYearsIdx];
            double C = vhtParamC[tradYearsIdx];
            double D = vhtParamD[tradYearsIdx];
            
            // Skew functional
            double skewFunc;
            if(Maths::isZero(A)) {
                skewFunc = 0.0;
            } else {
                skewFunc = A * tanh(C * x);
            }

            // Convexity functional
            double convxFunc;
            if(Maths::isZero(B)) {
                convxFunc = 0.0;
            } else {
                convxFunc = B * (1.0 - 1.0 / cosh(D * x));
            }
            
            // Local Var
            double smileFactor = 1.0 + skewFunc + convxFunc;
            if(Maths::isNegative(smileFactor)) {
                // throw ModelException("Negative Local var encountered");
                smileFactor = 0.0;
            }
            double localVar = atmFwdVar * smileFactor;
            return localVar;
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }

private:

    static IObject* defaultCtor(){
        return new VHTLocalVolCalculator();
    }

    /** Default constructor */
    VHTLocalVolCalculator(): CObject(TYPE) {}

    static void load(CClassSP& clazz){
        clazz->setPublic(); // make visible to EAS/spreadsheet
        REGISTER(VHTLocalVolCalculator, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(tradYears, "Trading time for backbone");
        FIELD(atmFwdVars, "Backbone forward var");
        FIELD(vhtParamA, "VHT Parameter A");
        FIELD(vhtParamB, "VHT Parameter B");
        FIELD(vhtParamC, "VHT Parameter C");
        FIELD(vhtParamD, "VHT Parameter D");
    }
    
    DoubleArray                             tradYears;              //!< Benchmark dates trading time
    DoubleArray                             atmFwdVars;             //!< Backbone forward var
    DoubleArray                             vhtParamA;              //!< Parameter A
    DoubleArray                             vhtParamB;              //!< Parameter B
    DoubleArray                             vhtParamC;              //!< Parameter C
    DoubleArray                             vhtParamD;              //!< Parameter D
};

typedef smartPtr<VHTLocalVolCalculator> VHTLocalVolCalculatorSP;

CClassConstSP const VHTLocalVolCalculator::TYPE =
CClass::registerClassLoadMethod("VHTLocalVolCalculator", typeid(VHTLocalVolCalculator), VHTLocalVolCalculator::load);


//////////////////////////////////////////////////////////////////////////

/** Interpolates vol at a given strike maturity */
class VHTLocalVolGrid: virtual public ILocalVolGrid {
public:
    VHTLocalVolGrid(const DateTimeArray& simDates,
                    const DateTime& valueDate,
                    TimeMetricConstSP timeMetric,
                    VHTLocalVolCalculatorSP calc):
    simDates(simDates), calc(calc) {
        int nbDates = simDates.size();
        simDatesTradYears.resize(nbDates);

        for(int i = 0 ; i < nbDates; i++) {
            simDatesTradYears[i] = timeMetric->yearFrac(valueDate, simDates[i]);
        }
    }
    
    /** calculate local volatility on a lattice of strikes */
    virtual void CalcLocVol(const CSliceDouble& strikes, int iDate, CSliceDouble& locVol) {
        for(int i = 0; i < strikes.size(); i++) {
            locVol[i] = interpLocVol(iDate, strikes[i]);
        }
    }
    
    /** Linear interpolation of local vol from grid */ 
    virtual double interpLocVar(int step, double logSpot) const {
        double tradYear = simDatesTradYears[step];
        double locVolSq = calc->calcLocalVolSquared(logSpot, tradYear);
        double sqrt_LocVarDt = sqrt(locVolSq * (simDatesTradYears[step + 1] - simDatesTradYears[step]));
        return sqrt_LocVarDt;
    }

    /** Linear interpolation of local vol from grid */ 
    double interpLocVol(int step, double logSpot) const {
        double tradYear = simDatesTradYears[step];
        double locVolSq = calc->calcLocalVolSquared(logSpot, tradYear);
        double locVol = sqrt(locVolSq);
        return locVol;
    }

    virtual double maxDrift() const {
        return 0.0;
    }
    
    /** Rolls LocalVolGrid forward */
    virtual void rollTime(const DateTime& valueDate) {

    }

private:
    DateTimeArray               simDates;
    DoubleArray                 simDatesTradYears;
    VHTLocalVolCalculatorSP     calc;
};

typedef refCountPtr<VHTLocalVolGrid> VHTLocalVolGridSP;

//////////////////////////////////////////////////////////////////////////


/** Type of CVolProcessedDVF. Returned following a LocalVolRequest */
class MARKET_DLL VolProcessedVHTLocalVol: public CVolProcessedDVF{
public:
    static CClassConstSP const TYPE;

    /** Full constructor */
    VolProcessedVHTLocalVol(
        const string& name, const DateTime& baseDate, TimeMetricConstSP timeMetric, VHTLocalVolCalculatorSP calc): 
    CVolProcessedDVF(TYPE), name(name), baseDate(baseDate), timeMetric(timeMetric), calc(calc) {};

    VolProcessedVHTLocalVol(): CVolProcessedDVF(TYPE) {}

    virtual double computeImpVol(const DateTime& maturity,
                                 double          strike) const {
        static const string method = "VolProcessedVHTLocalVol::computeImpVol";
        throw ModelException(method, "Cannot support Implied Vol functionality");
    }
    
    /** create a volatility calculator with a caching of forward values
        not implemented here */
    virtual CVolProcessedDVF::IVolCalculator* CreateVolCalculator(
        const DateTimeArray& maturity, bool isIntraDayInterp = true) const{
        // VHTLocalVolGrid* gridCalc = myVol->CreateVolCalculator(maturity);
        VHTLocalVolGrid* gridCalc = new VHTLocalVolGrid(maturity, baseDate, timeMetric, calc);
        return gridCalc;
    }
    
    virtual ~VolProcessedVHTLocalVol() {};

    /** Calculates LocalVol for a lattice of trikes */
    virtual void CalcLocVol(CLatticeDouble*     strikes,
                            DateTimeArray&      maturities, // the time axis of the lattice above
                            CLatticeDouble*	    locVol,
                            bool                isIntraDayInterp=true) const {
        static const string method = "VolProcessedVHTLocalVol::CalcLocVol";
        // Loop over dates
        for(int iDate = 0; iDate < maturities.size(); iDate++) {
            const CSliceDouble& tmpStrikes = (*strikes)[iDate];
            CSliceDouble& tmpLocVols = (*locVol)[iDate];
            double tradYear = timeMetric->yearFrac(baseDate, maturities[iDate]);
            // Loop over strikes
            for(int iStrike = 0; iStrike < tmpStrikes.size(); iStrike++) {
                double locVolSq = calc->calcLocalVolSquared(tmpStrikes[iStrike], tradYear);
                double locVol = sqrt(locVolSq);
                tmpLocVols[iStrike] = locVol;
            }
        }
    }
    
    /** Calculates LocalVol for a given strike / maturity */
    virtual double CalcLocVol(const DateTime& maturity,
                              double          strike,
                              bool            isIntraDayInterp=true) const {
        static const string method = "VolProcessedVHTLocalVol::CalcLocVol";
        
        double tradYear = timeMetric->yearFrac(baseDate, maturity);
        double locVolSq = calc->calcLocalVolSquared(strike, tradYear);
        double locVol = sqrt(locVolSq);
        return locVol;
    }
    
    /** Calculates variances */
    virtual void CalcLocVar(CLatticeDouble* 	strikes,
                            DateTimeArray& 	    maturities, // the time axis of the lattice above
                            CLatticeDouble*	    locVar,
                            bool                isIntraDayInterp=true) const {
        static const string method = "VolProcessedVHTLocalVol::CalcLocVar";
        throw ModelException(method, "Use VolCalculator instead");
    }
    
    virtual string getName() const {
        return name;
    }
    
    /** calculates the trading time between two dates */
    virtual double calcTradingTime(const DateTime &date1, 
                                   const DateTime &date2) const {
        return timeMetric->yearFrac(date1, date2);
    }
    
    /** retieve time measure for the vol */
    virtual TimeMetricConstSP GetTimeMetric()const {
        return timeMetric;
    }

private:
    static void load(CClassSP& clazz) {
        REGISTER(VolProcessedVHTLocalVol, clazz);
        SUPERCLASS(CVolProcessedDVF);
        EMPTY_SHELL_METHOD(defaultVolProcessedVHTLocalVol);
        
        FIELD(name, "Name");
        FIELD_MAKE_TRANSIENT(name);
        
        FIELD(baseDate, "Base date");
        FIELD_MAKE_TRANSIENT(baseDate);
        
        FIELD(timeMetric, "Time metric");
        FIELD_MAKE_TRANSIENT(timeMetric);

        FIELD(calc, "VHT Calculator");
        FIELD_MAKE_TRANSIENT(calc);
    }
    
    static IObject* defaultVolProcessedVHTLocalVol(){
        return new VolProcessedVHTLocalVol();
    }
    
    string                      name;
    DateTime                    baseDate;
    TimeMetricConstSP           timeMetric;
    VHTLocalVolCalculatorSP     calc;
};

typedef smartConstPtr<VolProcessedVHTLocalVol> VolProcessedVHTLocalVolConstSP;
typedef smartPtr<VolProcessedVHTLocalVol> VolProcessedVHTLocalVolSP;


//////////////////////////////////////////////////////////////////////////


/* Addin StrikeRefRatioMakerAddin */
class MARKET_DLL SimDatesGenerator: public CObject{
public:
    static CClassConstSP const TYPE;

    /** Interpolates Parameters implied vol */
    DateTimeArraySP createSimDates() {
        static const string method = "SimDatesGenerator::interpolateParameters";
        try {
            int dsec = DateTime::TIME_IN_DAY / (double)(nbSamplesPerDay + 1);
            if(dsec < secTolerance) {
                throw ModelException(method, "Too many samples per day relative to tolerance");
            }
            
            // First create daily tineline
            DateTimeArray dailyDates;
            DateTime thisDate = endDate;
            while (baseDate < thisDate){
                dailyDates.push_back(thisDate);
                thisDate = hols->addBusinessDays(thisDate, -1);
            }
            dailyDates.push_back(baseDate);
            sort(dailyDates.begin(), dailyDates.end());

            // Then insert intraday samples
            DateTimeArray intradayDates;
            
            for(int i = 0; i < dailyDates.size(); i++) {
                if(DateTime::START_OF_DAY_TIME + dsec < DateTime::END_OF_DAY_TIME &&
                   DateTime::START_OF_DAY_TIME < DateTime::END_OF_DAY_TIME - dsec) {
                
                    DateTime firstDate(dailyDates[i].getDate(), DateTime::START_OF_DAY_TIME + dsec);
                    DateTime lastDate(dailyDates[i].getDate(), DateTime::END_OF_DAY_TIME - secTolerance);
                    thisDate = firstDate;
                    while (thisDate <= lastDate){
                        intradayDates.push_back(thisDate);
                        int time = thisDate.getTime() + dsec;
                        if(time >= DateTime::TIME_IN_DAY) {
                            break;
                        }
                        thisDate = DateTime(thisDate.getDate(), time);
                    }
                }
            }
            sort(intradayDates.begin(), intradayDates.end());
            DateTimeArray merged = DateTime::merge(dailyDates, intradayDates);
            DateTimeArraySP output(new DateTimeArray(DateTime::doSortUniq(merged)));

            return output;
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }
    
private:
    
    /** For reflection */
    SimDatesGenerator(): CObject(TYPE) {}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(SimDatesGenerator, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(hols, "Holidays");
        FIELD(baseDate, "Base date");
        FIELD(endDate, "End date");
        FIELD(nbSamplesPerDay, "Nubers of points per day");
        FIELD(secTolerance, "Tolerance in seconds");

        Addin::registerObjectMethod("CREATE_SIM_DATES",
                                    Addin::MARKET,
                                    "Create simulation dates",
                                    false,
                                    Addin::returnHandle,
                                    &SimDatesGenerator::createSimDates);
    }
    
    static IObject* defaultCtor(){
        return new SimDatesGenerator();
    }
    
    // MANDATORY FIELDS
    HolidaySP           hols;                   //!< Holidays
    DateTime            baseDate;               //!< Base date
    DateTime            endDate;                //!< End date
    int                 nbSamplesPerDay;        //!< Number of simulation points per day
    int                 secTolerance;           //!< Tolerance in seconds
};


CClassConstSP const SimDatesGenerator::TYPE = CClass::registerClassLoadMethod(
    "SimDatesGenerator", typeid(SimDatesGenerator), load);


//////////////////////////////////////////////////////////////////////////


/* Addin StrikeRefRatioMakerAddin */
class MARKET_DLL LVGenerator: public CObject{
public:
    static CClassConstSP const TYPE;

    /** Interpolates Parameters implied vol */
    CDoubleMatrixSP computeLV() {
        static const string method = "LVGenerator::interpolateParameters";
        try {
            MarketDataFetcherSP mdf(new MarketDataFetcherLN("VHTLocalVol"));
            NonPricingModel model(mdf);
            vhtLV->getMarket(&model, market.get());
            asset.getData(&model, market.get());
            
            DateTime baseDate = market->GetReferenceDate();
            CDoubleMatrixSP output(new CDoubleMatrix(fwdMoneynessStrikes->size(), maturities->size()));
            
            CVolRequestSP request(new LocVolRequest(
                baseDate,
                false,
                false,
                true,
                true,
                1.0,
                1.0,
                1.0));
            
            IVolProcessedSP volProc(vhtLV->getProcessedVol(request.get(), asset.getSP().get()));
            CVolProcessedDVFSP dvfProc = CVolProcessedDVFSP::dynamicCast(volProc);
            CVolProcessedDVF::IVolCalculator* calc = dvfProc->CreateVolCalculator(*maturities);
            VHTLocalVolGridSP vhtCalc(dynamic_cast<VHTLocalVolGrid*>(calc));
            if(!vhtCalc) {
                throw ModelException(method, "VolCalculator must be of type VHTLocalVolGrid");
            }
            
            for(int i = 0; i < fwdMoneynessStrikes->size(); i++) {
                for(int j = 0; j < maturities->size(); j++) {
                    double yearFrac = vhtLV->calcTradingTime(baseDate, (*maturities)[j]);
                    if(vol) {
                        (*output)[i][j] = vhtCalc->interpLocVol(j, (*fwdMoneynessStrikes)[i]);
                    } else {
                        (*output)[i][j] = vhtCalc->interpLocVar(j, (*fwdMoneynessStrikes)[i]);
                    }
                }
            }
        
            return output;
        } catch(exception& e) {
            throw ModelException(e, method);
        }
    }
    
private:
    
    /** For reflection */
    LVGenerator(): CObject(TYPE) {}
    
    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz){
        REGISTER(LVGenerator, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultCtor);
        FIELD(vhtLV, "VHT Local Vol");
        FIELD(market, "Market");
        FIELD(asset, "Asset");
        FIELD(fwdMoneynessStrikes, "Strikes");
        FIELD(maturities, "Maturities");
        FIELD(vol, "True: vol, False: var");

        Addin::registerObjectMethod("VHTLV_LOCAL_VOL",
                                    Addin::MARKET,
                                    "Compute Local Vol",
                                    false,
                                    Addin::returnHandle,
                                    &LVGenerator::computeLV);
    }
    
    static IObject* defaultCtor(){
        return new LVGenerator();
    }
    
    // MANDATORY FIELDS
    VHTLocalVolSP       vhtLV;
    MarketDataSP        market;
    CAssetWrapper       asset;
    DoubleArraySP       fwdMoneynessStrikes;
    DateTimeArraySP     maturities;
    bool                vol;
};


CClassConstSP const LVGenerator::TYPE = CClass::registerClassLoadMethod(
    "LVGenerator", typeid(LVGenerator), load);


//////////////////////////////////////////////////////////////////////////


VHTLocalVol::VHTLocalVolParam::VHTLocalVolParam(): CVolParam(TYPE){}

VHTLocalVol::VHTLocalVolParam::~VHTLocalVolParam() {}

VHTLocalVolGrid* VHTLocalVol::CreateVolCalculator(
    const DateTimeArray& maturity, bool isIntraDayInterp) {
    
    VHTLocalVolGrid* lvCalc = new VHTLocalVolGrid(maturity,
                                                  baseDate,
                                                  timeMetric,
                                                  calc);
    return lvCalc;
}

void VHTLocalVol::VHTLocalVolParam::ComputeImpVol(const CVolBase*          vol,
                                                  const CLatticeDouble&    strikes,
                                                  const DateTimeArray&     maturities,
                                                  CLatticeDouble&          impV) const {
    // turn the vol into what we must have
    const VHTLocalVol* myVol = static_cast<const VHTLocalVol*>(vol);
    // then just pass through the parameterised vol
    myVol->computeImpVol(strikes, maturities, impV);
}


VolSurface* VHTLocalVol::VHTLocalVolParam::spotVolSurfaceFromStrikes(
    const CVolBase* vol, const CDoubleArray&   strikes) const{
    // turn the vol into what we must have
    const VHTLocalVol* myVol = static_cast<const VHTLocalVol*>(vol);
    // then just pass through the parameterised vol
    return myVol->spotVolSurfaceFromStrikes(strikes);
}


void VHTLocalVol::VHTLocalVolParam::load(CClassSP& clazz){
    REGISTER(VHTLocalVolParam, clazz);
    SUPERCLASS(CVolParam);
    EMPTY_SHELL_METHOD(defaultCtor);
}
        
IObject* VHTLocalVol::VHTLocalVolParam::defaultCtor(){
    return new VHTLocalVolParam();
}

VHTLocalVol::~VHTLocalVol() {}

void VHTLocalVol::validatePop2Object(){
    static const string method("VHTLocalVol::validatePop2Object");
    try{
        // Validate parameters are in range
        Calibrator::IAdjustable::checkRange(this);
    } catch(exception& e){
        throw ModelException(e, method);
    }
}


CVolParam* VHTLocalVol::createVolParam() const{
    return new VHTLocalVolParam();
}


string VHTLocalVol::getName() const{
    return CVolBaseParamSurface::getName();
}

IVolProcessed* VHTLocalVol::getProcessedVol(const CVolRequest* volRequest,
                                            const CAsset*      asset) const {
    static const string method = "VHTLocalVol::getProcessedVol";
    
    try{
        const LocVolRequest* lvRequest = dynamic_cast<const LocVolRequest*>(volRequest);
        if(lvRequest) {
            return new VolProcessedVHTLocalVol(getName(), baseDate, timeMetric, calc);
        }

        const VolRequestTime* tRequest = dynamic_cast<const VolRequestTime*>(volRequest);
        if(tRequest) {
            return VolRequestTime::createVolProcessed(getName());
        }
            
        throw ModelException(method, "Only LocVolRequest & VolRequestTime supported");

       // return VolBaseParam::getProcessedVol(volRequest, asset);
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


void VHTLocalVol::getMarket(const IModel* model, const MarketData* market){
    static const string method = "VHTLocalVol::getMarket";

    try {
        VolBaseParam::getMarket(model, market);
        baseDate = market->GetReferenceDate();

        int nbDates = atmVols.size();
        DoubleArray tradYears(nbDates);
        DoubleArray atmFwdVars(nbDates);

        // Convert from ATMVols to fwdVars
        for(int i = 0; i < nbDates; i++) {
            tradYears[i] = timeMetric->yearFrac(baseDate, benchmarkDates[i]);
            if(i) {
                atmFwdVars[i] = 
                    (Maths::square(atmVols[i]) * tradYears[i] - Maths::square(atmVols[i-1]) * tradYears[i-1]) / 
                    (tradYears[i] - tradYears[i-1]);
            } else {
                atmFwdVars[i] = Maths::square(atmVols[i]);
            }
        }

        // Create VHTCalculator
        calc = VHTLocalVolCalculatorSP(new VHTLocalVolCalculator(
            tradYears,
            atmFwdVars,
            vhtParamA,
            vhtParamB,
            vhtParamC,
            vhtParamD));
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}


void VHTLocalVol::computeImpVol(const CLatticeDouble&      strikes,
                                const DateTimeArray&       maturities,
                                CLatticeDouble&            impV) const {
    static const string method("VHTLocalVol::computeImpVol");
    throw ModelException(method, "Implied Vol computations not supported");
}


double VHTLocalVol::calcVol(double strike, double tradYear) const{
    static const string method("VHTLocalVol::calcVol");
    throw ModelException(method, "Implied Vol computations not supported");
}


VolSurface* VHTLocalVol::spotVolSurfaceFromStrikes(const CDoubleArray& strikes) const {
    static const string method("VHTLocalVol::spotVolSurfaceFromStrikes");
    throw ModelException(method, "Implied Vol computations not supported");
}


void VHTLocalVol::load(CClassSP& clazz){
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(VHTLocalVol, clazz);
    SUPERCLASS(VolBaseParam);
    IMPLEMENTS(IVolatilityBS);
    IMPLEMENTS(IVolatilityDVF);
    IMPLEMENTS(Calibrator::IAdjustable);
    EMPTY_SHELL_METHOD(defaultCtor);

    // MANDATORY
    FIELD(benchmarkDates, "Benchmark dates");
    FIELD(atmVols, "ATM vols");
    FIELD(vhtParamA, "VHT Parameter A");
    FIELD(vhtParamB, "VHT Parameter B");
    FIELD(vhtParamC, "VHT Parameter C");
    FIELD(vhtParamD, "VHT Parameter D");

    // TRANSIENT
    FIELD(baseDate, "Base date");
    FIELD_MAKE_TRANSIENT(baseDate);
    FIELD(calc, "LocalVol calculator");
    FIELD_MAKE_TRANSIENT(calc);

    // add our fields and their ranges to central list
    Calibrator::IAdjustable::registerField(
        clazz, "atmVols",
        new InfiniteRange());
    Calibrator::IAdjustable::registerField(
        clazz, "vhtParamA",
        new InfiniteRange());
    Calibrator::IAdjustable::registerField(
        clazz, "vhtParamB",
        new Range(ClosedBoundary(0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "vhtParamC",
        new Range(OpenBoundary(0), Infinity(Infinity::Plus)));
    Calibrator::IAdjustable::registerField(
        clazz, "vhtParamD",
        new Range(OpenBoundary(0), Infinity(Infinity::Plus)));
}


VHTLocalVol::VHTLocalVol(): VolBaseParam(TYPE) {}


IObject* VHTLocalVol::defaultCtor(){
    return new VHTLocalVol();
}


void VHTLocalVol::buildCache() {
    update();
}


void VHTLocalVol::fieldsUpdated(const CFieldArray& fields){
    update();
}


void VHTLocalVol::update(){
    static const string method = "VHTLocalVol::update";
    try {
        
    } catch(exception& e) {
        throw ModelException(e, method);
    }
}

CClassConstSP const VHTLocalVol::TYPE =
CClass::registerClassLoadMethod("VHTLocalVol", typeid(VHTLocalVol), load);

CClassConstSP const VolProcessedVHTLocalVol::TYPE =
CClass::registerClassLoadMethod("VolProcessedVHTLocalVol", typeid(VolProcessedVHTLocalVol), load);

CClassConstSP const VHTLocalVol::VHTLocalVolParam::TYPE =
CClass::registerClassLoadMethod("VHTLocalVol::VHTLocalVolParam", typeid(VHTLocalVolParam), load);


//////////////////////////////////////////////////////////////////////////////////////////////


// external symbol to allow class to be forced to be linked in
bool VHTLocalVolLinkIn(){
    return (VHTLocalVol::TYPE != 0);
}

DRLIB_END_NAMESPACE
