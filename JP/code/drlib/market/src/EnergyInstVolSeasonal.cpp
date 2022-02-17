trument for the specified maturity
InstrumentSP BootstrappedYieldCurve::getParInstrument(const ExpirySP maturity) const
{
    static const string method = "BootstrappedYieldCurve::getParInstrument";

    try
    {
        //find maturity in BootstrappedYieldCurve expiries
        //-throws exception if not found
        int i = maturity->search(allExpiries.get());

        //store maturity date of this par instrument;
        DateTime matDate = benchmarks[i]->getBenchmarkDate(valueDate);

        CashFlowArray cfl(0);
        if (benchmarks[i]->isCash())
        {
            if (!floatDcc)
            {
                throw ModelException(method, "swap floating leg day count convention not supplied");
            }

            //initial notional exchange at valueDate
            CashFlow cfs(valueDate, -1);
            cfl.push_back(cfs);

            //final notional exchange 1 + coupon
            CashFlow cfm(matDate, 1 + (benchmarks[i]->getRate() * floatDcc->years(valueDate, matDate)));
            cfl.push_back(cfm);
        }
        else if (benchmarks[i]->isSwap())
        {
            //multiple cashflows occurring at swapFrequency intervals
            CashFlow  cf;
            DateTime  cfDate;
            DateTime  prevCfDate;    //store previous date for accrual period calculations
            double    cfAmount;

            if (!fixedDcc)
            {
                throw ModelException(method, "swap fixed leg day count convention not supplied");
            }

            //initial notional exchange at valueDate
            cf.date = valueDate;
            cf.amount = -1;
            cfl.push_back(cf);
            prevCfDate = valueDate;

            //regular cashflows
            bool done = false;
            int count = 1;
            while (!done)
            {
                //cashflow date is regular from valueDate
                cfDate = MaturityPeriod::toDate(count * 12 / fixedIvl->approxAnnualFrequency(), "M", valueDate);

                //up to but not including the final cashflow date
                if (cfDate < matDate)
                {
                    //rate * accrual period year fraction * notional
                    cfAmount = benchmarks[i]->getRate() * fixedDcc->years(prevCfDate, cfDate); //notional of 1
                    cf.date = cfDate;
                    cf.amount = cfAmount;

                    //add to the list & continue loop
                    cfl.push_back(cf);
                    prevCfDate = cfDate;
                    count++;
                }
                else
                {
                    done = true;
                }
            }

            //final notional exchange 1 + coupon
            cfAmount = 1 + (benchmarks[i]->getRate() * fixedDcc->years(prevCfDate, matDate)); //notional of 1
            cf.date = matDate;
            cf.amount = cfAmount;
            cfl.push_back(cf);
        }

        SimpleCashFlowStreamSP cfs(new SimpleCashFlowStream(&cfl, name));
        return cfs;
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}

bool
BootstrappedYieldCurve::isIndexCurve() const
{
    return isIndexCurve_;
}

/*
 * Reflection and addin support.
 */

/** Invoked when Class is 'loaded' */
/* static */
void BootstrappedYieldCurve::load(CClassSP& clazz)
{
    REGISTER(BootstrappedYieldCurve, clazz);
    SUPERCLASS(YieldCurve);
    IMPLEMENTS(IDeterministicYieldCurve);
    IMPLEMENTS(IPrivateObject);
    IMPLEMENTS(IGetMarket);
    IMPLEMENTS(IRestorableWithRespectTo<RateParallel>);
    IMPLEMENTS(IRestorableWithRespectTo<RatePointwise>);
    IMPLEMENTS(IRestorableWithRespectTo<IRRatePointwise>);
    IMPLEMENTS(Theta::Shift);
    IMPLEMENTS(IRiskyCurve);
    IMPLEMENTS(Duration::IParHandlerWithClosedForm);
    IMPLEMENTS(Duration::IParHandlerWithoutClosedForm);
    IMPLEMENTS(YCWeightedShift::IShift);
    EMPTY_SHELL_METHOD(defaultBootstrappedYieldCurve);
    FIELD(ccy,                  "currency name");
    FIELD(name,                 "rate index");
    FIELD(today,                "today");
    FIELD(spotOffset,           "spot offset");
    FIELD(hols,                 "holidays");
    FIELD       (moneyMarketDayCount,  "money market day count");
    FIELD(benchmarks,           "benchmarks");
    FIELD       (zcMethod,             "zero curve methodology");
    FIELD(ccyBasis,             "currency basis");
    FIELD(irVol,                "Interest rate volatility");
    FIELD       (futMaturity,          "length of future");
    FIELD       (badDayConvention,     "swap date adjustment convention");
    FIELD       (fixedDcc,             "swap fixed leg day count convention");
    FIELD       (fixedIvl,             "swap fixed leg period");
    FIELD       (floatDcc,             "swap floating leg day count convention");
    FIELD       (floatIvl,             "swap floating leg period");
    FIELD       (fixDates,             "fixing dates");
    FIELD       (fixRates,             "fixing rates");
    FIELD       (basisDcc,             "swap basis leg day count convention");
    FIELD       (basisIvl,             "swap basis leg period");
    FIELD       (isIndexCurve_,         "act as index curve");

    FIELD_MAKE_OPTIONAL(today);
    FIELD_MAKE_OPTIONAL(ccyBasis);
    FIELD_MAKE_OPTIONAL(irVol);
    FIELD_MAKE_OPTIONAL(futMaturity);
    FIELD_MAKE_OPTIONAL(badDayConvention);
    FIELD_MAKE_OPTIONAL(floatDcc);
    FIELD_MAKE_OPTIONAL(floatIvl);
    FIELD_MAKE_OPTIONAL(basisDcc);
    FIELD_MAKE_OPTIONAL(basisIvl);
    FIELD_MAKE_OPTIONAL(isIndexCurve_);

    // transient fields
    FIELD_NO_DESC       (allExpiries);
    FIELD_NO_DESC(useProjectionCurve);
    FIELD_NO_DESC(nosort);

    FIELD_MAKE_TRANSIENT(allExpiries);
    FIELD_MAKE_TRANSIENT(useProjectionCurve);
    FIELD_MAKE_TRANSIENT(nosort);

    ClassSetAcceptMethod(BootstrappedYieldCurve::acceptValueDateCollector);
    ClassSetAcceptMethod(BootstrappedYieldCurve::acceptWrapperNameCollector);
    ClassSetAcceptMethod(BootstrappedYieldCurve::acceptYieldNameCollector);
    clazz->setPrivate(); // hide this class
    // record how to build it (since implements IPrivateObject)
    clazz->addConstructorClass(BootstrappedCurveAddin::TYPE);
}

/* static */
IObject* BootstrappedYieldCurve::defaultBootstrappedYieldCurve()
{
    return new BootstrappedYieldCurve();
}

CClassConstSP const BootstrappedYieldCurve::TYPE = CClass::registerClassLoadMethod(
    "BootstrappedYieldCurve", typeid(BootstrappedYieldCurve), BootstrappedYieldCurve::load);


// Addin for bootstrapping yield curves

/* static */
IObject* BootstrappedCurveAddin::defaultBootstrappedCurveAddin()
{
    return new BootstrappedCurveAddin();
}


/** Invoked when Class is 'loaded' */
/* static */
void BootstrappedCurveAddin::load(CClassSP& clazz)
{
    clazz->setDRIProxyType(BootstrappedYieldCurve::TYPE); //use BootstrappedYieldCurve for dri
    REGISTER(BootstrappedCurveAddin, clazz);
    SUPERCLASS(CObject);
    IMPLEMENTS(IPublicObject);
    EMPTY_SHELL_METHOD(defaultBootstrappedCurveAddin);
    FIELD(ccy,              "currency name");
    FIELD(name,             "yield curve name");
    FIELD(today,            "today");
    FIELD(spotOffset,       "spot offset");
    FIELD(hols,             "holidays");
    FIELD(benchmarks,       "benchmarks");
    FIELD       (factory,          "bootstrapping methodology");
    FIELD(moneyMarketDcc,   "money market day count");
    FIELD       (futuresMaturity,  "maturity of futures rates, eg. 3M or 1I");
    FIELD(irVol,            "Interest rate volatility");
    FIELD(badDayConvention, "swap bad day convention");
    FIELD(fixedDcc,         "swap fixed leg day count convention");
    FIELD       (fixedIvl,         "swap fixed leg period");
    FIELD(floatDcc,         "swap floating leg day count convention");
    FIELD       (floatIvl,         "swap floating leg period");
    FIELD       (fixDates,         "fixing dates");
    FIELD       (fixRates,         "fixing rates");
    FIELD(basisDcc,         "swap basis leg day count convention");
    FIELD       (basisIvl,         "swap basis leg period");
    FIELD(ccyBasis,         "currency basis");
    FIELD(isIndexCurve,     "whether to act as index curve");

    FIELD_MAKE_OPTIONAL(today);
    FIELD_MAKE_OPTIONAL(futuresMaturity);
    FIELD_MAKE_OPTIONAL(irVol);
    FIELD_MAKE_OPTIONAL(badDayConvention);
    FIELD_MAKE_OPTIONAL(floatDcc);
    FIELD_MAKE_OPTIONAL(floatIvl);
    FIELD_MAKE_OPTIONAL(fixDates);
    FIELD_MAKE_OPTIONAL(fixRates);
    FIELD_MAKE_OPTIONAL(basisDcc);
    FIELD_MAKE_OPTIONAL(basisIvl);
    FIELD_MAKE_OPTIONAL(ccyBasis);
    FIELD_MAKE_OPTIONAL(isIndexCurve);

    Addin::registerConstructor("BOOTSTRAP_YIELD_CURVE",
                               Addin::MARKET,
                               "Bootstrap a zero curve",
                               BootstrappedCurveAddin::TYPE);
}


BootstrappedCurveAddin::BootstrappedCurveAddin() : CObject(TYPE), isIndexCurve(false)
{
}


BootstrappedCurveAddin::BootstrappedCurveAddin(
    const string&                  ccy,
    const string&                  name,
    const DateTime&                today,
    int                            spotOffset,
    const HolidayWrapper&          hols,
    const ZeroCurveBenchmarkArray& benchmarks,
    const IZeroCurveFactory&       factory,
    const string&                  moneyMarketDcc,   // may be ""
    const MaturityPeriod*          futuresMaturity,
    const IRVolBaseWrapper&        irVol,
    const string&                  badDayConvention, // may be ""
    const string&                  fixedDcc,         // may be ""
    const MaturityPeriod*          fixedIvl,
    const string&                  floatDcc,         // may be ""
    const MaturityPeriod*          floatIvl,
    const ExpiryArray*             fixDates,
    const DoubleArray*             fixRates,
    const string&                  basisDcc,         // may be ""
    const MaturityPeriod*          basisIvl,
    const CurrencyBasisWrapper&    ccyBasis,
    bool                           isIndexCurve)
    : CObject(TYPE),
    ccy(ccy), name(name), today(today), spotOffset(spotOffset), hols(hols),
    benchmarks(benchmarks),
    factory(const_cast<IZeroCurveFactory*>(&factory)),
    moneyMarketDcc(moneyMarketDcc),
    futuresMaturity(futuresMaturity ? dynamic_cast<MaturityPeriod*>(futuresMaturity->clone()) : NULL),
    irVol(irVol),
    badDayConvention(badDayConvention),
    fixedDcc(fixedDcc),
    fixedIvl(fixedIvl ? dynamic_cast<MaturityPeriod*>(fixedIvl->clone()) : NULL),
    floatDcc(floatDcc),
    floatIvl(floatIvl ? dynamic_cast<MaturityPeriod*>(floatIvl->clone()) : NULL),
    fixDates(fixDates ? dynamic_cast<ExpiryArray*>(fixDates->clone()) : NULL),
    fixRates(fixRates ? dynamic_cast<DoubleArray*>(fixRates->clone()) : NULL),
    basisDcc(basisDcc),
    basisIvl(basisIvl ? dynamic_cast<MaturityPeriod*>(basisIvl->clone()) : NULL),
    ccyBasis(ccyBasis),
    isIndexCurve(isIndexCurve)
{
}


IPrivateObject* BootstrappedCurveAddin::toPrivateObject() const
{
    static const string method = "BootstrappedCurveAddin::toPrivateObject";

    try
    {
        if (!factory)
        {
            throw ModelException(method, "bootstrapping methodology not defined");
        }

        if (moneyMarketDcc.empty())
        {
            throw ModelException(method, "Money market day count convention must be defined");
        }

        if (fixedDcc.empty())
        {
            throw ModelException(method, "Swap fixed day count convention must be defined");
        }

        if (!fixedIvl.get())
        {
            throw ModelException(method, "Swap fixed interval must be defined");
        }

        return new BootstrappedYieldCurve(ccy,
                                 name,
                                 today,
                                 spotOffset,
                                 hols,
                                 benchmarks,
                                 *factory,
                                 *getDcc(moneyMarketDcc),
                                 futuresMaturity.get(),
                                 irVol,
                                 getBdc(badDayConvention).get(),
                                 *getDcc(fixedDcc),
                                 *fixedIvl,
                                 getDcc(floatDcc).get(),
                                 floatIvl.get(),
                                 fixDates.get(),
                                 fixRates.get(),
                                 getDcc(basisDcc).get(),
                                 basisIvl.get(),
                                 ccyBasis,
                                 isIndexCurve);
    }
    catch (exception &e)
    {
        throw ModelException(e, method, e.what());
    }
}


BadDayConventionSP  BootstrappedCurveAddin::getBdc(const string& bdc) const
{
    return BadDayConventionSP(bdc.empty() ? NULL : BadDayConventionFactory::make(bdc));
}


DayCountConventionSP  BootstrappedCurveAddin::getDcc(const string& dcc) const
{
    return DayCountConventionSP(dcc.empty() ? NULL : DayCountConventionFactory::make(dcc));
}


CClassConstSP const BootstrappedCurveAddin::TYPE =
CClass::registerClassLoadMethod("BootstrappedCurveAddin", typeid(BootstrappedCurveAddin),
                                BootstrappedCurveAddin::load);


// other addins

class ForwardCurveAddin: public CObject {
    static CClassConstSP const TYPE;

    // input parameters
    DateTime            forwardDate;
    BootstrappedYieldCurveSP yieldCurve;

    static IObjectSP getForwardSwapCurve(ForwardCurveAddin* params) {
        static const string routine = "ForwardCurveAddin::getForwardSwapCurve";
        try {
            IYieldCurveSP fwdCurve = IYieldCurveSP(params->yieldCurve->createForwardCurve(params->forwardDate));
            return fwdCurve;

        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    ForwardCurveAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(ForwardCurveAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultForwardCurveAddin);
        FIELD(forwardDate, "forward spot date of the yield curve");
        FIELD(yieldCurve,         "spot yield curve");
        Addin::registerClassObjectMethod("GET_FORWARD_YIELD_CURVE",
                                         Addin::MARKET,
                                         "Returns the forward yield curve for a given date",
                                         TYPE,
                                         true,
                                         Addin::returnHandle,
                                         (Addin::ObjMethod*)getForwardSwapCurve);

    }

    static IObject* defaultForwardCurveAddin(){
        return new ForwardCurveAddin();
    }

};

CClassConstSP const ForwardCurveAddin::TYPE = CClass::registerClassLoadMethod(
    "ForwardCurveAddin", typeid(ForwardCurveAddin), load);

class ZeroCurveAddin: public CObject {
    static CClassConstSP const TYPE;

    // input parameters
    BootstrappedYieldCurveSP     yieldCurve;
    bool                extractGrowth;

    static IObjectSP getZeroCurve(ZeroCurveAddin* params) {
        static const string routine = "ZeroCurveAddin::getZeroCurve";
        try {
            ZeroCurveSP     zc        = params->yieldCurve->get(params->extractGrowth);
            CashFlowArraySP cashFlows = zc->getRatesAndDates();


            DateTime::DateArraySP dateArray(new DateTime::DateArray(0));
            DoubleArraySP         doubleArray(new DoubleArray(0));

            for (int i = 0; i < cashFlows->size(); i++) {
                dateArray->push_back(DateTime::DateSP(new DateTime::Date((*cashFlows)[i].date.getDate())));
                doubleArray->push_back((*cashFlows)[i].amount);
            }

            ObjectArraySP objArray(new ObjectArray(0));

            objArray->push_back(dateArray);
            objArray->push_back(doubleArray);

            return objArray;
        }
        catch (exception& e){
            throw ModelException(e, routine);
        }
    }

    /** for reflection */
    ZeroCurveAddin():  CObject(TYPE), extractGrowth(true) {}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(ZeroCurveAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultZeroCurveAddin);
        FIELD(yieldCurve,         "spot yield curve");
        FIELD(extractGrowth, "return growth curve or discount");
        FIELD_MAKE_OPTIONAL(extractGrowth);
        Addin::registerClassObjectMethod("ZERO_CURVE_EXTRACT",
                                         Addin::MARKET,
                                         "Returns the zero curve for a yield curve",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)getZeroCurve);

    }

    static IObject* defaultZeroCurveAddin(){
        return new ZeroCurveAddin();
    }

};

CClassConstSP const ZeroCurveAddin::TYPE = CClass::registerClassLoadMethod(
    "ZeroCurveAddin", typeid(ZeroCurveAddin), load);


DRLIB_END_NAMESPACE
ze() ; i++)
        {
            const Expiry* expiry = benchmarks[i]->isFuture() ? benchmarks[i]->getStart() : benchmarks[i]->getEnd();
            (*allExpiries)[i] = ExpirySP(const_cast<Expiry*>(expiry));
        }
    }
    catch (exception& e)
    {
        throw ModelException(e, method);
    }
}

void BootstrappedYieldCurve::acceptValueDateCollector(
    BootstrappedYieldCurve*       yieldCurve,
    CValueDateCollector* collector)
{
    collector->valueDateValidate(yieldCurve->today,
                                 yieldCurve->getName());
}

void BootstrappedYieldCurve::acceptWrapperNameCollector(BootstrappedYieldCurve* yc,
                                               WrapperNameCollector* collector)
{
    collector->addName(yc->getName());
}

void BootstrappedYieldCurve::acceptYieldNameCollector(BootstrappedYieldCurve* yc,
                                             YieldNameCollector* collector)
{
    collector->addName(yc->getName());
}

/** Records name of iso code against yc name in market data object.
    It is invoked ONCE only
    - immediately after this object is placed in the cache. */
void BootstrappedYieldCurve::initialise(MarketData* market)
{
    market->setYieldCurveISOCode(name, ccy);
}

/** populate from market cache */
void BootstrappedYieldCurve::getMarket(const IModel* model, const MarketData* market)
{
    if (today.empty())
    {
        today = market->GetReferenceDate();
    }

    hols.getData(model, market);
    if (!ccyBasis.isEmpty())
    {
        ccyBasis.getData(model, market);
    }

    // if specified, ask for the vol
    if (!irVol.isEmpty())
    {
        // NB This might result in a null vol if the model reckons we don't
        // need it
        irVol.getData(model, market);
    }

    // calculate transient fields
    valueDate = hols->addBusinessDays(today, spotOffset);
    expiryCache();
}

/** return either the growth or discount curve  */
ZeroCurveSP BootstrappedYieldCurve::get(bool growthCurve) const
{
    static const string method = "BootstrappedYieldCurve::get";

    try
    {
        //TODO : could use ZCAccessor here to have "local" caching
        ZeroPairConstSP zp = CashSwapCurveGlobalCache::get(this);

        // clone the ZeroCurve to avoid modifying data in cache !
        const ZeroCurve* zc = zp->get(growthCurve);
        return ZeroCurveSP(dynamic_cast<ZeroCurve*>(zc->clone()));
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}


CashFlowArraySP BootstrappedYieldCurve::getRatesAndDates() const
{
    static const string method = "BootstrappedYieldCurve::getRatesAndDates";

    try
    {
        const ZeroCurve& zcurve = zc.get(*this);
        return zcurve.getRatesAndDates();
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}


//-----------------------------
//Duration::IParhandler methods
//-----------------------------

//return benchmarks for par instruments
//these will determine the durations to calculate unless specified
ExpiryArrayConstSP BootstrappedYieldCurve::getParBenchmarks() const
{
    ExpiryArraySP result(new ExpiryArray());
    for (int i = 0 ; i < benchmarks.size() ; i++)
    {
        if (benchmarks[i]->isCash() || benchmarks[i]->isSwap())
        {
            result->push_back(ExpirySP(const_cast<Expiry*>(benchmarks[i]->getEnd())));
        }
    }

    return ExpiryArrayConstSP(result.release());
}

//return closed form solution
//must be implemented if supportsDurationClosedForm = true
ExpiryResultArraySP BootstrappedYieldCurve::getDuration(const Duration* durationObj) const
{
    static const string method = "BootstrappedYieldCurve::getDuration";

    try
    {
        ExpiryArrayConstSP benchmarks = durationObj->getBenchmarks();
        if (!benchmarks) {
            benchmarks = getParBenchmarks();
        }

        //build storage for the duration results
        ExpiryResultArraySP durCalcs(new ExpiryResultArray(0));

        //for each tweakPoint
        for(int j=0;j<benchmarks->size();j++) {
            //get duration for this benchmark
            ExpirySP maturity((*benchmarks)[j]);

            //find maturity in BootstrappedYieldCurve expiries
            //-throws exception if not found
            int i = maturity->search(allExpiries.get());

            //store maturity date of this par instrument;
            DateTime matDate = this->benchmarks[i]->getBenchmarkDate(valueDate);
            //holiday adjust TODO

            double duration;

            if (this->benchmarks[i]->isCash())
            {
                //equivalent to ALIB_MM_SENS
                double yearFraction = floatDcc->years(valueDate, matDate);
                duration = yearFraction / (1 + (this->benchmarks[i]->getRate() * yearFraction));
            }
            else if (this->benchmarks[i]->isFuture())
            {
                continue;
            }
            else if (this->benchmarks[i]->isTurn())
            {
                throw ModelException(method, "Curve cannot have adjustments");
            }
            else if (this->benchmarks[i]->isSwap())
            {
                //clone this untweaked curve
                BootstrappedYieldCurveSP aCloneCSC(cloneYieldCurve(true));
                aCloneCSC->setProjectionCurve(true);

                //shift the appropriate rate by 1 bp
                aCloneCSC->benchmarks[i]->setRate(aCloneCSC->benchmarks[i]->getRate() + 0.0001);

                //equivalent to ALIB_BOND_DUR_EFF
                //with annual frequency and tweaked zero rate
                double zeroRate = aCloneCSC->zero(matDate);
                //always act/365F year fraction
                double yearFraction = matDate.daysDiff(valueDate)/365.0;
                duration = (1-pow(1+zeroRate,-yearFraction)) / aCloneCSC->benchmarks[i]->getRate();
            }

            //create result element
            ExpiryResult bmarkDuration(maturity, duration);
            //and store
            durCalcs->push_back(bmarkDuration);
        }

        return durCalcs;
    }
    catch (exception &e)
    {
        throw ModelException(e, method);
    }
}

//return a means of tweaking the MarketObject in a pointwise manner
//assumption is that the results of this twaeker is an ExpiryResultArray.....
VectorRiskPropertySensitivityConstSP BootstrappedYieldCurve::getPointwiseTweaker() const
{
    //1 bp shift explictly required
    return VectorRiskPropertySensitivityConstSP(new RhoPointwise(0.0001));
}

//return a par ins  5{                    ôX üù x "  z  Ø üúˇˇ íÿHíÿ     xÉÿíÿ    "         ïb ôX 	¥ x  3ê  ˆ ôWˇˇ íÿêíÿX     	¥ÉêíÿX    x         î  ïb ¸  	¥  ûç  ù ïaˇˇ íÿÿíÿ†     ¸ ÉHíÿ†    	¥         î  î  Ë  ¸          î ˇˇíŸ íÿË     Ë É íÿË    ¸          ë= î  Ãó Ë   ûÃ  ú ì˛ˇˇ íŸhíŸ0     ÃóÉ∏íŸ0    Ë          âP ë= ≥ Ãó  3à  ˜ ë;ˇˇ íŸ∞íŸx     ≥ÉpíŸx    Ãó         |∆ âP õˆ ≥  z  Ø âOˇˇ íŸ¯íŸ¿     õˆÉ(íŸ¿    ≥         l+ |∆ á’ õˆ  e  3 |≈ˇˇ í⁄@í⁄     á’É‡í⁄    õˆ         X
 l+ w: á’  3  e l*ˇˇ í⁄àí⁄P     w:Éòí⁄P    á’         @ X
 j∞ w:  Ø  z X	ˇˇ í⁄–í⁄ò     j∞ÉPí⁄ò    w:         'i @ b√ j∞  ˜  3à @Ôˇˇ í€í⁄‡     b√Éí⁄‡    j∞           'i `  b√  ú  ûÃ 'hˇˇ í€`í€(     ` É¿í€(    b√                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          @      Ä                                        @      Ä                                          @      Ä                                          @                                                                                                                                                                                                                                                                