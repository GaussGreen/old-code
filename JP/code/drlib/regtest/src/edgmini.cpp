/*************************************************************
 * Equity Derivatives Research: Example of use of Models Library Interface
 *
 *
*************************************************************/
/* this file is not needed - just used for precompiled headers */
#include "edginc/config.hpp"
#include "edginc/EDRInterface.h"

static void handleError(const char* routine){
    if (!EdrErrorCleared()){
        char* errorStack = EdrErrorStackTrace();
        if (errorStack) {
            printf(errorStack);
            free(errorStack);
        }
    }
    printf("%s: Failed\n", routine);
}

static EDRBool displayName(EDROutputName outputName){
    EDRBool     ok = EDR_FALSE;
    int numIDs;
    int i;
    if (!EdrOutputNameIDCount(outputName, &numIDs)){
        goto done;
    }
    printf("Name is '");
    for (i = 0; i < numIDs; i++){
        const char* componentID = EdrOutputNameIdGet(outputName, i);
        if (!componentID){
            goto done;
        }
        printf("%s%s", i > 0? ".": "", componentID);
    }
    printf("'\n");
    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("displayName");
    }
    return ok;
}

 /* print an ExpiryResultArray on the screen */
static EDRBool displayVector(const char* label, EDRObject vector) {
    EDRBool     ok = EDR_FALSE;
    EDRDataDict dd = 0;
    EDRObject   point = 0;
    EDRObject   value = 0;
    EDRObject   bm = 0;
    int         length;
    int         i;
    double      d;
    char*       bmp = 0;

    printf("%s\n", label);

    if (!EdrArrayLength(vector, &length)) {
        goto done;
    }

    for (i = 0; i < length; i++) {
        if (!(point = EdrArrayItem(vector, i))) {
            goto done;
        }
        if (!(dd = EdrObjectToDataDict(point))) {
            goto done;
        }
        if (!(value = EdrDataDictGet(dd, "result"))) {
            goto done;
        }
        if (!(bm = EdrDataDictGet(dd, "expiry"))) {
            goto done;
        }

        if (!EdrDoubleGet(value, &d)) {
            goto done;
        }

        if (!(bmp = EdrExpiryGet(bm))) {
            goto done;
        }
        printf("%3s %9f\n", bmp, d);

        EdrObjectFree(point);   
        EdrObjectFree(dd);   
        EdrObjectFree(value);   
        EdrObjectFree(bm);   
        free(bmp);
        point = 0;
        dd = 0;
        value = 0;
        bm = 0;
        bmp = 0;
    }

    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("displayVector");
    }
    EdrObjectFree(point);   
    EdrObjectFree(dd);   
    EdrObjectFree(value);   
    EdrObjectFree(bm);   
    free(bmp);
    return ok;
}

/** Display the results returned from the pricing call */
static EDRBool displayResults(EDRObject results){
    EDRBool        ok = EDR_FALSE;
    EDROutputName  name = 0;
    EDRObject      object = 0;
    const char*    resultsType = EdrObjectGetType(results);
    const char*    packet;
    EDRBool        iterStatus;

    EDRResultsIterator iter = 0;
    if (!results){
        goto done;
    }
    if (strcmp(resultsType, "Results")){
        printf("Can't cope with result object of type %s\n", resultsType);
        goto done;
    }
    
    if (!(iter = EdrResultsIteratorGet(results))){
        goto done;
    }

    while ((iterStatus = EdrResultsIteratorNext(iter, &packet,
                                                &name,
                                                &object)) &&
           packet){
        /* find out what sort of object the result is */
        const char* objType = EdrObjectGetType(object);
        if (!objType){
            goto done;
        }
        printf("%s:\n", packet);
        /* print out name  */
        if (!displayName(name)){
            goto done;
        }
        /* then display result (switch on type) */
        if (strcmp(objType, "Double") == 0){
            double value;
            if (!EdrDoubleGet(object, &value)) {
                goto done;
            }
            printf("Value is (%f)\n", value);
        } else if (strcmp(objType, "ExpiryResultArray") == 0){
            displayVector("Value is:", object);
        } else {
            printf("Don't know how to display output of type %s\n",
                   objType);
        }
    }
    ok = iterStatus;
done:
    EdrObjectFree(iter);
    if (!ok) {
        handleError("displayResults");
    }
    return ok;
}
                                 
/* build the benchmark dates for a yield curve */
static EDRObject benchmarks(void) {
    EDRBool ok = EDR_FALSE;
    static char* periods[] = {"ON",
                              "1M",
                              "2M",
                              "3M",
                              "6M",
                              "1Y",
                              "2Y",
                              "3Y",
                              "5Y",
                              "7Y",
                              "10Y",
                              "20Y",
                              "30Y"};

    int size = sizeof(periods)/sizeof(char*);
    int i;

    EDRObject bm = EdrArrayNew("Expiry", size);
    EDRObject expiry = 0;
    if (!bm){
        goto done;
    }

    for (i = 0; i < size; i++) {
        if (!(expiry = EdrMaturityPeriodNew(periods[i]))){
            goto done;
        }
        if (!EdrArraySet(bm, i, expiry)){
            goto done;
        }
        EdrObjectFree(expiry);
        expiry = 0;
    }


    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("benchmarks");
        EdrObjectFree(expiry);
    }
    return bm;
}

/* build the rates for a yield curve */
static EDRObject rates(void) {
    EDRBool ok = EDR_FALSE;
    static double rts[] = {0.05,
                           0.051,
                           0.052,
                           0.053,
                           0.054,
                           0.055,
                           0.056,
                           0.057,
                           0.058,
                           0.059,
                           0.06,
                           0.061,
                           0.062};

    int size = sizeof(rts)/sizeof(double);
    int i;

    EDRObject db = 0;
    EDRObject rt = EdrArrayNew("Double", size);
    if (!rt){
        goto done;
    }
    
    for (i = 0; i < size; i++) {
        if (!(db = EdrDoubleNew(rts[i]))) {
            goto done;
        }
        if (!EdrArraySet(rt, i, db)){
            goto done;
        }
        EdrObjectFree(db);
        db = 0;
    }

    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("rates");
        EdrObjectFree(db);
    }
    return rt;
}

/* build the list of cash/swap instruments for a yeild curve */
static EDRObject instruments(void) {
    EDRBool ok = EDR_FALSE;
    static char* inst[] = {"M",
                           "M",
                           "M",
                           "M",
                           "M",
                           "M",
                           "S",
                           "S",
                           "S",
                           "S",
                           "S",
                           "S",
                           "S"};

    int size = sizeof(inst)/sizeof(char*);
    int i;

    EDRObject st = 0;
    EDRObject instr = EdrArrayNew("String", size);

    for (i = 0; i < size; i++) {
        if (!(st = EdrStringNew(inst[i]))){
            goto done;
        }
        if (!EdrArraySet(instr, i, st)){
            goto done;
        }
        EdrObjectFree(st);
        st = 0;
    }

    ok = EDR_TRUE;
done:
    if (!ok) {
       handleError("instruments");
       EdrObjectFree(st);
    }
    return instr;
}

/* build a 4+i zero curve interpolator for a yield curve */
static EDRObject interp(void) {
    EDRBool ok = EDR_FALSE;
    EDRDataDict dd = 0;
    EDRObject   zci = 0;

    if (!(dd = EdrDataDictNew("FourPlusI"))) {
        goto done;
    }
      
    if (!(zci = EdrDataDictToObject(dd))) {
        goto done;
    }

    ok = EDR_TRUE;
done:
    if (!ok) {
       handleError("interp");
    }
    EdrObjectFree(dd);
    return zci;
}

/* build a "NO_WEEKENDS" style holiday */
static EDRObject noweekends(void) {
    EDRBool ok = EDR_FALSE;
    EDRDataDict dd = 0;
    EDRObject   hol = 0;
    EDRObject   wknds = 0;
    EDRObject   hols = 0;
    EDRObject   name = 0;

    if (!(name = EdrStringNew("noweekends"))) {
        goto done;
    }
    if (!(wknds = EdrBoolNew(EDR_FALSE))) {
        goto done;
    }
    if (!(hols = EdrArrayNew("DateTime",
                             0))) {
        goto done;
    }

    if (!(dd = EdrDataDictNew("Holiday"))) {
        goto done;
    }

    if (!EdrDataDictAdd(dd, "useWeekends", wknds)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "name", name)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "holidays", hols)) {
        goto done;
    }
    
    if (!(hol = EdrDataDictToObject(dd))) {
        goto done;
    }

    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("noweekends");
    }

    EdrObjectFree(wknds);
    EdrObjectFree(hols);
    EdrObjectFree(name);
    EdrObjectFree(dd);
    return hol;
}

/* build a cash/swap style yield curve */
static EDRObject yieldCurve(void)
{
    EDRBool ok = EDR_FALSE;

    EDRObject   yc = 0;
    EDRDataDict dd = 0;
    EDRObject   ycbm = 0;
    EDRObject   ycrt = 0;
    EDRObject   ycinst = 0;
    EDRObject   ccy = 0;
    EDRObject   name = 0;
    EDRObject   today = 0;
    EDRObject   spotOffset = 0;
    EDRObject   moneyMarketDenom = 0;
    EDRObject   swapDayCount = 0;
    EDRObject   swapFrequency = 0;
    EDRObject   zci = 0;
    EDRObject   hols = 0;

    printf("Building yield curve\n");

    if (!(ycbm = benchmarks())) {
        goto done;
    }
    if (!(hols = noweekends())) {
        goto done;
    }
    
    if (!(ycrt = rates())) {
        goto done;
    }

    if (!(ycinst = instruments())) {
        goto done;
    }

    if (!(ccy = EdrStringNew("GBP"))) {
        goto done;
    }
    if (!(name = EdrStringNew("GBP-LIBOR"))) {
        goto done;
    }
    if (!(today = EdrDateTimeNew("05-Apr-2001", "SOD"))) {
        goto done;
    }
    if (!(spotOffset = EdrIntNew(0))) {
        goto done;
    }
    if (!(moneyMarketDenom = EdrIntNew(365))) {
        goto done;
    }
    if (!(swapFrequency = EdrIntNew(2))) {
        goto done;
    }
    if (!(swapDayCount = EdrStringNew("Act/360"))) {
        goto done;
    }
    if (!(zci =interp())) {
        goto done;
    }

    if (!(dd = EdrDataDictNew("CashSwapCurve::Interface"))) {
        goto done;
    }
      
    if (!EdrDataDictAdd(dd, "ccy", ccy)) {
        goto done;
    }

    if (!EdrDataDictAdd(dd, "name", name)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "today", today)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "spotOffset", spotOffset)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "moneyMarketDenom", moneyMarketDenom)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "swapFrequency", swapFrequency)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "swapDayCount", swapDayCount)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "hols", hols)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "expiries", ycbm)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "instruments", ycinst)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "rates", ycrt)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "interp", zci)) {
        goto done;
    }
    
    if (!(yc = EdrDataDictToObject(dd))) {
        goto done;
    }
 
    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("yieldCurve");
    }
    EdrObjectFree(ccy);
    EdrObjectFree(name);
    EdrObjectFree(today);
    EdrObjectFree(spotOffset);
    EdrObjectFree(moneyMarketDenom);
    EdrObjectFree(swapDayCount);
    EdrObjectFree(swapFrequency);
    EdrObjectFree(dd);
    EdrObjectFree(ycbm);
    EdrObjectFree(ycrt);
    EdrObjectFree(ycinst);
    EdrObjectFree(zci);
    EdrObjectFree(hols);
    return yc;
}

/* build a cashflow */
static EDRObject cashflow(void) {
    EDRBool ok = EDR_FALSE;
    EDRDataDict dd = 0;
    EDRObject   date = 0;
    EDRObject   amount = 0;
    EDRObject   cf = 0;

    if (!(date = EdrDateTimeNew("15-Sep-2003", "EOD"))) {
        goto done;
    }
    if (!(amount = EdrDoubleNew(100.0))) {
        goto done;
    }

    if (!(dd = EdrDataDictNew("CashFlow"))) {
        goto done;
    }

    if (!EdrDataDictAdd(dd, "date", date)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "amount", amount)) {
        goto done;
    }
    
    if (!(cf = EdrDataDictToObject(dd))) {
        goto done;
    }

    ok = EDR_TRUE;
done:
    if (!ok) {
       handleError("cashflow");
    }
    EdrObjectFree(date);
    EdrObjectFree(amount);
    EdrObjectFree(dd);
    return cf;
}

/* build a CashFlowStream instrument */
static EDRObject cfinst(void) {
    EDRBool ok = EDR_FALSE;
    EDRDataDict dd = 0;
    EDRObject   cfl = 0;
    EDRObject   discount = 0;
    EDRObject   inst = 0;
    EDRObject   cf = 0;

    printf("Building CashFlowStream\n");

    if (!(discount = EdrStringNew("GBP-LIBOR"))) {
        goto done;
    }
    if (!(cf = cashflow())) {
        goto done;
    }

    if (!(cfl = EdrArrayNew("CashFlow",
                             1))) {
        goto done;
    }
    if (!EdrArraySet(cfl, 0, cf)){
        goto done;
    }

    if (!(dd = EdrDataDictNew("CashFlowStream"))) {
        goto done;
    }

    if (!EdrDataDictAdd(dd, "cfl", cfl)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "discount", discount)) {
        goto done;
    }
    
    if (!(inst = EdrDataDictToObject(dd))) {
        goto done;
    }

    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("cfinst");
    }

    EdrObjectFree(cf);
    EdrObjectFree(cfl);
    EdrObjectFree(discount);
    EdrObjectFree(dd);
    return inst;
}

/* create a simple scalar greek SensControl */
static EDRObject scalarGreek(const char* greek, double shift) {
    EDRBool ok = EDR_FALSE;
    EDRDataDict dd = 0;
    EDRObject   shiftSize = 0;
    EDRObject   sensctrl = 0;

    printf("Building %s\n", greek);

    if (!(shiftSize = EdrDoubleNew(shift))) {
        goto done;
    }

    if (!(dd = EdrDataDictNew(greek))) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "shiftSize", shiftSize)) {
        goto done;
    }
      
    if (!(sensctrl = EdrDataDictToObject(dd))) {
        goto done;
    }

    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("scalarGreek");
    }
    EdrObjectFree(dd);
    EdrObjectFree(shiftSize);
    return sensctrl;
}

/* build a Control */
static EDRObject control(void) {
    EDRBool ok = EDR_FALSE;
    EDRDataDict dd = 0;
    EDRObject   sens = 0;
    EDRObject   writeToFile = 0;
    EDRObject   ctrl = 0;
    EDRObject   rhopll = 0;
    EDRObject   rhopw = 0;

    printf("Building Control\n");

    if (!(writeToFile = EdrBoolNew(EDR_FALSE))) {
        goto done;
    }

    if (!(rhopll = scalarGreek("RhoParallel", 0.0001))) {
        goto done;
    }
    if (!(rhopw = scalarGreek("RhoPointwise", 0.0001))) {
        goto done;
    }

    /* take your pick on how to build arrays */
#if 1
    if (!(sens = EdrArrayNew("Sensitivity",
                             2))) {
        goto done;
    }
    if (!EdrArraySet(sens, 0, rhopll) || !EdrArraySet(sens, 1, rhopw)){
        goto done;
    }
#else
    if (!(sens = EdrArrayNew("Sensitivity", 0))) {
        goto done;
    }
    if (!EdrArrayAppend(sens, rhopll) || !EdrArrayAppend(sens, rhopw)){
        goto done;
    }
#endif
    if (!(dd = EdrDataDictNew("Control"))) {
        goto done;
    }

    if (!EdrDataDictAdd(dd, "sens", sens)) {
        goto done;
    }
    if (!EdrDataDictAdd(dd, "writeToFile", writeToFile)) {
        goto done;
    }
    
    if (!(ctrl = EdrDataDictToObject(dd))) {
        goto done;
    }

    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("control");
    }

    EdrObjectFree(sens);
    EdrObjectFree(rhopll);
    EdrObjectFree(rhopw);
    EdrObjectFree(writeToFile);
    EdrObjectFree(dd);
    return ctrl;
}

/* build a ClosedForm model */
static EDRObject model(void) {
    EDRBool ok = EDR_FALSE;
    EDRDataDict dd = 0;
    EDRObject   mdl = 0;

    printf("Building ClosedForm\n");

    if (!(dd = EdrDataDictNew("ClosedForm"))) {
        goto done;
    }
      
    if (!(mdl = EdrDataDictToObject(dd))) {
        goto done;
    }

    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("model");
    }
    EdrObjectFree(dd);
    return mdl;
}

/* demo how to build a simple cashflow instrument, get price & greeks,
   & display on screen. Rename this to main() to get a self contained
   executable */
REGTEST_DLL int edgminiMain (void)
{
    EDRBool ok = EDR_FALSE;

    int    components;
    int    i;
    char*  typekey;

    char* ccy = 0;

    EDRObject   today = 0;
    EDRObject   mdl = 0;
    EDRObject   inst = 0;
    EDRObject   ctrl = 0;
    EDRObject   discount = 0;

    EDRMarketDataCache mkt = 0;
    EDRObject          results = 0;

    EDRObject rho = 0;
    EDRObject rhopw = 0;
    double    multiplier = 10.0;
    double    weight = 1.0;

    /* turn library on */
    if (!EdrStartup()) {
        goto done;
    }

    /* what's in a type then ? */
    /*typekey = "CashSwapCurve::Interface"; */

    /* this bit just shows how to query a type to see what it's made of */
    typekey = "Equity"; 

    if (!EdrTypeNumItems(typekey, &components)) {
        goto done;
    }

    printf("A %s consists of:\n", typekey);
    for (i = 0; i < components; i++) {

        printf("%2d: field (%20s) type (%20s) description (%30s)\n",
               i,
               EdrTypeFieldName(typekey, i),
               EdrTypeFieldType(typekey, i),
               EdrTypeFieldDescription(typekey, i));
    }
    
    /* now build up the inputs to the RiskManager call */
    if (!(today = EdrDateTimeNew("05-Apr-2001", "SOD"))) {
        goto done;
    }
           
    if (!(mdl = model())) {
        goto done;
    }

    if (!(inst = cfinst())) {
        goto done;
    }
    if (!(ctrl = control())) {
        goto done;
    }

    printf("Building MarketDataCache\n");

    if (!(mkt = EdrMarketDataCacheNew(today))) {
        goto done;
    }
    /* create and then add yield curve to market data */
    if (!(discount = yieldCurve())) {
        goto done;
    }
    if (!EdrMarketDataCacheAdd(mkt, discount)){
        goto done;
    }
    
    printf("Running EdrMain\n");

    if (!(results = EdrMain(1,
                            &mdl,
                            &inst,
                            &ctrl,
                            &weight,
                            &multiplier,
                            0,    /* no scenarios */
                            mkt))) {
        goto done;
    }

    if (!displayResults(results)){
        goto done;
    }
    ok = EDR_TRUE;
done:
    EdrObjectFree(discount);
    EdrObjectFree(today);
    EdrObjectFree(mkt);
    EdrObjectFree(results);
    EdrObjectFree(mdl);
    EdrObjectFree(inst);
    EdrObjectFree(ctrl);
    EdrObjectFree(rho);
    EdrObjectFree(rhopw);
    free(ccy);

    if (!ok) {
        handleError("main");
    }        
    /* turn library off */
    EdrShutdown();
    return 0;
}
