/*************************************************************
 * Equity Derivatives Research: Example of use of DR Interface
 *
 *
 *************************************************************/
/* this file is not needed - just used for precompiled headers */
#include "edginc/config.hpp"
#include "edginc/EDGServices.h"
#include "edginc/DRUtil.hpp"

// Invoke only if (x) is a function that returns DRError or is a DRError type.
#define CHECK(x)    \
    DRI_CHECK(x, method, "", handleError, 0, ;, ;, ;)

// Invoke only if (x) is a function that returns DRError or is a DRError type.
#define FLUNK_IF(x) \
    DRI_CHECK(x, method, "", handleError, 0, ;, goto done, ;)
    
#define ERROR_EXIT(err)                         \
    {                                           \
        handleError(method, err, 0);            \
        goto done;                              \
    }

/* to save passing this around everywhere */
static DRService* g_svc = 0;

/************************* Utility Functions ***************************/

static void handleError(const char* routine, const char *err, 
                        void *cbParam = 0) {
    if (err) {
        printf("%s() failed: %s\n", routine, (char*) err);
    }
    else {
        printf("%s() failed!\n", routine);
    }
}

/** Make a string */
static DRValue mkStr(const char* str){
    static const char method[] = "mkStr";
    DRValue x;
    CHECK(EDGSDRValue(&x, DR_STRING, str));
    return x;
}

/** Make an int */
static DRValue mkInt(int val){
    static const char method[] = "mkInt";
    DRValue x;
    CHECK(EDGSDRValue(&x, DR_INT, val));
    return x;
}

/** Make a double */
static DRValue mkDb(double val){
    static const char method[] = "mkDb";
    DRValue x;
    CHECK(EDGSDRValue(&x, DR_DOUBLE, val));
    return x;
}

/** Make a bool */
static DRValue mkBl(DRBool val){
    static const char method[] = "mkBl";
    DRValue x;
    CHECK(EDGSDRValue(&x, DR_BOOL, val));
    return x;
}

/** Make a DoubleMatrix */
static DRValue mkDbMat(int cols, int rows, ...){
    static const char method[] = "mkDbMat";
    va_list        args;
    int            i, j;
    DRMatrix       obj = 0;
    DRValue        x;
    DRError        err = 0;

    CHECK(g_svc->fptr->matrixNew(g_svc, cols, rows, &obj));

    va_start(args, rows);
    for (i = 0; i < cols; i++){
        for (j = 0; j < rows; j++){
            double val = va_arg(args, double);
            if ((err = g_svc->fptr->matrixSet(g_svc, obj, i, j, val))) {
                DRI_FREE(g_svc, obj);
                DRI_ERROR_FREE(g_svc, err);
                x.type = DR_UNDEFINED;
                return x;
            }
        }
    }
    va_end(args);
    x.type = DR_OBJECT;
    x.value.object = obj;
    return x;
}


/** Create an array of strings */
static DRValue mkStrArray(int length, ...){
    static const char method[] = "mkStrArray";
    va_list     args;
    int         i;
    DRArray     array = 0;
    DRValue     x;

    FLUNK_IF(g_svc->fptr->arrayNew(g_svc, length, DRI_TYPE_STRING, &array));
    va_start(args, length);
    for (i = 0; i < length; i++){
        const char* str = va_arg(args, const char*);
        DRValue obj = mkStr(str);
        FLUNK_IF(g_svc->fptr->arraySet(g_svc, array, i, &obj));
        CHECK(g_svc->fptr->valueClear(g_svc, &obj));
    }
    va_end(args);
    x.type = DR_OBJECT;
    x.value.object = array;
    return x;
done:
    DRI_FREE(g_svc, array);
    x.type = DR_UNDEFINED;
    return x;
}

/** Create an array of ints */
#if 0
/* not currently used */
static EDRObject mkIntArray(int length, ...){
    static const char method[] = "mkIntArray";
    va_list     args;
    int         i;
    EDRObject array = EdrArrayNew(DRI_TYPE_INT, length);
    if (!array) {
        ERROR_EXIT(0);
    }
    va_start(args, length);
    for (i = 0; i < length; i++){
        int val = va_arg(args, int);
        EDRObject obj = EdrIntNew(val);
        if (!obj || !EdrArraySet(array, i, obj)){
            ERROR_EXIT(0);
        }
        EdrObjectFree(obj);
    }
    va_end(args);
    return array;
done:
    return 0;
}
#endif

/** Create an array of doubles */
static DRValue mkDbArray(int length, ...){
    static const char method[] = "mkDbArray";
    va_list     args;
    int         i;
    DRArray     array = 0;
    DRValue     x;

    FLUNK_IF(g_svc->fptr->arrayNew(g_svc, length,  DRI_TYPE_DOUBLE, &array));
    va_start(args, length);
    for (i = 0; i < length; i++){
        double val = va_arg(args, double);
        DRValue obj = mkDb(val);
        FLUNK_IF(g_svc->fptr->arraySet(g_svc, array, i, &obj));
        CHECK(g_svc->fptr->valueClear(g_svc, &obj));
    }
    va_end(args);
    x.type = DR_OBJECT;
    x.value.object = array;
    return x;
done:
    DRI_FREE(g_svc, array);
    x.type = DR_UNDEFINED;
    return x;
}

/** Create an array of specified type, inputs terminated by 0. Note inputs
    are also freed */
static DRValue mkObjArray(const char* type, int length, ...){
    static const char method[] = "mkObjArray";
    va_list     args;
    DRArray     array = 0;
    DRValue     x;
    int         i;

    FLUNK_IF(g_svc->fptr->arrayNew(g_svc, 0, type, &array));

    va_start(args, length);
    DRValue val;
    val.type = DR_UNDEFINED;
    for (i = 0; i < length; i++) {
        val = va_arg(args, DRValue);
        FLUNK_IF(EDGSArrayAppend(g_svc, array, &val));
        CHECK(g_svc->fptr->valueClear(g_svc, &val));
    } 

    va_end(args);
    x.type = DR_OBJECT;
    x.value.object = array;
    return x;
done:
    x.type = DR_UNDEFINED;
    return x;
}

/** Create a DR Wrapper of a specified type, inputs terminated by 0. Note inputs
    are also freed */
static DRValue mkWrapper2(const char* type, va_list args){
    static char  method[] = "mkWrapper2";
    const char*  field;
    DRValue      x;
    DRMap        drw = 0;
    DRObject     dd  = 0;
    DRValue      val;

    FLUNK_IF(g_svc->fptr->mapNew(g_svc, "DRWrapper", &drw));
    FLUNK_IF(EDGSDRValue(&val, DR_STRING, type));
    FLUNK_IF(g_svc->fptr->mapAddItem(g_svc, drw, "type", &val));
    FLUNK_IF(g_svc->fptr->mapToObject(g_svc, drw, &dd));
    do {
        field = va_arg(args, const char*);
        if (field){
            val = va_arg(args, DRValue);
            if (val.type == DR_UNDEFINED) {
                char tmp[100];
                sprintf(tmp, "value undefined when setting %s in wrapper %s",
                        field, type);
                ERROR_EXIT(tmp);
            }
            FLUNK_IF(g_svc->fptr->mapAddItem(g_svc, dd, field, &val));
            CHECK(g_svc->fptr->valueClear(g_svc, &val));
        }
    } while (field);
    x.type = DR_OBJECT;
    x.value.object = dd;
    return x;  
done:
    DRI_FREE(g_svc, drw);
    x.type = DR_UNDEFINED;
    return x;
}

static DRValue mkWrapper(const char* type, ...){
    va_list        args;
    va_start(args, type);
    return mkWrapper2(type, args);
}

/** Create an object of a specified type, inputs terminated by 0. Note inputs
    are also freed */
static DRValue mkObj2(const char* type, va_list args){
    static char   method[] = "mkObj2";
    const char*   field;
    DRMap         dd = 0;
    DRObject      obj = 0;
    DRValue       val;
    DRValue       x;

    FLUNK_IF(g_svc->fptr->mapNew(g_svc, type, &dd));
    do {
        field = va_arg(args, const char*);
        if (field){
            val = va_arg(args, DRValue);
            if (val.type == DR_UNDEFINED) {
                char tmp[100];
                sprintf(tmp, "value undefined when setting %s in map %s",
                        field, type);
                ERROR_EXIT(tmp);
            }
            FLUNK_IF(g_svc->fptr->mapAddItem(g_svc, dd, field, &val));
            CHECK(g_svc->fptr->valueClear(g_svc, &val));
        }
    } while (field);
    //cout << "Building " << type << endl;
    FLUNK_IF(g_svc->fptr->mapToObject(g_svc, dd, &obj));
    DRI_FREE(g_svc, dd);
    x.type = DR_OBJECT;
    x.value.object = obj;
    //{ DRString myString = EDGSValueXMLWrite(g_svc, &x);
    //cout << myString << endl; }
    return x;    
done:
    DRI_FREE(g_svc, dd);
    x.type = DR_UNDEFINED;
    return x;
}

/** Create an object of a specified type, inputs terminated by 0. Note inputs
    are also freed */
static DRValue mkObj(const char* type, ...){
    va_list        args;
    va_start(args, type);
    return mkObj2(type, args);
}

/** Creates object and stores it in supplied cache */
static DRError addObj(DRObject cache, const char* type, ...){
    static const char method[] = "addObj";
    va_list        args;
    va_start(args, type);
    DRValue obj = mkObj2(type, args);
    DRError err = 0;

    if (obj.type == DR_UNDEFINED) {
        char tmp[100];
        sprintf(tmp, "value undefined when adding to map %s", type);
        ERROR_EXIT(tmp);
    }
    FLUNK_IF(g_svc->fptr->mapAddItem(g_svc, cache, "", &obj));
    CHECK(g_svc->fptr->valueClear(g_svc, &obj));
    va_end(args);
    return 0;
done:
    //handleError("addObj", line, err);
    return err;
}

/** Writes object to stream and then reads it back again */
static DRObject testXMLReadWrite(DRValue* inObject){
    // read existing xml file into char*
    static const char method[] = "testXMLReadWrite";
    DRValue           object;

    DRString myString = EDGSValueXMLWrite(g_svc, inObject);
    if (!myString) {
        ERROR_EXIT("failure to stream object to XML");
    }
    FLUNK_IF(EDGSValueXMLRead(g_svc, myString, &object));
    DRI_STRING_FREE(g_svc, myString);
    return object.value.object;
done:
    return 0;
}

static DRValue mkDt(const char* val1, const char* val2){
    DRValue dt = mkObj("DateTime",
                       "date", mkStr(val1),
                       "time", mkStr(val2),
                       0);
    return dt;
}

/** Make a MaturityPeriod */
static DRValue mkMat(const char* val){
    DRValue x = mkObj("MaturityPeriod",
                      "period", mkStr(val),
                      0);
    return x;
}

/********************** End of Utility Functions *************************/

static const char* todayDate = "05-Apr-2001";
static const char* todayTime = "SOD";

// make a MarketData
static DRValue mkMarket(){
    DRValue x = mkObj("MarketData",
                      "today", mkDt(todayDate, todayTime),
                      0);
    return x;
}

/* these 'defines' highlight which names must match between the
   different market data objects and instruments. The format of the
   name is unimportant but certain characters are best avoided (eg '&'
   will cause problems in the xml). Note that the key requirement is
   that names used for the same type of data MUST be different Eg
   putting two cash swap curves both called "GBP-LIBOR" into the cache
   will result in the first one being removed.  */
#define FTSE "FTSE"
#define SPX "SPX"
#define DAX "DAX"
#define GBP_LIBOR "GBP-LIBOR"
#define USD_LIBOR "USD-LIBOR"
#define EUR_LIBOR "EUR-LIBOR"
#define FTSE_VOL  "FTSE-Vol"
#define DAX_VOL   "DAX-VolSurf"
#define SPX_VOL   "SPX-Vol"
#define FX_USD_GBP_VOL "GBP-USD_Vol"
#define FX_GBP_USD "USD/GBP"
#define FX_USD_EUR_VOL "EUR-USD_Vol"
#define FX_EUR_USD "USD/EUR"

/* Make a CashSwapCurve object */
static DRValue mkCashSwapCurve1() {
    return mkObj("CashSwapCurve::Interface",
                 "ccy", mkStr("GBP"),
                 "name", mkStr(GBP_LIBOR),
                 "today", mkDt(todayDate, todayTime),
                 "spotOffset", mkInt(0),
                 "moneyMarketDenom", mkInt(365),
                 "swapFrequency", mkInt(2),
                 "swapDayCount", mkStr("Act/360"),
                 "hols", mkObj("Holiday", 
                               "name", mkStr("noweekends"),
                               "useWeekends", mkBl(DR_FALSE),
                               "holidays", mkObjArray("DateTime", 0),
                               0),
                 "expiries", mkObjArray("Expiry", 13,
                                        mkMat("ON"), mkMat("1M"), 
                                        mkMat("2M"), mkMat("3M"),  
                                        mkMat("6M"), mkMat("1Y"), 
                                        mkMat("2Y"), mkMat("3Y"),
                                        mkMat("5Y"), mkMat("7Y"), 
                                        mkMat("10Y"), mkMat("20Y"), 
                                        mkMat("30Y"), 0),
                 "instruments", mkStrArray(13, "M", "M", "M", "M", "M", "M",
                                           "S", "S", "S", "S", "S", "S", "S"),
                 "rates", mkDbArray(13, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,
                                    0.05,  0.05, 0.05, 0.05, 0.05, 0.05,
                                    0.05),
                 "interp", mkObj("FourPlusI", 0),
                 0);
}

/* Just inserts a GBP yield curve into the supplied cache */
static DRBool initialiseMarketDataCache1(DRObject cache){
    // create GBP Yield Curve and put in cache
    static const char method[] = "initialiseMarketDataCache1";
    DRValue csCurve = mkCashSwapCurve1();
    FLUNK_IF(g_svc->fptr->mapAddItem(g_svc, cache, "", &csCurve));
    CHECK(g_svc->fptr->valueClear(g_svc, &csCurve));
    return DR_TRUE;
done:
    return 0;
}

/* Just inserts a GBP yield curve (with slightly different rates than
   initialiseMarketDataCache1) into the supplied cache */
static DRBool initialiseMarketDataCache2(DRObject cache){
    // create GBP Yield Curve and put in cache
    static const char method[] = "initialiseMarketDataCache2";
    FLUNK_IF(addObj(
                 cache, "CashSwapCurve::Interface",
                 "ccy", mkStr("GBP"),
                 "name", mkStr(GBP_LIBOR),
                 "today", mkDt(todayDate, todayTime),
                 "spotOffset", mkInt(0),
                 "moneyMarketDenom", mkInt(365),
                 "swapFrequency", mkInt(2),
                 "swapDayCount", mkStr("Act/360"),
                 "hols", mkObj("Holiday", 
                               "name", mkStr("noweekends"),
                               "useWeekends", mkBl(DR_FALSE),
                               "holidays", mkObjArray("DateTime", 0),
                               0),
                 "expiries", mkObjArray("Expiry", 13, 
                                        mkMat("ON"), mkMat("1M"), 
                                        mkMat("2M"), mkMat("3M"),  
                                        mkMat("6M"), mkMat("1Y"), 
                                        mkMat("2Y"), mkMat("3Y"),
                                        mkMat("5Y"), mkMat("7Y"), 
                                        mkMat("10Y"), mkMat("20Y"), 
                                        mkMat("30Y"), 0),
                 "instruments", mkStrArray(13, "M", "M", "M", "M", "M", "M",
                                           "S", "S", "S", "S", "S", "S", "S"),
                 "rates", mkDbArray(13, 0.05, 0.051, 0.052, 0.053, 0.054,
                                    0.055, 0.056,  0.057, 0.058, 0.059, 0.06,
                                    0.061, 0.062),
                 "interp", mkObj("FourPlusI", 0),
                 0));
    return DR_TRUE;
done:
    return 0;
}

/* populates the cache with the rest of the remaining data */
static DRBool initialiseMarketDataCache3(DRObject cache) {
    // create USD Yield Curve and put in cache
    static const char method[] = "initialiseMarketDataCache3";
    FLUNK_IF(addObj(
                 cache, "CashSwapCurve::Interface",
                 "ccy", mkStr("USD"),
                 "name", mkStr(USD_LIBOR),
                 "today", mkDt(todayDate, todayTime),
                 "spotOffset", mkInt(0),
                 "moneyMarketDenom", mkInt(365),
                 "swapFrequency", mkInt(2),
                 "swapDayCount", mkStr("Act/360"),
                 "hols", mkObj("Holiday", 
                               "name", mkStr("noweekends"),
                               "useWeekends", mkBl(DR_FALSE),
                               "holidays", mkObjArray("DateTime", 0),
                               0),
                 "expiries", mkObjArray("Expiry", 13,
                                        mkMat("ON"), mkMat("1M"), 
                                        mkMat("2M"), mkMat("3M"),  
                                        mkMat("6M"), mkMat("1Y"), 
                                        mkMat("2Y"), mkMat("3Y"),
                                        mkMat("5Y"), mkMat("7Y"), 
                                        mkMat("10Y"), mkMat("20Y"), 
                                        mkMat("30Y"), 0),
                 "instruments", mkStrArray(13, "M", "M", "M", "M", "M", "M",
                                           "S", "S", "S", "S", "S", "S", "S"),
                 "rates", mkDbArray(13, 0.04, 0.041, 0.042, 0.043, 0.044, 0.045,
                                    0.046,  0.047, 0.048, 0.049, 0.05, 0.051,
                                    0.052),
                 "interp", mkObj("FourPlusI", 0),
                 0));

    // create USD Yield Curve and put in cache
    FLUNK_IF(addObj(
                 cache, "CashSwapCurve::Interface",
                 "ccy", mkStr("EUR"),
                 "name", mkStr(EUR_LIBOR),
                 "today", mkDt(todayDate, todayTime),
                 "spotOffset", mkInt(0),
                 "moneyMarketDenom", mkInt(365),
                 "swapFrequency", mkInt(2),
                 "swapDayCount", mkStr("Act/360"),
                 "hols", mkObj("Holiday", 
                               "name", mkStr("noweekends"),
                               "useWeekends", mkBl(DR_FALSE),
                               "holidays", mkObjArray("DateTime", 0),
                               0),
                 "expiries", mkObjArray("Expiry", 13,
                                        mkMat("ON"), mkMat("1M"), 
                                        mkMat("2M"), mkMat("3M"),  
                                        mkMat("6M"), mkMat("1Y"), 
                                        mkMat("2Y"), mkMat("3Y"),
                                        mkMat("5Y"), mkMat("7Y"), 
                                        mkMat("10Y"), mkMat("20Y"), 
                                        mkMat("30Y"), 0),
                 "instruments", mkStrArray(13, "M", "M", "M", "M", "M", "M",
                                           "S", "S", "S", "S", "S", "S", "S"),
                 "rates", mkDbArray(13, 0.03, 0.031, 0.032, 0.033, 0.034, 0.035,
                                    0.036,  0.037, 0.038, 0.039, 0.04, 0.041,
                                    0.042),
                 "interp", mkObj("FourPlusI", 0),
                 0));

    FLUNK_IF(addObj(
                 cache, "VolSurface", 
                 "name", mkStr(FTSE_VOL),
                 "metric", mkObj("TimeMetric",
                                 "nonTradTimeFrac", mkDb(1.0),
                                 "marketHols", 
                                 mkObj("Holiday",
                                       "name", mkStr("noweekends"),
                                       "useWeekends", mkBl(DR_FALSE),
                                       "holidays", mkObjArray("DateTime",
                                                              0),
                                       0),
                                 0),
                 "strikes", mkDbArray(3, 5000.0, 6000.0, 7000.0),
                 "vol", mkDbMat(3, 4, 
                                0.5, 0.55, 0.6, 0.65,
                                0.4, 0.45, 0.5, 0.55,
                                0.3, 0.35, 0.4, 0.45),
                 "expiries", mkObjArray("Expiry", 4, 
                                        mkMat("1W"), mkMat("1M"), 
                                        mkMat("1Y"), mkMat("2Y"),  
                                        0),
                 "baseDate", mkDt(todayDate, todayTime),
                 0));

    FLUNK_IF(addObj(
                 cache, "VolSurface", 
                 "name", mkStr(DAX_VOL),
                 "metric", mkObj("TimeMetric",
                                 "nonTradTimeFrac", mkDb(1.0),
                                 "marketHols", 
                                 mkObj("Holiday",
                                       "name", mkStr("noweekends"),
                                       "useWeekends", mkBl(DR_FALSE),
                                       "holidays", mkObjArray("DateTime",
                                                              0),
                                       0),
                                 0),
                 "strikes", mkDbArray(3, 3000.0, 5000.0, 6000.0),
                 "vol", mkDbMat(3, 4, 
                                0.2, 0.25, 0.3, 0.35,
                                0.25, 0.35, 0.4, 0.45,
                                0.3, 0.35, 0.4, 0.45),
                 "expiries", mkObjArray("Expiry", 4, 
                                        mkMat("1W"), mkMat("1M"), 
                                        mkMat("1Y"), mkMat("2Y"),  
                                        0),
                 "baseDate", mkDt(todayDate, todayTime),
                 0));

    FLUNK_IF(addObj(
                 cache, "VolSurface", 
                 "name", mkStr(SPX_VOL),
                 "metric", mkObj("TimeMetric",
                                 "nonTradTimeFrac", mkDb(1.0),
                                 "marketHols", 
                                 mkObj("Holiday",
                                       "name", mkStr("noweekends"),
                                       "useWeekends", mkBl(DR_FALSE),
                                       "holidays", mkObjArray("DateTime",
                                                              0),
                                       0),
                                 0),
                 "strikes", mkDbArray(3, 2000.0, 4000.0, 7000.0),
                 "vol", mkDbMat(3, 4, 
                                0.5, 0.55, 0.6, 0.65,
                                0.4, 0.45, 0.5, 0.55,
                                0.3, 0.35, 0.4, 0.45),
                 "expiries", mkObjArray("Expiry", 4, 
                                        mkMat("1W"), mkMat("1M"), 
                                        mkMat("1Y"), mkMat("2Y"),  
                                        0),
                 "baseDate", mkDt(todayDate, todayTime),
                 0));

    FLUNK_IF(addObj(
                 cache, "SimpleEquity", 
                 "equity", mkObj(
                     "Equity",
                     "name", mkStr(FTSE),
                     "stockPrice", mkDb(6000.0),
                     "settlement", 
                     mkObj("RollingSettlement", 
                           "period", mkInt(3),
                           "hols", 
                           mkObj("Holiday",
                                 "name", mkStr("noweekends"),
                                 "useWeekends", 
                                 mkBl(DR_FALSE),
                                 "holidays",
                                 mkObjArray("DateTime", 0),
                                 0),
                           0),
                     "yc", mkStr(GBP_LIBOR),
                     "valueDate", mkDt(todayDate, todayTime),
                     "stockDate", mkDt(todayDate, todayTime),
                     "divList", 
                     mkObj("DividendList",
                           "divArray",
                           mkObjArray("Dividend", 3, 
                                      mkObj("Dividend",
                                            "exDivDate", 
                                            mkDt("01-Apr-2001", "SOD"),
                                            "divAmount", mkDb(90.0),
                                            "divType", mkInt(0),
                                            "payDivDate", 
                                            mkDt("15-Apr-2001", "SOD"),
                                            0),
                                      mkObj("Dividend",
                                            "exDivDate", 
                                            mkDt("01-Jul-2001",
                                                 "SOD"),
                                            "divAmount", mkDb(10.0),
                                            "divType", mkInt(0),
                                            "payDivDate", 
                                            mkDt("01-Jul-2001", "SOD"),
                                            0),
                                      mkObj("Dividend",
                                            "exDivDate", 
                                            mkDt("01-Oct-2001", "SOD"),
                                            "divAmount", mkDb(0.01),
                                            "divType", mkInt(1),
                                            "payDivDate", 
                                            mkDt("01-Oct-2001", "SOD"),
                                            0),
                                     
                                      0),
                           0),
                     "borrowCurve", 
                     mkObj("BorrowCurve::Interface",
                           "baseDate", mkDt(todayDate, todayTime),
                           "expiries", 
                           mkObjArray("Expiry", 4, 
                                      mkMat("1W"), mkMat("1M"), 
                                      mkMat("2M"), mkMat("1Y"),  
                                      0),
                           "rates", mkDbArray(4, 0.01, 0.02, 0.03, 0.03),
                           "basis", mkInt(5000),
                           "name", mkStr("FTSE-BORROW"),
                           "dayCount",  mkStr("Act/365F"),
                           0),
                     0),
                 "vol", mkStr(FTSE_VOL),
                 0));

    FLUNK_IF(addObj(
                 cache, "SimpleEquity", 
                 "equity", mkObj(
                     "Equity",
                     "name", mkStr(DAX),
                     "stockPrice", mkDb(5000.0),
                     "settlement", 
                     mkObj("RollingSettlement", 
                           "period", mkInt(5),
                           "hols", 
                           mkObj("Holiday",
                                 "name", mkStr("noweekends"),
                                 "useWeekends", 
                                 mkBl(DR_FALSE),
                                 "holidays",
                                 mkObjArray("DateTime", 
                                            0),
                                 0),
                           0),
                     "yc", mkStr(EUR_LIBOR),
                     "valueDate", mkDt(todayDate, todayTime),
                     "stockDate", mkDt(todayDate, todayTime),
                     "divList", 
                     mkObj("DividendList",
                           "divArray",
                           mkObjArray("Dividend", 0),
                           0),
                     "borrowCurve", 
                     mkObj("BorrowCurve::Interface",
                           "baseDate", mkDt(todayDate, todayTime),
                           "expiries", 
                           mkObjArray("Expiry", 4, 
                                      mkMat("1W"), mkMat("1M"), 
                                      mkMat("2M"), mkMat("1Y"),  
                                      0),
                           "rates", mkDbArray(4, 0.00, 0.00, 0.00, 0.00),
                           "basis", mkInt(5000),
                           "name", mkStr("DAX-BORROWCOSTS"),
                           "dayCount",  mkStr("Act/365F"),
                           0),
                     0),
                 "vol", mkStr(DAX_VOL),
                 0));

    FLUNK_IF(addObj(
                 cache, 
                 "SimpleEquity", 
                 "equity", mkObj(
                     "Equity",
                     "name", mkStr(SPX),
                     "stockPrice", mkDb(4000.0),
                     "settlement", 
                     mkObj("RollingSettlement", 
                           "period", mkInt(1),
                           "hols", 
                           mkObj("Holiday",
                                 "name", mkStr("noweekends"),
                                 "useWeekends", 
                                 mkBl(DR_FALSE),
                                 "holidays",
                                 mkObjArray("DateTime", 
                                            0),
                                 0),
                           0),
                     "yc", mkStr(USD_LIBOR),
                     "valueDate", mkDt(todayDate, todayTime),
                     "stockDate", mkDt(todayDate, todayTime),
                     "divList", 
                     mkObj("DividendList",
                           "divArray",
                           mkObjArray("Dividend", 1, 
                                      mkObj("Dividend",
                                            "exDivDate", 
                                            mkDt("01-Apr-2001", "SOD"),
                                            "divAmount", mkDb(0.01),
                                            "divType", mkInt(2),
                                            "payDivDate", 
                                            mkDt("01-Apr-2001", "SOD"),
                                            0),
                                      0),
                           0),
                     "borrowCurve", 
                     mkObj("BorrowCurve::Interface",
                           "baseDate", mkDt(todayDate, todayTime),
                           "expiries", 
                           mkObjArray("Expiry", 4, 
                                      mkMat("1W"), mkMat("1M"), 
                                      mkMat("2M"), mkMat("1Y"),  
                                      0),
                           "rates", mkDbArray(4, 0.01, 0.01, 0.01, 0.01),
                           "basis", mkInt(5000),
                           "name", mkStr(SPX),
                           "dayCount",  mkStr("Act/365F"),
                           0),
                     0),
                 "vol", mkStr(SPX_VOL),
                 0));
    
    FLUNK_IF(addObj(
                 cache, "FlatFXVol", 
                 "name", mkStr(FX_USD_GBP_VOL),
                 "flatFXVol", mkDb(0.3),
                 "baseDate", mkDt(todayDate, todayTime),
                 "timeMetric", 
                 mkObj("TimeMetric",
                       "nonTradTimeFrac", mkDb(1.0),
                       "marketHols", mkObj("Holiday",
                                           "name", mkStr("noweekends"),
                                           "useWeekends", mkBl(DR_FALSE),
                                           "holidays", mkObjArray("DateTime", 0),
                                           0),
                       0),
                 0));

    FLUNK_IF(addObj(
                 cache, "FXAsset", 
                 "name", mkStr(FX_GBP_USD),
                 "today", mkDt(todayDate, todayTime),
                 "holidays", mkObj("Holiday",
                                   "name", mkStr("noweekends"),
                                   "useWeekends", mkBl(DR_FALSE),
                                   "holidays", mkObjArray("DateTime", 0),
                                   0),
                 "riskCcy", mkStr(GBP_LIBOR),
                 "baseCcy", mkStr(USD_LIBOR),
                 "spotFX", mkDb(1.42),
                 "fxVol", mkStr(FX_USD_GBP_VOL),
                 0));

    FLUNK_IF(addObj(
                 cache, "FlatFXVol", 
                 "name", mkStr(FX_USD_EUR_VOL),
                 "flatFXVol", mkDb(0.2),
                 "baseDate", mkDt(todayDate, todayTime),
                 "timeMetric", 
                 mkObj("TimeMetric",
                       "nonTradTimeFrac", mkDb(1.0),
                       "marketHols", mkObj("Holiday",
                                           "name", mkStr("noweekends"),
                                           "useWeekends", mkBl(DR_FALSE),
                                           "holidays", mkObjArray("DateTime", 0),
                                           0),
                       0),
                 0));

    FLUNK_IF(addObj(
                 cache, "FXAsset", 
                 "name", mkStr(FX_EUR_USD),
                 "today", mkDt(todayDate, todayTime),
                 "holidays", mkObj("Holiday",
                                   "name", mkStr("noweekends"),
                                   "useWeekends", mkBl(DR_FALSE),
                                   "holidays", mkObjArray("DateTime", 0),
                                   0),
                 "riskCcy", mkStr(EUR_LIBOR),
                 "baseCcy", mkStr(USD_LIBOR),
                 "spotFX", mkDb(1.1),
                 "fxVol", mkStr(FX_USD_EUR_VOL),
                 0));

    FLUNK_IF(addObj(cache, "Correlation", 
                    "name", mkStr("Corr-FTSE_USD"),
                    "asset1", mkStr(FTSE),
                    "asset2", mkStr(FX_GBP_USD),
                    "correlation", mkDb(0.3),
                    0));

    FLUNK_IF(addObj(cache, "Correlation", 
                    "name", mkStr("Corr-DAX_USD"),
                    "asset1", mkStr(DAX),
                    "asset2", mkStr(FX_EUR_USD),
                    "correlation", mkDb(0.3),
                    0));

    FLUNK_IF(addObj(cache, "Correlation", 
                    "name", mkStr("Corr-FTSE-DAX"),
                    "asset1", mkStr(FTSE),
                    "asset2", mkStr(DAX),
                    "correlation", mkDb(0.1),
                    0));

    FLUNK_IF(addObj(cache, "Correlation", 
                    "name", mkStr("Corr-FTSE-SPX"),
                    "asset1", mkStr(FTSE),
                    "asset2", mkStr(SPX),
                    "correlation", mkDb(0.2),
                    0));

    FLUNK_IF(addObj(cache, "Correlation", 
                    "name", mkStr("Corr-DAX-SPX"),
                    "asset1", mkStr(DAX),
                    "asset2", mkStr(SPX),
                    "correlation", mkDb(0.3),
                    0));

    FLUNK_IF(addObj(
                 cache, "XCB", 
                 "name", mkStr("INDEX_TRIO"),
                 "assets", mkStrArray(3, FTSE, DAX, SPX),
                 "ccyTreatments", mkStrArray(3, "S", "P", "N"),
                 "basketYCName", mkStr(USD_LIBOR),
                 "unitWeights", mkBl(DR_FALSE),
                 "pubWeights", mkDbArray(3, 1.0/3, 1.0/3, 1.0/3),
                 "marketHols", mkObj("Holiday",
                                     "name", mkStr("noweekends"),
                                     "useWeekends", mkBl(DR_FALSE),
                                     "holidays", mkObjArray("DateTime", 0),
                                     0),
                 "startDate", mkDt("01-Jan-2001", todayTime),
                 "spotsAtStart", mkDbArray(3, 6500.0, 5500.0, 4500.0),
                 "smileType", mkStr("E"),
                 "timeMetric", mkObj("TimeMetric",
                                     "nonTradTimeFrac", mkDb(1.0),
                                     "marketHols", 
                                     mkObj("Holiday",
                                           "name", mkStr("noweekends"),
                                           "useWeekends", mkBl(DR_FALSE),
                                           "holidays", mkObjArray("DateTime", 0),
                                           0),
                                     0),
                 0));

    return DR_TRUE;
done:
    return 0;
}

/* Prints out the strings making up an EDROutputName */
static DRBool displayName(DRValue* outputName, char* buffer){
    static const char method[] = "displayName";
    DRBool     ok = DR_FALSE;
    DRMap      map = 0;
    DRValue    names;
    DRValue    id;
    int        numIDs = 0;
    int        i;

    FLUNK_IF(g_svc->fptr->objectToMap(g_svc, outputName->value.object, &map));
    FLUNK_IF(g_svc->fptr->mapGetItem(g_svc, map, "names", &names));
    FLUNK_IF(g_svc->fptr->arrayLength(g_svc, names.value.object, &numIDs));

    for (i = 0; i < numIDs; i++){
        const char* nm;
        FLUNK_IF(g_svc->fptr->arrayGet(g_svc, names.value.object, i, &id));
        nm = id.value.string;
        if (*nm != '\0'){
            char nameBuffer[200];
            sprintf(nameBuffer, "%s%s", (i > 0? ".": "_"), nm);
            strcat(buffer, nameBuffer);
        }
        CHECK(g_svc->fptr->valueClear(g_svc, &id));
    }
    ok = DR_TRUE;
done:
    DRI_FREE(g_svc, map);
    CHECK(g_svc->fptr->valueClear(g_svc, &names));
    CHECK(g_svc->fptr->valueClear(g_svc, &id));
    return ok;
}

static DRBool displayObject(DRValue* object, char* buffer, FILE* file);

static void top(FILE* file) {
    fprintf(file, "<OutputFile>\n");
    fprintf(file, "<Output>\n<Summary>\nPermutation 1:\n");
}
    
static void tail(FILE* file) {
    fprintf(file, "</Summary>\n</Output>\n");
    fprintf(file, "</OutputFile>\n");
}

/** list all available types */
static void typeList(const char* filename) {
    static const char method[] = "typeList";
    DRValue  x;
    DRArray  types;
    DRString elemtype = 0;
    FILE*    file = 0;
    DRError  err = 0;

    if (filename) {
        file = fopen(filename, "w");
    }

    if ((err = g_svc->fptr->typeList(g_svc, &types))) {
        handleError("typeList", err);
        DRI_ERROR_FREE(g_svc, err);
    }
    else {
        CHECK(EDGSDRValue(&x, DR_OBJECT, types));
        if ((err = g_svc->fptr->arrayElementType(g_svc, types, &elemtype))) {
            DRI_ERROR_FREE(g_svc, err);
            handleError("typeList", err);
        }
        else {
            printf("Type list is an array of: %s\n", elemtype);
        }
        printf("All available types in EDG:\n");
        top(file ? file : stdout);
        displayObject(&x, "", file ? file : stdout);
        tail(file ? file : stdout);
    }

    DRI_STRING_FREE(g_svc, elemtype);
    DRI_FREE(g_svc, types);

    if (file) {
        fclose(file);
    }
}

/** Display the results returned from the pricing call */
static DRBool displayResults(DRObject results, char* buffer, FILE* file){
    static const char method[] = "displayResults";
    DRBool        ok = DR_FALSE;
    DRString      resultsType = 0;
    DRString      objType = 0;
    DRMap         map = 0;
    DRValue       packet;
    DRValue       id;
    DRValue       value;
    int           numResults;
    int           i;

    if (!results){
        ERROR_EXIT("null results cannot be displayed");
    }

    FLUNK_IF(g_svc->fptr->objectGetType(g_svc, results, &resultsType));

    if (strcmp(resultsType, "Results")){
        char tmp[100];
        sprintf(tmp, "Can't cope with result object of type %s\n", 
                resultsType);
        ERROR_EXIT(tmp);
    }

    FLUNK_IF(g_svc->fptr->objectToMap(g_svc, results, &map));
    
    /** pull out packet, id and value arrays */
    FLUNK_IF(g_svc->fptr->mapGetItem(g_svc, map, "packet", &packet));

    FLUNK_IF(g_svc->fptr->mapGetItem(g_svc, map, "identifier", &id));

    FLUNK_IF(g_svc->fptr->mapGetItem(g_svc, map, "value", &value));

    FLUNK_IF(g_svc->fptr->arrayLength(g_svc, packet.value.object,
                                      &numResults));

    for (i = 0; i < numResults; i++) {
        DRValue packetElem;
        DRValue idElem;
        DRValue valueElem;
        char localBuffer[500];
        strcpy(localBuffer, buffer);
        if (*localBuffer != '\0'){
            strcat(localBuffer, "_");
        }

        FLUNK_IF(g_svc->fptr->arrayGet(g_svc, packet.value.object, i,
                                       &packetElem));

        FLUNK_IF(g_svc->fptr->arrayGet(g_svc,
                                       id.value.object,
                                       i,
                                       &idElem));

        FLUNK_IF(g_svc->fptr->arrayGet(g_svc, value.value.object, i,
                                       &valueElem));

        strcat(localBuffer,
               packetElem.value.string);
        /* print out name  */
        if (!displayName(&idElem, localBuffer)){
            ERROR_EXIT("name cannot be displayed");
        }

        /* then display object (switch on type) */
        if (!displayObject(&valueElem, localBuffer, file)){
            ERROR_EXIT("object cannot be displayed");
        }

        CHECK(g_svc->fptr->valueClear(g_svc, &packetElem));  
        CHECK(g_svc->fptr->valueClear(g_svc, &idElem));  
        CHECK(g_svc->fptr->valueClear(g_svc, &valueElem));  
    }
    ok = DR_TRUE;
done:
    CHECK(g_svc->fptr->valueClear(g_svc, &packet));  
    CHECK(g_svc->fptr->valueClear(g_svc, &id));  
    CHECK(g_svc->fptr->valueClear(g_svc, &value));  
    DRI_STRING_FREE(g_svc, resultsType);
    DRI_STRING_FREE(g_svc, objType);
    DRI_FREE(g_svc, map);
    return ok;
}
                                 

/* Recurse through an object printing its components */
static DRBool displayObject(DRValue* object, char* buffer, FILE* file){
    static const char method[] = "displayObject";
    DRBool        ok = DR_FALSE;
    DRString      objType = 0;
    DRBool        isArray;
    DRBool        isMatrix;
    DRMapIterator iter = 0;

    if (object->type == DR_OBJECT) {
        /* find out what sort of object the result is */
        FLUNK_IF(g_svc->fptr->objectIsArray(g_svc, object->value.object,
                                            &isArray));
        FLUNK_IF(g_svc->fptr->objectIsMatrix(g_svc, object->value.object,
                                             &isMatrix));

        DRI_ERROR_FREE(g_svc, 
                       g_svc->fptr->mapIteratorGet(g_svc, object->value.object,
                                                   &iter));

        FLUNK_IF(g_svc->fptr->objectGetType(g_svc, object->value.object,
                                            &objType));
    }
    
    /* display result (switch on type) */
    if (object->type == DR_DOUBLE){
        char tmpBuffer[500];
        sprintf(tmpBuffer, ": %f\n", object->value.real);
        strcat(buffer, tmpBuffer);
        fprintf(file, "Output %s", buffer);
    }
    else if (object->type == DR_INT){
        char tmpBuffer[500];
        sprintf(tmpBuffer, ": %d\n", object->value.integer);
        strcat(buffer, tmpBuffer);
        fprintf(file, "Output %s", buffer);
    }
    else if (object->type == DR_BOOL){
        char tmpBuffer[500];
        sprintf(tmpBuffer, ": %s\n", object->value.boolean ? "Y" : "N");
        strcat(buffer, tmpBuffer);
        fprintf(file, "Output %s", buffer);
    }
    else if (object->type == DR_STRING){
        fprintf(file, "Output %s: %s\n", buffer, object->value.string);
    }
    else if (objType && strcmp(objType, "NotApplicable") == 0){
        fprintf(file, "Output %s: NotApplicable\n", buffer);
    } 
    else if (objType && strcmp(objType, "Untweakable") == 0){
        fprintf(file, "Output %s: Untweakable\n", buffer);
    } 
    else if (objType && strcmp(objType, "Results") == 0){
        char localBuffer[500];
        strcpy(localBuffer, buffer);
        displayResults(object->value.object, localBuffer, file);
    }
    else if (isArray){
        int i, length;
        FLUNK_IF(g_svc->fptr->arrayLength(g_svc, object->value.object,
                                          &length));
        for (i = 0; i < length; i++){
            DRValue elt;
            FLUNK_IF(g_svc->fptr->arrayGet(g_svc, object->value.object,
                                           i, &elt));

            char tmpBuffer[10];
            sprintf(tmpBuffer, "[%d]", i);
            char localBuffer[500];
            strcpy(localBuffer, buffer);
            strcat(localBuffer, tmpBuffer);
            if (!displayObject(&elt, localBuffer, file)){
                ERROR_EXIT("failed to display object");
            }
            CHECK(g_svc->fptr->valueClear(g_svc, &elt));
        }
    } 
    else if (isMatrix) {
        int numCol;
        int numRow;
        int x;
        int y;
        double val;
        char tmpBuffer[500];

        FLUNK_IF(g_svc->fptr->matrixSize(g_svc, object->value.object,
                                         &numCol, &numRow));

        for (x = 0; x < numCol; x++) {
            fprintf(file, "Output %s [%d]:", buffer, x);
            for (y = 0; y < numRow; y++) {
                FLUNK_IF(g_svc->fptr->matrixGet(
                             g_svc, object->value.object,x,y,&val));

                sprintf(tmpBuffer,
                        "%f\t",
                        val);
                fprintf(file, "%s", tmpBuffer);
            }
            fprintf(file, "\n");
        }
        fprintf(file, "\n");
    }
    else if (iter) {
        DRString      field = 0;
        DRValue       obj;
        DRError       err = 0;
 
        while (!(err = g_svc->fptr->mapIteratorNext(g_svc,
                                                    iter,
                                                    &field,
                                                    &obj)) && field) {
            char localBuffer[500];
            strcpy(localBuffer, buffer);
            if (*localBuffer != '\0'){
                strcat(localBuffer, "_");
            }
            strcat(localBuffer, field);
            displayObject(&obj, localBuffer, file);
            DRI_STRING_FREE(g_svc, field);
            CHECK(g_svc->fptr->valueClear(g_svc, &obj));
        }

        if (err){
            DRI_ERROR_FREE(g_svc, err);
            goto done;
        }
        DRI_FREE(g_svc, iter);
    }
    else {
        DRValue v;
        DRMap   dd = 0;
        FLUNK_IF(g_svc->fptr->objectToMap(g_svc, object->value.object, &dd));
        v.type = DR_OBJECT;
        v.value.object = dd;
        if (!dd || !displayObject(&v, buffer, file)){
            ERROR_EXIT("Either DataDictionary is null or object \
cannot be displayed.");
        }
        DRI_FREE(g_svc, dd);
    }

    ok = DR_TRUE;
  done:
    DRI_STRING_FREE(g_svc, objType)
    return ok;
}


/* build a CashFlowStream instrument */
static DRValue cfinst(void) {
    DRValue inst;
    inst = mkObj("CashFlowStream", 
                 "cfl", mkObjArray("CashFlow", 1, 
                                   mkObj("CashFlow",
                                         "date", mkDt("15-Sep-2003", "EOD"),
                                         "amount", mkDb(100.0),
                                         0),
                                   0),
                 "discount", mkStr(GBP_LIBOR),
                 0);
    if (inst.type == DR_UNDEFINED){
        printf("cfinst failed\n");
    }
    return inst;
}

/* build a CashFlowStream DR Wrapper */
static DRValue cfwrapper(void) {
    DRValue inst;
    inst = mkWrapper("CashFlowStream", 
                     "cfl", mkObjArray("CashFlow", 1, 
                                       mkObj("CashFlow",
                                             "date",mkDt("15-Sep-2003","EOD"),
                                             "amount", mkDb(100.0),
                                             0),
                                       0),
                     "discount", mkStr(GBP_LIBOR),
                     0);
    if (inst.type == DR_UNDEFINED){
        printf("cfwrapper failed\n");
    }
    return inst;
}


/* build an average spot instument - vanilla ccy treatment */
static DRValue avgSpotInst1(void) {
    DRValue inst;
    inst = mkObj("AverageSpot", 
                 "isCall", mkBl(DR_TRUE),
                 "maturity", mkDt("15-Apr-2002", "EOD"),
                 "strike", mkDb(6000),
                 "avgOut", mkObj("SampleList",
                                 "dates", 
                                 mkObjArray("DateTime", 1, 
                                            mkDt("15-Apr-2002", "EOD"),
                                            0),
                                 "values", mkDbArray(1, 0.0),
                                 "weights", mkDbArray(1, 1.0),
                                 0),
                 "fwdStarting", mkBl(DR_FALSE),
                 "startDate", mkDt(todayDate, todayTime),
                 "oneContract", mkBl(DR_TRUE),
                 /* for one contract can skip spot at start and notional */
                 "instSettle", mkObj("CashSettleDate", 
                                     "settleDate", mkDt("15-Jun-2002", "EOD"),
                                     0),
                 "premiumSettle", mkObj("CashSettleDate", 
                                        "settleDate", mkDt(todayDate, todayTime),
                                     0),
                 "asset", mkStr(FTSE),
                 "ccyTreatment", mkStr("N"),
                 "discount", mkStr(GBP_LIBOR),
                 0);
    if (inst.type == DR_UNDEFINED){
        printf("avgSpotInst1 failed\n");
    }
    return inst;
}

/* build an average spot instument - protected ccy treatment (FTSE into $)*/
static DRValue avgSpotInst2(void) {
    DRValue inst;
    inst = mkObj("AverageSpot", 
                 "isCall", mkBl(DR_TRUE),
                 "maturity", mkDt("15-Apr-2002", "EOD"),
                 "strike", mkDb(6000),
                 "avgOut", mkObj("SampleList",
                                 "dates", 
                                 mkObjArray("DateTime", 1, 
                                            mkDt("15-Apr-2002", "EOD"),
                                            0),
                                 "values", mkDbArray(1, 0.0),
                                 "weights", mkDbArray(1, 1.0),
                                 0),
                 "fwdStarting", mkBl(DR_FALSE),
                 "startDate", mkDt(todayDate, todayTime),
                 "oneContract", mkBl(DR_TRUE),
                 /* for one contract can skip spot at start and notional */
                 "instSettle", mkObj("CashSettleDate", 
                                     "settleDate", mkDt("15-Jun-2002", "EOD"),
                                     0),
                 "premiumSettle", mkObj("CashSettleDate", 
                                        "settleDate", mkDt(todayDate, todayTime),
                                     0),
                 "asset", mkStr(FTSE),
                 "ccyTreatment", mkStr("P"),
                 "discount", mkStr(USD_LIBOR),
                 0);
    if (inst.type == DR_UNDEFINED){
        printf("avgSpotInst2 failed\n");
    }
    return inst;
}

/* build an average spot instument - struck ccy treatment (FTSE into $)*/
static DRValue avgSpotInst3(void) {
    DRValue inst;
    inst = mkObj("AverageSpot", 
                 "isCall", mkBl(DR_TRUE),
                 "maturity", mkDt("15-Apr-2002", "EOD"),
                 "strike", mkDb(6000),
                 "avgOut", mkObj("SampleList",
                                 "dates", 
                                 mkObjArray("DateTime", 1, 
                                            mkDt("15-Apr-2002", "EOD"),
                                            0),
                                 "values", mkDbArray(1, 0.0),
                                 "weights", mkDbArray(1, 1.0),
                                 0),
                 "fwdStarting", mkBl(DR_FALSE),
                 "startDate", mkDt(todayDate, todayTime),
                 "oneContract", mkBl(DR_TRUE),
                 /* for one contract can skip spot at start and notional */
                 "instSettle", mkObj("CashSettleDate", 
                                     "settleDate", mkDt("15-Jun-2002", "EOD"),
                                     0),
                 "premiumSettle", mkObj("CashSettleDate", 
                                        "settleDate", mkDt(todayDate, todayTime),
                                     0),
                 "asset", mkStr(FTSE),
                 "ccyTreatment", mkStr("S"),
                 "discount", mkStr(USD_LIBOR),
                 0);
    if (inst.type == DR_UNDEFINED){
        printf("avgSpotInst3 failed\n");
    }
    return inst;
}

/* build an average spot instument - on an XCB ($ denominated basket) */
static DRValue avgSpotInst4(void) {
    DRValue inst;
    inst = mkObj("AverageSpot", 
                 "isCall", mkBl(DR_TRUE),
                 "maturity", mkDt("15-Apr-2002", "EOD"),
                 "strike", mkDb(100),
                 "avgOut", mkObj("SampleList",
                                 "dates", 
                                 mkObjArray("DateTime", 1, 
                                            mkDt("15-Apr-2002", "EOD"),
                                            0),
                                 "values", mkDbArray(1, 0.0),
                                 "weights", mkDbArray(1, 1.0),
                                 0),
                 "fwdStarting", mkBl(DR_FALSE),
                 "startDate", mkDt(todayDate, todayTime),
                 "oneContract", mkBl(DR_TRUE),
                 /* for one contract can skip spot at start and notional */
                 "instSettle", mkObj("CashSettleDate", 
                                     "settleDate", mkDt("15-Jun-2002", "EOD"),
                                     0),
                 "premiumSettle", mkObj("CashSettleDate", 
                                        "settleDate", mkDt(todayDate, todayTime),
                                     0),
                 "asset", mkStr("INDEX_TRIO"),
                 "ccyTreatment", mkStr("N"),
                 "discount", mkStr(USD_LIBOR),
                 0);
    if (inst.type == DR_UNDEFINED){
        printf("avgSpotInst4 failed\n");
    }
    return inst;
}

/* build a set of MU_S bucket dates */
static DRValue mudates(void) {
    static const char method[] = "mudates";
    DRValue action = mkObj("MuDateBuilder",
                             "today", mkDt(todayDate, todayTime),
                             "mutype", mkStr("SPECIAL"),
                             0);

    DRValue mudates;
    CHECK(g_svc->fptr->execute(g_svc, action.value.object, &mudates));
    CHECK(g_svc->fptr->valueClear(g_svc, &action));
    return mudates;
}


/* build a Control */
static DRValue control(void) {
    DRValue  ctrl;
    printf("Building Control\n");
    ctrl = mkObj("Control", 
                 /* dump to file here only refers to single instrument without
                    scenarios etc */
                 "writeToFile", mkBl(DR_FALSE),
                 "outputRequests",
                 /* OutputRequest's ask for additional outputs that aren't
                    sensitivities */
                 mkObjArray("OutputRequest", 1, 
                            mkObj("OutputRequest", 
                                  "requestName", mkStr("FWD_AT_MAT"),
                                  0),
                            0),
                 "sens", mkObjArray("Sensitivity", 18, 
                                    mkObj("RhoBorrowPointwise",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("RhoBorrowParallel",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("FXPhi",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("RhoPointwise",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("RhoParallel",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("VegaSkewParallel",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("RootTimeVega",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("VegaMatrix",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("VegaPointwise",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("VegaParallel",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("Phi",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("VegaSkewPointwise",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("FXVega",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("SpotPrice",
                                          0),
                                    mkObj("MuParallel",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("MuSpecial",
                                          "shiftSize", mkDb(0.0001),
                                          "expiries", mudates(),
                                          0),
                                    mkObj("Delta",
                                          "shiftSize", mkDb(0.0001),
                                          0),
                                    mkObj("CrossGamma",
                                          "shiftSize", mkDb(0.0001),
                                          0),
#if 0 /* add in if you want theta */
                                    mkObj("Theta",
                                          "offset", mkInt(1),
                                          "hols", 
                                          mkObj("Holiday", 
                                                "name", mkStr("noweekends"),
                                                "useWeekends", mkBl(DR_FALSE),
                                                "holidays", 
                                                mkObjArray("DateTime", 0),
                                                0)),

#endif
                                    0),
                 0);
    if (ctrl.type == DR_UNDEFINED){
        printf("control failed\n");
    }
    return ctrl;
}

// build a CompositeInstrument
static DRObject composite(DRService* svc, 
                          int        numInsts,
                          DRValue*   model,
                          DRValue*   inst,
                          DRValue*   ctrl,
                          double*    weight,
                          double*    multiplier,
                          DRValue    scenario,
                          DRObject   mkt) {
    static const char method[] = "composite";
    DRBool     ok = DR_FALSE;
    DRMap      map = 0;
    DRObject   compo = 0;
    DRValue    models;
    DRValue    insts;
    DRValue    ctrls;
    DRValue    m;
    DRValue    w;
    DRValue    mktV;

    DRValue    weights;
    DRValue    mults;
    int        i;

    mktV.type = DR_OBJECT;
    mktV.value.object = mkt;

    printf("Building CompositeInstrument\n");

    models.type = DR_OBJECT;
    FLUNK_IF(svc->fptr->arrayNew(svc, numInsts, DRI_TYPE_VARIANT,
                                 &models.value.object));

    insts.type = DR_OBJECT;
    FLUNK_IF(svc->fptr->arrayNew(svc, numInsts, DRI_TYPE_VARIANT,
                                 &insts.value.object));

    ctrls.type = DR_OBJECT;
    FLUNK_IF(svc->fptr->arrayNew(svc, numInsts, DRI_TYPE_VARIANT,
                                 &ctrls.value.object));

    weights.type = DR_OBJECT;
    FLUNK_IF(svc->fptr->arrayNew(svc, numInsts, DRI_TYPE_DOUBLE,
                                 &weights.value.object));

    mults.type = DR_OBJECT;
    FLUNK_IF(svc->fptr->arrayNew(svc, numInsts, DRI_TYPE_DOUBLE,
                                 &mults.value.object));

    for (i = 0; i < numInsts; i++) {
        FLUNK_IF(svc->fptr->arraySet(svc, models.value.object, i, &model[i]));
        FLUNK_IF(svc->fptr->arraySet(svc, insts.value.object, i, &inst[i]));
        FLUNK_IF(svc->fptr->arraySet(svc, ctrls.value.object, i, &ctrl[i]));

        FLUNK_IF(EDGSDRValue(&w, DR_DOUBLE, weight[i]));
        FLUNK_IF(svc->fptr->arraySet(svc, weights.value.object, i, &w));

        FLUNK_IF(EDGSDRValue(&m, DR_DOUBLE, multiplier[i]));
        FLUNK_IF(svc->fptr->arraySet(svc, mults.value.object, i, &m));
    }
    
    FLUNK_IF(svc->fptr->mapNew(svc, "CompositeInstrument", &map));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "model", &models));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "inst", &insts));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "ctrl", &ctrls));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "multiplier", &mults));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "weight", &weights));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "scenario", &scenario));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "market", &mktV));
    FLUNK_IF(svc->fptr->mapToObject(svc, map, &compo));

    ok = DR_TRUE;
done:
    CHECK(svc->fptr->valueClear(g_svc, &models));
    CHECK(svc->fptr->valueClear(g_svc, &insts));
    CHECK(svc->fptr->valueClear(g_svc, &ctrls));
    CHECK(svc->fptr->valueClear(g_svc, &mults));
    CHECK(svc->fptr->valueClear(g_svc, &weights));
    DRI_FREE(svc, map);
    return compo;
}

/* demo building an object with a matrix, demo pulling them both apart */
static int mapMaker(const char* filename) {
    static const char method[] = "mapMaker";
    DRBool   ok = DR_FALSE;
    DRString mapType = 0;
    DRMap    dd = 0;
    DRValue  elem;

    char* getme = "name";
    FILE* file = 0;
    if (filename) {
        file = fopen(filename, "w");
    }

    DRValue obj = mkObj("VolSurface", 
                        "name", mkStr(FTSE_VOL),
                        "metric", mkObj("TimeMetric",
                                        "nonTradTimeFrac", mkDb(1.0),
                                        "marketHols", 
                                        mkObj("Holiday",
                                              "name", mkStr("noweekends"),
                                              "useWeekends", mkBl(DR_FALSE),
                                              "holidays", mkObjArray("DateTime",
                                                                     0),
                                              0),
                                        0),
                        "strikes", mkDbArray(3, 5000.0, 6000.0, 7000.0),
                        "vol", mkDbMat(3, 4, 
                                       0.5, 0.55, 0.6, 0.65,
                                       0.4, 0.45, 0.5, 0.55,
                                       0.3, 0.35, 0.4, 0.45),
                        "expiries", mkObjArray("Expiry", 4, 
                                               mkMat("1W"), mkMat("1M"), 
                                               mkMat("1Y"), mkMat("2Y"),  
                                               0),
                        "baseDate", mkDt(todayDate, todayTime),
                        0);

    /* flip this back to a map */
    FLUNK_IF(g_svc->fptr->objectToMap(g_svc, obj.value.object, &dd));

    FLUNK_IF(g_svc->fptr->mapGetTypeName(g_svc, dd, &mapType));
    
    printf("Map is of type: %s\n", mapType);

    /* pull an element out of the map */
    CHECK(EDGSDRValue(&elem, DR_UNDEFINED));
    FLUNK_IF(g_svc->fptr->mapGetItem(g_svc, dd, getme, &elem));

    top(file ? file : stdout);
    if (!displayObject(&elem, getme, file ? file : stdout)) {
        ERROR_EXIT("failed to display object");
    }

    /* display the whole thing (covers pulling matrices apart) */
    if (!displayObject(&obj, "VolSurface", file ? file : stdout)) {
        ERROR_EXIT("failed to display VolSurface object");
    }
    tail(file ? file : stdout);
    ok = DR_TRUE;
done:
    if (file) {
        fclose(file);
    }
    DRI_FREE(g_svc, dd);
    CHECK(g_svc->fptr->valueClear(g_svc, &obj));
    CHECK(g_svc->fptr->valueClear(g_svc, &elem));
    DRI_STRING_FREE(g_svc, mapType);
    return ok;
}

/* demo building an object with enums in it, demo pulling it apart */
static int mapMaker3(const char* filename) {
    static const char method[] = "mapMaker3";
    DRBool   ok = DR_FALSE;
    DRString mapType = 0;
    DRMap    dd = 0;
    DRValue  elem;

    char* getme = "myEnum1";
    FILE* file = 0;
    if (filename) {
        file = fopen(filename, "w");
    }

    DRValue obj = mkObj("EnumTester", 
                        /* enums go through as string for DRI */
                        "myEnum1", mkStr("FIRST_VALUE"),
                        "myEnum2", mkStr("THIRD_VALUE"),
                        "myEnums", mkStrArray(3,
                                              "THIRD_VALUE",
                                              "SECOND_VALUE",
                                              "FIRST_VALUE"),
                        0);

    /* flip this back to a map */
    FLUNK_IF(g_svc->fptr->objectToMap(g_svc, obj.value.object, &dd));

    FLUNK_IF(g_svc->fptr->mapGetTypeName(g_svc, dd, &mapType));
    
    printf("Map is of type: %s\n", mapType);

    /* pull an element out of the map */
    CHECK(EDGSDRValue(&elem, DR_UNDEFINED));
    FLUNK_IF(g_svc->fptr->mapGetItem(g_svc, dd, getme, &elem));

    top(file ? file : stdout);
    if (!displayObject(&elem, getme, file ? file : stdout)) {
        ERROR_EXIT("failed to display object");
    }

    /* display the whole thing (covers pulling matrices apart) */
    if (!displayObject(&obj, "EnumTester", file ? file : stdout)) {
        ERROR_EXIT("failed to display VolSurface object");
    }
    tail(file ? file : stdout);
    ok = DR_TRUE;
done:
    if (file) {
        fclose(file);
    }
    DRI_FREE(g_svc, dd);
    CHECK(g_svc->fptr->valueClear(g_svc, &obj));
    CHECK(g_svc->fptr->valueClear(g_svc, &elem));
    DRI_STRING_FREE(g_svc, mapType);
    return ok;
}

/** build an exported type (this is the documentation service) */
static DRObject docMaker(const char* filename, const char* classname){
    static const char method[] = "docMaker";
    DRValue  vfile;
    DRValue  vclass;
    DRValue  params[2];
    DRMap    dd;
    void*    id;
    DRObject o;
    DRBool   ok;

    CHECK(EDGSDRValue(&vfile, DR_STRING, filename));
    CHECK(EDGSDRValue(&vclass, DR_STRING, classname));

    params[0] = vfile;
    params[1] = vclass;

    CHECK(g_svc->fptr->idFromName(g_svc, "ClassDocumentation", &id));
    CHECK(g_svc->fptr->mapNewFromIdAndValues(g_svc,
                                             id,
                                             params,
                                             2,
                                             &dd));
    CHECK(g_svc->fptr->mapToObject(g_svc, dd, &o));
    DRI_FREE(g_svc, dd);

    /** illustrate call to objectIsExecutable */
    CHECK(g_svc->fptr->objectIsExecutable(g_svc, o, &ok));

    if (!ok) {
        handleError("docMaker - something's gone very wrong here", 0);
    }

    return o;
}


/* other ways to pull maps apart */
static int mapValues(const char* filename) {
    static const char method[] = "mapValues";
    DRBool   ok = DR_FALSE;
    DRObject o = 0;
    DRMap    dd = 0;
    DRArray  elem = 0;
    DRValue  elemV;

    FILE* file = 0;
    if (filename) {
        file = fopen(filename, "w");
    }

    /* make an exported type */
    if (!(o = docMaker("ClassDocumentation.html", "ClassDocumentation"))) {
        ERROR_EXIT("failed to make ClassDocumentation.html");
    }

    /* flip this back to a map */
    FLUNK_IF(g_svc->fptr->objectToMap(g_svc, o, &dd));

    /* rip it too bits */
    FLUNK_IF(g_svc->fptr->mapGetValues(g_svc, dd, &elem));

    CHECK(EDGSDRValue(&elemV, DR_OBJECT, elem));

    top(file ? file : stdout);
    if (!displayObject(&elemV, "ClassDocumentation", file ? file : stdout)) {
        ERROR_EXIT("failed to display ClassDocumentation object");
    }
    tail(file ? file : stdout);

    ok = DR_TRUE;
done:
    if (file) {
        fclose(file);
    }
    DRI_FREE(g_svc, o);
    DRI_FREE(g_svc, dd);
    DRI_FREE(g_svc, elem);
    return ok;
}

/* cause a deliberate error */
static void makeItFail(const char* filename) {
    DRError error = 0;
    DRMap   dd = 0;
    FILE*   file = 0;
    if (filename) {
        file = fopen(filename, "w");
    }

    error = g_svc->fptr->mapNew(g_svc, "TypeThatDoesNotExist", &dd);

    if (error) {
        top(file ? file : stdout);
        fprintf(file ? file : stdout, "%s\n", error);
        tail(file ? file : stdout);
        DRI_ERROR_FREE(g_svc, error);
    }

    if (file) {
        fclose(file);
    }
}

static void typeInfo(char* clazz, const char* filename) {
    static const char method[] = "typeInfo";
    DRObject desc = 0;
    FILE* file = 0;
    DRError err = 0;

    if (filename) {
        file = fopen(filename, "w");
    }

    if ((err = g_svc->fptr->typeDescription(g_svc, clazz, &desc))) {
        handleError("typeInfo", err);
        DRI_ERROR_FREE(g_svc, err);
    }
    else {
        DRValue x;
        CHECK(EDGSDRValue(&x, DR_OBJECT, desc));
        displayObject(&x, clazz, file ? file : stdout);
    }

    DRI_FREE(g_svc, desc);

    if (file) {
        fclose(file);
    }
}

/* demo how to build and price a variety of instruments, get price &
   greeks, & display on screen. Rename this to main() to get a self
   contained executable */
REGTEST_DLL int drimaxiMain (int argc, char* argv[])
{
    static const char method[] = "drimaxiMain";
    DRError err = 0;

    DRBool ok = DR_FALSE;

    DRValue   today = {0};
    DRValue   mdlCF = {0};
    DRValue   mdlCFLN = {0};
    DRValue   inst[2];
    DRValue   ctrl[2];
    DRValue   scenario = {0};
    DRValue   model[2];

    DRObject  compImnt = 0;
    DRValue   mkt1;
    DRObject  mkt2 = 0;
    DRValue   results = {0};
    DRValue   results2 = {0};

    double    multiplier[2] = {10.0, 1.0};
    double    weight[2] = {1.0, 1.0};
    int       i;
    int       numInst = 1; /* default */

    /* turn library on */
    DRServiceInitArgs args;
    EDGSInitServiceCreateArgs(&args);

    if (!DRCreateService(&args, &g_svc)) {
        printf("Failed to create EDR Service!\n");
        EDGSGetServiceInvocationError(&err);
        // Because svc = 0, err won't be freed in handleError().
        ERROR_EXIT(err);
    }

    /* init these DRValues */
    CHECK(EDGSDRValue(&inst[0], DR_UNDEFINED));
    CHECK(EDGSDRValue(&inst[1], DR_UNDEFINED));
    CHECK(EDGSDRValue(&ctrl[0], DR_UNDEFINED));
    CHECK(EDGSDRValue(&ctrl[1], DR_UNDEFINED));
    CHECK(EDGSDRValue(&model[0], DR_UNDEFINED));
    CHECK(EDGSDRValue(&model[1], DR_UNDEFINED));

    if (argc == 1){
        DRString version = 0;
        CHECK(g_svc->fptr->serviceDescription(g_svc, &version));
        const char* pg = argv[0];
        printf("Version: %s\n", version);
        printf("Usage: \n");
        printf("1. %s -typeList [optional file name]\n", pg);
        printf("2. %s -map1 [optional file name]\n", pg);
        printf("3. %s -map2 [optional file name]\n", pg);
        printf("4. %s -fail\n", pg);
        printf("5. %s [instrumentNumber] [outFile1] [outFile2]\n", pg);
        printf("\t where instrumentNumber is\n");
        printf("\t 0 for cashflow\n");
        printf("\t 1 for average on plain equity\n");
        printf("\t 2 for average on protected equity\n");
        printf("\t 3 for average on struck equity\n");
        printf("\t 4 for average on XCB\n");
        printf("\t 5 for composite\n");
        printf("\t 6 for DR Wrapper\n");
        printf("6. %s -typeInfo [class] [optional file name]\n", pg);
        printf("7. %s -addmkt", pg);
        DRI_STRING_FREE(g_svc, version);
        return 0;
    }
    if (argc > 1 && strcmp(argv[1], "-typeList") == 0){
        typeList(argc == 3 ? argv[2]: 0);
    }  
    else if (argc > 1 && strcmp(argv[1], "-map1") == 0){
        mapMaker(argc == 3 ? argv[2]: 0);
    }  
    else if (argc > 1 && strcmp(argv[1], "-map2") == 0){
        mapValues(argc == 3 ? argv[2]: 0);
    }  
    else if (argc > 1 && strcmp(argv[1], "-map3") == 0){
        mapMaker3(argc == 3 ? argv[2]: 0);
    }  
    else if (argc > 1 && strcmp(argv[1], "-fail") == 0){
        makeItFail(argc == 3 ? argv[2]: 0);
    } 
    else if (argc > 1 && strcmp(argv[1], "-typeInfo") == 0){
        typeInfo(argv[2], argc == 4 ? argv[3]: 0);
    } 
    else {
        /* now build up the inputs to the RiskManager call */
        int instID = 0;
        char* file1 = 0;
        char* file2 = 0;

        DRString description = 0;

        CHECK(g_svc->fptr->serviceDescription(g_svc, &description));
        printf("Using DR service: %s\n", description);
        DRI_STRING_FREE(g_svc, description);

        int instIndex = 1;
        if (argc > 1 && strcmp(argv[1], "-addmkt") == 0){
            ++instIndex;
        }
        if (argc > 2){
            sscanf(argv[instIndex], "%i", &instID);
        }
        if (argc >= 3 && strcmp(argv[2],"-") != 0){
            file1 = argv[instIndex + 1];
        }
        if (argc >= 4 && strcmp(argv[3],"-") != 0){
            file2 = argv[instIndex + 2];
        }

        mdlCF = mkObj("ClosedForm", 0);
        if (mdlCF.type == DR_UNDEFINED) {
            ERROR_EXIT("failed to make ClosedForm");
        }

        mdlCFLN = mkObj("ClosedFormLN", 
                        "volType", mkStr("IVolatilityBS"),
                        0);

        if (mdlCFLN.type == DR_UNDEFINED) {
            ERROR_EXIT("failed to make ClosedFormLN");
        }
        
        ctrl[0] = control();
        if (ctrl[0].type == DR_UNDEFINED) {
            ERROR_EXIT("failed to make control");
        }
        ctrl[1] = ctrl[0]; /* reuse the control for simplicity */

        // Preparing markets
        printf("Building MarketDataCache\n");

        int numMkts = 1;
      
        /* create two market data caches - both the same except for GBP
           yield curve */
        mkt1 = mkMarket();
        if (mkt1.type == DR_UNDEFINED) {
            ERROR_EXIT("failed to build MarketDataCache");
        }
 
        /* add data to market data cache */
        if (!initialiseMarketDataCache3(mkt1.value.object)){
            ERROR_EXIT("failed to initialize MarketDataCache3");
        }

        if (argc > 1 && strcmp(argv[1], "-addmkt") == 0){
            printf("Testing EXTEND_MARKET addin\n");
            /* add [different] GBP yield curves to mkt1 */

#ifdef _ADDMARKET_ADDIN_TEST
            /* This AddMarketAddin test is only useful if AddMarketAddin is
               set public and made to implement ClientRunnable. */
            DRMap addmkt;
            FLUNK_IF(g_svc->fptr->mapNew(g_svc, "AddMarketAddin", &addmkt));
            FLUNK_IF(g_svc->fptr->mapAddItem(g_svc, addmkt, "market", &mkt1));

            DRArray array;
            FLUNK_IF(g_svc->fptr->arrayNew(g_svc, 1, "MarketObject", &array));
            DRValue csCurve = mkCashSwapCurve1();
            FLUNK_IF(g_svc->fptr->arraySet(g_svc, array, 0, &csCurve));
            CHECK(g_CHECK(svc->fptr->valueClear(g_svc, &csCurve)));

            DRValue objs;
            objs.type = DR_OBJECT;
            objs.value.object = array;
            
            FLUNK_IF(g_svc->fptr->mapAddItem(g_svc, addmkt, "objects", &objs));

            DRObject addMktObj;
            FLUNK_IF(g_svc->fptr->mapToObject(g_svc, addmkt, &addMktObj));

            DRValue addResult;
            FLUNK_IF(g_svc->fptr->execute(g_svc, addMktObj, &addResult));

#else                                             
            if (!initialiseMarketDataCache1(mkt1.value.object)) {
                ERROR_EXIT("failed to initialize MarketDataCache1");
            }
#endif
        }
        else {
            /* create copy of mk1 */
            if (!(mkt2 = EDGSObjectClone(g_svc, mkt1.value.object))){
                goto done;
            }

            ++numMkts;

            /* add [different] GBP yield curves to mkt1 and mkt2 */
            if (!initialiseMarketDataCache1(mkt1.value.object)) {
                ERROR_EXIT("failed to initialize MarketDataCache1");
            }
            if (!initialiseMarketDataCache2(mkt2)){
                ERROR_EXIT("failed to initialize MarketDataCache2");
            }
        }

#ifdef _ADDMARKET_ADDIN_TEST
        { DRString myString = EDGSValueXMLWrite(g_svc, &mkt1);
            cout << myString << endl; }
#endif

        scenario = mkObj("Scenario", 
                         "shifts", 
                         mkObjArray("ScenarioShift", 2, 
                                                 mkObj("ScenarioShift",
                                                       "sensCtrl",
                                                       mkObj("SpotLevel", 
                                                             "shiftSize", mkDb(7000),
                                                             0),
                                                       "marketDataName", mkStr(FTSE),
                                                       0),
                                    mkObj("ScenarioShift",
                                          "sensCtrl",
                                          mkObj("SpotLevel", 
                                                "shiftSize", mkDb(2000),
                                                             0),
                                          "marketDataName", mkStr(DAX),
                                          0),
                                    0),
                         0);
        if (scenario.type == DR_UNDEFINED){
            ERROR_EXIT("failed to make Scenario");
        }

        switch (instID){
        case 0:
            inst[0]  = cfinst();
            model[0] = mdlCF;
            break;
        case 1:
            inst[0]  = avgSpotInst1();
            model[0] = mdlCFLN;
            break;
        case 2:
            inst[0]  = avgSpotInst2();
            model[0] = mdlCFLN;
            break;
        case 3:
            inst[0]  = avgSpotInst3();
            model[0] = mdlCFLN;
            break;
        case 4:
            inst[0]  = avgSpotInst4();
            model[0] = mdlCFLN;
            break;
        case 5:
            inst[0]  = cfinst();
            model[0] = mdlCF;
            inst[1]  = avgSpotInst4();
            model[1] = mdlCFLN;
            numInst = 2;
            break;
        case 6:
            inst[0]  = cfwrapper();
            model[0] = mdlCF;
            break;
        default:
            ERROR_EXIT("Instrument choice must lie in [0,6]");
        }
        /* Test out xml read/write from char stream */
        if (instID == 0){
            printf("Testing out xml read/write methods\n");
            DRObject newInst = testXMLReadWrite(&inst[0]);
            if (!newInst){
                ERROR_EXIT("failed to xml read/write new instrument");
            }
            CHECK(g_svc->fptr->valueClear(g_svc, &inst[0]));
            inst[0].type = DR_OBJECT;
            inst[0].value.object = newInst;
        }
        printf("Running\n");
        for (i = 0; i < numMkts; i++){
            printf("Running for market data cache %d\n", i+1);
            
            if (!(compImnt = composite(g_svc,
                                       numInst,
                                       model,
                                       inst,
                                       ctrl,
                                       weight,
                                       multiplier,
                                       scenario,
                                       i == 0? mkt1.value.object: mkt2))) {
                ERROR_EXIT("failed to make composite");
            }
            
            FLUNK_IF(g_svc->fptr->execute(g_svc, compImnt, &results));

            /* demonstrate ability to create regression file */
            if ((i == 0 || i == 6) && file1){
                DRValue c;
                char* inpFile = (char*)malloc(strlen(file1)+5);
                if (!inpFile){
                    ERROR_EXIT(0);
                }
                strcpy(inpFile, file1);
                strcat(inpFile, ".xml");

                CHECK(EDGSDRValue(&c, DR_OBJECT, compImnt));

                FLUNK_IF(EDGSValueWrite(g_svc, inpFile, &c));

                /* then see if we can read the file back in */
                {
                    DRValue theObject;
                    FLUNK_IF(EDGSValueRead(g_svc, inpFile, &theObject));
                    CHECK(g_svc->fptr->valueClear(g_svc, &theObject));
                }

                free(inpFile);
                FLUNK_IF(EDGSWriteResultsToFile(g_svc, file1,
                                                results.value.object));
                
                /* if it's a wrapper read it in, then write that out
                   to show that it's associative
                */
                results2.type = DR_OBJECT;
                if (!(results2.value.object = EDGSDRWrapperRead(g_svc, file1))) {
                    ERROR_EXIT(0);
                }

                FLUNK_IF(EDGSWriteResultsToFile(g_svc, file1,
                                                results2.value.object));
            }
            char buffer[500] = "";
            FILE* file = 0;
            if (file2){
                if (!(file = fopen(file2, i == 0? "w": "a"))){
                    ERROR_EXIT(0);
                }
                if (i == 0){
                    fprintf(file, "<OutputFile>\n");
                }
                fprintf(file, "<Output>\n<Summary>\nPermutation %d:\n",
                        i+1);
            }
            if (!displayObject(&results, buffer, file2? file: stdout)){
                ERROR_EXIT(0);
            }
            if (file2){
                fprintf(file, "</Summary>\n</Output>\n");
                if (i == 1){
                    fprintf(file, "</OutputFile>\n");
                }
                fclose(file);
            }
            CHECK(g_svc->fptr->valueClear(g_svc, &results));
            DRI_FREE(g_svc, compImnt);
        }
    }
    ok = DR_TRUE;
done:
    if (g_svc){
        CHECK(g_svc->fptr->valueClear(g_svc, &today));
        CHECK(g_svc->fptr->valueClear(g_svc, &scenario));
        CHECK(g_svc->fptr->valueClear(g_svc, &mkt1));
        DRI_FREE(g_svc, mkt2);
        CHECK(g_svc->fptr->valueClear(g_svc, &mdlCF));
        CHECK(g_svc->fptr->valueClear(g_svc, &mdlCFLN));
        CHECK(g_svc->fptr->valueClear(g_svc, &inst[0]));
        CHECK(g_svc->fptr->valueClear(g_svc, &inst[1]));
        CHECK(g_svc->fptr->valueClear(g_svc, &ctrl[0]));
        CHECK(g_svc->fptr->valueClear(g_svc, &results2));
        CHECK(g_svc->fptr->valueClear(g_svc, &results));

        /* turn library off */
        CHECK(g_svc->fptr->serviceFree(g_svc));
    }
    return 0;
}
