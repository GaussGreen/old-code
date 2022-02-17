/*************************************************************
 * Equity Derivatives Research: Example of use of DR Interface
 *
 *
 *
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

static void handleError(const char* routine, const char *err, 
                        void *cbParam = 0) {
    if (err) {
        printf("%s() failed: %s\n", routine, (char*) err);
    }
    else {
        printf("%s() failed!\n", routine);
    }
}

static DRValue makeDate(DRService* svc, const char* dt, const char* tm) {
    static const char method[] = "makeDate";
    DRValue vdt;
    DRValue vtm;
    DRValue x;
    DRMap   dd;
    DRObject o;

    CHECK(EDGSDRValue(&vdt, DR_STRING, dt));
    CHECK(EDGSDRValue(&vtm, DR_STRING, tm));
        
    CHECK(svc->fptr->mapNew(svc, "DateTime", &dd));
    CHECK(svc->fptr->mapAddItem(svc, dd, "date", &vdt));
    CHECK(svc->fptr->mapAddItem(svc, dd, "time", &vtm));
    CHECK(svc->fptr->mapToObject(svc, dd, &o));
    DRI_FREE(svc, dd);

    CHECK(EDGSDRValue(&x, DR_OBJECT, o));
    return x;
}    

static DRValue makeMat(DRService* svc, const char* period) {
    static const char method[] = "makeMat";
    DRValue p;
    DRValue x;
    DRMap   dd;
    DRObject o;

    CHECK(EDGSDRValue(&p, DR_STRING, period));
        
    CHECK(svc->fptr->mapNew(svc, "MaturityPeriod", &dd));
    CHECK(svc->fptr->mapAddItem(svc, dd, "period", &p));
    CHECK(svc->fptr->mapToObject(svc, dd, &o));
    DRI_FREE(svc, dd);

    CHECK(EDGSDRValue(&x, DR_OBJECT, o));
    return x;
}    

static DRBool displayName(DRService* svc, DRValue* outputName){
    static const char method[] = "displayName";
    DRBool      ok = DR_FALSE;
    DRMap       map = 0;
    DRValue     names;
    DRValue     id;
    int         numIDs = 0;
    int         i;

    FLUNK_IF(svc->fptr->objectToMap(svc, outputName->value.object, &map));
    FLUNK_IF(svc->fptr->mapGetItem(svc, map, "names", &names));
    FLUNK_IF(svc->fptr->arrayLength(svc, names.value.object, &numIDs));

    printf("Name is '");
    for (i = 0; i < numIDs; i++){
        FLUNK_IF(svc->fptr->arrayGet(svc, names.value.object, i, &id));
        printf("%s%s", i > 0? ".": "", id.value.string);
        CHECK(svc->fptr->valueClear(svc, &id));
    }
    printf("'\n");
    ok = DR_TRUE;
done:
    DRI_FREE(svc, map);
    CHECK(svc->fptr->valueClear(svc, &names));
    CHECK(svc->fptr->valueClear(svc, &id));
    return ok;
}

 /* print an ExpiryResultArray on the screen */
static DRBool displayVector(DRService* svc, const char* label, DRArray vector) {
    static const char method[] = "displayVector";
    DRBool      ok = DR_FALSE;
    DRMap       map = 0;
    DRValue     point;
    DRValue     value;
    DRValue     bm;
    int         length;
    int         i;
    double      d;
 
    printf("%s\n", label);

    FLUNK_IF(svc->fptr->arrayLength(svc, vector, &length));

    for (i = 0; i < length; i++) {
        FLUNK_IF(svc->fptr->arrayGet(svc, vector, i, &point));
        FLUNK_IF(svc->fptr->objectToMap(svc, point.value.object, &map));
        FLUNK_IF(svc->fptr->mapGetItem(svc, map, "result", &value));
        FLUNK_IF(svc->fptr->mapGetItem(svc, map, "expiry", &bm));

        d = value.value.real;

        printf("%d %9f\n", i, d);

        CHECK(svc->fptr->valueClear(svc, &point));
        CHECK(svc->fptr->valueClear(svc, &value));
        CHECK(svc->fptr->valueClear(svc, &bm));
        DRI_FREE(svc, map);
        map = 0;
    }

    ok = DR_TRUE;
done:
    CHECK(svc->fptr->valueClear(svc, &point));
    CHECK(svc->fptr->valueClear(svc, &value));
    CHECK(svc->fptr->valueClear(svc, &bm));
    DRI_FREE(svc, map);
    return ok;
}

/** Display the results returned from the pricing call */
static DRBool displayResults(DRService* svc, DRObject results){
    static const char method[] = "displayResults";
    DRBool      ok = DR_FALSE;
    DRString    resultsType = 0;
    DRString    objType = 0;
    DRMap       map = 0;
    DRValue     packet;
    DRValue     id;
    DRValue     value;
    int         numResults;
    int         i;

    if (!results){
        ERROR_EXIT(0);
    }

    /* In DRI 1.2, array and matrix will give errors. */
    CHECK(svc->fptr->objectGetType(svc, results, &resultsType));

    if (resultsType && strcmp(resultsType, "Results")){
        printf("Can't cope with result object of type %s\n", resultsType);
        ERROR_EXIT(0);
    }

    FLUNK_IF(svc->fptr->objectToMap(svc, results, &map));
    
    /** pull out packet, id and value arrays */
    FLUNK_IF(svc->fptr->mapGetItem(svc, map, "packet", &packet));
    FLUNK_IF(svc->fptr->mapGetItem(svc, map, "identifier", &id));
    FLUNK_IF(svc->fptr->mapGetItem(svc, map, "value", &value));
    FLUNK_IF(svc->fptr->arrayLength(svc, packet.value.object, &numResults));

    for (i = 0; i < numResults; i++) {
        DRValue packetElem;
        DRValue idElem;
        DRValue valueElem;

        FLUNK_IF(svc->fptr->arrayGet(svc, packet.value.object, i, 
                                     &packetElem));
        FLUNK_IF(svc->fptr->arrayGet(svc, id.value.object, i, &idElem));
        FLUNK_IF(svc->fptr->arrayGet(svc, value.value.object, i, &valueElem));

        printf("%s:\n", packetElem.value.string);
        /* print out name  */
        if (!displayName(svc, &idElem)){
            ERROR_EXIT(0);
        }

        /* then display result (switch on type) */
        if (valueElem.type == DR_DOUBLE){
            printf("Value is (%f)\n", value.value.real);
        }
        else if (valueElem.type == DR_OBJECT) {
            /* Find out what sort of object the result is. 
               DRI 1.2: Can no longer use objectGetType to get the type of the
               array directly. */
            DRBool isArray;
            FLUNK_IF(svc->fptr->objectIsArray(svc, valueElem.value.object,
                                              &isArray));

            if (isArray) {
                FLUNK_IF(svc->fptr->arrayElementType(svc,
                                                     valueElem.value.object,
                                                     &objType));

                if (strcmp(objType, "ExpiryResult") == 0){
                    displayVector(svc, "Value is:", valueElem.value.object);
                } 
                else {
                    printf("Don't know how to display output of type %s\n",
                           objType);
                }
            }
            else {
                printf("Don't know how to display output.\n");
            }
        }
        else {
            printf("Don't know how to display output of this type\n");
        }  
        CHECK(svc->fptr->valueClear(svc, &packetElem)); 
        CHECK(svc->fptr->valueClear(svc, &idElem));  
        CHECK(svc->fptr->valueClear(svc, &valueElem));  
    }
    ok = DR_TRUE;
done:
    CHECK(svc->fptr->valueClear(svc, &packet));
    CHECK(svc->fptr->valueClear(svc, &id));
    CHECK(svc->fptr->valueClear(svc, &value));
    DRI_STRING_FREE(svc, resultsType);
    DRI_STRING_FREE(svc, objType);
    DRI_FREE(svc, map);
    return ok;
}
                                 
/* build the benchmark dates for a yield curve */
static DRObject benchmarks(DRService* svc) {
    static const char method[] = "benchmarks";
    DRBool ok = DR_FALSE;
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

    DRArray bm = 0;

    FLUNK_IF(svc->fptr->arrayNew(svc, size, "Expiry", &bm));

    DRValue expiry;

    for (i = 0; i < size; i++) {
        expiry = makeMat(svc, periods[i]);
        FLUNK_IF(svc->fptr->arraySet(svc, bm, i, &expiry));
        CHECK(svc->fptr->valueClear(svc, &expiry));
    }

    ok = DR_TRUE;
done:
    if (!ok) {
        CHECK(svc->fptr->valueClear(svc, &expiry));
        DRI_FREE(svc, bm);
        bm = 0;
    }
    return bm;
}

/* build the rates for a yield curve */
static DRObject rates(DRService* svc) {
    static const char method[] = "rates";
    DRBool ok = DR_FALSE;
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

    DRValue db;
    DRArray rt = 0;

    FLUNK_IF(svc->fptr->arrayNew(svc, size, DRI_TYPE_DOUBLE, &rt));

    for (i = 0; i < size; i++) {
        CHECK(EDGSDRValue(&db, DR_DOUBLE, rts[i]));
        FLUNK_IF(svc->fptr->arraySet(svc, rt, i, &db));
    }

    ok = DR_TRUE;
done:
    if (!ok) {
        DRI_FREE(svc, rt);
        rt = 0;
    }
    return rt;
}

/* build the list of cash/swap instruments for a yeild curve */
static DRObject instruments(DRService* svc) {
    static char  method[] = "instruments";
    DRBool       ok = DR_FALSE;
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

    DRValue st;
    DRArray instr = 0; 
    FLUNK_IF(svc->fptr->arrayNew(svc, size, DRI_TYPE_STRING, &instr));

    for (i = 0; i < size; i++) {
        FLUNK_IF(EDGSDRValue(&st, DR_STRING, inst[i]));
        FLUNK_IF(svc->fptr->arraySet(svc, instr, i, &st));
        CHECK(svc->fptr->valueClear(svc, &st));
    }

    ok = DR_TRUE;
done:
    if (!ok) {
        CHECK(svc->fptr->valueClear(svc, &st));
        DRI_FREE(svc, instr);
        instr = 0;
    }
    return instr;
}

/* build a 4+i zero curve interpolator for a yield curve */
static DRObject interp(DRService* svc) {
    static const char method[] = "interp";
    DRBool      ok = DR_FALSE;
    DRMap       map = 0;
    DRObject    zci = 0;

    FLUNK_IF(svc->fptr->mapNew(svc, "FourPlusI", &map));
    FLUNK_IF(svc->fptr->mapToObject(svc, map, &zci));
    ok = DR_TRUE;
done:
    DRI_FREE(svc, map);
    return zci;
}

/* build a "NO_WEEKENDS" style holiday */
static DRObject noweekends(DRService* svc) {
    static const char method[] = "noweekends";
    DRBool      ok = DR_FALSE;
    DRMap       map = 0;
    DRArray     hol = 0;
    DRValue     wknds;
    DRObject    hols = 0;
    DRValue     name;
    DRValue     holsV;

    FLUNK_IF(EDGSDRValue(&name, DR_STRING,  "noweekends"));

    FLUNK_IF(EDGSDRValue(&wknds, DR_BOOL, DR_FALSE));

    FLUNK_IF(svc->fptr->arrayNew(svc, 0, "DateTime", &hols));

    FLUNK_IF(svc->fptr->mapNew(svc, "Holiday", &map));

    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "useWeekends", &wknds));

    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "name", &name));

    holsV.type = DR_OBJECT;
    holsV.value.object = hols;

    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "holidays", &holsV));

    FLUNK_IF(svc->fptr->mapToObject(svc, map, &hol));

    ok = DR_TRUE;
done:
    CHECK(svc->fptr->valueClear(svc, &wknds));
    CHECK(svc->fptr->valueClear(svc, &name));
    DRI_FREE(svc, hols);
    DRI_FREE(svc, map);
    return hols;
}

/* build a cash/swap style yield curve */
static DRObject yieldCurve(DRService* svc)
{
    static const char method[] = "yieldCurve";
    DRBool      ok = DR_FALSE;
    DRObject    yc = 0;
    DRMap       map = 0;
    DRValue     ycbm;
    DRValue     ycrt;
    DRValue     ycinst;
    DRValue     ccy;
    DRValue     name;
    DRValue     today;
    DRValue     spotOffset;
    DRValue     moneyMarketDenom;
    DRValue     swapDayCount;
    DRValue     swapFrequency;
    DRValue     zci;
    DRValue     hols;

    printf("Building yield curve\n");

    ycbm.type = DR_OBJECT;
    if (!(ycbm.value.object = benchmarks(svc))) {
        ERROR_EXIT(0);
    }

    hols.type = DR_OBJECT;
    if (!(hols.value.object = noweekends(svc))) {
        ERROR_EXIT(0);
    }
    
    ycrt.type = DR_OBJECT;
    if (!(ycrt.value.object = rates(svc))) {
        ERROR_EXIT(0);
    }

    ycinst.type = DR_OBJECT;
    if (!(ycinst.value.object = instruments(svc))) {
        ERROR_EXIT(0);
    }

    FLUNK_IF(EDGSDRValue(&ccy, DR_STRING, "GBP"));

    FLUNK_IF(EDGSDRValue(&name, DR_STRING, "GBP-LIBOR"));

    today = makeDate(svc, "05-Apr-2001", "SOD");

    FLUNK_IF(EDGSDRValue(&spotOffset, DR_INT, 0));

    FLUNK_IF(EDGSDRValue(&moneyMarketDenom, DR_INT, 365));

    FLUNK_IF(EDGSDRValue(&swapFrequency, DR_INT, 2));

    FLUNK_IF(EDGSDRValue(&swapDayCount, DR_STRING, "Act/360"));

    zci.type = DR_OBJECT;
    if (!(zci.value.object = interp(svc))) {
        ERROR_EXIT(0);
    }

    FLUNK_IF(svc->fptr->mapNew(svc, "CashSwapCurve::Interface", &map));

    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "ccy", &ccy));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "name", &name));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "today", &today));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "spotOffset", &spotOffset));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "moneyMarketDenom",
                                   &moneyMarketDenom));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "swapFrequency",
                                   &swapFrequency));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "swapDayCount", &swapDayCount));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "hols", &hols));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "expiries", &ycbm));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "instruments", &ycinst));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "rates", &ycrt));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "interp", &zci));

    FLUNK_IF(svc->fptr->mapToObject(svc, map, &yc));

    ok = DR_TRUE;
done:
    CHECK(svc->fptr->valueClear(svc, &ccy));
    CHECK(svc->fptr->valueClear(svc, &name));
    CHECK(svc->fptr->valueClear(svc, &today));
    CHECK(svc->fptr->valueClear(svc, &spotOffset));
    CHECK(svc->fptr->valueClear(svc, &moneyMarketDenom));
    CHECK(svc->fptr->valueClear(svc, &swapFrequency));
    CHECK(svc->fptr->valueClear(svc, &swapDayCount));
    CHECK(svc->fptr->valueClear(svc, &hols));
    CHECK(svc->fptr->valueClear(svc, &ycbm));
    CHECK(svc->fptr->valueClear(svc, &ycinst));
    CHECK(svc->fptr->valueClear(svc, &ycrt));
    CHECK(svc->fptr->valueClear(svc, &zci));
    DRI_FREE(svc, map);
    return yc;
}

/* build a cashflow */
static DRObject cashflow(DRService* svc) {
    static const char method[] = "cashflow";
    DRBool      ok = DR_FALSE;
    DRMap       map = 0;
    DRValue     date;
    DRValue     amount;
    DRObject    cf = 0;
    
    date = makeDate(svc, "15-Sep-2003", "EOD");

    FLUNK_IF(EDGSDRValue(&amount, DR_DOUBLE, 100.0));
    
    FLUNK_IF(svc->fptr->mapNew(svc, "CashFlow", &map));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "date", &date));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "amount", &amount));
    FLUNK_IF(svc->fptr->mapToObject(svc, map, &cf));

    ok = DR_TRUE;
done:
    CHECK(svc->fptr->valueClear(svc, &date));
    CHECK(svc->fptr->valueClear(svc, &amount));
    DRI_FREE(svc, map);
    return cf;
}

/* build a CashFlowStream instrument */
static DRObject cfinst(DRService* svc) {
    static const char method[] = "cfinst";
    DRBool      ok = DR_FALSE;
    DRMap       map = 0;
    DRArray     cfl = 0;
    DRValue     cflV;
    DRValue     discount;
    DRObject    inst = 0;
    DRValue     cf;

    printf("Building CashFlowStream\n");

    FLUNK_IF(EDGSDRValue(&discount, DR_STRING, "GBP-LIBOR"));

    cf.type = DR_OBJECT;
    if (!(cf.value.object = cashflow(svc))) {
        ERROR_EXIT(0);
    }

    FLUNK_IF(svc->fptr->arrayNew(svc, 1, "CashFlow", &cfl));
    FLUNK_IF(svc->fptr->arraySet(svc, cfl, 0, &cf));

    cflV.type = DR_OBJECT;
    cflV.value.object = cfl;

    FLUNK_IF(svc->fptr->mapNew(svc, "CashFlowStream", &map));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "cfl", &cflV));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "discount", &discount));
    FLUNK_IF(svc->fptr->mapToObject(svc, map, &inst));

    ok = DR_TRUE;
done:
    CHECK(svc->fptr->valueClear(svc, &cf));
    CHECK(svc->fptr->valueClear(svc, &discount));
    DRI_FREE(svc, cfl);
    DRI_FREE(svc, map);
    return inst;
}

/* create a simple scalar greek SensControl */
static DRObject scalarGreek(DRService* svc, const char* greek, double shift) {
    static const char method[] = "scalarGreek";
    DRBool      ok = DR_FALSE;
    DRMap       map = 0;
    DRValue     shiftSize;
    DRObject    sensctrl = 0;
    
    printf("Building %s\n", greek);

    FLUNK_IF(EDGSDRValue(&shiftSize, DR_DOUBLE, shift));
    FLUNK_IF(svc->fptr->mapNew(svc, greek, &map));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "shiftSize", &shiftSize));
    FLUNK_IF(svc->fptr->mapToObject(svc, map, &sensctrl));

    ok = DR_TRUE;
done:
    CHECK(svc->fptr->valueClear(svc, &shiftSize));
    DRI_FREE(svc, map);
    return sensctrl;
}

/* build a Control */
static DRObject control(DRService* svc) {
    static const char method[] = "control";
    DRBool      ok = DR_FALSE;
    DRMap       map = 0;
    DRArray     sens = 0;
    DRValue     writeToFile;
    DRObject    ctrl = 0;
    DRValue     rhopll;
    DRValue     rhopw;
    DRValue     sensV;

    printf("Building Control\n");
    FLUNK_IF(EDGSDRValue(&writeToFile, DR_BOOL, DR_FALSE));

    rhopll.type = DR_OBJECT;
    if (!(rhopll.value.object = scalarGreek(svc, "RhoParallel", 0.0001))) {
        ERROR_EXIT(0);
    }
    rhopw.type = DR_OBJECT;
    if (!(rhopw.value.object = scalarGreek(svc, "RhoPointwise", 0.0001))) {
        ERROR_EXIT(0);
    }

    /* take your pick on how to build arrays */
#if 1
    FLUNK_IF(svc->fptr->arrayNew(svc, 2, "Sensitivity", &sens));
    FLUNK_IF(svc->fptr->arraySet(svc, sens, 0, &rhopll));
    FLUNK_IF(svc->fptr->arraySet(svc, sens, 1, &rhopw));
#else
    FLUNK_IF(svc->fptr->arrayNew(svc, 0, "Sensitivity", &sens));
    FLUNK_IF(EDGSArrayAppend(svc, sens, &rhopll));
    FLUNK_IF(EDGSArrayAppend(svc, sens, &rhopw));
#endif

    sensV.type = DR_OBJECT;
    sensV.value.object = sens;

    FLUNK_IF(svc->fptr->mapNew(svc, "Control", &map));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "sens", &sensV));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "writeToFile", &writeToFile));
    FLUNK_IF(svc->fptr->mapToObject(svc, map, &ctrl));

    ok = DR_TRUE;
done:
    CHECK(svc->fptr->valueClear(svc, &writeToFile));
    CHECK(svc->fptr->valueClear(svc, &rhopll));
    CHECK(svc->fptr->valueClear(svc, &rhopw));
    DRI_FREE(svc, sens);
    DRI_FREE(svc, map);
    return ctrl;
}

/* build a ClosedForm model */
static DRObject model(DRService* svc) {
    static const char method[] = "model";
    DRBool      ok = DR_FALSE;
    DRMap       map = 0;
    DRObject    mdl = 0;

    printf("Building ClosedForm\n");
    FLUNK_IF(svc->fptr->mapNew(svc, "ClosedForm", &map));
    FLUNK_IF(svc->fptr->mapToObject(svc, map, &mdl));

    ok = DR_TRUE;
done:
    DRI_FREE(svc, map);
    return mdl;
}

// build a CompositeInstrument
static DRObject composite(DRService* svc, 
                          DRValue    model,
                          DRValue    inst,
                          DRValue    ctrl,
                          DRObject   mkt) {
    static const char method[] = "model";
    DRBool      ok = DR_FALSE;
    DRMap       map = 0;
    DRObject    compo = 0;
    DRValue     models;
    DRValue     insts;
    DRValue     ctrls;
    DRValue     mults;
    DRValue     weights;
    DRValue     mktV;

    DRValue     weight;
    DRValue     multiplier;
    
    mktV.type = DR_OBJECT;
    mktV.value.object = mkt;

    printf("Building CompositeInstrument\n");

    models.type = DR_OBJECT;
    FLUNK_IF(svc->fptr->arrayNew(svc, 1, DRI_TYPE_VARIANT,
                                 &models.value.object));
    FLUNK_IF(svc->fptr->arraySet(svc, models.value.object, 0, &model));
    
    insts.type = DR_OBJECT;
    FLUNK_IF(svc->fptr->arrayNew(svc, 1, DRI_TYPE_VARIANT,
                                 &insts.value.object));
    FLUNK_IF(svc->fptr->arraySet(svc, insts.value.object, 0, &inst));

    ctrls.type = DR_OBJECT;
    FLUNK_IF(svc->fptr->arrayNew(svc, 1, DRI_TYPE_VARIANT,
                                 &ctrls.value.object));
    FLUNK_IF(svc->fptr->arraySet(svc, ctrls.value.object, 0, &ctrl));

    FLUNK_IF(EDGSDRValue(&weight, DR_DOUBLE, 1.0));

    weights.type = DR_OBJECT;
    FLUNK_IF(svc->fptr->arrayNew(svc, 1, DRI_TYPE_DOUBLE,
                                 &weights.value.object));
    FLUNK_IF(svc->fptr->arraySet(svc, weights.value.object, 0, &weight));

    FLUNK_IF(EDGSDRValue(&multiplier, DR_DOUBLE, 10.0));

    mults.type = DR_OBJECT;
    FLUNK_IF(svc->fptr->arrayNew(svc, 1, DRI_TYPE_DOUBLE,
                                 &mults.value.object));
    FLUNK_IF(svc->fptr->arraySet(svc, mults.value.object, 0, &multiplier));
    FLUNK_IF(svc->fptr->mapNew(svc, "CompositeInstrument", &map));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "model", &models));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "inst", &insts));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "ctrl", &ctrls));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "multiplier", &mults));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "weight", &weights));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "market", &mktV));
    FLUNK_IF(svc->fptr->mapToObject(svc, map, &compo));
    ok = DR_TRUE;
done:
    CHECK(svc->fptr->valueClear(svc, &models));
    CHECK(svc->fptr->valueClear(svc, &insts));
    CHECK(svc->fptr->valueClear(svc, &ctrls));
    CHECK(svc->fptr->valueClear(svc, &mults));
    CHECK(svc->fptr->valueClear(svc, &weights));
    DRI_FREE(svc, map);
    return compo;
}

static DRObject market(DRService* svc, DRValue* today) {
    static const char method[] = "market";
    DRBool      ok = DR_FALSE;
    DRMap       map = 0;
    DRObject    mkt = 0;
    
    FLUNK_IF(svc->fptr->mapNew(svc, "MarketData", &map));
    FLUNK_IF(svc->fptr->mapAddItem(svc, map, "today", today));
    FLUNK_IF(svc->fptr->mapToObject(svc, map, &mkt));
    ok = DR_TRUE;
done:
    DRI_FREE(svc, map);
    return mkt;
}


/* demo how to build a simple cashflow instrument, get price & greeks,
   & display on screen. Rename this to main() to get a self contained
   executable */
REGTEST_DLL int driminiMain (void)
{
    static const char method[] = "main";
    DRBool      ok = DR_FALSE;
    DRString    err = 0;
    DRValue     today;
    DRValue     mdl;
    DRValue     inst;
    DRValue     ctrl;
    DRValue     discount;
    
    DRObject    mkt = 0;
    DRValue     results;
    DRObject    compInst = 0;
    
    /* turn library on */
    DRServiceInitArgs args;
    EDGSInitServiceCreateArgs(&args);

    DRService *svc = 0;

    if (!DRCreateService(&args, &svc)) {
        printf("Failed to create EDR Service!\n");
        EDGSGetServiceInvocationError(&err);
        // Because svc = 0, err won't be freed in handleError().
        ERROR_EXIT(0);
    }

    /* now build up the inputs to the RiskManager call */
    today = makeDate(svc, "05-Apr-2001", "SOD");
    
    mdl.type = DR_OBJECT;
    if (!(mdl.value.object = model(svc))) {
        ERROR_EXIT(0);
    }

    inst.type = DR_OBJECT;
    if (!(inst.value.object = cfinst(svc))) {
        ERROR_EXIT(0);
    }
    ctrl.type = DR_OBJECT;
    if (!(ctrl.value.object = control(svc))) {
        ERROR_EXIT(0);
    }

    printf("Building MarketDataCache\n");
    if (!(mkt = market(svc, &today))) {
        ERROR_EXIT(0);
    }
    /* create and then add yield curve to market data */
    discount.type = DR_OBJECT;
    if (!(discount.value.object = yieldCurve(svc))) {
        ERROR_EXIT(0);
    }
    /* can overload mapAddItem for a MarketData */
    FLUNK_IF(svc->fptr->mapAddItem(svc, mkt, "", &discount));

    printf("Building CompositeInstrument\n");
    
    if (!(compInst = composite(svc, mdl, inst, ctrl, mkt))) {
        ERROR_EXIT(0);
    }
    
    printf("Running CompositeInstrument\n");

    FLUNK_IF(svc->fptr->execute(svc, compInst, &results));
    
    if (!displayResults(svc, results.value.object)){
        ERROR_EXIT(0);
    }
    ok = DR_TRUE;
done:
    CHECK(svc->fptr->valueClear(svc, &discount));
    DRI_FREE(svc, mkt);
    CHECK(svc->fptr->valueClear(svc, &results));
    DRI_FREE(svc, compInst);
    CHECK(svc->fptr->valueClear(svc, &today));
    CHECK(svc->fptr->valueClear(svc, &mdl));
    CHECK(svc->fptr->valueClear(svc, &inst));
    CHECK(svc->fptr->valueClear(svc, &ctrl));

    /* turn library off */
    CHECK(svc->fptr->serviceFree(svc));
    return 0;
}
