/*************************************************************
 * Equity Derivatives Research: Example of use of Models Library Interface
 *
 *
*************************************************************/
/* this file is not needed - just used for precompiled headers */
#include "edginc/config.hpp"
#include "edginc/EDRInterface.h"

/************************* Utility Functions ***************************/

static EDRBool displayObject(EDRObject object, char* buffer, FILE* file);

static void top(FILE* file) {
    fprintf(file, "<OutputFile>\n");
    fprintf(file, "<Output>\n<Summary>\nPermutation 1:\n");
}
    
static void tail(FILE* file) {
    fprintf(file, "</Summary>\n</Output>\n");
    fprintf(file, "</OutputFile>\n");
}
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

/** Make a string */
static EDRObject mkStr(const char* str){
    EDRObject obj = EdrStringNew(str);
    if (!obj){
        handleError("mkStr");
    }
    return obj;
}

/** Make an int */
static EDRObject mkInt(int val){
    EDRObject obj = EdrIntNew(val);
    if (!obj){
        handleError("mkInt");
    }
    return obj;
}
/** Make a double */
static EDRObject mkDb(double val){
    EDRObject obj = EdrDoubleNew(val);
    if (!obj){
        handleError("mkDb");
    }
    return obj;
}
/** Make a bool */
static EDRObject mkBl(EDRBool val){
    EDRObject obj = EdrBoolNew(val);
    if (!obj){
        handleError("mkBl");
    }
    return obj;
}

/** Make an enum */
static EDRObject mkEnum(const char* enumType, const char* enumValue){
    EDRObject obj = EdrEnumNew(enumType, enumValue);
    if (!obj){
        handleError("mkEnum");
    }
    return obj;
}

/** Make a DateTime */
static EDRObject mkDt(const char* val1, const char* val2){
    EDRObject obj = EdrDateTimeNew(val1, val2);
    if (!obj){
        handleError("mkDt");
    }
    return obj;
}
/** Make a MaturityPeriod */
static EDRObject mkMat(const char* val){
    EDRObject obj = EdrMaturityPeriodNew(val);
    if (!obj){
        handleError("mkMat");
    }
    return obj;
}
/** Make a DoubleMatrix */
static EDRObject mkDbMat(int cols, int rows, ...){
    va_list        args;
    int            i, j;
    EDRObject      obj = EdrDoubleMatrixEmptyNew(cols, rows);
    if (!obj){
        handleError("mkMat");
    }
    va_start(args, rows);
    for (i = 0; i < cols; i++){
        for (j = 0; j < rows; j++){
            double val = va_arg(args, double);
            if (!EdrDoubleMatrixSet(obj, i, j, val)){
                handleError("mkDbMat");
                return 0;
            }
        }
    }
    return obj;
}

/** Make a MaturityTimePeriod */
#if 0
/* not currently used */
static EDRObject mkMatTime(const char* val1, const char* val2){
    EDRObject obj = EdrMaturityTimePeriodNew(val1, val2);
    if (!obj){
        handleError("mkMatTime");
    }
    return obj;
}
#endif

/** Create an array of strings */
static EDRObject mkStrArray(int length, ...){
    va_list        args;
    int            i;
    EDRObject      array = EdrArrayNew("String", length);
    if (!array){
        goto done;
    }
    va_start(args, length);
    for (i = 0; i < length; i++){
        const char* str = va_arg(args, const char*);
        EDRObject obj = EdrStringNew(str);
        if (!obj || !EdrArraySet(array, i, obj)){
            goto done;
        }
        EdrObjectFree(obj);
    }
    return array;
done:
    handleError("mkStrArray");
    return 0;
}

/** Create an array of ints */
#if 0
/* not currently used */
static EDRObject mkIntArray(int length, ...){
    va_list        args;
    int            i;
    EDRObject      array = EdrArrayNew("Int", length);
    if (!array){
        goto done;
    }
    va_start(args, length);
    for (i = 0; i < length; i++){
        int val = va_arg(args, int);
        EDRObject obj = EdrIntNew(val);
        if (!obj || !EdrArraySet(array, i, obj)){
            goto done;
        }
        EdrObjectFree(obj);
    }
    return array;
done:
    handleError("mkIntArray");
    return 0;
}
#endif

/** Create an array of doubles */
static EDRObject mkDbArray(int length, ...){
    va_list        args;
    int            i;
    EDRObject      array = EdrArrayNew("Double", length);
    if (!array){
        goto done;
    }
    va_start(args, length);
    for (i = 0; i < length; i++){
        double val = va_arg(args, double);
        EDRObject obj = EdrDoubleNew(val);
        if (!obj || !EdrArraySet(array, i, obj)){
            goto done;
        }
        EdrObjectFree(obj);
    }
    return array;
done:
    handleError("mkDbArray");
    return 0;
}

/** Create an array of specified type, inputs terminated by 0. Note inputs
    are also freed */
static EDRObject mkObjArray(const char* type, ...){
    va_list        args;
    EDRObject      array = EdrArrayNew(type, 0);
    if (!array){
        goto done;
    }
    va_start(args, type);
    EDRObject val;
    do {
        val = va_arg(args, EDRObject);
        if (val && !EdrArrayAppend(array, val)){
            goto done;
        }
        EdrObjectFree(val);
    } while (val);
    return array;
done:
    handleError("ObjArray");
    return 0;
}

/** Create a DR Wrapper of a specified type, inputs terminated by 0. Note inputs
    are also freed */
static EDRObject mkWrapper2(const char* type, va_list args){
    const char*    field;
    EDRDataDict    dd = EdrDRWrapperNew(type);
    if (!dd){
        goto done;
    }
    do {
        field = va_arg(args, const char*);
        if (field){
            EDRObject val = va_arg(args, EDRObject);
            if (!val || !EdrDataDictAdd(dd, field, val)){
                goto done;
            }
            EdrObjectFree(val);
        }
    } while (field);
    return dd;
done:
    handleError("mkObj2");
    return 0;
}

static EDRObject mkWrapper(const char* type, ...){
    va_list        args;
    va_start(args, type);
    return mkWrapper2(type, args);
}

/** Create an object of a specified type, inputs terminated by 0. Note inputs
    are also freed */
static EDRObject mkObj2(const char* type, va_list args){
    const char*    field;
    EDRObject      dd = EdrDataDictNew(type);
    EDRObject      obj;
    if (!dd){
        goto done;
    }
    do {
        field = va_arg(args, const char*);
        if (field){
            EDRObject val = va_arg(args, EDRObject);
            if (!val || !EdrDataDictAdd(dd, field, val)){
                goto done;
            }
            EdrObjectFree(val);
        }
    } while (field);
    if (!(obj = EdrDataDictToObject(dd))) {
        goto done;
    }
    EdrObjectFree(dd);
    return obj;
done:
    handleError("mkObj2");
    return 0;
}

/** Create an object of a specified type, inputs terminated by 0. Note inputs
    are also freed */
static EDRObject mkObj(const char* type, ...){
    va_list        args;
    va_start(args, type);
    return mkObj2(type, args);
}

/** Creates object and stores it in supplied cache */
static EDRBool addObj(EDRMarketDataCache cache,  const char* type, ...){
    va_list        args;
    va_start(args, type);
    EDRObject obj = mkObj2(type, args);
    if (!obj){
        goto done;
    }
    if (!EdrMarketDataCacheAdd(cache, obj)){
        goto done;
    }
    EdrObjectFree(obj);
    return EDR_TRUE;
done:
    handleError("addObj");
    return EDR_FALSE;
}

/********* Utility Functions for querying type information ************/

/* returns the method needed to create an array of a particular type. Returns
   0 if the generic route should be taken */
const char* arrayMethodName(const char* type){
    const char* lookUp[] = {"Double", "mkDbArray()",
                            "Int",    "mkIntArray()",
                            "String", "mkStrArray()"};
    static int size = sizeof(lookUp)/(2 * sizeof(char*));
    int i;
    for (i = 0; i < size; i++){
        if (strcmp(lookUp[2*i], type) == 0){
            return lookUp[2*i+1];
        }
    }
    return 0;
}

/* returns the method needed to create an object of a particular type. Returns
   0 if the generic route should be taken  */
const char* objMethodName(const char* type){
    static const char* lookUp[] = {"Double", "mkDb()",
                                   "Int",    "mkInt()", 
                                   "String", "mkStr()",
                                   "Bool",   "mkBl()",
                                   "DateTime", "mkDt(, )",
                                   "MaturityPeriod", "mkMat()",
                                   "DoubleMatrix", "mkDbMat(, , ,)",
                                   "MaturityTimePeriod", "mkMatTime()"};
    static int size = sizeof(lookUp)/(2 * sizeof(char*));
    int i;
    for (i = 0; i < size; i++){
        if (strcmp(lookUp[2*i], type) == 0){
            return lookUp[2*i+1];
        }
    }
    return 0;
}

static void typeTableDemoMain(const char* type);
                    
/** Illustrates use of type information to build a 'recipe' for
    building an object up. It generates the templates used for the 'mkObj'
    calls used for creating the objects in the rest of this program */
void typeTableDemo(const char* type){
    int    components;
    int    i;
    const char* method;
    if (!EdrTypeNumItems(type, &components)) {
        goto done;
    }
    if ((method = objMethodName(type))){
        printf("%s,\n", method);
    } else {
        for (i = 0; i < components; i++) {
            EDRBool isArray;
            EDRBool isWrapper;
            EDRBool isEnum;
            const char* name = EdrTypeFieldName(type, i);
            const char* fieldType = EdrTypeFieldType(type, i);
            const char* desc = EdrTypeFieldDescription(type, i);
            if (!desc || !fieldType || !name ||
                !EdrTypeIsArray(fieldType, &isArray) ||
                !EdrTypeIsEnum(fieldType, &isEnum)){
                goto done;
            }
            printf("\"%s\", ", name);
            if (!EdrTypeIsAssignableFrom("MarketObjectWrapper", fieldType, 
                                         &isWrapper)){
                goto done;
            }
            /* override wrappers with just the string holding the
               object's name */
            if (isWrapper){
                printf("mkStr(),\n");
            } else if (isArray){
                const char* eltType = EdrArrayElemType(fieldType);
                if (!eltType){
                    goto done;
                }
                if (!EdrTypeIsAssignableFrom("MarketObjectWrapper", eltType, 
                                             &isWrapper)){
                    goto done;
                }
                if (isWrapper){
                    printf("mkStrArray(),\n");
                } else {
                    method = arrayMethodName(eltType);
                    if (method){
                        printf("%s,\n", method);
                    } else {
                        printf("mkObjArray(\"%s\", \n", eltType);
                        typeTableDemoMain(eltType);
                        printf("0),\n");
                    }
                }
            } else if (isEnum){
                /* find out what the possible values are */
                const char** enumValues = EdrEnumValues(fieldType, EDR_FALSE);
                if (!enumValues){
                    goto done;
                }
                printf("mkEnum(\"%s\", {", fieldType);
                const char** p = enumValues;
                while (*p){
                    if (p != enumValues){
                        printf("|");
                    }
                    printf(*p);
                    p++;
                }
                printf("}),\n");
                free(enumValues);
            } else {
                method = objMethodName(fieldType);
                if (method){
                    printf("%s,\n", method);
                } else {
                    printf("mkObj(\"%s\",\n", fieldType);
                    typeTableDemo(fieldType);
                    printf("0),\n");
                }
            }
        }
    }
    return;
done:
    handleError("typeTableDemo");
}

static void typeTableDemoMain(const char* type){
    printf("mkObj(\"%s\", \n", type);
    typeTableDemo(type);
    printf("0)\n");
}

/** Writes object to stream and then reads it back again */
static EDRObject testXMLReadWrite(EDRObject inObject){
    // read existing xml file into char*
    EDRObject object = 0;
    char* myString = EdrXMLWrite(inObject);
    if (!myString){
        goto done;
    }
    if (!(object = EdrXMLRead(myString))){
        goto done;
    }
    free(myString);
    return object;
done:
    handleError("testXMLReadWrite");
    return object;
}

/** Prints out all the known types which are derived from the supplied
    type. If the supplied type is 0, all types are printed */
static void typeDemo2(const char* type){
    const char** allTypes = EdrTypeList();
    const char** pos = allTypes;
    while(*pos){
        EDRBool printOut = EDR_TRUE;
        if (type){
            if (!EdrTypeIsAssignableFrom(type, *pos, 
                                         &printOut)){
                goto done;
            }
        }
        if (printOut){
            int modifiers;
            if (!EdrTypeGetModifiers(*pos, &modifiers)){
                goto done;
            }
            printf("%s\t( ", *pos);
            if (modifiers & EDR_PUBLIC){
                printf("public ");
            }
            if (modifiers & EDR_PROTECTED){
                printf("protected ");
            }
            if (modifiers & EDR_ABSTRACT){
                printf("abstract ");
            }
            if (modifiers & EDR_INTERFACE){
                printf("interface ");
            }
            printf(")\n");
        }
        pos++;
    }
done:
    free(allTypes);
}

/** Prints out all the known 'constructor' types which can be used to build
    an instance of type or a derived instance of type */
static void typeDemo3(const char* type){
    const char** allTypes = EdrTypeListConstructorTypes(type);
    const char** pos = allTypes;
    if (!allTypes){
        goto done;
    }
    while(*pos){
        printf("%s\n", *pos);
        pos++;
    }
done:
    free(allTypes);
}

/* demo building an object with enums in it, demo pulling it apart. Note
   that clients should use reflection to find out possible values for enums.
   The example here of FIRST_VALUE etc are just the names chosen for the
   examples here. Each enum will typically have different values. */
static int mapMaker3(const char* filename, EDRBool useStringsForEnums) {
    EDRBool     ok = EDR_FALSE;
    EDRDataDict dd = 0;
    EDRObject   cmpt = 0;
    char* getme = "myEnum1";
    FILE* file = 0;
    if (filename) {
        file = fopen(filename, "w");
    }

    EDRObject obj = useStringsForEnums?
        mkObj("EnumTester", 
              /* enums can go through either as strings
                 or as explicit types. Here we demo using strings */
              "myEnum1", mkStr("FIRST_VALUE"),
              "myEnum2", mkStr("THIRD_VALUE"),
              "myEnums", mkStrArray(3,
                                    "THIRD_VALUE",
                                    "SECOND_VALUE",
                                    "FIRST_VALUE"),
              0):
        mkObj("EnumTester", 
              /* enums can go through either as strings
                 or as explicit types. Here we demo using explicit types */
              "myEnum1", mkEnum("EnumTester::MyValues", "FIRST_VALUE"),
              "myEnum2", mkEnum("EnumTester::MyValues", "THIRD_VALUE"),
              "myEnums", mkObjArray("EnumTester::MyValues",
                                    mkEnum("EnumTester::MyValues", 
                                           "THIRD_VALUE"),
                                    mkEnum("EnumTester::MyValues", 
                                           "SECOND_VALUE"),
                                    mkEnum("EnumTester::MyValues", 
                                           "FIRST_VALUE"),
                                    0),
              0);


    /* flip this back to a map */
    dd = EdrObjectToDataDict(obj);
    if (!dd){
        goto done;
    }

    /* pull an element out of the map */
    if (!(cmpt = EdrDataDictGet(dd, getme))){
        goto done;
    }

    top(file ? file: stdout);
    if (!displayObject(cmpt, getme, file ? file : stdout)) {
        goto done;
    }

    /* display the whole thing (covers pulling matrices apart) */
    if (!displayObject(obj, "EnumTester", file ? file : stdout)) {
        goto done;
    }
    tail(file ? file: stdout);
    ok = EDR_TRUE;
done:
    if (file) {
        fclose(file);
    }
    EdrObjectFree(obj);
    EdrObjectFree(dd);
    EdrObjectFree(cmpt);
    return ok;
}

/********************** End of Utility Functions *************************/

static const char* todayDate = "05-Apr-2001";
static const char* todayTime = "SOD";

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

/* Just inserts a GBP yield curve into the supplied cache */
static EDRBool initialiseMarketDataCache1(EDRMarketDataCache cache){
    // create GBP Yield Curve and put in cache
    if (!addObj(cache, "CashSwapCurve::Interface",
                "ccy", mkStr("GBP"),
                "name", mkStr(GBP_LIBOR),
                "today", mkDt(todayDate, todayTime),
                "spotOffset", mkInt(0),
                "moneyMarketDenom", mkInt(365),
                "swapFrequency", mkInt(2),
                "swapDayCount", mkStr("Act/360"),
                "hols", mkObj("Holiday", 
                              "name", mkStr("noweekends"),
                              "useWeekends", mkBl(EDR_FALSE),
                              "holidays", mkObjArray("DateTime", 0),
                              0),
                "expiries", mkObjArray("Expiry",
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
                0)){
        goto done;
    }

    return EDR_TRUE;
done:
    handleError("initialiseMarketDataCache1");
    return 0;

}

/* Just inserts a GBP yield curve (with slightly different rates than
   initialiseMarketDataCache1) into the supplied cache */
static EDRBool initialiseMarketDataCache2(EDRMarketDataCache cache){
    // create GBP Yield Curve and put in cache
    if (!addObj(cache, "CashSwapCurve::Interface",
                "ccy", mkStr("GBP"),
                "name", mkStr(GBP_LIBOR),
                "today", mkDt(todayDate, todayTime),
                "spotOffset", mkInt(0),
                "moneyMarketDenom", mkInt(365),
                "swapFrequency", mkInt(2),
                "swapDayCount", mkStr("Act/360"),
                "hols", mkObj("Holiday", 
                              "name", mkStr("noweekends"),
                              "useWeekends", mkBl(EDR_FALSE),
                              "holidays", mkObjArray("DateTime", 0),
                              0),
                "expiries", mkObjArray("Expiry",
                                       mkMat("ON"), mkMat("1M"), 
                                       mkMat("2M"), mkMat("3M"),  
                                       mkMat("6M"), mkMat("1Y"), 
                                       mkMat("2Y"), mkMat("3Y"),
                                       mkMat("5Y"), mkMat("7Y"), 
                                       mkMat("10Y"), mkMat("20Y"), 
                                       mkMat("30Y"), 0),
                "instruments", mkStrArray(13, "M", "M", "M", "M", "M", "M",
                                          "S", "S", "S", "S", "S", "S", "S"),
                "rates", mkDbArray(13, 0.05, 0.051, 0.052, 0.053, 0.054, 0.055,
                                   0.056,  0.057, 0.058, 0.059, 0.06, 0.061,
                                   0.062),
                "interp", mkObj("FourPlusI", 0),
                0)){
        goto done;
    }

    return EDR_TRUE;
done:
    handleError("initialiseMarketDataCache2");
    return 0;

}

/* populates the cache with the rest of the remaining data */
static EDRBool initialiseMarketDataCache3(EDRMarketDataCache cache){
    // create USD Yield Curve and put in cache
    if (!addObj(cache, "CashSwapCurve::Interface",
                "ccy", mkStr("USD"),
                "name", mkStr(USD_LIBOR),
                "today", mkDt(todayDate, todayTime),
                "spotOffset", mkInt(0),
                "moneyMarketDenom", mkInt(365),
                "swapFrequency", mkInt(2),
                "swapDayCount", mkStr("Act/360"),
                "hols", mkObj("Holiday", 
                              "name", mkStr("noweekends"),
                              "useWeekends", mkBl(EDR_FALSE),
                              "holidays", mkObjArray("DateTime", 0),
                              0),
                "expiries", mkObjArray("Expiry",
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
                0)){
        goto done;
    }

    // create USD Yield Curve and put in cache
    if (!addObj(cache, "CashSwapCurve::Interface",
                "ccy", mkStr("EUR"),
                "name", mkStr(EUR_LIBOR),
                "today", mkDt(todayDate, todayTime),
                "spotOffset", mkInt(0),
                "moneyMarketDenom", mkInt(365),
                "swapFrequency", mkInt(2),
                "swapDayCount", mkStr("Act/360"),
                "hols", mkObj("Holiday", 
                              "name", mkStr("noweekends"),
                              "useWeekends", mkBl(EDR_FALSE),
                              "holidays", mkObjArray("DateTime", 0),
                              0),
                "expiries", mkObjArray("Expiry",
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
                0)){
        goto done;
    }

    if (!addObj(cache, "VolSurface", 
                "name", mkStr(FTSE_VOL),
                "metric", mkObj("TimeMetric",
                                "nonTradTimeFrac", mkDb(1.0),
                                "marketHols", 
                                mkObj("Holiday",
                                      "name", mkStr("noweekends"),
                                      "useWeekends", mkBl(EDR_FALSE),
                                      "holidays", mkObjArray("DateTime",
                                                             0),
                                      0),
                                0),
                "strikes", mkDbArray(3, 5000.0, 6000.0, 7000.0),
                "vol", mkDbMat(3, 4, 
                               0.5, 0.55, 0.6, 0.65,
                               0.4, 0.45, 0.5, 0.55,
                               0.3, 0.35, 0.4, 0.45),
                "expiries", mkObjArray("Expiry", 
                                       mkMat("1W"), mkMat("1M"), 
                                       mkMat("1Y"), mkMat("2Y"),  
                                       0),
                "baseDate", mkDt(todayDate, todayTime),
                0)){
        goto done;
    }

    if (!addObj(cache, "VolSurface", 
                "name", mkStr(DAX_VOL),
                "metric", mkObj("TimeMetric",
                                "nonTradTimeFrac", mkDb(1.0),
                                "marketHols", 
                                mkObj("Holiday",
                                      "name", mkStr("noweekends"),
                                      "useWeekends", mkBl(EDR_FALSE),
                                      "holidays", mkObjArray("DateTime",
                                                             0),
                                      0),
                                0),
                "strikes", mkDbArray(3, 3000.0, 5000.0, 6000.0),
                "vol", mkDbMat(3, 4, 
                               0.2, 0.25, 0.3, 0.35,
                               0.25, 0.35, 0.4, 0.45,
                               0.3, 0.35, 0.4, 0.45),
                "expiries", mkObjArray("Expiry", 
                                       mkMat("1W"), mkMat("1M"), 
                                       mkMat("1Y"), mkMat("2Y"),  
                                       0),
                "baseDate", mkDt(todayDate, todayTime),
                0)){
        goto done;
    }

    if (!addObj(cache, "VolSurface", 
                "name", mkStr(SPX_VOL),
                "metric", mkObj("TimeMetric",
                                "nonTradTimeFrac", mkDb(1.0),
                                "marketHols", 
                                mkObj("Holiday",
                                      "name", mkStr("noweekends"),
                                      "useWeekends", mkBl(EDR_FALSE),
                                      "holidays", mkObjArray("DateTime",
                                                             0),
                                      0),
                                0),
                "strikes", mkDbArray(3, 2000.0, 4000.0, 7000.0),
                "vol", mkDbMat(3, 4, 
                               0.5, 0.55, 0.6, 0.65,
                               0.4, 0.45, 0.5, 0.55,
                               0.3, 0.35, 0.4, 0.45),
                "expiries", mkObjArray("Expiry", 
                                       mkMat("1W"), mkMat("1M"), 
                                       mkMat("1Y"), mkMat("2Y"),  
                                       0),
                "baseDate", mkDt(todayDate, todayTime),
                0)){
        goto done;
    }

    if (!addObj(cache, 
                "SimpleEquity", 
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
                                mkBl(EDR_FALSE),
                                "holidays",
                                mkObjArray("DateTime", 
                                           0),
                                0),
                          0),
                    "yc", mkStr(GBP_LIBOR),
                    "valueDate", mkDt(todayDate, todayTime),
                    "stockDate", mkDt(todayDate, todayTime),
                    "divList", 
                    mkObj("DividendList",
                          "divArray",
                          mkObjArray("Dividend",
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
                          mkObjArray("Expiry", 
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
                0)){
        goto done;
    }

    if (!addObj(cache, 
                "SimpleEquity", 
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
                                mkBl(EDR_FALSE),
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
                          mkObjArray("Expiry", 
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
                0)){
        goto done;
    }

    if (!addObj(cache, 
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
                                mkBl(EDR_FALSE),
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
                          mkObjArray("Dividend",
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
                          mkObjArray("Expiry", 
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
                0)){
        goto done;
    }

    if (!addObj(cache, "FlatFXVol", 
                "name", mkStr(FX_USD_GBP_VOL),
                "flatFXVol", mkDb(0.3),
                "baseDate", mkDt(todayDate, todayTime),
                "timeMetric", 
                mkObj("TimeMetric",
                      "nonTradTimeFrac", mkDb(1.0),
                      "marketHols", mkObj("Holiday",
                                          "name", mkStr("noweekends"),
                                          "useWeekends", mkBl(EDR_FALSE),
                                          "holidays", mkObjArray("DateTime", 
                                                                 0),
                                          0),
                     0),
                0)){
        goto done;
    }
    if (!addObj(cache, "FXAsset", 
                "name", mkStr(FX_GBP_USD),
                "today", mkDt(todayDate, todayTime),
                "holidays", mkObj("Holiday",
                                  "name", mkStr("noweekends"),
                                  "useWeekends", mkBl(EDR_FALSE),
                                  "holidays", mkObjArray("DateTime", 
                                                         0),
                                  0),
                "riskCcy", mkStr(GBP_LIBOR),
                "baseCcy", mkStr(USD_LIBOR),
                "spotFX", mkDb(1.42),
                "fxVol", mkStr(FX_USD_GBP_VOL),
                0)){
        goto done;
    }
    if (!addObj(cache, "FlatFXVol", 
                "name", mkStr(FX_USD_EUR_VOL),
                "flatFXVol", mkDb(0.2),
                "baseDate", mkDt(todayDate, todayTime),
                "timeMetric", 
                mkObj("TimeMetric",
                      "nonTradTimeFrac", mkDb(1.0),
                      "marketHols", mkObj("Holiday",
                                          "name", mkStr("noweekends"),
                                          "useWeekends", mkBl(EDR_FALSE),
                                          "holidays", mkObjArray("DateTime", 
                                                                 0),
                                          0),
                     0),
                0)){
        goto done;
    }
    if (!addObj(cache, "FXAsset", 
                "name", mkStr(FX_EUR_USD),
                "today", mkDt(todayDate, todayTime),
                "holidays", mkObj("Holiday",
                                  "name", mkStr("noweekends"),
                                  "useWeekends", mkBl(EDR_FALSE),
                                  "holidays", mkObjArray("DateTime", 
                                                         0),
                                  0),
                "riskCcy", mkStr(EUR_LIBOR),
                "baseCcy", mkStr(USD_LIBOR),
                "spotFX", mkDb(1.1),
                "fxVol", mkStr(FX_USD_EUR_VOL),
                0)){
        goto done;
    }
    if (!addObj(cache, "Correlation", 
                "name", mkStr("Corr-FTSE_USD"),
                "asset1", mkStr(FTSE),
                "asset2", mkStr(FX_GBP_USD),
                "correlation", mkDb(0.3),
                0)){
        goto done;
    }
    if (!addObj(cache, "Correlation", 
                "name", mkStr("Corr-DAX_USD"),
                "asset1", mkStr(DAX),
                "asset2", mkStr(FX_EUR_USD),
                "correlation", mkDb(0.3),
                0)){
        goto done;
    }
    if (!addObj(cache, "Correlation", 
                "name", mkStr("Corr-FTSE-DAX"),
                "asset1", mkStr(FTSE),
                "asset2", mkStr(DAX),
                "correlation", mkDb(0.1),
                0)){
        goto done;
    }
    if (!addObj(cache, "Correlation", 
                "name", mkStr("Corr-FTSE-SPX"),
                "asset1", mkStr(FTSE),
                "asset2", mkStr(SPX),
                "correlation", mkDb(0.2),
                0)){
        goto done;
    }
    if (!addObj(cache, "Correlation", 
                "name", mkStr("Corr-DAX-SPX"),
                "asset1", mkStr(DAX),
                "asset2", mkStr(SPX),
                "correlation", mkDb(0.3),
                0)){
        goto done;
    }

    if (!addObj(cache, "XCB", 
                "name", mkStr("INDEX_TRIO"),
                "assets", mkStrArray(3, FTSE, DAX, SPX),
                "ccyTreatments", mkStrArray(3, "S", "P", "N"),
                "basketYCName", mkStr(USD_LIBOR),
                "unitWeights", mkBl(EDR_FALSE),
                "pubWeights", mkDbArray(3, 1.0/3, 1.0/3, 1.0/3),
                "marketHols", mkObj("Holiday",
                                    "name", mkStr("noweekends"),
                                    "useWeekends", mkBl(EDR_FALSE),
                                    "holidays", mkObjArray("DateTime", 
                                                           0),
                                    0),
                "startDate", mkDt("01-Jan-2001", todayTime),
                "spotsAtStart", mkDbArray(3, 6500.0, 5500.0, 4500.0),
                "smileType", mkStr("E"),
                "timeMetric", mkObj("TimeMetric",
                                    "nonTradTimeFrac", mkDb(1.0),
                                    "marketHols", 
                                    mkObj("Holiday",
                                          "name", mkStr("noweekends"),
                                          "useWeekends", mkBl(EDR_FALSE),
                                          "holidays", mkObjArray("DateTime",
                                                                 0),
                                          0),
                                    0),
                0)){
        goto done;
    }

    return EDR_TRUE;
done:
    handleError("initialiseMarketDataCache3");
    return 0;

}

/* Prints out the strings making up an EDROutputName */
static EDRBool displayName(EDROutputName outputName, char* buffer){
    EDRBool     ok = EDR_FALSE;
    int numIDs;
    int i;
    if (!EdrOutputNameIDCount(outputName, &numIDs)){
        goto done;
    }
    for (i = 0; i < numIDs; i++){
        const char* componentID = EdrOutputNameIdGet(outputName, i);
        if (!componentID){
            goto done;
        }
        if (*componentID != '\0'){
            char nameBuffer[200];
            sprintf(nameBuffer, "%s%s", (i > 0? ".": "_"), componentID);
            strcat(buffer, nameBuffer);
        }
    }
    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("displayName");
    }
    return ok;
}

/** Display the results returned from the pricing call */
static EDRBool displayResults(EDRObject results, char* buffer, FILE* file){
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
        char localBuffer[500];
        strcpy(localBuffer, buffer);
        if (*localBuffer != '\0'){
            strcat(localBuffer, "_");
        }
        strcat(localBuffer, packet);
        /* print out name  */
        if (!displayName(name, localBuffer)){
            goto done;
        }
        /* then display object (switch on type) */
        if (!displayObject(object, localBuffer, file)){
            goto done;
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
                                 

/* Recurse through an object printing its components */
static EDRBool displayObject(EDRObject object, char* buffer, FILE* file){
    EDRBool        ok = EDR_FALSE;
    EDRBool        isExpiry;
    /* find out what sort of object the result is */
    const char* objType = EdrObjectGetType(object);
    if (!objType){
        goto done;
    }
    /* display result (switch on type) */
    if (!EdrTypeIsAssignableFrom("Expiry", objType, &isExpiry)){
        goto done;
    }
    if (isExpiry){
        char* expiry = EdrExpiryGet(object);
        if (!expiry){
            goto done;
        }
        fprintf(file, "Output %s: %s\n", buffer, expiry);
        free(expiry);
    } else if (strcmp(objType, "Double") == 0){
        double value;
        if (!EdrDoubleGet(object, &value)) {
            goto done;
        }
        char tmpBuffer[500];
        sprintf(tmpBuffer, ": %f\n", value);
        strcat(buffer, tmpBuffer);
        fprintf(file, "Output %s", buffer);
        /* to do - check for other primitive types including int etc ...*/
    }else if (strcmp(objType, "String") == 0){
        char* s = EdrStringGet(object);
        if (!s){
            goto done;
        }
        fprintf(file, "Output %s: %s\n", buffer, s);
        free(s);
    } else if (strcmp(objType, "NotApplicable") == 0){
        fprintf(file, "Output %s: NotApplicable\n", buffer);
    } else if (strcmp(objType, "Untweakable") == 0){
        fprintf(file, "Output %s: Untweakable\n", buffer);
    } else if (strcmp(objType, "DataDictionary") == 0){
        // printf("Type: %s\n", objType);
        EDRBool iterStatus;
        const char* field;
        EDRObject   cmpt;
        EDRDataDictIterator iter = EdrDataDictIteratorGet(object);
        if (!iter){
            goto done;
        }
        while ((iterStatus = EdrDataDictIteratorNext(iter, &field, &cmpt)) &&
               field){
            char localBuffer[500];
            strcpy(localBuffer, buffer);
            if (*localBuffer != '\0'){
                strcat(localBuffer, "_");
            }
            strcat(localBuffer, field);
            displayObject(cmpt, localBuffer, file);
        }
        if (!iterStatus){
            goto done;
        }
        EdrObjectFree(iter);
    } else if (strcmp(objType, "Results") == 0){
        char localBuffer[500];
        strcpy(localBuffer, buffer);
        displayResults(object, localBuffer, file);
    } else {
        EDRBool isArray;
        if (!EdrTypeIsArray(objType, &isArray)){
            goto done;
        }
        if (isArray){
            int i, length;
            if (!EdrArrayLength(object, &length)){
                goto done;
            }
            for (i = 0; i < length; i++){
                EDRObject elt = EdrArrayItem(object, i);
                if (!elt){
                    goto done;
                }
                char tmpBuffer[10];
                sprintf(tmpBuffer, "[%d]", i);
                char localBuffer[500];
                strcpy(localBuffer, buffer);
                strcat(localBuffer, tmpBuffer);
                if (!displayObject(elt, localBuffer, file)){
                    goto done;
                }
                EdrObjectFree(elt);
            }
        } else {
            EDRBool isEnum;
            if (!EdrTypeIsEnum(objType, &isEnum)){
                goto done;
            }
            if (isEnum){
                const char* enumValue = EdrEnumGet(object);
                if (!enumValue){
                    goto done;
                }
                fprintf(file, "Output %s: %s\n", buffer, enumValue);
            } else {
                EDRDataDict dd = EdrObjectToDataDict(object);
                if (!dd || !displayObject(dd, buffer, file)){
                    goto done;
                }
                EdrObjectFree(dd);
            }
        }
    }
    ok = EDR_TRUE;
done:
    if (!ok) {
        handleError("displayObject");
    }
    return ok;
}


/* build a CashFlowStream instrument */
static EDRObject cfinst(void) {
    EDRObject inst;
    inst = mkObj("CashFlowStream", 
                 "cfl", mkObjArray("CashFlow", 
                                   mkObj("CashFlow",
                                         "date", mkDt("15-Sep-2003", "EOD"),
                                         "amount", mkDb(100.0),
                                         0),
                                   0),
                 "discount", mkStr(GBP_LIBOR),
                 0);
    if (!inst){
        handleError("cfinst");
    }
    return inst;
}

/* build a CashFlowStream DR Wrapper */
static EDRObject cfwrapper(void) {
    EDRObject inst;
    inst = mkWrapper("CashFlowStream", 
                     "cfl", mkObjArray("CashFlow", 
                                       mkObj("CashFlow",
                                             "date",mkDt("15-Sep-2003","EOD"),
                                             "amount", mkDb(100.0),
                                             0),
                                       0),
                     "discount", mkStr(GBP_LIBOR),
                     0);
    if (!inst){
        handleError("cfinst");
    }
    return inst;
}


/* build an average spot instument - vanilla ccy treatment */
static EDRObject avgSpotInst1(void) {
    EDRObject inst;
    inst = mkObj("AverageSpot", 
                 "isCall", mkBl(EDR_TRUE),
                 "maturity", mkDt("15-Apr-2002", "EOD"),
                 "strike", mkDb(6000),
                 "avgOut", mkObj("SampleList",
                                 "dates", 
                                 mkObjArray("DateTime", 
                                            mkDt("15-Apr-2002", "EOD"),
                                            0),
                                 "values", mkDbArray(1, 0.0),
                                 "weights", mkDbArray(1, 1.0),
                                 0),
                 "fwdStarting", mkBl(EDR_FALSE),
                 "startDate", mkDt(todayDate, todayTime),
                 "oneContract", mkBl(EDR_TRUE),
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
    if (!inst){
        handleError("cfinst");
    }
    return inst;
}

/* build an average spot instument - protected ccy treatment (FTSE into $)*/
static EDRObject avgSpotInst2(void) {
    EDRObject inst;
    inst = mkObj("AverageSpot", 
                 "isCall", mkBl(EDR_TRUE),
                 "maturity", mkDt("15-Apr-2002", "EOD"),
                 "strike", mkDb(6000),
                 "avgOut", mkObj("SampleList",
                                 "dates", 
                                 mkObjArray("DateTime", 
                                            mkDt("15-Apr-2002", "EOD"),
                                            0),
                                 "values", mkDbArray(1, 0.0),
                                 "weights", mkDbArray(1, 1.0),
                                 0),
                 "fwdStarting", mkBl(EDR_FALSE),
                 "startDate", mkDt(todayDate, todayTime),
                 "oneContract", mkBl(EDR_TRUE),
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
    if (!inst){
        handleError("avgSpotInst2");
    }
    return inst;
}

/* build an average spot instument - struck ccy treatment (FTSE into $)*/
static EDRObject avgSpotInst3(void) {
    EDRObject inst;
    inst = mkObj("AverageSpot", 
                 "isCall", mkBl(EDR_TRUE),
                 "maturity", mkDt("15-Apr-2002", "EOD"),
                 "strike", mkDb(6000),
                 "avgOut", mkObj("SampleList",
                                 "dates", 
                                 mkObjArray("DateTime", 
                                            mkDt("15-Apr-2002", "EOD"),
                                            0),
                                 "values", mkDbArray(1, 0.0),
                                 "weights", mkDbArray(1, 1.0),
                                 0),
                 "fwdStarting", mkBl(EDR_FALSE),
                 "startDate", mkDt(todayDate, todayTime),
                 "oneContract", mkBl(EDR_TRUE),
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
    if (!inst){
        handleError("avgSpotInst3");
    }
    return inst;
}

/* build an average spot instument - on an XCB ($ denominated basket) */
static EDRObject avgSpotInst4(void) {
    EDRObject inst;
    inst = mkObj("AverageSpot", 
                 "isCall", mkBl(EDR_TRUE),
                 "maturity", mkDt("15-Apr-2002", "EOD"),
                 "strike", mkDb(100),
                 "avgOut", mkObj("SampleList",
                                 "dates", 
                                 mkObjArray("DateTime", 
                                            mkDt("15-Apr-2002", "EOD"),
                                            0),
                                 "values", mkDbArray(1, 0.0),
                                 "weights", mkDbArray(1, 1.0),
                                 0),
                 "fwdStarting", mkBl(EDR_FALSE),
                 "startDate", mkDt(todayDate, todayTime),
                 "oneContract", mkBl(EDR_TRUE),
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
    if (!inst){
        handleError("avgSpotInst4");
    }
    return inst;
}

/* build a set of MU_S bucket dates */
static EDRObject mudates(void) {

    EDRObject action = mkObj("MuDateBuilder",
                             "today", mkDt(todayDate, todayTime),
                             "mutype", mkStr("SPECIAL"),
                             0);

    EDRObject mudates = EdrAction(action);

    if (!mudates){
        handleError("mudates");
    }
    EdrObjectFree(action);
    return mudates;
}


/* build a Control */
static EDRObject control(void) {
    EDRObject  ctrl;
    printf("Building Control\n");
    ctrl = mkObj("Control", 
                 /* dump to file here only refers to single instrument without
                    scenarios etc */
                 "writeToFile", mkBl(EDR_FALSE),
                 "outputRequests",
                 /* OutputRequest's ask for additional outputs that aren't
                    sensitivities */
                 mkObjArray("OutputRequest", 
                            mkObj("OutputRequest", 
                                  "requestName", mkStr("FWD_AT_MAT"),
                                  0),
                            0),
                 "sens", mkObjArray("Sensitivity",
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
                                                "useWeekends", mkBl(EDR_FALSE),
                                                "holidays", 
                                                mkObjArray("DateTime", 0),
                                                0)),

#endif
                                    0),
                 0);
    if (!ctrl){
        handleError("control");
    }
    return ctrl;
}


/* demo how to build and price a variety of instruments, get price &
   greeks, & display on screen. Rename this to main() to get a self
   contained executable */
REGTEST_DLL int edgmaxiMain(int argc, char* argv[])
{
    EDRBool ok = EDR_FALSE;

    EDRObject   today = 0;
    EDRObject   mdlCF = 0;
    EDRObject   mdlCFLN = 0;
    EDRObject   inst[2] = {0, 0};
    EDRObject   ctrl[2] = {0, 0};
    EDRObject   scenario = 0;
    EDRObject   model[2] = {0, 0};
    EDRMarketDataCache mkt1 = 0;
    EDRMarketDataCache mkt2 = 0;
    EDRObject          results = 0;
    EDRObject          results2 = 0;

    double    multiplier[2] = {10.0, 1.0};
    double    weight[2] = {1.0, 1.0};
    int       i;
    int       numInst = 1; /* default */

    /* turn library on */
    if (!EdrStartup()) {
        goto done;
    }
    if (argc == 1){
        char *version = EdrVersion();
        const char* pg = argv[0];
        printf("Version: %s\n", version);
        printf("Usage: \n");
        printf("1. %s -typeDemo1 typeKey\n", pg);
        printf("2. %s -typeDemo2 [typeKeyToQuery]\n", pg);
        printf("3. %s -typeDemo3 typeKeyToQuery\n", pg);
        printf("4. %s -map3 [outFile] for demo using strings as enums\n", pg);
        printf("5. %s -map4 [outFile] for demo using explicit enums\n", pg);
        printf("6. %s [instrumentNumber] [outFile1] [outFile2]\n", pg);
        printf("\t where instrumentNumber is\n");
        printf("\t 0 for cashflow\n");
        printf("\t 1 for average on plain equity\n");
        printf("\t 2 for average on protected equity\n");
        printf("\t 3 for average on struck equity\n");
        printf("\t 4 for average on XCB\n");
        printf("\t 5 for composite\n");
        printf("\t 6 for DR Wrapper\n");
        free(version);
        return 0;
    }
    if (argc > 2 && strcmp(argv[1], "-typeDemo1") == 0){
        typeTableDemoMain(argv[2]);
        //testXMLReadWrite(argv[2]);
    } else if (argc > 1 && strcmp(argv[1], "-typeDemo2") == 0){
        typeDemo2(argc > 2? argv[2]: 0);
    } else if (argc > 1 && strcmp(argv[1], "-typeDemo3") == 0){
        typeDemo3(argc > 2? argv[2]: 0);
    } else if (argc > 1 && strcmp(argv[1], "-map3") == 0){
        mapMaker3(argc > 2? argv[2]: 0, true);
    } else if (argc > 1 && strcmp(argv[1], "-map4") == 0){
        mapMaker3(argc > 2? argv[2]: 0, false);
    } else {
        /* now build up the inputs to the RiskManager call */
        int instID = 0;
        char* file1 = 0;
        char* file2 = 0;
        char *version = EdrVersion();
        printf("Using version: %s\n", version);
        free(version);

        if (argc > 1){
            sscanf(argv[1], "%i", &instID);
        }
        if (argc >= 3 && strcmp(argv[2],"-") != 0){
            file1 = argv[2];
        }
        if (argc >= 4 && strcmp(argv[3],"-") != 0){
            file2 = argv[3];
        }
        if (!(mdlCF = mkObj("ClosedForm", 0))){
            goto done;
        }
        if (!(mdlCFLN = mkObj("ClosedFormLN", 
                              "volType", mkStr("IVolatilityBS"),
                              0))){
            goto done;
        }
        
        if (!(ctrl[0] = control())) {
            goto done;
        }
        ctrl[1] = ctrl[0]; /* reuse the control for simplicity */

        printf("Building MarketDataCache\n");
        if (!(today = EdrDateTimeNew(todayDate, todayTime))) {
            goto done;
        }

        /* create two market data caches - both the same except for GBP
           yield curve */
        if (!(mkt1 = EdrMarketDataCacheNew(today))) {
            goto done;
        }
        /* add data to market data cache */
        if (!initialiseMarketDataCache3(mkt1)){
            goto done;
        }
        /* create copy of mk1 */
        if (!(mkt2 = EdrMarketDataCacheClone(mkt1))){
            goto done;
        }
        /* add [different] GBP yield curves to mkt1 and mkt2 */
        if (!initialiseMarketDataCache1(mkt1) ||
            !initialiseMarketDataCache2(mkt2)){
            goto done;
        }

        scenario = mkObj("Scenario", 
                         "shifts", 
                         mkObjArray("ScenarioShift",
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
        if (!scenario){
            goto done;
        }
        switch (instID){
        case 0:
            inst[0] = cfinst();
            model[0] = mdlCF;
            break;
        case 1:
            inst[0] = avgSpotInst1();
            model[0] = mdlCFLN;
            break;
        case 2:
            inst[0] = avgSpotInst2();
            model[0] = mdlCFLN;
            break;
        case 3:
            inst[0] = avgSpotInst3();
            model[0] = mdlCFLN;
            break;
        case 4:
            inst[0] = avgSpotInst4();
            model[0] = mdlCFLN;
            break;
        case 5:
            inst[0] = cfinst();
            model[0] = mdlCF;
            inst[1] = avgSpotInst4();
            model[1] = mdlCFLN;
            numInst = 2;
            break;
        case 6:
            inst[0] = cfwrapper();
            model[0] = mdlCF;
            break;
        default:
            printf ("Instrument choice must lie in [0,6]");
            goto done;
        }
        /* Test out xml read/write from char stream */
        if (instID == 0){
            printf("Testing out xml read/write methods\n");
            EDRObject newInst = testXMLReadWrite(inst[0]);
            if (!newInst){
                goto done;
            }
            EdrObjectFree(inst[0]);
            inst[0] = newInst;
        }
        printf("Running EdrMain\n");
        for (i = 0; i < 2; i++){
            printf("Running EdrMain for market data cache %d\n", i+1);
            if (!(results = EdrMain(numInst,
                                    model,
                                    inst,
                                    ctrl,
                                    weight,
                                    multiplier,
                                    scenario,
                                    i == 0? mkt1: mkt2))) {
                goto done;
            }
            
            /* demonstrate ability to create regression file */
            if ((i == 0 || i == 6) && file1){
                char* inpFile = (char*)malloc(strlen(file1)+5);
                if (!inpFile){
                    goto done;
                }
                strcpy(inpFile, file1);
                strcat(inpFile, ".xml");
                if (!EdrWriteInputsToFile(inpFile,
                                          numInst,
                                          model,
                                          inst,
                                          ctrl,
                                          weight,
                                          multiplier,
                                          scenario,
                                          mkt1)) {
                    goto done; /* leaks on failure */
                }
                /* then see if we can read the file back in */
                {
                    EDRObject theObject = EdrReadObjectFromFile(inpFile);
                    if (!theObject){
                        goto done;
                    }
                    EdrObjectFree(theObject);
                }

                free(inpFile);
                if (!EdrWriteResultsToFile(file1,
                                           results)){
                    goto done;
                }
                
                /* if it's a wrapper read it in, then write that out
                   to show that it's associative
                */
                if (!(results2 = EdrDRWrapperRead(file1))) {
                    goto done;
                }

                if (!EdrWriteResultsToFile(file1, results2)) {
                    goto done;
                }
            }
            char buffer[500] = "";
            FILE* file = 0;
            if (file2){
                if (!(file = fopen(file2, i == 0? "w": "a"))){
                    goto done;
                }
                if (i == 0){
                    fprintf(file, "<OutputFile>\n");
                }
                fprintf(file, "<Output>\n<Summary>\nPermutation %d:\n",
                        i+1);
            }
            if (!displayObject(results, buffer, file2? file: stdout)){
                goto done;
            }
            if (file2){
                fprintf(file, "</Summary>\n</Output>\n");
                if (i == 1){
                    fprintf(file, "</OutputFile>\n");
                }
                fclose(file);
            }
            EdrObjectFree(results);
            results = 0;
        }
    }
    ok = EDR_TRUE;
done:
    EdrObjectFree(today);
    EdrObjectFree(scenario);
    EdrObjectFree(mkt1);
    EdrObjectFree(mkt2);
    EdrObjectFree(mdlCF);
    EdrObjectFree(mdlCFLN);
    EdrObjectFree(inst[0]);
    EdrObjectFree(inst[1]);
    EdrObjectFree(ctrl[0]);
    EdrObjectFree(results2);

    if (!ok) {
        handleError("main");
        EdrObjectFree(results);
    }        

    /* turn library off */
    EdrShutdown();
    return 0;
}
