/****************************************************************
 * Equity Derivatives Research: Models Library Interface
 ****************************************************************/

#ifndef EDR_EASINTERFACE_H
#define EDR_EASINTERFACE_H

/* set up DLL_EXPORT correctly. Define EDR_NO_DLL if you want to link
   against the static libraries */
#ifdef EDR_NO_DLL
#define DLL_EXPORT
#else
#ifdef _MSC_VER
/* only for NT of course */
#ifdef EDG_DR_PRIVATE
/* EDG_DR_PRIVATE should only be defined when the source code for EDR is 
   being compiled */
#define DLL_EXPORT __declspec(dllexport)
#else
#define DLL_EXPORT __declspec(dllimport)
#endif
#else
/* on unix this define expands to nothing */
#define DLL_EXPORT
#endif
#endif

#ifdef __cplusplus
extern "C" {
#endif

/** and true/false */
#define EDR_TRUE  1
#define EDR_FALSE 0

    /** opaque structure for holding all objects */
    typedef struct _EDRObject_  {int opaque;} * EDRObject;
    /* more specialised objects: */
    typedef EDRObject EDRDataDict;
    typedef EDRObject EDRDataDictIterator;
    typedef EDRObject EDRResultsIterator;
    typedef EDRObject EDRMarketDataCache;
    typedef EDRObject EDROutputName;

    typedef int EDRBool;

    /** VITAL to call this before DOING ANYTHING AT ALL - return TRUE/FALSE*/
    DLL_EXPORT EDRBool EdrStartup(void);

    DLL_EXPORT void EdrShutdown(void);

    /** return version as a dynamic C style string - 
     use free() to cleanup */
    DLL_EXPORT char* EdrVersion();

    /* Note: Internally refence counting is used. Therefore once an object
       has been created and stored either in a data dictionary or in the
       market data cache, that object can safely be freed using the
       EdrObjectFree method below */

    /* Note: This interface is suitable for calling from C. No
       exceptions are thrown - all error codes are given by the return
       value. THEREFORE ALL RETURN CODES MUST BE CHECKED */

    /************* Creation of Atomic Objects  *********************/

    /** Creates a wrapper to an integer. Returns NULL on failure */
    DLL_EXPORT EDRObject EdrIntNew(int value);

    /** Creates a wrapper to a double. Returns NULL on failure */
    DLL_EXPORT EDRObject EdrDoubleNew(double value);

    /** Creates a wrapper to a C style null terminated string (will
        take a copy of the supplied string) . Returns NULL on failure */
    DLL_EXPORT EDRObject EdrStringNew(const char* value);

    /** Creates a wrapper to a bool. value must be 0 or 1. Returns NULL on
        failure */
    DLL_EXPORT EDRObject EdrBoolNew(EDRBool value);

    /** Creates a wrapper to a double matrix. The matrix is created
        with all elements set to 0. Returns NULL on failure */
    DLL_EXPORT EDRObject EdrDoubleMatrixEmptyNew(int             numCols,
                                      int             numRows);

    /** Sets a value in a double matrix created with via
        EdrDoubleMatrixEmptyNew */
    DLL_EXPORT EDRBool EdrDoubleMatrixSet(EDRObject       matrix,
                                          int             col,
                                          int             row,
                                          double          value);
    
    /** Creates a datetime object. Date is in dd-mmm-yyyy format.
        Time is hh:mm:ss or "SOD" or "EOD". Returns NULL on failure */
    DLL_EXPORT EDRObject EdrDateTimeNew(const char* date, const char* time);

    /** Creates an expiry (e.g. "6M"). Returns NULL on failure */
    DLL_EXPORT EDRObject EdrMaturityPeriodNew(const char* period);
    
    /** Creates an expiry (e.g. "6M") with a fixed time 
        (either hh-mm-ss or "SOD" or "EOD"). Returns NULL on failure */
    DLL_EXPORT EDRObject EdrMaturityTimePeriodNew(const char* period, 
                                                  const char* time);

    /** Creates an expiry which corresponds to a fixed date. Date is
        in dd-mmm-yyyy format.  Time is hh:mm:ss or "SOD" or "EOD".
        Returns NULL on failure */
    DLL_EXPORT EDRObject EdrBenchmarkDateNew(const char* date, 
                                             const char* time);

    /** Creates a wrapper to an enum which is specifed by the type of the
        enum and a string (null terminated) representing the value to
        use. (No references are kept to either strings) . Returns NULL on
        failure.  A list of possible values can be obtained via
        EdrEnumValues. */
    DLL_EXPORT EDRObject EdrEnumNew(const char* enumType,
                                    const char* value);

    /************ Decomposing Atomic Objects *******************/
    /** Populates d with the value of the supplied object which must
        be a wrapped double (see EdrDoubleNew) . The return value indicates
        SUCCESS/FAILURE. */
    DLL_EXPORT EDRBool EdrDoubleGet(const EDRObject object,
                                    double*         d);
    
    /** returns dynamic string (use free) representing the value of
        the supplied object which must be a wrapped expiry (see
        EdrMaturityPeriodNew and EdrMaturityTimePeriodNew). Returns
        NULL on failure */
    DLL_EXPORT char* EdrExpiryGet(const EDRObject object);

    /** returns dynamic string (use free) representing the value of
        the supplied object which must be of type DateTime. Returns
        NULL on failure */
    DLL_EXPORT char* EdrDateTimeGet(const EDRObject object);

    /** returns dynamic string (use free) representing the value of the
        supplied object which must be a wrapped string (see
        EdrStringNew). Returns NULL on failure */
    DLL_EXPORT char* EdrStringGet(const EDRObject object);

    /** Populates b with the value of the supplied object which must
        be a wrapped bool (see EdrBoolNew). The return value indicates
        SUCCESS/FAILURE. */
    DLL_EXPORT EDRBool EdrBoolGet(const EDRObject object,
                                  EDRBool*        b);
    
    /** Populates i with the value of the supplied object which must
        be a wrapped int (see EdrIntNew) . The return value indicates
        SUCCESS/FAILURE. */
    DLL_EXPORT EDRBool EdrIntGet(const EDRObject object,
                                 int*            i);

    /** Returns the number of rows and columns of the matrix. The
        supplied object must be a wrapped double matrix (see
        EdrDoubleMatrixEmptyNew). Returns EDR_FALSE on failure*/
    DLL_EXPORT EDRBool EdrDoubleMatrixGetSize(const EDRObject object, /* (I) */
                                              int*            numCols,/* (O) */
                                              int*            numRows);/* (O)*/

    /** Returns an element of a double matrix together. The supplied
        object must be a wrapped double matrix (see
        EdrDoubleMatrixNew). Returns EDR_FALSE on failure */
    DLL_EXPORT EDRBool EdrDoubleMatrixGet(const EDRObject object,  /* (I) */
                                          int             col,     /* (I) */
                                          int             row,     /* (I) */
                                          double*         value);  /* (O) */


    /** returns a const char* style string (do NOT use free) representing
        the value of the supplied enum which must be a wrapped enum
        (see EdrEnumNew). Returns NULL on failure */
    DLL_EXPORT const char* EdrEnumGet(const EDRObject object);

    /************* Creation of Array Objects  *********************/

    /** Creates an empty array of object whose components are of
        componentType. The length of the array created is determined by
        numElements. Returns NULL on failure. */
    DLL_EXPORT EDRObject EdrArrayNew(const char* componentType,
                                     int         numElements);

    /** Sets the supplied object in the supplied array at the
        specified index location (index must lie in [0,numElements-1]).
        The supplied object is not deep copied when it is inserted.
        Returns EDR_FALSE on failure */
    DLL_EXPORT EDRBool EdrArraySet(EDRObject array,    /* (M) */
                                   int       index,    /* (I) */
                                   EDRObject object);  /* (I) */

    /** Appends the supplied object to the end of supplied array. The
        length of the array is increased by 1. (aka push_back). 
        The supplied object is not deep copied when it is appended.
        Returns EDR_FALSE on failure*/
    DLL_EXPORT EDRBool EdrArrayAppend(EDRObject array,
                                      EDRObject object);

    /** how long is an array ? */
    DLL_EXPORT EDRBool EdrArrayLength(const EDRObject object,
                                      int*            length);

    /** get the i'th item from an array (starts at 0) */
    DLL_EXPORT EDRObject EdrArrayItem(EDRObject       object, /* (I) */
                                      int             index); /* (I) */
    

    /************* DataDictionary Support *********************/

    /** Creates a data dictionary corresponding to the object of the
        supplied type. Returns NULL on failure */
    DLL_EXPORT EDRDataDict EdrDataDictNew(const char* type);

    /** Adds the specified object to the supplied data dictionary. Returns
        TRUE/FALSE. */
    DLL_EXPORT EDRBool EdrDataDictAdd(EDRDataDict  dataDict,
                                      const char*  fieldName,
                                      EDRObject    objectToAdd);
    
    /** Converts supplied data dictionary into its object
        equivalent. Each of the objects within the data dictionary are
        deep copied in this process */
    DLL_EXPORT EDRObject EdrDataDictToObject(const EDRDataDict dataDict);

    /** Converts supplied data dictionary into its object
        equivalent. No deep copies made */
    DLL_EXPORT EDRObject EdrDataDictToObjectShallow(const EDRDataDict dataDict);

    /** and the other way round. Each of the objects within the data
        dictionary are deep copies. */
    DLL_EXPORT EDRDataDict EdrObjectToDataDict(const EDRObject object);

    /** Get an object out of a data dictionary. The returned object is
        not a deep copy */
    DLL_EXPORT EDRObject EdrDataDictGet(EDRDataDict    dataDict,
                                        const char*    fieldName);

    /** Creates an EDRDataDictIterator object (which must eventually be freed)
        This object can then be used to explore the contents of a data
        dictonary */
    DLL_EXPORT EDRDataDictIterator EdrDataDictIteratorGet(
        EDRDataDict  dataDict);

    /** Returns next <field, object> pair in iteration - field and
        object will be null when there are no more elements in the
        iteration. Neither the field nor the object should be freed. The
        object returned is not a deep copy */
    DLL_EXPORT EDRBool EdrDataDictIteratorNext(
        EDRDataDictIterator iterator, /* (M) */
        const char**        field,    /* (O) */
        EDRObject*          object);  /* (O) */

    /************* Type Driven Interface Support *********************/

    /** Methods to Return information regarding the components needed to build
        an object etc. */

    /** Returns a non dynamic string giving the runtime class of the
        object */
    DLL_EXPORT const char* EdrObjectGetType(const EDRObject object);

    /** Returns a null terminated array of const char * containing the
        types within the library - each element in the array gives the
        name of a type. Use free to release the memory used by the
        array. Do not free the elements of the array */
    DLL_EXPORT const char** EdrTypeList(void);

    /** List all types that can be used to build (via the data
        dictionary route) an object of type typeKey or an object
        derived from type typeKey */
    DLL_EXPORT const char** EdrTypeListConstructorTypes(const char* typeKey);

    /** Returns, via the modifiers parameter, the modifiers for this
        class or interface, encoded in an integer. The values
        EDR_ABSTRACT, EDR_INTERFACE, EDR_PROTECTED and EDR_PUBLIC
        should be used to decode the integer using 'bitwise and'.
        Returns TRUE/FALSE for success/failure */
    DLL_EXPORT EDRBool EdrTypeGetModifiers(const char* type,       /* (I) */
                                           int*        modifiers); /* (O) */

#define EDR_ABSTRACT  0x0001 /* value representing the abstract
                                modifier. Abstract classes cannot be
                                instantiated */

#define EDR_INTERFACE 0x0004 /* value representing the interface
                                modifier. An interface is a class with no
                                data fields and which defines a set of
                                pure virtual methods */

#define EDR_PROTECTED 0x0010 /* value representing the protected
                                modifier. A protected class cannot be
                                constructed via the data dictionary
                                route. Instead a 'helper' type must be
                                used to construct the object. The list of
                                helper types for a protected type may be
                                obtained by calling
                                EdrTypeListConstructorTypes */

#define EDR_PUBLIC    0x0020 /* value representing the public modifier. A
                                public class may be constructed via the
                                data dictionary route.  */
    
    /** find out how many items are in a type description */
    DLL_EXPORT EDRBool EdrTypeNumItems(
        const char* type,
        int*        numComponents);     /* (O) */
        
    /** return the field name of the i'th item in a type description 
        do NOT free the output of this function */
    DLL_EXPORT const char* EdrTypeFieldName(
        const char* type,
        int         index);

    /** return the type name of the i'th item in a type description 
        do NOT free the output of this function */
    DLL_EXPORT const char* EdrTypeFieldType(
        const char* type,
        int         index);

    /** is a field in a type optional ? */
    DLL_EXPORT EDRBool EdrTypeFieldIsOptional(
        const char* type,
        int         index,
        EDRBool*    optional);     /* (O) */

    /** return the description of the i'th item in a type description 
        do NOT free the output of this function */
    DLL_EXPORT const char* EdrTypeFieldDescription(
        const char* type,
        int         index);

    /** is a type an array ? The return value indicates SUCCESS/FAILURE */
    DLL_EXPORT EDRBool EdrTypeIsArray(
        const char* type,
        EDRBool*    isArray);     /* (O) */

    /** what is the type of element in an array ? 
        do NOT free the output of this function */
    DLL_EXPORT const char* EdrArrayElemType(const char* type);

    /** Does a type represent an enum? The return value indicates
        SUCCESS/FAILURE */
    DLL_EXPORT EDRBool EdrTypeIsEnum(
        const char* type,
        EDRBool*    isEnum);     /* (O) */

    /** Returns a null terminated array of const char * containing either
        the possible values for the specified enum type (getDescriptions =
        FALSE) or the descriptions for each of these values
        (getDescriptions = TRUE). When getDescriptions = FALSE, each
        element in the array gives a possible value with which an enum of
        the specified type can be built. 
        For either value of getDescriptions, use free to release the memory
        used by the array but do not free the elements of the array.  Returns
        NULL on failure. */
    DLL_EXPORT const char** EdrEnumValues(
        const char* enumType,          /* (I) */
        EDRBool     getDescriptions);  /* (I) */

    /** Determines if an object whose type is given by typeOfObject is
        derived from typeOfClass. A java style isInstance() method could be
        implemented via clazz->isAssignableFrom(object->getClass()) */
    DLL_EXPORT EDRBool EdrTypeIsAssignableFrom(const char* typeOfClass,
                                               const char* typeOfObject,
                                               EDRBool*    assignable);

    /************* Market Data Cache support ******************/

    /** Creates a new cache of market data. Returns NULL on failure */
    DLL_EXPORT EDRMarketDataCache EdrMarketDataCacheNew(
        const EDRObject valueDate);

    /** Adds the specified object to the supplied market data
        cache. Returns TRUE/FALSE. Supplied object must be derived
        from MarketObject. The supplied object is not deep copied when
        it is added to the cache */
    DLL_EXPORT EDRBool EdrMarketDataCacheAdd(EDRMarketDataCache  cache,
                                             EDRObject           objectToAdd);

    /** Creates a copy of a market data cache. Performance wise: the
        contents of the cache are not deep copied, ie references are
        taken to the data within the cache */
    DLL_EXPORT EDRMarketDataCache EdrMarketDataCacheClone(
        EDRMarketDataCache cache);

    /************* DR Wrapper/Regression Test Creaton support **************/
    /** Creates a DR Wrapper corresponding to the object of the
        supplied type. Returns NULL on failure */
    DLL_EXPORT EDRDataDict EdrDRWrapperNew(const char* type);

    /** Write inputs to file */
    DLL_EXPORT EDRBool EdrWriteInputsToFile(
        const char*          filename,
        int                  lengthOfArrays,
        EDRObject*           models,
        EDRObject*           instruments,
        EDRObject*           controls,
        double*              multipliers,    
        double*              weights,
        EDRObject            scenario,
        EDRMarketDataCache   cache);
    
    /** Writes the supplied results object to file (includes a 'summary' of the
        results at the top of the file) */
    DLL_EXPORT EDRBool EdrWriteResultsToFile(const char*          filename,
                                             EDRObject            results);

    /** Writes any object to file */
    DLL_EXPORT EDRBool EdrWriteObjectToFile(const char* filename,
                                            EDRObject   object);

    /** returns object written to file using EdrWriteObjectToFile. 
        Also works for file created using EdrWriteInputsToFile - in this
        case the object returned will be of type "CompositeInstrument". Use
        EDR_TYPE_COMPONENTS("CompositeInstrument") in Excel to see
        what fields this object contains */
    DLL_EXPORT EDRObject EdrReadObjectFromFile(const char* filename);

    /** read wrapper results back from file. Same output as EdrMain */
    DLL_EXPORT EDRObject EdrDRWrapperRead(const char*  filename);    /* (I) */

    /** return wrapper errors from file as a dynamic C style string - 
     use free() to cleanup */
    DLL_EXPORT char* EdrDRWrapperError(const char* filename);  /* (I) */

    /** Parses the supplied string and build the object that the xml describes.
        The xml must conform to the EDR representation */
    DLL_EXPORT EDRObject EdrXMLRead(const char* xmlString);

    /** Takes the supplied object and writes it, in xml format, to a
        string. The returned string is dynamic and must be freed using
        free() */
    DLL_EXPORT char* EdrXMLWrite(const EDRObject object);


    /************* Pricing/Sensitivities/Scenario support ******************/
    
    /** Main Pricing Call. lengthOfArrays reflects the number of
        instruments to price. The object returned will be of different
        types depending on whether lengthOfArrays == 1 or not. In the
        former case, an object of type "Results" will be returned. In the
        latter, a typed object will be returned. It is expected that
        this will be a "ResultsSet" object which will contain the
        individual results as well as the composite (ie combined)
        result.  Note that the weights and multipliers are still used
        even if lengthOfArrays == 1. */
    DLL_EXPORT EDRObject EdrMain(int                  lengthOfArrays, /* (I) */
                                 EDRObject*           models,         /* (I) */
                                 EDRObject*           instruments,    /* (I) */
                                 EDRObject*           controls,       /* (I) */
                                 double*              multipliers,    /* (I) */
                                 double*              weights,        /* (I) */
                                 EDRObject            scenario,       /* (I) */
                                 EDRMarketDataCache   cache);         /* (I) */

    /** Generic function for all other actions (e.g. event handling).
        Create object of appropriate type, pass it in, and get results
        back
    */
    DLL_EXPORT EDRObject EdrAction(EDRObject input);  /* (I) */
        

    /***************** Results support ****************************
     * Results are stored within packets with each packet typically
     * holding the results for a different sensitivity. As such
     * results are identified by which packet they belong to and their
     * OutputName which typically holds the name(s) of what was
     * tweaked to get this result. 
     *************************************************************/
    
    /** Creates an EDRResultsIterator object (which must eventually be freed)
        This object can then be used to explore the contents of the results
        object. The supplied object must be of type "Results" */
    DLL_EXPORT EDRResultsIterator EdrResultsIteratorGet(EDRObject  results);

    /** Returns next <packet name, identifier, object> triple in
        iteration - packet, name and object will be null when there
        are no more elements in the iteration. Neither the packet nor
        the name nor the object should be freed. Both the name and
        object returned are not deep copies */
    DLL_EXPORT EDRBool EdrResultsIteratorNext(
        EDRResultsIterator   iterator, /* (M) */
        const char**         packet,   /* (O) */
        EDROutputName*       name,     /* (O) */
        EDRObject*           object);  /* (O) */
    
    /** Returns the number of strings which make up an output name */
    DLL_EXPORT EDRBool EdrOutputNameIDCount(
        const EDROutputName outputName,  /* (I) */
        int*                numStrings); /* (O) */

    /** Returns the component string, identified by index, which
        makes up the supplied EDROutputName */
    DLL_EXPORT const char* EdrOutputNameIdGet(const EDROutputName outputName,  /* (I) */
                                   int                 index);      /* (I) */
        
    /************* Memory Support *********************/
    
    /** Frees the memory contained within an EDRObject. */
    DLL_EXPORT void EdrObjectFree(EDRObject object);

    /** Frees the memory of a C string pointer (char*) returned by EDR.
        this is provided to avoid the need of explict call to free() by client which may 
        give rise problems due to lib versions. */
    DLL_EXPORT void EdrCStrFree(char* str);

    /************* Error Handling Support *********************/
    /** Returns true if there are no errors outstanding */
    DLL_EXPORT EDRBool EdrErrorCleared();

    /** return stack trace as a dynamic C style string - 
     use free() to cleanup */
    DLL_EXPORT char* EdrErrorStackTrace();

#ifdef __cplusplus
}
#endif

#endif
