/* -----------------------------------------------------------------------
 * This proprietary software has been developed strictly for 
 * internal use.  Any use or misuse, intentional or otherwise, which
 * contradicts or places this policy in jeopardy is strictly forbidden.
 *
 * Copyright 2002 J.P. Morgan Chase & Co. All rights reserved.
 * -----------------------------------------------------------------------
 *
 * $Header$
 *
 * $Log$
 * Revision 1.35.2.9  2005/02/09 20:07:12  dzhuravy
 * no message
 *
 * Revision 1.35.2.8  2005/01/18 15:35:27  dzhuravy
 * no message
 *
 * Revision 1.35.2.7  2005/01/14 16:51:00  dzhuravy
 * no message
 *
 * Revision 1.35.2.6  2005/01/11 10:53:42  cbauer
 * Added interfaceVersion and devkitID to DRServiceInterface_
 *
 * Revision 1.35.2.5  2005/01/11 00:03:06  dzhuravy
 * no message
 *
 * Revision 1.35.2.4  2005/01/04 19:43:21  dzhuravy
 * no message
 *
 * Revision 1.35.2.3  2005/01/04 14:55:53  dzhuravy
 * no message
 *
 * Revision 1.35.2.2  2004/11/24 23:29:02  smcguire
 * Fixed typo
 *
 * Revision 1.35.2.1  2004/11/09 21:07:01  dzhuravy
 * - Modified error handling
 * - Standardized primitive type names
 *
 * Revision 1.35  2004/03/10 17:01:35  vpaskhav
 * DRI_VERSION macro
 *
 * Revision 1.34  2004/03/09 16:28:29  aswain
 * move const char*s for array typenames to end of structure
 * gives a degree of backwards compatibility with older headers
 *
 * Revision 1.33  2004/02/27 11:06:48  aswain
 * update arrayNew comment
 *
 * Revision 1.32  2004/02/26 16:08:33  abayroot
 * Changed DRI_TYPE_ARRAY to DRI_TYPE_ARRAY_POSTFIX
 * Dropped DRIArrayType_typeName from reflection information
 * as it is redundant once we have the postfix policy on arrays.
 *
 * Revision 1.31  2004/02/25 20:32:33  abayroot
 * Added const char* DRI_TYPE_ARRAY.
 *
 * Revision 1.30  2004/02/25 17:36:46  aswain
 * added #defines for ALIB and EDG (for objectVersionGet)
 * added const char* to hold typenames for primitives (e.g. DRI_TYPE_INT)
 * added typeName to Array reflection element
 *
 *
 */

#ifndef DR_ANALYTICS_INTERFACE
#define DR_ANALYTICS_INTERFACE

#ifdef _MSC_VER
/* only for NT of course */
#ifdef DR_COMPILATION
/* DR_COMPILATION should only be defined when the source code is 
   being compiled */
#define DRI_DLL_EXPORT __declspec(dllexport)
#else
#define DRI_DLL_EXPORT __declspec(dllimport)
#endif
#else
/* on unix this define expands to nothing */
#define DRI_DLL_EXPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

struct DRServiceInterface_;

typedef struct DRService_ {
    const struct DRServiceInterface_ *fptr;
} DRService;

typedef struct DRObjectInterface_ {
    /*service pointer is available as long as object handle is valid */
    DRService *svc; 
} DRObjectInterface;


#define DR_TRUE  1
#define DR_FALSE 0

typedef int                 DRBool; 
typedef const char        * DRString;
typedef DRObjectInterface * DRObject;
typedef DRObject            DRMap;
typedef DRObject            DRArray;
typedef DRObject            DRMatrix;
typedef DRObject            DRMapIterator;
typedef DRString            DRError;

typedef enum DR_TYPE_
{
    DR_UNDEFINED = 0, 
    DR_BOOL,
    DR_INT,
    DR_DOUBLE,
    DR_STRING,
    DR_OBJECT,
    DR_DATE,
    DR_MAX_TYPE = DR_DATE
} DR_TYPE;

typedef struct DRDate_ {
  short year;
  char month;
  char day;
} DRDate;

typedef struct DRValue_ {
    int type; /* one of the constants from DR_TYPE */
    union {
        DRBool      boolean;
        int         integer;
        double      real;
        DRString    string;
        DRObject    object;
        DRDate      date;
    } value;
} DRValue;

/** How to start up a DRService */

/** The following macro converts a 4-tuple to an integer
    It is used to convert versions such as 3.6.0.0 to 
    an integer.
    @sa DRServiceInitArgs
*/
#define DRLIB_VERSION(a,b,c,d) (a<<24 | b<<16 | c<<8 | d)
#define DRI_VERSION DRLIB_VERSION(1,2,0,0)

/** these define the identity of development kits currently implementing DRI */
#define DRI_LIBRARY_ID_UNDEFINED 0
#define DRI_LIBRARY_ID_ALIB      1
#define DRI_LIBRARY_ID_EDG       6
#define DRI_LIBRARY_ID_QLIB      7
#define DRI_LIBRARY_ID_TS_QLIB   8
#define DRI_LIBRARY_ID_MAGNET    13

/** these define the "typeName" to use for arrays of primitive types */
#define DRI_TYPE_BOOL          "Bool"
#define DRI_TYPE_INT           "Int"
#define DRI_TYPE_DOUBLE        "Double"
#define DRI_TYPE_STRING        "String"
#define DRI_TYPE_DATE          "GDRDate"
#define DRI_TYPE_ARRAY_POSTFIX "Array"
#define DRI_TYPE_MATRIX        "DoubleMatrix"
#define DRI_TYPE_VARIANT       "IObject"

typedef enum DR_ERROR_
{
    DR_SUCCEEDED = 0,
    DR_FAILED,
    DR_INVALID_PARAMETER,
    DR_UNKNOWN_SERVICE_NAME,
    DR_INVALID_SERVICE_VERSION,
    DR_INVALID_INTERFACE_VERSION,
    DR_UNKNOWN_OPTION,
    DR_INVALID_OPTION,
    DR_MISSING_OPTION,
    DR_MAX_ERROR = DR_MISSING_OPTION
} DR_ERROR;

/** @brief Error handler function.
 *  @param userData    Arbitrary user data.
 *  @param errorString Error description.
 *  @param errorCode   Error code.
 *
 *  This function is only being called once during a call to DRCreateService.
 *
 *  Note that @p errorString can only be used during this call and should not be stored,
 *  as it may be deleted by the caller at any time.
 * 
 *  @sa DRCreateService
 */
typedef void (* DRErrorHandler)(
    void * userData,
    const char * errorString,
    DR_ERROR errorCode );

/** every option is a name/value pair 
    @sa DRServiceInitArgs */
typedef struct _DRServiceOption {
    char    *optionName;
    DRValue optionValue;
} DRServiceOption;

typedef struct _DRServiceInitArgs {
    /** serviceName is the name of the library to load. eg "fxCom" */
    const char*      serviceName; 
    /** serviceVersion is the version of the desired library.  
       Since library versions are in the form 3.6.0.0, use the 
       DRLIB_VERSION macro in this header to convert to int.  */
    int              serviceVersion;  
    /** interfaceVersion is the version of the DRServiceInterface. */
    int              interfaceVersion;

    /** Array of options (key/value pairs).  This allows additional
       information to be passed when loading specific DR Libraries */
    int              nOptions;
    DRServiceOption* options;

    /** true means unrecognized options are ignored
        false means unrecognized options cause failure */
    DRBool           ignoreUnrecognized;

    /** error handling function, can be NULL */
    DRErrorHandler errorHandler;
    /** arbitrary user data, which will be passed in errorHandler */
    void * userData;

} DRServiceInitArgs;

/** @brief Create service function. */
DRI_DLL_EXPORT DRBool DRCreateService(
    DRServiceInitArgs* args,
    DRService** service );

typedef DRBool (*DR_CREATE_PROC)(
    DRServiceInitArgs* args,
    DRService** service );

#define DR_CREATE_SERVICE "DRCreateService"

typedef void (* DRCreateServiceArgs)(
    void * userData,
    DRServiceInitArgs* args );

/**
   @defgroup interface DR Service Interface Functions 
   @{
 all values of type DRString have to be memory managed through DRService
 all values of type DRObject, DRArray, DRMap, DRMapIterator must be freed using objectFree
 order of freeing is independent of order of construction
 service must be alive to free objects created by it
 freeing a service does not free objects it created

 If succeeded all functions return NULL, If failed functions return an error description,
 which should be freed later with stringFree.

 Actual return values (if there are any) are passed through the last argument

*/
struct DRServiceInterface_ {

    /** @defgroup data Service Data 
        @{ */
    /** service data */        

    /** Version of the interface service supports */
    int              interfaceVersion;

    /** Identity of the development kit used to build service */
    int              devkitId;

    /** @} */

    /** @defgroup service Service Functions 
        @{ */
    /** service functions */        

    /** Returns service name and version */

    DRError (*serviceCreateArgs)(DRService* service,
                                 DRCreateServiceArgs createServiceArgs,
                                 void * userData );

    /** Returns description string of a service */
    DRError (*serviceDescription)(DRService* service, 
                                  DRString*  description);
    
    /** Shut down a DR service */
    DRError (*serviceFree)(DRService* service);
    
    /** Run an executable DRObject
     @sa serviceList */
    DRError (*execute)(DRService* service, 
                       DRObject   input,
                       DRValue*   output);

    /** @} */

    /** @defgroup object Object Functions
        @{ */
    /** object functions */

    /** Get the data type of an object */
    DRError (*objectGetType)(DRService* service, 
                             DRObject   object,
                             DRString*  typeName);

    /** Retrieves object unique identifier for this service instance */
    DRError (*objectGetId)(DRService* service,
                           DRObject   object,
                           void**     id);

    /** Free a DRObject */
    DRError (*objectFree)(DRService* service, 
                          DRObject   object);

    /** Does an object represent an array ?
        @sa arrayNew arrayLength arraySet arrayGet */
    DRError (*objectIsArray)(DRService* service,
                             DRObject object,
                             DRBool*  isArray);

    /** Does an object represent a matrix od doubles ? 
        @sa matrixNew matrixSet matrixGet matrixSize */
    DRError (*objectIsMatrix)(DRService* service,
                              DRObject object,
                              DRBool*  isMatrix);

    /** Can execute() be called on the object ?
        @sa execute */
    DRError (*objectIsExecutable)(DRService* service,
                                  DRObject object,
                                  DRBool*  isExecutable);
    
    /** Convert an object into a map where its component parts are
        accessible as key/value pairs
        @sa mapNew mapToObject mapAddItem mapGetItem mapIteratorGet mapIteratorNext
    */
    DRError (*objectToMap)(DRService* service,
                          DRObject   object, 
                          DRMap*     map );

    /** Convert a map of key/value pairs into an object
        @sa objectToMap objectGetType
    */
    DRError (*mapToObject)(DRService* service,
                          DRMap      map,
                          DRObject*  object);

    /** Returns version of the DR interface implementation that created the DRObject
        This function is internal and is supposed to be called from DR libraries only.

        When the library gets object handle, it checks whether the object is created with 
        the same DRService instance. If it's not, the library can call this function to
        determine the devkit that created the object. For example if it determines that
        object was created with ALIB devkit, it can access the object using ALIB object
        manager.

        devKitId identifies the library that implements DRObject (ALIB, EDG, CMLib etc.)
        devKitVersion is the version of that implementation.

        For libraries that do not implement multiversion object handles this function
        returns devKitId = 0 and devKitVersion = 0
    */
    DRError (*objectVersionGet)(DRService  *service,
                                DRObject   object,
                                int        *devKitId,
                                int        *devKitVersion); 
    /** @} */

    /** @defgroup string String Functions
        @{ */
    /** string functions */
    /** Create a new DRString from a C char* 
        @sa stringFree
    */
    DRError (*stringNew)(DRService* service,
                         const char* inString,
                         DRString*  item);
    
    /** Free a DRString
        @sa stringNew
    */
    DRError (*stringFree)(DRService* service, 
                          DRString   string);
    
    /** Free any DRString or DRObject in a DRValue */
    DRError (*valueClear)(DRService* service, 
                         DRValue    *value);
    /** @} */

    /** @defgroup array Array Functions
        @{ */
    /** array functions */

    /** Create an array of type 'typeName' and size numElements.
        @sa arraySet arrayGet arrayLength objectIsArray */
    DRError (*arrayNew)(DRService*  service, 
                        int         numElements,
                        const char* typeName,
                        DRArray*    array);

    /** get element type name
        The result is the same type name that was passed to arrayNew.
        For variant arrays result is 0.
        @sa arrayNew
    */
    DRError (*arrayElementType)(DRService* service,
                                DRArray    array,
                                DRString*  typeName);
    
    /** array access functions check for valid bounds */
    /** set the i'th item from an array (starts at 0)
        @sa arrayNew arrayGet arrayLength objectIsArray */
    DRError (*arraySet)(DRService*     service,
                        DRArray        array,
                        int            index,
                        const DRValue *item); 

    /** get the i'th item from an array (starts at 0)
        @sa arrayNew arraySet arrayLength objectIsArray */
    DRError (*arrayGet)(DRService* service,
                        DRArray    object,
                        int        index,
                        DRValue*   item);

    /** get the number of elements in the array 
        @sa arrayNew arraySet arrayGet objectIsArray */
    DRError (*arrayLength)(DRService* service,
                           DRArray    array, 
                           int*       length);
    /** @} */

    /** @defgroup matrix Matrix Functions
        @{ */
    /** matrix functions */

    /** create a new, empty matrix of doubles
        @sa matrixSet matrixGet matrixSize
    */
    DRError (*matrixNew)(DRService* service, 
                        int        numCols, 
                        int        numRows, 
                        DRMatrix*  grid);

    /** matrix access functions check for valid bounds */

    /** set the [i][j]th element of a matrix
        @sa matrixNew matrixGet matrixSize
    */       
    DRError (*matrixSet)(DRService* service, 
                        DRMatrix   matrix, 
                        int        col,
                        int        row,
                        double     value);

    /** get the [i][j]th element of a matrix
        @sa matrixNew matrixSet matrixSize
    */ 
    DRError (*matrixGet)(DRService* service, 
                        DRMatrix   matrix, 
                        int        col,
                        int        row,
                        double     *value);
    
    /** get the dimensions of a matrix 
        @sa matrixNew matrixSet matrixGet
    */        
    DRError (*matrixSize)(DRService* service, 
                         DRMatrix   matrix, 
                         int        *numCol,
                         int        *numRow);
    /** @} */

    /** map functions */
    /** 
        @defgroup map Map Functions
        @{
        A map contains key-value pairs of data used to represent the
        elements of classes. Maps can be converted into objects 
        and vice-versa 
    */

    /** Create a new (empty) map to represent the "typeName" class
        @param service  the DR service
        @param typeName name of the class that the new map represents
        @param map      the new map
        @sa objectToMap mapToObject mapAddItem mapGetItem mapGetTypeName 
        mapIteratorGet mapIteratorNext
    */
    DRError (*mapNew)(DRService*  service,
                      const char* typeName,
                      DRMap*      map);

    /* Calling mapNewFromValues has the same effect as
       calling mapNew and then calling mapAddItem, where 
       "fieldName" argument  takes first "numItems" field names
       of the type that is being created.
       This function is supported only by types marked as exported.
       To omit an optional value set type of corresponding DRValue
       to DR_UNDEFINED.
       @sa mapAddItem mapNew
    */

    DRError (*mapNewFromValues)(DRService  *service,
                                const char *typeName,
                                DRValue    *items,
                                int        numItems,
                                DRMap      *map);

    /* get a type id that can be used instead of a typename*/
    DRError (*idFromName)(DRService  *service,
                          const char *typeName,
                          void**     id);

    /* make a map from a type id instead of a typename */
    DRError (*mapNewFromId)(DRService  *service,
                            void*      typeId,
                            DRMap      *map);

    /*This function is supported only by types marked as exported*/
    DRError (*mapNewFromIdAndValues)(DRService  *service,
                                     void*      typeId,
                                     DRValue    *items,
                                     int        numItems,
                                     DRMap      *map);

    /** add a value to a map
     @sa mapNew mapGetItem mapGetTypeName mapIteratorGet mapIteratorNext
    */
    DRError (*mapAddItem)(DRService*     service,
                         DRMap          map,
                         const char*    fieldName,
                         const DRValue* valueToAdd);
    
    /** get a value out of a  map
     @sa mapNew mapAddItem mapGetTypeName mapIteratorGet mapIteratorNext
    */        
    DRError (*mapGetItem)(DRService*  service,
                         DRMap       map,
                         const char* fieldName,
                         DRValue*    item);

    /** get type of the object that is represented by the map
    */
    DRError (*mapGetTypeName)(DRService  *service,  
                             DRMap         map,
                             DRString   *typeName);


    /** get array of values from the map.
        Values are returned in the same order as corresponding 
        data members in the type represnted by the map.
        Omitted values have type DR_UNDEFINED.
        This function is supported only by types marked as exported.
        Calling this function for other types will set *array to 0.

    */
    DRError (*mapGetValues)(DRService  *service,  
                                       DRMap         map,
                                       DRArray    *array);

    /** Creates a map iterator (which must eventually be freed)
        This object can then be used to explore the contents of a map
     @sa objectToMap mapNew mapAddItem mapGetTypeName mapIteratorNext
    */
    DRError (*mapIteratorGet)(DRService*     service,
                             DRMap          map, 
                             DRMapIterator* iterator);
    
    /** Returns next <field, object> pair in iteration 
        field and object will be null when there are no more elements 
        in the iteration. 
     @sa objectToMap mapNew mapAddItem mapGetTypeName mapIteratorGet
    */
    DRError (*mapIteratorNext)(
        DRService*    service,
        DRMapIterator iter, 
        DRString*     field,   
        DRValue*      object);
    /** @} */

    /** @defgroup type Type Information
        @{ */ 
    /** type information */
    /** Return a list of all data types that can be created via the service
        i.e. classes built via maps
        @sa typeByName serviceList serviceByName mapToObject
    */
    DRError (*typeList)(DRService* service,
                       DRArray*   types);

   /** This returns an object representing the 
       type description of the "typename" class.
        @sa typeList serviceList serviceByName mapToObject
    */
    DRError (*typeDescription)( DRService*  service,
                               const char* typeName,
                               DRObject*   type);

    /** @} */
};
/** @} */

#ifdef __cplusplus
}
#endif

/* abstract type DRIType */

/*
    struct DRIScalarType : DRIType {
        int typeId; //DR_TYPE constant where 0 means Variant
    }
*/
#define DRIScalarType__name "DRIScalarType" 
#define DRIScalarType_typeId "typeId" 

/*
    struct DRIArrayType : DRIType {
        DRIType *elementType;
        }
*/
#define DRIArrayType__name "DRIArrayType" 
#define DRIArrayType_elementType "elementType" 

/*
    struct DRIMatrixType : DRIType {
    }
*/
#define DRIMatrixType__name "DRIMatrixType"
 
/*
    struct DRIEnumerationType {
        string typeName;
        string description;
        vector<string> identifiers;
        vector<string> identifierDescriptions;
    }
*/
#define DRIEnumerationType__name "DRIEnumerationType"
#define DRIEnumerationType_typeName "typeName"
#define DRIEnumerationType_description "description"
#define DRIEnumerationType_identifiers "identifiers"
#define DRIEnumerationType_identifierDescriptions "identifierDescriptions"

/*
    struct DRIAbstractType : DRIType 
    {
        string typeName;
        vector<DRIType*> baseTypes; // direct base classes (one level up)
        string description;
    }
*/
#define DRIAbstractType__name "DRIAbstractType" 
#define DRIAbstractType_typeName "typeName"
#define DRIAbstractType_baseTypes "baseTypes"
#define DRIAbstractType_description "description"

/*
    struct DRIConcreteType : DRIAbstractType 
    {
        vector<DRIArgumentDescriptor*> members;
        bool isExported;
    }
*/
#define DRIConcreteType__name "DRIConcreteType" 
#define DRIConcreteType_typeName "typeName"
#define DRIConcreteType_baseTypes "baseTypes"
#define DRIConcreteType_description "description"
#define DRIConcreteType_members "members" 
#define DRIConcreteType_isExported "isExported" 

/*
    struct DRIExecutableType : DRIConcreteType 
    {
        DRIType* resultType;
        vector<DRIResultDescriptor*> outArguments;
    }
*/
#define DRIExecutableType__name "DRIExecutableType" 
#define DRIExecutableType_typeName "typeName"
#define DRIExecutableType_baseTypes "baseTypes"
#define DRIExecutableType_description "description"
#define DRIExecutableType_members "members" 
#define DRIExecutableType_isExported "isExported" 
#define DRIExecutableType_resultType "resultType"
#define DRIExecutableType_outArguments "outArguments"

/*
    struct DRIArgumentDescriptor {
        string name;
        DRIType *type;
        string description;
        bool isOptional;
  
        int groupIndex;         // ==0 unless it is FX library
        int alternationIndex;   // ==0 unless it is FX library
        int indexInAlternation; // ==index of the argument unless it is FX library
    }

*/
#define DRIArgumentDescriptor__name "DRIArgumentDescriptor"
#define DRIArgumentDescriptor_name "name"
#define DRIArgumentDescriptor_type "type"
#define DRIArgumentDescriptor_description "description"
#define DRIArgumentDescriptor_isOptional "isOptional"
#define DRIArgumentDescriptor_groupIndex "groupIndex"
#define DRIArgumentDescriptor_alternationIndex "alternationIndex"
#define DRIArgumentDescriptor_indexInAlternation "indexInAlternation"

/*    struct DRIResultDescriptor {
        string name;
        DRIType *type;
        string description;
    }
*/
#define DRIResultDescriptor__name "DRIResultDescriptor"
#define DRIResultDescriptor_name "name"
#define DRIResultDescriptor_type "type"
#define DRIResultDescriptor_description "description"

#endif
