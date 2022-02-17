//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Handle.hpp
//
//   Description : Class for holding handles - which are string references
//                 to objects
//
//   Author      : Mark A Robson
//
//   Date        : 21 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_HANDLE_HPP
#define EDG_HANDLE_HPP

#include "edginc/Hashtable.hpp"
#include "edginc/Array.hpp"

DRLIB_BEGIN_NAMESPACE

/** Holds all handles in use. Handles are essentially string
    references to objects. They are primarily used so that pointers
    can be 'returned' to the spreadsheet - and provide for type safety */
class TOOLKIT_DLL Handle: public CObject {
public:
    friend class HandleHelper;
    static CClassConstSP const TYPE;

    /** max handle id. The id part is the digits at the start of a handle */
    static const int MAX_ID; 

    /** number of character needed for MAX_ID */
    static const int MAX_ID_LENGTH; 

    /** Method which returns a string uniquely identifying the cell
        on the spreadsheet. (In VB, this just returns the same string) */
    typedef string (SpreadSheetCellNameMethod)();

    /** Sets the method used for finding the current cell in the spreadsheet.
        Method can be null (default), in which case "NULL" is used */
    static void setSpreadSheetCellNameMethod(
        SpreadSheetCellNameMethod* sSheetMethod);

    /** create and store a handle to an object. Returns complete handle name */
    static string create(const string&     userName,
                         const IObjectSP&  object);

    /** Returns handle with given name. objectClass, if non null, is used
        to validate the type of the object */
    static IObjectSP fetch(const string&  name,
                           CClassConstSP  objectClass); // can be null 

    /** Returns true if given string has the right format to be a handle. */
    static bool valid(const string& name);

    /** Returns true if given string has the right format to be a handle. */
    static bool valid(const char* countedStringL);

    /** Returns true if the handle with the given name exists */
    static bool exists(const string& name);

    /** Returns the user part of a full handle name eg for 0001Divs.[Sheet1!A1]
        returns Divs */
    static string getUserHandleName(const string& fullHandle);

    /** Removes the specified handle from storage */
    static void destroy(const string&  name);

    /** Removes all handles from storage */
    static void destroyAll();

    /** Addin function classes can implement this interface in order to
        supply a default handle name. The method is invoked if the addin
        function returns a handle but has no handle name parameter or if
        the addin parameter (from the spreadsheet) is empty */
    class TOOLKIT_DLL IName{
    public:
        /** Return a default handle name when this class is used for an
            addin function */
        virtual string defaultHandleName() const = 0;
        virtual ~IName();
    };

    /** Class for accessing raw strings (ie accessing handle names without
        the object being pulled from the handle) */
    class TOOLKIT_DLL RawString: public CObject{
    public:
        static CClassConstSP const TYPE;
        /** returns the string - which could be a handle name */
        const string& getString() const;
        /** constructor */
        RawString(const string& rawString);
    private:
        // field
        string rawString;

        static void load(CClassSP& clazz);
        static IObject* defaultCreate();
        RawString();
    };
    typedef smartPtr<RawString> RawStringSP;
    typedef array<RawStringSP, RawString> RawStringArray;
    typedef smartPtr<RawStringArray> RawStringArraySP;
    
private:
    class FromObject;
    class ToObject;
    static const int MIN_HANDLE_LENGTH; // min length for a handle
    static const int MAX_HANDLE_LENGTH; // max length of a handle
    static SpreadSheetCellNameMethod* sheetNameMethod;
    static HashtableSP handles;     // contains all the handles

    int        id;       // The last unique identifier - avoid stale handles $unregistered
    IObjectSP  object;   // what the handle contains $unregistered
    char       idAsString[5]; /* holds id as a string - compiler problems // $unregistered
                                 if MAX_ID_LENGTH +1 is used */
    
    Handle(const IObjectSP& object);
};

typedef smartPtr<Handle> HandleSP;

#ifndef QLIB_HANDLE_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<Handle>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<Handle>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<Handle::RawString>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<Handle::RawStringSP _COMMA_ Handle::RawString>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<Handle::RawStringArray>);
#endif

DRLIB_END_NAMESPACE

#endif
