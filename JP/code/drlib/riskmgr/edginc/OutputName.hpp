

#ifndef EDG_OUTPUT_NAME_H
#define EDG_OUTPUT_NAME_H
#include <string>
#include "edginc/Object.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/Array.hpp"

DRLIB_BEGIN_NAMESPACE

class OutputName;

typedef smartConstPtr<OutputName> OutputNameConstSP;
typedef smartPtr<OutputName> OutputNameSP;
typedef array<OutputNameSP, OutputName> OutputNameArray;
typedef smartPtr<OutputNameArray> OutputNameArraySP;
typedef smartConstPtr<OutputNameArray> OutputNameArrayConstSP;

#ifndef QLIB_OUTPUTNAME_CPP
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<OutputName>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<OutputName>);
EXTERN_TEMPLATE(class RISKMGR_DLL array<OutputNameSP _COMMA_ OutputName>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartConstPtr<OutputNameArray>);
EXTERN_TEMPLATE(class RISKMGR_DLL_SP smartPtr<OutputNameArray>);
EXTERN_TEMPLATE(IObjectSP RISKMGR_DLL FieldGetSmartPtr<OutputNameSP>(OutputNameSP* t));
EXTERN_TEMPLATE(void RISKMGR_DLL FieldSetSmartPtr<OutputNameSP>(OutputNameSP* t,
                                                    IObjectSP o));
#else
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<OutputName>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<OutputName>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL array<OutputNameSP _COMMA_ OutputName>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartConstPtr<OutputNameArray>);
INSTANTIATE_TEMPLATE(class RISKMGR_DLL smartPtr<OutputNameArray>);
INSTANTIATE_TEMPLATE(IObjectSP RISKMGR_DLL FieldGetSmartPtr<OutputNameSP>(OutputNameSP* t));
INSTANTIATE_TEMPLATE(void RISKMGR_DLL FieldSetSmartPtr<OutputNameSP>(OutputNameSP* t,
                                                         IObjectSP o));
#endif

/** Captures the name needed to identify an output - can depend upon more 
    than one identifier eg equity v equity correlations. Output names are
    immutable for performance reasons and ease of use */
class RISKMGR_DLL OutputName: public CObject{
public:
    friend class OutputNameHelper;
    static CClassConstSP const TYPE;

    /** Creates a new output name using the given string */
    OutputName(const string& name1);

    /** Creates a new output name using the given pair of strings */
    OutputName(const string& name1, const string& name2);
    
    /** Creates a new output name using the two output names. */
    OutputName(const OutputName* output1, const OutputName* output2);

    /** Creates a new output name using the two output names. */
    static OutputNameSP SP(OutputNameConstSP output1,
                           OutputNameConstSP output2);

    /** Creates a new output name using the given string. */
    static OutputNameSP SP(const string& name);

    /** Destructor */
    ~OutputName();

    /** returns a hashcode for the object based upon the strings making up
        the name */
    virtual int hashCode() const;

    /** Returns true if the this and the given output name are identical */
    bool equals(const OutputName* name) const;

    /** Returns true if the this name identity is given by a single string
        which matches the given one */
    bool equals(const string& val) const;

    /** Returns true if the this name identity is given by two strings
        which matche the one provided */
    bool equals(const string& val1, const string val2) const;

    /** formatted version of outputname - ignores any alias */
    string toString() const;

    /** Returns true if output name is just "" */
    bool isEmpty() const;

    /** Returns the number of strings which identifies this output name */
    int idCount() const;

    /** Returns the component string, identified by index, which
        identifies this output name */
    const string& idGet(int index) const;

    /** boil an array of name pairs down into a list of single names */
    static OutputNameArraySP singleNameArray(OutputNameArrayConstSP names);

    /** is an array of names the same name again & again ? */
    static bool uniqueNameArray(OutputNameArrayConstSP names);
    
    /** remove duplicates and empty names from an array. Order is otherwise
        unchanged */
    static OutputNameArraySP trim(OutputNameArrayConstSP names);

    /** Returns the intersection of the two arrays ie the names that appear
        in both lists */
    static OutputNameArraySP intersection(OutputNameArrayConstSP set1,
                                          OutputNameArrayConstSP set2);

    /** Returns the 'difference' of the two arrays ie the names that appear
        in the first list but not in the second */
    static OutputNameArraySP difference(OutputNameArrayConstSP namesToInclude,
                                        OutputNameArrayConstSP namesToExclude);

    /** Appends to an array */
    static void append(OutputNameArray& names,
                       OutputNameArrayConstSP namesToAppend);

    /* class needed by hash_map template */
    class RISKMGR_DLL HashUtil{
    public:
        bool operator()(const OutputNameConstSP& name1, 
                        const OutputNameConstSP& name2) const;
        int operator()(const OutputNameConstSP& name) const;
    };
        

private:
    // could make more general by using string lists
    string name1;
    // following field is optional
    string* name2;
    /* used for formatting */
    static string const DELIMITER;
    /** for reflection purposes */
    OutputName();

};

// for those who like a C at the front
typedef OutputName COutputName;
typedef OutputNameArray COutputNameArray;
typedef OutputNameArraySP COutputNameArraySP;
typedef OutputNameArrayConstSP COutputNameArrayConstSP;

DRLIB_END_NAMESPACE

#endif

