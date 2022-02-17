
#ifndef EDG_ATOMICARRAY_H
#define EDG_ATOMICARRAY_H
#include "edginc/Array.hpp"

DRLIB_BEGIN_NAMESPACE
/** specialisation of arrayObjectCast for arrays of doubles */
template <> class TOOLKIT_DLL arrayObjectCast<double>{
public:
    /** Casts array element to an IObject */
    static IObjectSP toIObject(double value);

    /** Sets the value of the indexed component of the specified
        array object to the specified new value. */
    static double fromIObject(IObjectSP value);
};

template <> class TOOLKIT_DLL arrayClone<double>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};
   
template <> class TOOLKIT_DLL arrayCompare<double>{
public:
    /** Overridden for performance */
    static int hashCode(const CArray* arrayToHash);
    static bool equalTo(const CArray* arrayToCompare, const IObject* obj);
};
   
/** specialisation of arrayObjectCast for arrays of ints */
template <> class TOOLKIT_DLL arrayObjectCast<int>{
public:
    /** Casts array element to an IObject */
    static IObjectSP toIObject(int value);

    /** Sets the value of the indexed component of the specified
        array object to the specified new value. */
    static int fromIObject(IObjectSP value);
};

template <> class TOOLKIT_DLL arrayClone<int>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

template <> class TOOLKIT_DLL arrayCompare<int>{
public:
    /** Overridden for performance */
    static int hashCode(const CArray* arrayToHash);
    static bool equalTo(const CArray* arrayToCompare, const IObject* obj);
};
   
/** specialisation of arrayObjectCast for arrays of bools */
template <> class TOOLKIT_DLL arrayObjectCast<bool>{
public:
    /** Casts array element to an IObject */
    static IObjectSP toIObject(bool value);

    /** Sets the value of the indexed component of the specified
        array object to the specified new value. */
    static bool fromIObject(IObjectSP value);
};

template <> class TOOLKIT_DLL arrayCompare<bool>{
public:
    /** Overridden for performance */
    static int hashCode(const CArray* arrayToHash);
    static bool equalTo(const CArray* arrayToCompare, const IObject* obj);
};
   
/** specialisation of arrayObjectCast for arrays of strings */
template <> class TOOLKIT_DLL arrayObjectCast<string>{
public:
    /** Casts array element to an IObject */
    static IObjectSP toIObject(const string& value);

    /** Sets the value of the indexed component of the specified
        array object to the specified new value. */
    static string fromIObject(IObjectSP value);
};

template <> class TOOLKIT_DLL arrayClone<string>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

template <> class TOOLKIT_DLL arrayCompare<string>{
public:
    /** Overridden for performance */
    static int hashCode(const CArray* arrayToHash);
    static bool equalTo(const CArray* arrayToCompare, const IObject* obj);
};

/** template for wrapping stl vector to form class derived from CArray */
template <> class TOOLKIT_DLL array<bool, bool> : public CArray{
public:
#if defined(_MSC_VER) || \
    (defined(__GNUC__) && ((__GNUC__ == 3 && __GNUC_MINOR__ >= 2)) || __GNUC__ > 3)
    typedef unsigned char MY_BOOL;
#else
    typedef unsigned int MY_BOOL;
#endif
public:
    typedef arrayObjectCast<bool> eltCast;
    static CClassConstSP const TYPE;

    ~array();

    explicit array();

    explicit array(int length);
        
    explicit array(int n, bool t);

    /** override default 'datadict' based copy */
    virtual IObject* clone() const;

    /** Returns a hash code value for the object by XORing the hashCode for
        each element in the array */
    virtual int hashCode() const;

    /** Indicates whether some other object is "equal to" this one by comparing
        each element in the array */
    virtual bool equalTo(const IObject* obj) const;
        
    bool equals(const array<bool, bool>& rhs) const;

    /** Returns the value of the indexed component in the
        specified array object. */
    virtual IObjectConstSP get(int index) const;

    virtual IObjectSP get(int index);

    /** Returns the length of the specified array object, as an int. */
    virtual int getLength() const;

    /** Sets the value of the indexed component of the specified
        array object to the specified new value. */
    virtual void set(int index, IObjectSP value);

    /** Appends the supplied object to the end of supplied array. The
        length of the array is increased by 1. (aka push_back).  */
    virtual void append(IObjectSP value);

    // from stl
    bool& operator[](int index);

    const bool& operator[](int index) const;

    int size() const;

    bool empty() const;

    void push_back(const bool& val);

    void pop_back();

    vector<MY_BOOL>::iterator begin();

    vector<MY_BOOL>::iterator end();

    vector<MY_BOOL>::const_iterator begin() const;

    vector<MY_BOOL>::const_iterator end() const;

    bool& front();

    bool& back();

    const bool& front() const;

    const bool& back() const;

    void reserve(int n);

    int capacity() const;

    void clear();

    void resize(int n, bool t = bool());

private:
    // in stl vector<bool> is stored as an array of bits - which makes 
    // everything a complete nightmare. 
    vector<MY_BOOL> theArray; // $unregistered

    // type registration
    static IArray* createArray(int length);

    static void load(CClassSP& clazz);
};
#if !defined(DEBUG) || defined(QLIB_ATOMICARRAY_CPP)
// avoid cost of including it for debug
#include "edginc/AtomicArray.inl"
#endif

typedef array<double, double> CDoubleArray;
typedef array<int, int> CIntArray;
typedef array<bool, bool> CBoolArray;
typedef array<string, string> CStringArray;
// typedefs for a smart pointer to each of these arrays
typedef smartPtr<CDoubleArray> CDoubleArraySP;
typedef smartConstPtr<CDoubleArray> CDoubleArrayConstSP;
typedef smartPtr<CIntArray> CIntArraySP;
typedef smartPtr<CStringArray> CStringArraySP;
typedef smartConstPtr<CStringArray> CStringArrayConstSP;
typedef smartPtr<CBoolArray> CBoolArraySP;
typedef smartConstPtr<CBoolArray> CBoolArrayConstSP;

typedef CDoubleArray          DoubleArray;
typedef CDoubleArraySP        DoubleArraySP;
typedef CDoubleArrayConstSP   DoubleArrayConstSP;
typedef CIntArray             IntArray;
typedef CIntArraySP           IntArraySP;
typedef CBoolArray            BoolArray;
typedef CBoolArraySP          BoolArraySP;
typedef CBoolArrayConstSP     BoolArrayConstSP;
typedef CStringArray          StringArray;
typedef CStringArrayConstSP   StringArrayConstSP;


/* Computing of a sum of 2 arrays. Sizes must match. */
DoubleArraySP operator+(const DoubleArray & g1, const DoubleArray & g2);

/* Adding of g2 to g1. Sizes must match. */
void operator +=(DoubleArray & g1, const DoubleArray & g2);

/* Multiplying of gArray (element by element) by gDouble */
DoubleArraySP operator*(double gDouble, const DoubleArray & gArray);

/* Multiplying of gArray (element by element) by gDouble */
DoubleArraySP operator*(const DoubleArray & gArray, double gDouble);

/* Adding of g2 to g1. Sizes must match. */
void operator *=(DoubleArray & gArray, double gDouble);


/* Array of SPs to IntArrays */
typedef smartPtr<IntArray> IntArraySP;
typedef smartConstPtr<IntArray> IntArrayConstSP;

typedef array<IntArraySP, IntArray> IntArrayArray;

typedef smartPtr<IntArrayArray> IntArrayArraySP;
typedef smartConstPtr<IntArrayArray> IntArrayArrayConstSP;

// Array of strings and array of string arrays
typedef smartPtr<StringArray> StringArraySP;
typedef smartConstPtr<StringArray> StringArrayConstSP;

typedef array<StringArraySP, StringArray> StringArrayArray;

typedef smartPtr<StringArrayArray> StringArrayArraySP;
typedef smartConstPtr<StringArrayArray> StringArrayArrayConstSP;


// Array of double arrays
typedef array<DoubleArray> DoubleArrayArray;

// Array of DoubleArrayArray's (i.e. 3-d arrays).  Useful for debugging credit 
typedef array<DoubleArrayArray> DoubleArrayArrayArray;

/** specialisations of arrayObjectCast (needed as the array is not an array
    of pointers) */
template <> class TOOLKIT_DLL arrayObjectCast<DoubleArray>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const DoubleArray& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(DoubleArray& value);

    /** Turns the IObjectSP into a DateTime */
    static const DoubleArray& fromIObject(IObjectSP& value);
};

/** specialisations of arrayObjectCast for DoubleArrayArrayArray */
template <> class TOOLKIT_DLL arrayObjectCast<DoubleArrayArray>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const DoubleArrayArray& value);

    /** Casts array element to an IObject */
    static IObjectSP toIObject(DoubleArrayArray& value);

    /** Turns the IObjectSP into a DateTime */
    static const DoubleArrayArray& fromIObject(IObjectSP& value);
};

/** specialisation of arrayClone */
template <> class TOOLKIT_DLL arrayClone<DoubleArray>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

/** specialisation of arrayClone */
template <> class TOOLKIT_DLL arrayClone<DoubleArrayArray>{
public:
    /** Overridden for performance */
    static IObject* clone(const CArray* arrayToClone);
};

typedef smartPtr<DoubleArrayArray> DoubleArrayArraySP;
typedef smartConstPtr<DoubleArrayArray> DoubleArrayArrayConstSP;

typedef smartPtr<DoubleArrayArrayArray> DoubleArrayArrayArraySP;
typedef smartConstPtr<DoubleArrayArrayArray> DoubleArrayArrayArrayConstSP;

#ifndef QLIB_ATOMICARRAY_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL array<double _COMMA_ double>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<int _COMMA_ int>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<string _COMMA_ string>);

EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<DoubleArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<DoubleArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<IntArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<IntArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<StringArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<StringArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<BoolArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<BoolArray>);

EXTERN_TEMPLATE(class TOOLKIT_DLL array<DoubleArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<DoubleArrayArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<DoubleArrayArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<DoubleArrayArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<DoubleArrayArrayArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<DoubleArrayArrayArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<StringArraySP _COMMA_ StringArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<IntArraySP _COMMA_ IntArray>);

EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<IntArrayArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<IntArrayArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<StringArrayArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<StringArrayArray>);
// again, avoid code bloat by declaring common implementation of function
// templates as extern
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<DoubleArray>(DoubleArray* t));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<IntArray>(IntArray* t));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<StringArray>(StringArray* t));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetInLine<BoolArray>(BoolArray* t));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<DoubleArray>(DoubleArray* t, IObjectSP o));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<IntArray>(IntArray* t, IObjectSP o));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<StringArray>(StringArray* t, IObjectSP o));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetInLine<BoolArray>(BoolArray* t, IObjectSP o));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<CDoubleArraySP>(CDoubleArraySP* t));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<CIntArraySP>(CIntArraySP* t));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<CStringArraySP>(CStringArraySP* t));
EXTERN_TEMPLATE(IObjectSP TOOLKIT_DLL FieldGetSmartPtr<CBoolArraySP>(CBoolArraySP* t));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<CDoubleArraySP>(CDoubleArraySP* t,
                                                     IObjectSP o));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<CIntArraySP>(CIntArraySP* t,
                                                   IObjectSP o));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<CStringArraySP>(CStringArraySP* t,
                                                     IObjectSP o));
EXTERN_TEMPLATE(void TOOLKIT_DLL FieldSetSmartPtr<CBoolArraySP>(CBoolArraySP* t, 
                                                    IObjectSP o));
#endif

DRLIB_END_NAMESPACE
#endif

