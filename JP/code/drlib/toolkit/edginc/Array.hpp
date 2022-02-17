
#ifndef EDR_ARRAY_HPP
#define EDR_ARRAY_HPP
#include "edginc/Class.hpp"
DRLIB_BEGIN_NAMESPACE

/** Interfaces for arrays. Nearly all arrays derive from CArray but the
    separation of CArray from IArray allows other implementations (eg XArray) */
class TOOLKIT_DLL IArray: public virtual IObject{
public:
    /** attribute for recording an array's length */
    static const string ARRAY_LENGTH;

    static CClassConstSP const TYPE;

    //// constructors
    IArray();
    IArray(const IArray& rhs);

    /** Returns the value of the indexed component in the
        specified array object. */
    virtual IObjectConstSP get(int index) const = 0;

    /** Non const version of above */
    virtual IObjectSP get(int index) = 0;

    /** Returns the length of the specified array object, as an int. */
    virtual int getLength() const = 0;

    /** Sets the value of the indexed component of the specified
        array object to the specified new value */
    virtual void set(int index, IObjectSP value) = 0;

    /** Appends the supplied object to the end of supplied array. The
        length of the array is increased by 1. (aka push_back).  */
    virtual void append(IObjectSP value) = 0;

    virtual ~IArray();

private:
    static void load(CClassSP& classToLoad);
};

#if !defined(DEBUG) || defined(QLIB_ARRAY_CPP)
//// in line for performance
OPT_INLINE IArray::IArray(){}
OPT_INLINE IArray::IArray(const IArray& rhs) : IObject() {}
OPT_INLINE IArray::~IArray(){} 
#endif

//   Provides wrapper around stl vector template to give object like
//   functionality. The array is an array of non constant objects. An array
//   of constant objects might prove necessary

/** base class for array wrapper classes.  */
class TOOLKIT_DLL CArray: public CObject,
              public virtual IArray {
public:
    static CClassConstSP const TYPE;

    //// constructor
    CArray(const CArray& rhs);

    /** override default 'datadict' based copy */
    virtual IObject* clone() const;

    /** write object out to Writer */
    virtual void write(const string& tag, Writer* writer) const;

    /** populate an empty object from Reader */
    virtual void import(Reader::Node* elem, Reader* reader);

    /** write object out in 'output' format - ie suitable for comparing
        regression files with */
    virtual void outputWrite(const string& linePrefix,
                             const string& prefix, ostream& stream) const;

    /** Returns a hash code value for the object by XORing the hashCode for
        each element in the array */
    virtual int hashCode() const;

    /** Indicates whether some other object is "equal to" this one by comparing
        each element in the array */
    virtual bool equalTo(const IObject* obj) const;

    virtual ~CArray();

    /** same as clone but a static method - used in default clone
        template below */
    static IObject* copyArray(const CArray* arrayToClone);

    /** utility method */
    static int hashArray(const CArray* arrayToHash);
    /** utility method */
    static bool equalToArray(const CArray* arrayToCompare, const IObject* obj);
protected:
    CArray(const CClassConstSP& objClass);

private:
    static void load(CClassSP& classToLoad);
        
};

#if !defined(DEBUG) || defined(QLIB_ARRAY_CPP)
//// in line for performance
//// in line for performance
OPT_INLINE CArray::CArray(const CArray& rhs): IObject(), IArray(), CObject(rhs.getClass()) {}
OPT_INLINE CArray::~CArray(){} 
OPT_INLINE CArray::CArray(const CClassConstSP& objClass):
    CObject(objClass){}
#endif

/** Template class allows a common view onto elements of arrays
    regardless of the way in which the data is stored in the array.
    The default implementation is for array of smart pointers. Arrays where
    the data is held as an array of structures (eg DoubleArray) need to
    provide a specialised version of this template */
template <class X> class arrayObjectCast{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const X& value){
        IObjectConstSP objValue(value);
        return objValue;
    }
    
    /** Casts array element to an IObject */
    static IObjectSP toIObject(X& value){
        return value;
    }

    /** Sets the value of the indexed component of the specified
        array object to the specified new value. */
    static X fromIObject(IObjectSP& value){
        X  newVal(X::dynamicCast(value));
        return newVal;
    }
};

/** Template class for cloning arrays - allows array of structures (rather
    than arrays of pointers) to override this for performance reasons. ie
    can use stl copy constructor for cloning */
template <class X> class arrayClone{
public:
    /** Default clone method for arrays */
    static IObject* clone(const CArray* arrayToClone){
        return CArray::copyArray(arrayToClone);
    }
};

/** Template class for hashing/comparing arrays - allows specific arrays to
    override this for performance reasons. */
template <class X> class arrayCompare{
public:
    /** Default hashCode method for arrays */
    static int hashCode(const CArray* arrayToHash){
        return CArray::hashArray(arrayToHash);
    }
    /** Default equalTo method for arrays */
    static bool equalTo(const CArray* arrayToCompare, const IObject* obj){
        return CArray::equalToArray(arrayToCompare, obj);
    }
};

/** specialisation of arrayObjectCast for arrays of IObjectSP */
template <> class TOOLKIT_DLL arrayObjectCast<IObjectSP>{
public:
    /** Casts array element to an IObject */
    static IObjectConstSP toIObject(const IObjectSP& value);
    
    /** Casts array element to an IObject */
    static IObjectSP toIObject(IObjectSP& value);
    
    /** Sets the value of the indexed component of the specified
        array object to the specified new value. */
    static IObjectSP fromIObject(IObjectSP& value);
};

/** template for wrapping stl vector to form class derived from CArray */
template <class X, class component = X> class array: public CArray{
private:
    typedef arrayObjectCast<X> eltCast;
    
public:
    static CClassConstSP const TYPE;
    typedef	                X								       value_type;
    typedef	          const X&						          const_reference;
    typedef typename vector<X>::const_iterator                 const_iterator;
    typedef typename vector<X>::iterator                             iterator;
    typedef typename vector<X>::reverse_iterator             reverse_iterator;
    typedef typename vector<X>::const_reverse_iterator const_reverse_iterator;
    
    explicit array();

    explicit array(int length);
    
    explicit array(int n, const X& t);

    array(const array<X, component>& rhs);

    array<X, component>& operator=(const array<X, component>& rhs);

    static smartPtr<array> SP();

    static smartPtr<array> SP(int length);

    static smartPtr<array> SP(int n, const X& t);

    bool equals(const array<X, component>& rhs) const;

    /** creates vector from a range given by two iterators */
    explicit array(const_iterator start, 
                   const_iterator end);
    
    ~array();

    /** override default 'datadict' based copy */
    virtual IObject* clone() const;

    /** Returns a hash code value for the object by XORing the hashCode for
        each element in the array */
    virtual int hashCode() const;

    /** Indicates whether some other object is "equal to" this one by comparing
        each element in the array */
    virtual bool equalTo(const IObject* obj) const;

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
    X& operator[](int index);

    const X& operator[](int index) const;

    int size() const;

    bool empty() const;

    void push_back(const X& val);

    void pop_back();

    iterator begin();

    iterator end();

    const_iterator begin() const;

    const_iterator end() const;

    reverse_iterator rbegin();

    reverse_iterator rend();

    const_reverse_iterator rbegin() const;

    const_reverse_iterator rend() const;

    X& front();

    X& back();

    const X& front() const;

    const X& back() const;

    void reserve(int n);

    int capacity() const;

    void clear();

    void resize(int n, X t = X());

    iterator erase(iterator pos);

    iterator erase(iterator first, iterator last);

    iterator insert(iterator pos, const X& x);

    /* two definitions of this depending on whether member templates
       work or not - apparently they don't in MSVC */
    void insert(iterator       __position,
                const_iterator __first,
                const_iterator __last);

    /** getSubarrayGivenPlaces returns a subarray given the indexes of the elements
    of subarray in the array 'this'. An error will be thrown if at least one
    of the elements of gPlaces is out of range [ 0, this->size()-1 ] */
    smartPtr< array<X, component> > 
        getSubarrayGivenPlaces(const array<int> & gPlaces) const;

    //// just does what you'd get if you did '==' with 2 vectors
    bool operator==(const array<X, component>& rhs) const;

    // g++ barfs in Array.inl definition of back_inserter if you write out the
    // type in full

    typedef back_insert_iterator<vector<X> > back_insert_iterator_vector_X;

	back_insert_iterator<vector<X> > back_inserter();

    // to do - others

private:
    vector<X> theArray; // $unregistered
    // type registration

    static IArray* createArray(int length);
    DLL_FIX_FOR_TEMPLATE_TYPE; // work around for VC71 and dlls
    static void load(CClassSP& clazz);
    
};

// then include templated definitions
#include "edginc/Array.inl"

// typedefs for smartPtrs
typedef smartConstPtr<IArray> IArrayConstSP;
typedef smartPtr<IArray> IArraySP;

// typedefs for an array of IObjectSP's
typedef array<IObjectSP, IObject> ObjectArray;
typedef smartPtr<ObjectArray> ObjectArraySP;
typedef smartConstPtr<ObjectArray> ObjectArrayConstSP;
#ifndef QLIB_ARRAY_CPP
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartConstPtr<IArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<IArray>);
EXTERN_TEMPLATE(class TOOLKIT_DLL array<IObjectSP _COMMA_ IObject>);
EXTERN_TEMPLATE(class TOOLKIT_DLL_SP smartPtr<ObjectArray>);
#endif
DRLIB_END_NAMESPACE
#endif

