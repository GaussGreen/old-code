
#ifndef EDG_SMARTPTR_H
#define EDG_SMARTPTR_H

#include "edginc/config.hpp"
#include "edginc/ModelException.hpp"

//#include "edginc/refCountPtr.hpp"

DRLIB_BEGIN_NAMESPACE
class IObject;
class CClass;
typedef const CClass* CClassConstSP;

/** This implementation of a smart pointer has been borrowed from the stl
    implementation of the auto_ptr. All the additional classes and casting
    operators allow casting to/from smart pointers to derived classes */
class TOOLKIT_DLL sp_base {
public:
    void* p;
    void  __set(const void* p);
    void  __set(void* p);
};

template <class _Tp> class sp_ref {
public:
    sp_base&   _M_r;
    _Tp* const p;

    sp_ref (sp_base& __r, _Tp* __p);

    _Tp* release() const;

};

template <class _Tp> class sp_const_ref {
public:
    sp_base&         _M_r;
    const _Tp* const p;

    sp_const_ref (sp_base& __r, const _Tp* __p);

    const _Tp* release() const;

};

/** EDR_DELETE must be used in place of delete when freeing memory that has
    interacted with smart pointers */
#define EDR_DELETE(__ptr) {if ((__ptr) && (__ptr)->getRefCount() == 0){ \
    delete (__ptr);}}

/** 'Smart' pointer for const pointers - uses reference counting with the
    count being contained within the object. Can only be used for classes
    that have a getRefCount() method such as IObject. Use refCountPtr
    in these situation where there is no getRefCount() method.
    Note that the constructors which take pointers always take ownership of
    the memory.

    Technical note: It is possible to use a single smartPtr template
    for const and non const pointers. The problem is that (1) this
    requires the clone operator to be reworked, (2) causes problems
    when paremeters to a function call need to be casted from
    smartPtr<A> to smartPtr<const A> and the function call takes a non
    const ref */
template<class _Tp> class smartConstPtr{
protected:
    _Tp*  p;
    void  __set(const _Tp* p);
    void  __set(_Tp* p);
public:
    /* General note:
       getRefCount() method is used to hold count of objects. The sign of
       the number is used to determine if delete should be invoked when
       the count hits zero. (This is to enable release to work in a sensible
       manner then other smart pointers have a hold on the object.) If
       getRefCount() == 0 - this indicates a special value meaning that the
       smart pointer(s) does not own the memory and never have/never will */

    typedef _Tp            element_type;
    typedef smartConstPtr<_Tp>  _Self;

    //// clones the contained object (ie invokes clone and casts pointer
    //// back to this type)
    _Tp* clone() const;

    /** Releases the pointer from this smartPointer and any other smart
        pointers referencing this pointer. The pointer is returned */
    const _Tp* release();

    /** clears current pointer and takes ownership of supplied pointer */
    void reset(const _Tp* __px=0);

    const _Tp* get() const;

    const _Tp& obj() const;

    const _Tp* operator->() const;

    const _Tp& operator*() const;

    smartConstPtr();

    /** Constructor from 'raw' pointer. smart ptr ups ref count - for
        a object just created with new this means it takes ownership
        of memory. Memory management still works even if pointer has come
        from smartConstPtr.get() */
    explicit smartConstPtr(const _Tp* __px);

    // looks like stuff put in for older compilers messes up newer compilers
    // Inclusion/exclusion of these has been done on a trial and error basis
#if (!defined(__GNUC__) || __GNUC__ < 4) && (!defined(_MSC_VER) || _MSC_VER < 1400)
    /** Constructor from smart pointer. Ref count is
        incremented. Obscure technical point: ideally the parameter
        __r is not const (see stl auto_ptr) as this helps preserve
        const correctness. However if you do that then you can't use
        the class inside containers such as vector! */
    template<class _Tp1> smartConstPtr(const smartConstPtr<_Tp1>& __r) {
        const _Tp* __conversionCheck = __r.get();
        incRefCount(__conversionCheck);
        this->__set(__conversionCheck);
    }
#endif
    /** Copy Constructor from smart pointer. Ref count is incremented.
        See above comment regarding constness of __r */
    template<class _Tp1> smartConstPtr<_Tp>& operator=(
        const smartConstPtr<_Tp1>& __r){
        const _Tp* __conversionCheck = __r.get();
        update(__conversionCheck);
        return *this;
    }

    smartConstPtr(const smartConstPtr<_Tp>& __r);

    _Self& operator=(const _Self& __r);

    ~smartConstPtr();

    // looks like stuff put in for older compilers messes up newer compilers
    // Inclusion/exclusion of these has been done on a trial and error basis
#if (!defined(__GNUC__) || __GNUC__ < 4) && (!defined(_MSC_VER) || _MSC_VER < 1400)
    smartConstPtr(sp_const_ref<_Tp> __r);

    _Self& operator=(sp_const_ref<_Tp> __r);
#endif
#if 0
    /** If this is included then it breaks gcc 2.95 - unclear why. Nor is
        what the effect if missing this is */
    template<class _Tp1> operator sp_const_ref<_Tp1>() const {
        return sp_const_ref<_Tp1>(*this, this->get());
    }
#endif
#if !defined(__ICC)
    template<class _Tp1> operator smartConstPtr<_Tp1>() const{
        return smartConstPtr<_Tp1>::attachToRef(get());
    }
#endif
    // extra bits and bobs ...
    //// dynamic cast
    template<class _Tp1> static smartConstPtr<_Tp> dynamicCast(
        const smartConstPtr<_Tp1>& __r){
        const _Tp1* __object = __r.get();
        CClassConstSP __clazz = _Tp::TYPE;
        const _Tp* value = ((const _Tp*)classDynamicCast(__clazz,__object));
        return smartConstPtr<_Tp>::attachToRef(value);
    }

    //// operator!
    bool operator!() const;

    //// operator '==' (just compares pointers ie the same as if the SP really
    //// was just a pointer
    bool operator==(const _Self& __r) const;

    /** attach to ref - intelligent version. Looks at count, if non-zero
        ups count. */
    static smartConstPtr attachToRef(const _Tp* __r);
protected:
    /* updates smart pointer with new pointer - does not assume ownership
       of memory */
    void update(const _Tp* __px=0);

    static void incRefCount(const _Tp* __px);

    static void decRefCount(const _Tp* __px);

    /** become owner of pointer */
    static void becomeOwner(const _Tp* __px);
};

/** 'Smart' pointer for non-const pointers - uses reference counting with the
    count being contained within the object. Can only be used for classes
    that have a getRefCount() method such as IObject. Use refCountPtr
    in these situation where there is no getRefCount() method.
    Note that the constructors which take pointers always take ownership of
    the memory.
*/
template<class _Tp> class smartPtr: public smartConstPtr<_Tp>  {
public:
    typedef _Tp                 element_type;
    typedef smartPtr<_Tp>      _Self;
    typedef smartConstPtr<_Tp> _Parent;

    ~smartPtr();

    _Tp* release();

    /** clears current pointer and takes ownership of supplied pointer */
    void reset(_Tp* __px=0);

    _Tp* get() const;

    _Tp& obj() const;

    _Tp* operator->() const;

    _Tp& operator*() const;

    smartPtr();

    /** Constructor from 'raw' pointer. smart ptr ups ref count - for
        a object just created with new this means it takes ownership
        of memory. Memory management still works even if pointer has come
        from smartPtr.get() */
    explicit smartPtr(_Tp* __px);

    // looks like stuff put in for older compilers messes up newer compilers
    // Inclusion/exclusion of these has been done on a trial and error basis
#if (!defined(__GNUC__) || __GNUC__ < 4) && (!defined(_MSC_VER) || _MSC_VER < 1400)
    /** Constructor from smart pointer. Ref count is
        incremented. Obscure technical point: ideally the parameter
        __r is not const (see stl auto_ptr) as this helps preserve
        const correctness. However if you do that then you can't use
        the class inside containers such as vector! */
    template<class _Tp1> smartPtr(const smartPtr<_Tp1>& __r) {
        _Tp* __conversionCheck = __r.get();
        incRefCount(__conversionCheck);
        this->__set(__conversionCheck);
    }
#endif
    /** Copy Constructor from smart pointer. Ref count is incremented.
        See above comment regarding constness of __r */
    template<class _Tp1> smartPtr<_Tp>& operator=(const smartPtr<_Tp1>& __r){
        _Tp* __conversionCheck = __r.get();
        update(__conversionCheck);
        return *this;
    }

    smartPtr(const smartPtr<_Tp>& __r);

    _Self& operator=(const _Self& __r);

    // looks like stuff put in for older compilers messes up newer compilers
    // Inclusion/exclusion of these has been done on a trial and error basis
#if (!defined(__GNUC__) || __GNUC__ < 4) && (!defined(_MSC_VER) || _MSC_VER < 1400)
    smartPtr(sp_ref<_Tp> __r);

    _Self& operator=(sp_ref<_Tp> __r);
#endif

    template<class _Tp1> operator sp_ref<_Tp1>() const {
        return sp_ref<_Tp1>(*this, this->get());
    }
#if !defined(__ICC)
    template<class _Tp1> operator smartPtr<_Tp1>() const {
        return smartPtr<_Tp1>::attachToRef(get());
    }
#endif
    // extra bits and bobs ...
    // dynamic cast
    template<class _Tp1> static smartPtr<_Tp> dynamicCast(
        const smartPtr<_Tp1>& __r){
        _Tp1* __object = __r.get();
        CClassConstSP __clazz = _Tp::TYPE;
        _Tp* value = ((_Tp*)classDynamicCast(__clazz, __object));
        return smartPtr<_Tp>::attachToRef(value);
    }

    /** attach to ref - intelligent version. Looks at count, if non-zero
        ups count. */
    static smartPtr attachToRef(_Tp* __r);

    static smartPtr constCast(const smartConstPtr<_Tp>& __r);
};

// then include templated definitions
#include "edginc/smartPtr.inl"

#include "edginc/oldRefCountPtr.hpp"
//    const X& obj() const throw() { return *ptrConst; }
//    X& obj() const throw() { return *this->ptr; }

DRLIB_END_NAMESPACE

#include "edginc/refCountPtr.hpp"
#endif
