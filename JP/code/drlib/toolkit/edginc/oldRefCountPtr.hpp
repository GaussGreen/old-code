// Filename:  oldRefCountPtr.hpp
// Purpose:   legacy implementation of reference counded smart pointers that provides release() method.

// This header defines [drlib::]attic::refCountPtr<T> template -- the old implementation of reference counted smart pointers.
// It is expected that people don't use this file unless for legacy purposes (use  <refCountPtr.hpp> instead).
// The main differences between refCountPtr and the QlibAttic::refCountPtr are the following:
// The new one (based on Boost) has no  refCountConstPtr -- use refCountPtr<const T>
// The syntax for dynamic_cast and const_cast is slightly different, use DYNAMIC_POINT_CAST<NewType>(ptrSP) macro
// The new one doesn't have "attach()" method -- instead you can create shared pointer with null-deleter (i.e. ptrSP = typeSP(new type(), NullDeleter())). NullDeleter is defined in <refCountPtr.hpp>
// The new one doesn't have "release()" method and there is no easy way to emulate it. See Boost shared_pt FAQ on why it was decided to removed this feature (basically, if one uses it, memory management is suddenly an issue again). Practically, we keep the old implementation to support two places in Qlib that use release() functionality (only XObject and FRController have refCountedPtr and use release()).



#ifndef OLD_REFCOUNT_PTR
#define OLD_REFCOUNT_PTR

namespace attic { // non-standard ref. counted smart pointers

USING_CORE_NAMESPACE;

// DRLIB_BEGIN_NAMESPACE
/* structure holds reference count details - currently sizeof
   this structure is 5. We could reduce it to 4 by using signed int.
   However it is only created when an object is created so the
   impact of the extra byte will be small */
/** Structure holds reference count details for smart pointers */
typedef struct _TObjectRef
{
    unsigned int count;   // number of references outstanding
    bool delObj;  // true: delete object
}
TObjectRef;

/** old style version of smartConstPtr - works with
    objects not derived from IObject */
template <class X>
class refCountConstPtr
{
protected:
    union{
        const X* ptrConst;
        X* ptr;
    };
    mutable TObjectRef* ref; /* null implies no ownership else stores
            number of references to pointer etc */

    refCountConstPtr( const X* p, TObjectRef* ref ) throw() :
            ptrConst( p ), ref( ref )
    {
        if ( ref ) {
            ref->count++;
        }
    }

    /** release the memory from this smart pointer - this may not always
            be possible if the memory was never owned in the first place
            or it has already been released */
    X* releasePtr() const
    {
        if ( !ref ) {
            throw ModelException( "refCountPtr::release",
                                  "Smart pointer never owned memory" );
        }
        if ( !ref->delObj ) {
            throw ModelException( "refCountPtr::release",
                                  "Memory has already been released" );
        }
        ref->delObj = false;
        return ptr;
    }

public:
    // takes ownership of memory
    explicit refCountConstPtr( const X* p ) : ptrConst( p )
    {
        if ( p ) {
            ref = new TObjectRef;
            ref->count = 1;
            ref->delObj = true;
        } else {
            ref = 0; // set to null
        }
    }

    explicit refCountConstPtr() throw() :
            ptrConst( 0 ), ref( 0 )
    {}

    /** Constructor from smart pointer. Ref count is incremented.*/
    template <class _Tp1>
    refCountConstPtr( const refCountConstPtr<_Tp1>& __r ) :
            ptrConst( __r.get() ), ref( __r.refCount__() )
    {
        if ( ref ) {
            ref->count++;
        }
    }
    /** Copy Constructor from smart pointer. Ref count is incremented. */
    template <class _Tp1>
    refCountConstPtr<X>& operator=(
        const refCountConstPtr<_Tp1>& rhs )
    {
        if ( this != &rhs ) {
            this->~refCountConstPtr();
            ref = rhs.refCount__(); // common reference count
            if ( ref ) {
                ref->count++; // up reference count
            }
            ptrConst = rhs.ptrConst;
        }
    }

    refCountConstPtr( const refCountConstPtr& a ) throw() :
            ptrConst( a.ptrConst ), ref( a.ref )
    {
        if ( ref ) {
            ref->count++;
        }
    }

    refCountConstPtr& operator=( const refCountConstPtr& rhs ) throw()
    {
        if ( this != &rhs ) {
            this->~refCountConstPtr();
            ref = rhs.ref; // common reference count
            if ( ref ) {
                ref->count++; // up reference count
            }
            ptrConst = rhs.ptrConst;
        }
        return ( *this );
    }

    /** dynamic cast from type _Tp1 to type X */
    template <class _Tp1>
    static refCountConstPtr<X> dynamicCast(
        const refCountConstPtr<_Tp1>& __r )
    {
        const _Tp1 * __object = __r.get();
        const X* ptr = &dynamic_cast<const X&>( *__object );
        return refCountConstPtr<X>::attach__( ptr, __r.refCount__() );
    }

    ~refCountConstPtr()
    {
        if ( ref ) {
            ref->count--;
            if ( ref->count == 0 ) {
                if ( ref->delObj ) {
                    delete ptr;
                }
                delete ref;
            }
        }
    }

    // for creating an refCountConstPtr from a reference - will not free memory
    static refCountConstPtr attachToRef( const X* p ) throw()
    {
        return refCountConstPtr( p, 0 );
    }

    /** release the memory from this smart pointer - this may not always
                    be possible if the memory was never owned in the first place
                    or it has already been released */
    const X* release() const
    {
        return releasePtr(); // cast X* to const X*
    }

    const X& operator*() const throw()
    {
        return * ptrConst;
    }
    const X* operator->() const throw()
    {
        return ptrConst;
    }
    const X* get
        () const throw()
    {
        return ptrConst;
    }
    bool operator!() const throw()
    {
        return ptrConst == 0;
    }

    /** returns true if this smart pointer owns the object ie will delete it
                    when ref count hits zero */
    bool owns() const
    {
        return ref ? ref->delObj : false;
    }

    /** clears current pointer and takes ownership of supplied pointer */
    void reset( const X* p = 0 )
    {
        *this = refCountConstPtr<X>( p );
    }

    /** ideally would make accessible to all templates via friends - but
                    not supported by compiler - so therefore give obscure name */
    static refCountConstPtr attach__( const X* p, TObjectRef* ref )
    throw()
    {
        return refCountConstPtr( p, ref );
    }
    /** ditto for this function */
    TObjectRef* refCount__() const throw()
    {
        return ref;
    }

};

/** old style version of smartPtr - works with objects not derived
    from IObject */
template <class X>
class refCountPtr: public refCountConstPtr<X>
{
    typedef refCountConstPtr<X> Parent_;
public:
    explicit refCountPtr( X* p = 0 ) throw() :
            Parent_( p )
    {}

    /** Constructor from smart pointer. Ref count is incremented.*/
    template <class _Tp1>
    refCountPtr( const refCountPtr<_Tp1>& __r ) :
            Parent_( __r.get(), __r.refCount__() )
    {}

    /** Copy Constructor from smart pointer. Ref count is incremented. */
    template <class _Tp1>
    refCountPtr<X>& operator=(
        const refCountPtr<_Tp1>& rhs )
    {
        if ( get () != rhs.get() ) {
                this->~refCountPtr();
                this->ref = rhs.refCount__(); // common reference count
                if ( this->ref ) {
                    this->ref->count++; // up reference count
                }
                this->ptrConst = rhs.get();
            }
        return ( *this );
    }

    refCountPtr( const refCountPtr& a ) throw() :
            Parent_( a.ptr, a.ref )
    {}

    refCountPtr& operator=( const refCountPtr& rhs ) throw()
    {
        if ( this != &rhs ) {
            this->~refCountPtr();
            this->ref = rhs.ref;
            this->ptr = rhs.ptr;
            if ( this->ref ) {
                this->ref->count++;
            }
        }
        return ( *this );
    }

    // for creating an refCountPtr from a reference - will not free memory
    static refCountPtr attachToRef( X* p ) throw()
    {
        return refCountPtr( p, 0 );
    }


    /** allow const cast from constSP to this class */
    static refCountPtr constCast( const Parent_& a )
    {
        return refCountPtr( a );
    }

    /** dynamic cast from type _Tp1 to type _Tp*/
    template <class _Tp1>
    static refCountPtr<X> dynamicCast(
        const refCountPtr<_Tp1>& __r )
    {
        _Tp1 * __object = __r.get();
        X* ptr = &dynamic_cast<X&>( *__object );
        return refCountPtr<X>::attach__( ptr, __r.refCount__() );
    }

    /** release the memory from this smart pointer - this may not always
            be possible if the memory was never owned in the first place
            or it has already been released */
    X* release() const
    {
        return this->releasePtr();
    }

    bool operator==( const refCountPtr& rhs ) const
    {
        return this->ptr == rhs.ptr;
    }
    X& operator*() const throw()
    {
        return * this->ptr;
    }
    X* operator->() const throw()
    {
        return this->ptr;
    }
    X* get
        () const throw()
    {
        return this->ptr;
    }
    bool operator!() const throw()
    {
        return this->ptr == 0;
    }

    /** clears current pointer and takes ownership of supplied pointer */
    void reset( X* p = 0 )
    {
        *this = refCountPtr<X>( p );
    }

    /** ideally would make accessible to all templates via friends - but
            not supported by compiler - so therefore give obscure name */
    static refCountPtr attach__( X* p, TObjectRef* ref ) throw()
    {
        return refCountPtr( p, ref );
    }

private:
    refCountPtr( X* p, TObjectRef* ref ) throw() :
            Parent_( p, ref )
    {}
    /** allow const cast from constSP to SP */
    refCountPtr( const Parent_& a ) : Parent_( a )
    {}
}
;

// for using refCountPtr<T> within STL
template <class T>
bool
lessThan( const refCountPtr<T>& first, const refCountPtr<T>& second )
{
    if ( first.get() == NULL ) {
        return true;
} else if ( second.get() == NULL ) {
        return true;
} else {
        return *first < *second;
}
}

template <class T>
struct oldCountPtr {
        typedef refCountPtr<T> type;
};

} // namespace qlib::attic

#endif // OLD_REFCOUNT_PTR
