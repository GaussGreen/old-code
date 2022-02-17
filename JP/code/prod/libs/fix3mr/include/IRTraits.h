// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2000 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 10/27/2003 Afshin Bayrooti
//
// $Header: /nasdev/export2/home/drdev/cvsadmin/cvs/libs/fix3/include/IRTraits.h,v 1.1 2004/04/19 15:50:39 markss Exp $
//

// #pragma once

#if ! defined(_IR_TRAITS_)
#define _IR_TRAITS_

#include <Lego/MacroInterface.h>
#include "IRTraitTypes.h"


#ifdef _MSC_VER
#pragma warning(disable:4786)
#pragma warning(disable:4100)
#pragma warning(disable:4512)
#pragma warning(disable:4018)
#pragma warning(disable:4127)

#endif
namespace IR {

namespace detail {
    
    // Use function overloading to detect whether a type is a array, matrix or pointer.
    // If call a function and pass in pointer to type, the sizeof the return type
    // can be used to determine which overload was selected at compilation time.

    template<class T> 
    char isArrayHelper(std::vector<T>*);
    int  isArrayHelper(...);
    template<class T> 
    char isMatrixHelper(Matrix<T>*);
    int  isMatrixHelper(...);
    template<class T> 
    char isPointerHelper(shared_ptr<T>*);
    int  isPointerHelper(...);

    // maps for basic types
    // Use template specialization.  What is the DRI type what corresponds to the
    // IR / Hybrids library type.

    // If it maps to NoType, then it's not a basic type. 
    template <class T> struct ToLego { typedef lego::NoType Type; }; 

    // Basic types are specified below.  Types on the right are DRI types which can only be
    // bool, int, double, String, Date, Variant.  Template arguments correspond to the
    // types used in the IR / Hybrids libraries.
    
    template<> struct ToLego<bool> { typedef bool Type; };
    // What should char be?  Should it look like a string or like a number?
    template<> struct ToLego<char> { typedef lego::String Type; };
    template<> struct ToLego<int> { typedef int Type; };
    template<> struct ToLego<long> { typedef int Type; };
    template<> struct ToLego<double> { typedef double Type; };
    template<> struct ToLego<float> { typedef double Type; };
    template<> struct ToLego<std::string> { typedef lego::String Type; };
    template<> struct ToLego<Date> { typedef lego::Date Type; };

    // For automatic coersion.  When a function expects an object but a basic type 
    // is provided.  Here we specify how to map DRI basic type to IR / Hybrids
    // basic type.  Constructor is automatically called from the IR / Hybrids basic
    // type if it is provided.  Example coersions could be from strings ("Act/360"),
    // int (360), etc. Each type must be explicitely specified below or automatic
    // coersion is not attempted.  If you do not specify anything, there will
    // be no automatic coersion.

    template <class T> struct FromLego { typedef lego::NoType Type; };
    template<> struct FromLego<lego::String> { typedef std::string Type; };

}

// The following templates must be defined in the IR / Hybrids namespace.  Lego
// doesn't know what the libraries namespace will be so they are defined here.
// We may provide a macro which eases the pain on your eyes.

template <class T> struct ReflectionOf { enum { value = lego::TypeNotDeclared }; };
template <class T> struct PrefixOf { enum { value = lego::TypeNotDeclared }; };
template <class T> struct ClassTraitsBase {
    typedef ReflectionOf<T> ReflectionOfType;
    typedef PrefixOf<T> PrefixOfType;
};

// When building an object, Lego can be configured to call a function or member of the class
// to validate the data immediately after creating the object.  This ensures errors are caught
// close to where they occur.  For IR / Hybrids, we've configured Lego to call the following function.
// Overloading will assure the correct version is called for a specified type.

inline void Finalize(...){}

struct LibraryTraits 
{

    // Every library defines these traits in it's namespace.

    struct BasicTypeTraits 
    {

        // The following imports the basic types defined at the top of this file
        // to this object.  The following can be thought of as compile time functions.
        // For example, to tell the compiler to convert from IR / Hybrids Date to Lego::Date,
        // one would have the following in the code:
        // LibraryTraits::BasicStructure::MapToLego::_Result<IR::Date>::Type
        // The reslt of this expression is Lego::Date.
        // Users of Lego do not need to do any of this compile time function invocation.
        // Rather, all type conversions are done automatically by Lego at compile time.

        struct MapToLego   { template<class T> struct _Result : detail::ToLego<T> { }; };
        struct MapFromLego { template<class T> struct _Result : detail::FromLego<T> { }; };

        // If specify mapping, first Lego tries to use implicite conversion (eg char to int would
        // work but IR::Date to Lego::Date would not automatically work). The convert method
        // below allows this information to be explicitely provided (but is only invoked
        // after the implicite conversions have failed so one couldn't overwrite the char to int
        // conversion).


        static void convert( const lego::Date &src, Date &dst ) 
        {
            dst.day = src.day;
            dst.month = src.month;
            dst.year = src.year;
        }
        static void convert( const Date &src, lego::Date &dst) 
        {
            dst.day = src.day;
            dst.month = src.month;
            dst.year = src.year;
        }
        
        static void convert( const std::string &src, lego::String &result) 
        {
            result = lego::String(src.c_str());
        }
        
        static void convert( const lego::String &src, std::string &result) 
        {
            result = std::string(src.getString());
        }

        static void convert( char src, lego::String &result) 
        {
            char temp[2] = { src, '\0' };
            result = lego::String(temp);
        }
        
        static void convert( const lego::String &src, char& result) 
        {
            const char *str = src.getString();
            if(!str || str[0] =='\0' || str[1] != '\0') 
            {
                throw ::lego::Error( 
                    std::string("Cannot convert ")+ 
                    (str ? str : "null string") + " to a character");
            }
            result = str[0];
        }

    };

    
    struct IsType { 
        template<class T> struct _Result : lego::DefaultIsType<T>
        {
            // enums are compile-time integer constants.  Think of this as a compile time
            // function which takes a type and returns a boolean.  It allows Magnet to 
            // check whether an object is an array, matrix, or pointer at compile time.
            // If your library doesn't support one of these, just set the value to 
            // false. Note that isPointer really means isSmartPointer.
            enum 
            {
                isArray = sizeof(detail::isArrayHelper((T*)0))==sizeof(char),
                isMatrix = sizeof(detail::isMatrixHelper((T*)0))==sizeof(char),
                isPointer = sizeof(detail::isPointerHelper((T*)0))==sizeof(char)
            };
        }; 
    };

    //////////////////////////////////////////////////////////////////////////

    struct ArrayTraits 
    {
        template<class T> 
        struct _Result 
        {
            // Input here is T which is an array type.  IR / Hybrids library implementation
            // of array must provide the following.

            typedef  typename T::value_type ElementType;
            static int getSize(const T &array) 
            {
                return array.size();
            }
            static void  setSize(T &array, int sz) 
            {
                array.resize(sz);
            }
            static ElementType &getArrayElement(T &array, int index) 
            {
                return array[index];
            }
        }; 
    };

    //////////////////////////////////////////////////////////////////////////

    struct MatrixTraits 
    { 
        template<class T> 
        struct _Result 
        {
            // Input here is T which is an matrix type.  IR / Hybrids library implementation
            // of matrix must provide the following.
        
            typedef  typename T::value_type ElementType;
            static int getRowSize(const T &matrix) 
            {
                return matrix.rows();
            }
            static int getColSize(const T &matrix) 
            {
                return matrix.cols();
            }
            static void setSize(T &matrix, int rows, int cols ) 
            {
                matrix.resize(rows, cols);
            }
            static ElementType &getMatrixElement(T &matrix, int row,int col) 
            {
                return matrix.at(row, col);
            }

        };
    };

    //////////////////////////////////////////////////////////////////////////

    struct SmartPointerTraits 
    { 
        template<class T> 
        struct _Result 
        {

            // T is a smart pointer to an object (eg. SharedPointer<ZeroCurve> in CMLib).  
            // RawPointerType must be defined in the way that makes sense in the IR / Hybrids
            // library.

            // If your smart pointer is just a pointer (ie. no smart pointer), then
            // the following typedef would be "typedef T RawPointerType;"  However,
            // in this case, the other methods in this structure would have to be
            // modified to do something which makes sense (eg, get returns
            // returns ptr, setRawPointer and setManagedPointer should copy the object.

            typedef typename T::value_type* RawPointerType;
            static RawPointerType get(const T &ptr ) 
            {
                return ptr.get();
            }
            static void setRawPointer(T &ptr, RawPointerType val ) 
            {
                ptr = T(val);
            }
            // The following method allows one smart pointer to be set to another one.
            // For example, say b is a smartpointer to an object of type B and d is a 
            // smart pointer to an derived of type D, where D is derived from B. When 
            // assigning b = d, we will have:
            //      srcHolder = &d
            //      dstHolder = b
            //      newValue = (B*)d.get()
            // The point is you put new value in dstHolder and copy the pointer to reference
            // count in the srcHolder (since the reference count is stored outside of the object
            // itself in the implementation of shared_ptr).  Other implementations of reference
            // counted objects would require a different implemtation of this function.

            static void setManagedPointer(
                void *srcHolder, 
                T &dstHolder, 
                RawPointerType newValue) 
            {
                shared_ptr_base *dstPtr = static_cast<shared_ptr_base*>(&dstHolder);
                shared_ptr_base *srcPtr = static_cast<shared_ptr_base*>(srcHolder);
                if(srcPtr->m_count)
                    srcPtr->m_count->addRef();
                if(dstPtr->m_count)
                    dstPtr->m_count->release();
                dstPtr->m_count = srcPtr->m_count;
                dstPtr->m_ptr = newValue;
            }
        };
    };

    //////////////////////////////////////////////////////////////////////////

    struct ClassTraits { 
        
        template<class T> 
        struct _Result : ClassTraitsBase<T>, lego::DefaultClassTraits<T>
        {
            // Specify whether one could declare smart pointers which refer T.
            // For the IR / Hybrids libraries, since we do reference counting outside
            // of the objects, we can support smart pointers on all objects.  However,
            // in CMLib, only objects derived from a certain class support reference 
            // counting so there, the implementation looks like this:
            //  supportsSmartPointer = lego::IsConvertible<T*, SharedInterface*>::value
            // where the lego::IsConvertible "compiler" function tells you whether
            // one type can be converted into another at compilation time.

            // We need to know this because both the client of DRI and the IR / Hybrids
            // library can hold pointers to the same object. They need to share the
            // same reference count object or else one of them (the client or the IR / Hybrids
            // internal object) will become invalid.

            enum { supportsSmartPointer = true };
            typedef shared_ptr<T> SmartPointerType;

            // The following method is invoked right after an object is created.  It 
            // provides the IR / Hybrids library with a way to check that the object
            // is valid or return an error as early as possible.  
   
            static void finalize( T &that ) {
                Finalize( &that );
            }

            template <class PublicClass>
            static void toPublicInstance(
                T *instancePtr, 
                PublicClass *publicInstancePtr) 
            {
                ToPublicObject(instancePtr, publicInstancePtr);
            }

            template <class PublicClass>
            static void fromPublicInstance(
                T *instancePtr, 
                PublicClass *publicInstancePtr) 
            {
                FromPublicObject(instancePtr, publicInstancePtr);
            }

        }; 
    };

    // When specifying the external view of class, enum, and data members, NameConverters
    // can be used to strip out internal conventions.  For example, some libraries use
    // m_ before all data members but we do not wish to display this to external users.
    // Similarly, the "C" before class names can be automatically removed for a library
    // whose convention is the preceed any class with "C" (eg. CZeroCurve).

    struct NameConverters : lego::DefaultNameConverters
    {
        static std::string className(std::string name) 
        { 
            return lego::MaybeStripPrefix(name.c_str(), "_");
        }
    };

    //////////////////////////////////////////////////////////////////////////

    // The following provides the type of the exception and how to get the description 
    // from that type.

    struct ExceptionTraits 
    {
        typedef Error Type;
        static const char *what(const Type &e) 
        {
            return e.what().c_str();
        }
    };
};


} // IR

#endif
