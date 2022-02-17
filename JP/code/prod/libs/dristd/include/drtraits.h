// -----------------------------------------------------------------------
// This proprietary software has been developed strictly for J.P. Morgan's
// internal use.  Any use or misuse, intentional or otherwise, which
// contradicts or places this policy in jeopardy is strictly forbidden.
//
// Copyright 2000 J.P. Morgan & Co. Incorporated.  All rights reserved.
// -----------------------------------------------------------------------
//
// 1.0 04/23/2003 Victor Paskhaver
//
// $Header$
//

#pragma once

#ifndef _dr_traits_h
#define _dr_traits_h

#include <driMagnet/MacroInterface.h>
#include <driMagnet/NameConverters.h>

#include <driMagnet/DefaultTraits.h>

namespace drdri {
    
#include <stdexcept>

#include <driMagnet/UserTypes/SharedPointer.h>
#include <driMagnet/UserTypes/Matrix.h>


/** Common conversions - I would consider making them part of the traits file
*/
struct DateAsLong {
    typedef lego::Date ExportAsType;
    static void toPublic( long src, ExportAsType &value ) {
        value.year  =  src / 10000L;
        value.month = (src % 10000L) / 100L;
        value.day   = (src % 100L  );
    }
    static void fromPublic( long &dst, ExportAsType value ) {
        dst = value.year * 10000L + value.month * 100L + value.day;
    }
};





/*
 When building an object, Lego performs following steps:
 - Objects is created on allocated memory. This step can be customized using 
    ClassTraits::createInstance, ClassTraits::deleteInstance 
 - Data from a map is copied into object data members
 - Object is finalized by calling ClassTraits::finalize. In this example finalize is called 
    only on types that are derived from ISupportsFinalize

(SK) We should consider mandating that each exported class provides 'finalize' method.

*/
    
class ISupportsFinalize{
public:
    virtual void finalize() = 0;
};

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

// if the type is derived from ISupportsFinalize, finalize it
 void finalizeHelper(ISupportsFinalize *that) {
    that->finalize();
 }
 void  finalizeHelper(...) {}

// ToLego<T>::Type maps library basic types into basic types understood by lego
// Use template specialization.  What is the DRI type what corresponds to the library type.
//
// If it maps to NoType, then it's not a basic type. 

template <class T> struct ToLego { typedef lego::NoType Type; };

// Basic types are specified below.  Types on the right are DRI types which can only be
// bool, int, double, String, Date, Variant.  Template arguments correspond to the
// types used in the DR libraries.

template<> struct ToLego<bool> { typedef bool Type; };

template<> struct ToLego<int> { typedef int Type; };
template<> struct ToLego<long> { typedef int Type; };

template<> struct ToLego<double> { typedef double Type; };
template<> struct ToLego<float> { typedef double Type; };

template<> struct ToLego<std::string> { typedef lego::String Type; };
template<> struct ToLego<char> { typedef lego::String Type; };

// FromLego<T>::Type is used to coerce from basic types to user defined types
// When a function expects an object but a basic type 
// is provided.  Here we specify how to map DRI basic type to DR library
// basic type.  Constructor is automatically called from the DR library basic
// type if it is provided.  Example coersions could be from strings ("Act/360"),
// int (360), etc. Each type must be explicitely specified below or automatic
// coersion is not attempted.  If you do not specify anything, there will
// be no automatic coersion.
// The conversion works like this: 
//      lego basic type ==(convert using FromLego<T>)==> 
//      library basic type ==(constructor that takes library basic type)==> user defined type


template<class T> struct FromLego { typedef lego::NoType Type; };
template<> struct FromLego<lego::Date> { typedef lego::Date Type; };
template<> struct FromLego<lego::String> { typedef std::string Type; };



// The following templates must be defined in the DR library namespace.  Lego
// doesn't know what the libraries namespace will be so they are defined here.

template <class T> struct ReflectionOf { enum { value = lego::TypeNotDeclared }; };
template <class T> struct PrefixOf { enum { value = lego::TypeNotDeclared }; };
template <class T> struct ClassTraitsBase {
    typedef ReflectionOf<T> ReflectionOfType;
    typedef PrefixOf<T> PrefixOfType;
};

// Every library defines these traits in it's namespace.
struct LibraryTraits : lego::DefaultLibraryTraits 
{


    struct BasicTypeTraits {
// The following imports the basic types defined at the top of this file
// to this object.  The following can be thought of as compile time functions.
// For example, to tell the compiler to convert from library Date to Lego::Date,
// one would have the following in the code:
// LibraryTraits::BasicTypeTraits::MapToLego::_Result<IR::Date>::Type
// The reslt of this expression is Lego::Date.
// Users of Lego do not need to do any of this compile time function invocation.
// Rather, all type conversions are done automatically by Lego at compile time.
        struct MapToLego   { template<class T> struct _Result : ToLego<T> { }; };
        struct MapFromLego { template<class T> struct _Result : FromLego<T> { }; };

// If specify mapping, first Lego tries to use implicite conversion (eg char to int would
// work but IR::Date to Lego::Date would not automatically work). The convert method
// below allows this information to be explicitely provided (but is only invoked
// after the implicite conversions have failed so one couldn't overwrite the char to int
// conversion).

        static void convert( const std::string &src, lego::String &result) {
            result = lego::String(src.c_str());
        }
        
        static void convert( const lego::String &src, std::string &result) {
            result = std::string(src.getString());
        }
        
        static void convert( char src, lego::String &result) {
            char temp[2] = { src, '\0' };
            result = lego::String(temp);
        }
        
        static void convert( const lego::String &src, char& result) {
            const char *str = src.getString();
            if(!str || str[0] =='\0' || str[1] != '\0') {
                throw std::runtime_error( 
                    std::string("Cannot convert ")+ 
                    (str ? str : "null string") + " to a character");
            }
            result = str[0];
        }
    };


// enums are compile-time integer constants.  Think of this as a compile time
// function which takes a type and returns a boolean.  It allows Lego to 
// check whether an object is an array, matrix, or pointer at compile time.
// If your library doesn't support one of these, just set the value to 
// false. Note that isPointer really means isSmartPointer.
    struct IsType { 
        template<class T> struct _Result : lego::DefaultIsType<T>
        {
            enum { 
                isArray = sizeof(isArrayHelper((T*)0))==sizeof(char),
                isMatrix = sizeof(isMatrixHelper((T*)0))==sizeof(char),
                isPointer = sizeof(isPointerHelper((T*)0))==sizeof(char)
            };
        }; 
    };

//////////////////////////////////////////////////////////////////////////

    struct ArrayTraits { template<class T> struct _Result {
// Input here is T which is an array type. DR library implementation
// of array must provide the following.
        typedef  typename T::value_type ElementType;
        static int getSize(const T &array) {
            return array.size();
        }
        static void  setSize(T &array, int sz) {
            array.resize(sz);
        }
        static ElementType &getArrayElement(T &array, int index) {
            return array[index];
        }
    }; };

//////////////////////////////////////////////////////////////////////////

    struct MatrixTraits { template<class T> struct _Result {
// Input here is T which is an matrix type. DR library implementation
// of matrix must provide the following.
        typedef  typename T::value_type ElementType;
        static int getRowSize(const T &matrix) {
            return matrix.rows();
        }
        static int getColSize(const T &matrix) {
            return matrix.cols();
        }
        static void setSize(T &matrix, int rows, int cols ) {
            matrix.initSize(rows, cols);
        }
        static ElementType &getMatrixElement(T &matrix, int row,int col) {
            return matrix.at(row, col);
        }

    }; };

//////////////////////////////////////////////////////////////////////////

    struct SmartPointerTraits { template<class T> struct _Result {
// RawPointerType is pointer to type that smart pointer points to
        typedef typename T::value_type *RawPointerType;

// Extract pointer to object from smart pointer        
        static RawPointerType get(const T &ptr ) {
            return ptr.get();
        }
// Assign pointer to just allocated object. "ptr" is going to be the first reference to "val"
        static void setRawPointer(T &ptr, RawPointerType val ) {
            ptr = T(val);
        }

// Assign pointer to existing object. "dstHolder" is going to be a reference to "newValue".
// "scrHolder" is one of the existing references. If reference count is external, 
// this function should copy pointer to ref count from srcHolder to dstHolder
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
    }; };

//////////////////////////////////////////////////////////////////////////
// This type is used to describe memory allocation policy for LG_MEMBER_ARRAY
// where the type of the array is C pointer
    struct DynamicArrayTraits { template<class ElementType> struct _Result {
        static ElementType *allocate(int sz) {
            return (ElementType*)malloc(sizeof(ElementType)*sz);
        }
    }; };
    
//////////////////////////////////////////////////////////////////////////
    
    struct ClassTraits { 
        template<class T> struct _Result : ClassTraitsBase<T>, lego::DefaultClassTraits<T> 
        {
// Specify whether one could declare smart pointers which refer T.
// We need to know this because both the client of DRI and the DR
// library can hold pointers to the same object. They need to share the
// same reference count object or else one of them (the client or the DR
// internal object) will become invalid.
            enum { supportsSmartPointer = true };
            typedef shared_ptr<T> SmartPointerType;


// The following method is invoked right after an object is created.  It 
// provides the DR library with a way to check that the object
// is valid or return an error as early as possible.  

// require finalize method to be present in all exported types
            static void finalize( T &that ) {
                that.finalize();
            }

            template <class PublicClass>
            static void toPublicInstance(
                T *instancePtr, 
                PublicClass *publicInstancePtr) 
            {
                instancePtr->toPublicObject(publicInstancePtr);
            }

            template <class PublicClass>
            static void fromPublicInstance(
                T *instancePtr, 
                PublicClass *publicInstancePtr) 
            {
                instancePtr->fromPublicObject(publicInstancePtr);
            }

        }; 
    };

//////////////////////////////////////////////////////////////////////////
// When specifying the external view of class, enum, and data members, NameConverters
// can be used to strip out internal conventions.  For example, some libraries use
// m_ before all data members but we do not wish to display this to external users.
// Similarly, the "C" before class names can be automatically removed for a library
// whose convention is the preceed any class with "C" (eg. CZeroCurve).
    struct NameConverters : lego::DefaultNameConverters
    {
        static std::string dataMemberName(const char *name) { 
            return lego::MaybeStripPrefix(name, "m");
        }
    };

//////////////////////////////////////////////////////////////////////////
// The following provides the type of the exception and how to get the description 
// from that type.
// (SK) we should use standard exceptions instead of inventing our own

    struct ExceptionTraits {
        typedef std::exception Type;
        static const char *what(const Type &e) {
            return e.what();
        }
    };
};

#include <driMagnet/UserTypes/SharedPointer.inl>

} // drdri

#endif
