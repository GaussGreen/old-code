//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : TypeConvert.hpp
//
//   Description : Defines interface for the ability for one type to convert
//                 itself into another type. This can be useful for 'factory'
//                 type methods or for backward compatibility
//
//   Author      : Mark A Robson
//
//   Date        : 27 Mar 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_TYPE_CONVERT_HPP
#define EDG_TYPE_CONVERT_HPP

DRLIB_BEGIN_NAMESPACE

/** Defines interface for the ability for one type to convert itself
    into another type. This can be useful for 'factory' type
    approaches to building objects or for backwards compatibility eg
    converting a legacy class into its newer counterpart */
class TOOLKIT_DLL ITypeConvert {
public:
    static CClassConstSP const TYPE; // in Object.cpp

    ITypeConvert(); //  in Object.cpp

    virtual ~ITypeConvert(); //  in Object.cpp
    ITypeConvert(const ITypeConvert& rhs); // in Object.cpp
    ITypeConvert& operator=(const ITypeConvert& rhs); // in Object.cpp

    /** Converts this object to an instance of the requiredType. Throws an
        exception if a conversion to the required Type is not supported.
        Note that the supplied object is a smart pointer to this. The 
        converted object should be stored in the object parameter. */
    virtual void convert(IObjectSP&    object,
                         CClassConstSP requiredType) const = 0;
};


DRLIB_END_NAMESPACE

#endif
