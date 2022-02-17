//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLConvertFactory.hpp
//
//   Description : Provide the correct instance of an XLConvert class
//
//   Author      : Mark A Robson
//
//   Date        : 12 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_XLCONVERTFACTORY_HPP
#define EDG_XLCONVERTFACTORY_HPP

#include "edginc/XLConvert.hpp"
DRLIB_BEGIN_NAMESPACE


class ADDINS_DLL XLConvertFactory{
public:
    /** Create an instance of an XLConvert which corresponds to the
        given class */
    static const XLConvert* create(CClassConstSP clazz);

    /** Supply instance of XLConvert corresponding to given class
        which should not be an array.  NB XLConvertFactory takes
        ownership of memory */
    static bool registerXLObjectConvert(CClassConstSP clazz,
                                        XLConvert*    instance);

    /** Supply instance of XLConvert corresponding to arrays whose
        elements are of, or are derived from, given class.  NB
        XLConvertFactory takes ownership of memory */
    static bool registerXLArrayConvert(CClassConstSP  arrayComponentType,
                                       XLConvert*    instance);

private:
    /** Does look up for exact type. Returns null if no match */
    static const XLConvert* lookUp(
        CClassConstSP clazz,
        bool          clazzIsArrayComponentType);

};

DRLIB_END_NAMESPACE
#endif
