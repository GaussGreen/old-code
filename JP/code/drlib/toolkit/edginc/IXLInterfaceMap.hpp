//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : IXLInterfaceMap.hpp
//
//   Description : Allow objects to have a different Interface 
//                 on the spreadsheet to the internal class model
//
//   Author      : Stephen Hope
//
//   Date        : 5th March 2002
//
//
//----------------------------------------------------------------------------
#ifndef IXL_INTERFACE_MAP_HPP
#define IXL_INTERFACE_MAP_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

/** Interface that objects must inherit in order to map from an objects
    internal class data model to the spreadsheet representation OR back */

class TOOLKIT_DLL IXLInterfaceMap {
public:
    static CClassConstSP const TYPE;

    virtual ~IXLInterfaceMap();

    /** Map the object */
    virtual IObjectSP map()const = 0;

protected:
    IXLInterfaceMap();

};

DRLIB_END_NAMESPACE
#endif
