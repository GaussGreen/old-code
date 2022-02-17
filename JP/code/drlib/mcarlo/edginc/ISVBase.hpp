//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ISVBase.hpp
//
//   Description : An abstract base for all the SV interfaces that plan to 
//                 use path() access method 
//
//----------------------------------------------------------------------------

#ifndef ISVBase_HPP
#define ISVBase_HPP

#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/ElemStateVariableGen.hpp"
#include "edginc/SVPath.hpp"

DRLIB_BEGIN_NAMESPACE

//class SVPath;

/** An abstract base for all the SV interfaces that plan to use path() access method */
class ISVBase : public virtual IStateVariable
{
public:
    virtual const SVPath& path() const = 0;
};

DRLIB_END_NAMESPACE

#endif // ISVBase_HPP
