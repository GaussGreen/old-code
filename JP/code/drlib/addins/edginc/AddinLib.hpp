//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : AddinLib.cpp
//
//   Description : Ensures all classes get linked in (important if there are
//                 no other dependencies upon a class)
//
//   Author      : Mark A Robson
//
//   Date        : 28 Feb 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDG_ADDINLIB_HPP
#define EDG_ADDINLIB_HPP

DRLIB_BEGIN_NAMESPACE

/** class insures all classes in the addin directory are linked in. */
class ADDINS_DLL AddinLib{
public:
    // calling this method causes all classes in the addin library to be
    // linked in
    static bool linkInClasses();
    // calling this method causes all classes in the whole library to be
    // linked in
    static void linkInAllLibs();
private:
    // cannot be instantiated
    AddinLib();
};
        
DRLIB_END_NAMESPACE

#endif
