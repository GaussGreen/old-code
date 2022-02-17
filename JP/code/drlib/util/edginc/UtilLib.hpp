//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : UtilLib.hpp
//
//   Description : ensure all classes are loaded
//
//   Author      : Andrew J Swain
//
//   Date        : 14 June 2001
//
//
//----------------------------------------------------------------------------

#ifndef UTILLIB_HPP
#define UTILLIB_HPP

DRLIB_BEGIN_NAMESPACE

class UTIL_DLL UtilLib{
public:
    /** The sole purpose of this method is to ensure that the linker 
        includes all symbols out of the util directory. Many symbols
        are automatically linked because they are used by other classes
        which are already included.
        
        An example of symbols that could be dropped would be an entire class
        representing a product which was referenced by no other classes */
    static void linkInClasses();
private:
    UtilLib();
};

DRLIB_END_NAMESPACE

#endif
