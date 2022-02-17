//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : CreditLib.hpp
//
//   Description : Force linker to include all relevant files/symbols
//
//   Author      : Jay Blumenstein
//
//   Date        : 22 Aug 2002
//
//


#ifndef EDG_CREDIT_LIB_HPP
#define EDG_CREDIT_LIB_HPP

DRLIB_BEGIN_NAMESPACE

class CREDIT_DLL CCreditLib{
public:
    /** The sole purpose of this method is to ensure that the linker 
        includes all symbols out of the credit directory. Many symbols
        are automatically linked because they are used by other classes
        which are already included.
        
        An example of symbols that could be dropped would be an entire class
        representing a product which was referenced by no other classes */
    static void linkInClasses();
private:
    CCreditLib();
};

DRLIB_END_NAMESPACE

#endif
