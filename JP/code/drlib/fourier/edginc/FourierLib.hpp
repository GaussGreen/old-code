//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FourierLib.hpp
//
//   Description : Force linker to include all relevant files/symbols
//
//   Author      : Regis Guichard
//
//   Date        : 09 March 2004
//
//
//----------------------------------------------------------------------------

#ifndef EDG_FOURIER_LIB_HPP
#define EDG_FOURIER_LIB_HPP

DRLIB_BEGIN_NAMESPACE

class FOURIER_DLL FourierLib{
public:
    /** The sole purpose of this method is to ensure that the linker 
        includes all symbols out of the riskmr directory. Many symbols
        are automatically linked because they are used by other classes
        which are already included.
        
        An example of symbols that could be dropped would be an entire class
        representing a product which was referenced by no other classes */
    static void linkInClasses();

private:
    FourierLib();
};

DRLIB_END_NAMESPACE

#endif
