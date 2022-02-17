//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Modifier.cpp
//
//   Description : The Modifier class provides static methods and 
//                 constants to decode class and member access modifiers. 
//
//   Author      : Mark A Robson
//
//   Date        : 4th June 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Modifier.hpp"

DRLIB_BEGIN_NAMESPACE

/** The Modifier class provides static methods and constants to decode
    class and member access modifiers. Limited functionality at present */

//// The int value representing the abstract modifier.  
int const Modifier::ABSTRACT = 0x0001;

#if 0
//// The int value representing the final modifier.    
int const Modifier::FINAL = 0x0002;
#endif
               
//// The int value representing the interface modifier.
int const Modifier::INTERFACE = 0x0004;

//// The int value representing the private modifier.
int const Modifier::PRIVATE = 0x0008;

//// The int value representing the protected modifier.
int const Modifier::PROTECTED = 0x0010;

//// The int value representing the public modifier.
int const Modifier::PUBLIC = 0x0020;

//// The int value representing the modifier.
int const Modifier::STATIC = 0x0040;

//// The int value representing the transient modifier.
int const Modifier::TRANSIENT = 0x0080;

//// The int value representing the export modifier.
int const Modifier::EXPORT = 0x0200;

#if 0
//// The int value representing the volatile modifier.
int const Modifier::VOLATILE = 0x0100;
#endif

////Return true if the specifier integer includes the abstract modifier.
bool Modifier::isAbstract(int mod){
    return (mod & ABSTRACT)? true: false;
}

#if 0
//// Return true if the specified integer includes the final modifier.
bool Modifier::isFinal(int mod){
    return (mod & FINAL)? true: false;
}
#endif

//// Return true if the specifier integer includes the interface modifier.
bool Modifier::isInterface(int mod){
    return (mod & INTERFACE)? true: false;
}

//// Return true if the specifier integer includes the private modifier.
bool Modifier::isPrivate(int mod){
    return (mod & PRIVATE)? true: false;
}

//// Return true if the specifier integer includes the protected modifier.
bool Modifier::isProtected(int mod){
    return (mod & PROTECTED)? true: false;
}

//// Return true if the specified integer includes the public modifier.
bool Modifier::isPublic(int mod){
    return (mod & PUBLIC)? true: false;
}

//// Return true if the specifier integer includes the modifier.
bool Modifier::isStatic(int mod){
    return (mod & STATIC)? true: false;
}

//// Return true if the specifier integer includes the transient modifier.
bool Modifier::isTransient(int mod){
    return (mod & TRANSIENT)? true: false;
}

//// Return true if the specifier integer includes the export modifier.
bool Modifier::isExport(int mod){
    return (mod & EXPORT)? true: false;
}

#if 0
//// Return true if the specified integer includes the volatile modifier.
bool Modifier::isVolatile(int mod){
    return (mod & VOLATILE)? true: false;
}
#endif

/** Return a string describing the access modifier flags in the
    specified modifier. */
string Modifier::toString(int mod){
    // to add others when they get supported
    return (string(isPrivate(mod)? "private ": "")+
            string(isProtected(mod)? "protected ": "")+        
            string(isPublic(mod)? "public ": "")+        
            string(isAbstract(mod)? "abstract ": "")+
            string(isInterface(mod)? "interface ": "")+
            string(isStatic(mod)? "static ": "")+
            string(isTransient(mod)? "transient ": "") + 
            string(isExport(mod)? "export ": ""));
}


DRLIB_END_NAMESPACE
