//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Modifier.hpp
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

#ifndef EDR_MODIFIER_HPP
#define EDR_MODIFIER_HPP
#include <string>
using namespace std;    // string


DRLIB_BEGIN_NAMESPACE

/** The Modifier class provides static methods and constants to decode
    class TOOLKIT_DLL and member access modifiers. Limited functionality at present.
    In particular, fields support transient only.
    Classes support public, protected, private, interface, abstract, static.
    Protected classes are those that 
    implement IPrivateObject (which is a reasonable analogy). Private
    classes are those which call Class.setPrivate() */
class TOOLKIT_DLL Modifier{
public:
    //// The int value representing the abstract modifier.  
    static int const ABSTRACT;

#if 0
    //// The int value representing the final modifier.    
    static int const FINAL;
#endif
               
    //// The int value representing the interface modifier.
    static int const INTERFACE;

    //// The int value representing the private modifier.
    static int const PRIVATE;

    //// The int value representing the protected modifier.
    static int const PROTECTED;

    //// The int value representing the public modifier.
    static int const PUBLIC;

    //// The int value representing the static modifier.
    static int const STATIC;

    //// The int value representing the transient modifier.
    static int const TRANSIENT;

    //// The int value representing the export modifier.
    static int const EXPORT;

#if 0
    //// The int value representing the volatile modifier.
    static int const VOLATILE;

#endif
    ////Return true if the specifier integer includes the abstract modifier.
    static bool isAbstract(int mod);

#if 0
    //// Return true if the specified integer includes the final modifier.
    static bool isFinal(int mod);
#endif

    //// Return true if the specifier integer includes the interface modifier.
    static bool isInterface(int mod);

    //// Return true if the specifier integer includes the private modifier.
    static bool isPrivate(int mod);

    //// Return true if the specifier integer includes the protected modifier.
    static bool isProtected(int mod);

    //// Return true if the specified integer includes the public modifier.
    static bool isPublic(int mod);

    //// Return true if the specifier integer includes the static modifier.
    static bool isStatic(int mod);

    //// Return true if the specifier integer includes the transient modifier.
    static bool isTransient(int mod);

    //// Return true if the specifier integer includes the export modifier.
    static bool isExport(int mod);

#if 0
    //// Return true if the specified integer includes the volatile modifier.
    static bool isVolatile(int mod);
#endif
    /** Return a string describing the access modifier flags in the
        specified modifier. */
    static string toString(int mod);

private:
    Modifier();
};


DRLIB_END_NAMESPACE
#endif
