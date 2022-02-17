//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : SensitivityFactory.hpp
//
//   Description : Class for managing instantiation of Sensitivities
//
//   Author      : Andrew J Swain
//
//   Date        : 2 November 2001
//
//
//----------------------------------------------------------------------------
#ifndef _SENSITIVITYFACTORY_HPP
#define _SENSITIVITYFACTORY_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE
class Sensitivity;
class Control;

/** Support ability to build sensitivities based upon their name with
    either default parameters or with a given shift size */
class RISKMGR_DLL SensitivityFactory {
public:
    // list all the different types of construction we might support

    /** marker interface for all type of creation supported */
    class RISKMGR_DLL ICreation{
    public:
        ICreation();
        // empty except for virtual destructor (forces class to be polymorphic)
        virtual ~ICreation();
    };

    /** Create a default sensitivity */
    class RISKMGR_DLL IDefault: virtual public ICreation{
    public:
        IDefault();
        virtual ~IDefault();
        virtual Sensitivity* createDefault() = 0;
    };

    /** Create a sensitivity characterised by a scalar shift */
    class RISKMGR_DLL IScalar: virtual public ICreation{
    public:
        IScalar();
        virtual ~IScalar();
        virtual Sensitivity* createScalar(double shiftSize) = 0;
    };

    /** Register a sensitivity with the factory. Method takes ownership of 
        Creation* */
    static void addSens(
        const string& name,            // name for this sensitivity eg DELTA
        ICreation*    factory,         // capable of building Sensitivity
        Sensitivity*  example,         // any old instance of the sensitivity
        CClassConstSP shiftInterface); /* the associate shift interface for
                                          this tweak eg Delta::Shift */

    static Sensitivity* defaultSensitivity(const string& sensName);

    static Sensitivity* scalarSensitivity(const string& sensName,
                                          double        shiftSize);

    /** create a control with as much defaulted as possible */
    static Control* megaControl();
    static Control* megaControl(bool requests); // false=no defaulted requests

private:
    SensitivityFactory();
    SensitivityFactory(const SensitivityFactory &rhs);
    SensitivityFactory& operator=(const SensitivityFactory& rhs);
};

DRLIB_END_NAMESPACE

#endif
