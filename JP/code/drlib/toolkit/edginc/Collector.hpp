//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : Coolector.hpp
//
//   Description : Prototype collector
//
//   Author      : Andrew J Swain
//
//   Date        : 20 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef COLLECTOR_HPP
#define COLLECTOR_HPP
#include "edginc/smartPtr.hpp"
#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

/** Collector is a mechanism to be able to collect data from a variety
    of sources. The source of the data (ie which function is invoked to
    collect the data) is a function of the type of the collector and the
    the type of the object being queried
    There wil be different types of Collector for different aspects
    of data and classes
    
    This interface is just a marker interface for signalling the types of
    object that the accept method can take  */
class TOOLKIT_DLL ICollector: public virtual IObject {
public:
    /** the class representing the IPrivateObject interface */
    static CClassConstSP const TYPE; // in Object.cpp
};

DRLIB_END_NAMESPACE

#endif
