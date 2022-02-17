//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XLGetCorrelations.cpp
//
//   Description : Addin which retrieves FWD_AT_MAT from a Results handle
//
//   Author      : André Segger
//
//   Date        : 04 Oct 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Addin.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE

/* non static variable to force linkage of this file. Outside of class to
   avoid necessity of header file */
bool XLGetOutputPacketAddinsRegistered = true;

/** Addin to determine run-time type of a handle */
class XLGetOutputPacketAddin: public CObject{
    static CClassConstSP const TYPE;

    /**  parameters that our addin takes */
    ResultsSP     results;
    string        packetName;

    /** the 'addin function' - construct object from components */
    static IObjectSP getFwdAtMat(XLGetOutputPacketAddin* params){
        static const string routine = "XLGetOutputPacketAddin::getFwdAtMat";

        ObjectArraySP fwdAtMat(new ObjectArray(2)); // 2 columns
        ResultsSP     results    = params->results;
        string        packetName = params->packetName;

        OutputNameArraySP underlyingNames = results->packetContents(
                                                    packetName);

        CStringArraySP underlyingIdentifiers(new CStringArray(0));

        CDoubleArraySP  fwdValues(new CDoubleArray(0));

        int i;
        for (i=0; i<underlyingNames->size() ; ++i) {
            fwdValues->push_back(results->retrieveScalarGreek(
                                            packetName,
                                            (*underlyingNames)[i]));

            underlyingIdentifiers->push_back((*underlyingNames)[i]->toString());    
        }

        (*fwdAtMat)[0] = underlyingIdentifiers;
        (*fwdAtMat)[1] = fwdValues;

        return fwdAtMat;
    }

    /** for reflection */
    XLGetOutputPacketAddin():  CObject(TYPE){}

    /** Invoked when Class is 'loaded' */
    static void load(CClassSP& clazz) {
        REGISTER(XLGetOutputPacketAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultXLGetOutputPacketAddin);
        FIELD(results, "Results handle");
        FIELD(packetName,  "Output packet name, eg. DELTA, FWD_AT_MAT, etc.");

        Addin::registerClassObjectMethod("GET_OUTPUT_PACKET",
                                         Addin::UTILITIES,
                                         "retrieves outputs from an output handle",
                                         TYPE,
                                         false,
                                         Addin::expandMulti,
                                         (Addin::ObjMethod*)getFwdAtMat);
    }

    static IObject* defaultXLGetOutputPacketAddin() {
        return new XLGetOutputPacketAddin();
    }
};

CClassConstSP const XLGetOutputPacketAddin::TYPE = CClass::registerClassLoadMethod(
    "XLGetOutputPacketAddin", typeid(XLGetOutputPacketAddin), load);

DRLIB_END_NAMESPACE
