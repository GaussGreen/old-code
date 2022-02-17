//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : EDRTester.hpp
//
//   Description : Models Library Interface Tester
//
//   Author      : Andrew J Swain
//
//   Date        : 15 June 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDRTESTER_HPP
#define EDRTESTER_HPP

#include "edginc/DoubleMatrix.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/Atomic.hpp"

DRLIB_BEGIN_NAMESPACE

class ADDINS_DLL EDRTester: public CObject {

public:
    static CClassConstSP const TYPE;

    class ADDINS_DLL Atomic: public CObject {
    public:
        static CClassConstSP const TYPE;
        
        static IObjectSP run(Atomic* params);

    private:
        int    i;
        double d;
        string s;
        int    b;

        // dbl matrix
        DoubleMatrix dm;

        // datetime
        string dt;
        string tm;

        // expiries
        string matp;
        string matpdt;
        string matptm;
        string bmdt;
        string bmtm;
        CStringSP enumType;
        CStringSP enumValue;
        Atomic();
        friend class EDRAtomicHelper;
    };

    class ADDINS_DLL Array: public CObject {
    public:
        static CClassConstSP const TYPE;       
        static IObjectSP run(Array* params);

    private:
        // bit crap, only tests array of doubles
        DoubleArraySP dbls;
        double        extra;  // gets appended

        Array();
        friend class EDRArrayHelper;
    };

    class ADDINS_DLL TypeInfo: public CObject {
    public:
        static CClassConstSP const TYPE;       
        static IObjectSP run(TypeInfo* params);

    private:
        string typeName;

        TypeInfo();
        friend class EDRTypeInfoHelper;
    };


    class ClassDescription;
private:
    class PublicInterface;
    friend class EDRTesterHelper;
    EDRTester();
};

DRLIB_END_NAMESPACE

#endif
