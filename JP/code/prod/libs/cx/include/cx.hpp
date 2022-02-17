/*
***************************************************************************
** HEADER FILE: cx.hpp
**
** Provides mini C++ classes wrapping up the C data structures.
***************************************************************************
*/

#ifndef _CX_CX_HPP
#define _CX_CX_HPP

#include "cx.h"
#include "creditcurve.h"
#include "recovery.h"

#include <alib/convert.h>
#include <alib/dtivlo.h>
#include <alib/tcurve.h>
#include <cxutils/include/cxmacros.h>
#include <exception>
#include <vector>

using namespace std;

#ifndef __STRING_EXCEPTION__
#define __STRING_EXCEPTION__

    class StringException : public exception
    {
    private:
        const StringException& operator= (const StringException&);
    public:
        StringException() : errors(0) {}
        StringException(const string& msg) : errors(1,msg) {}
        StringException(exception &e, const string &routine) : 
            errors(1, routine + ": " + e.what()) {}
        StringException(const StringException& e) { errors = e.errors; }
        virtual ~StringException() throw() {}
        void addMsg(const string &msg) { errors.push_back(msg); }
        virtual const char* what() const throw() {
            if (errors.size() > 0) return errors[0].c_str();
            return "";
        }
    private:
        vector<string> errors;
    };

#endif /* __STRING_EXCEPTION__ */

class CxTRecoveryCurveSP {
public:
    CxTRecoveryCurveSP() : p(0) {}
    CxTRecoveryCurveSP(CxTRecoveryCurve* ptr) : p(ptr) {}
    CxTRecoveryCurveSP(const CxTRecoveryCurve* ptr) {copy(ptr);}
    CxTRecoveryCurveSP(const CxTRecoveryCurveSP &src) {copy(src.p);}
    const CxTRecoveryCurveSP& operator= (const CxTRecoveryCurveSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTRecoveryCurveSP() {deallocate();}
    CxTRecoveryCurve* get() {return p;}
//  CxTRecoveryCurve* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTRecoveryCurve *p;
    void deallocate() {CxRecoveryCurveFree(p); p=0; }
    void copy(const CxTRecoveryCurve*ptr) {
        if (ptr) {
            p = CxRecoveryCurveCopy((CxTRecoveryCurve*)ptr);
            if (!p) throw StringException ("Could not copy CxTRecoveryCurve");
        } else {
            p = 0;
        }
    }
};

class CxTCdsCurveDateAdjSP {
public:
    CxTCdsCurveDateAdjSP() : p(0) {}
    CxTCdsCurveDateAdjSP(CxTCdsCurveDateAdj* ptr) : p(ptr) {}
    CxTCdsCurveDateAdjSP(const CxTCdsCurveDateAdj* ptr) {copy(ptr);}
    CxTCdsCurveDateAdjSP(const CxTCdsCurveDateAdjSP &src) {copy(src.p);}
    const CxTCdsCurveDateAdjSP& operator= (const CxTCdsCurveDateAdjSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTCdsCurveDateAdjSP() {deallocate();}
    CxTCdsCurveDateAdj* get() {return p;}
//  CxTCdsCurveDateAdj* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTCdsCurveDateAdj *p;
    void deallocate() {CxCdsCurveDateAdjFree(p); p=0; }
    void copy(const CxTCdsCurveDateAdj*ptr) {
        if (ptr) {
            p = CxCdsCurveDateAdjCopy((CxTCdsCurveDateAdj*)ptr);
            if (!p) throw StringException ("Could not copy CxTCdsCurveDateAdj");
        } else {
            p = 0;
        }
    }
};

class CxTFeeLegSP {
public:
    CxTFeeLegSP() : p(0) {}
    CxTFeeLegSP(CxTFeeLeg* ptr) : p(ptr) {}
    CxTFeeLegSP(const CxTFeeLeg* ptr) {copy(ptr);}
    CxTFeeLegSP(const CxTFeeLegSP &src) {copy(src.p);}
    const CxTFeeLegSP& operator= (const CxTFeeLegSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTFeeLegSP() {deallocate();}
    CxTFeeLeg* get() {return p;}
//  CxTFeeLeg* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTFeeLeg *p;
    void deallocate() {CxFeeLegFree(p); p=0; }
    void copy(const CxTFeeLeg*ptr) {
        if (ptr) {
            p = CxFeeLegCopy((CxTFeeLeg*)ptr);
            if (!p) throw StringException ("Could not copy CxTFeeLeg");
        } else {
            p = 0;
        }
    }
};

class CxTCdsSP {
public:
    CxTCdsSP() : p(0) {}
    CxTCdsSP(CxTCds* ptr) : p(ptr) {}
    CxTCdsSP(const CxTCds* ptr) {copy(ptr);}
    CxTCdsSP(const CxTCdsSP &src) {copy(src.p);}
    const CxTCdsSP& operator= (const CxTCdsSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTCdsSP() {deallocate();}
    CxTCds* get() {return p;}
//  CxTCds* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTCds *p;
    void deallocate() {CxCdsFree(p); p=0; }
    void copy(const CxTCds*ptr) {
        if (ptr) {
            p = CxCdsCopy((CxTCds*)ptr);
            if (!p) throw StringException ("Could not copy CxTCds");
        } else {
            p = 0;
        }
    }
};

class CxTIrMarketDataSP {
public:
    CxTIrMarketDataSP() : p(0) {}
    CxTIrMarketDataSP(CxTIrMarketData* ptr) : p(ptr) {}
    CxTIrMarketDataSP(const CxTIrMarketData* ptr) {copy(ptr);}
    CxTIrMarketDataSP(const CxTIrMarketDataSP &src) {copy(src.p);}
    const CxTIrMarketDataSP& operator= (const CxTIrMarketDataSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTIrMarketDataSP() {deallocate();}
    CxTIrMarketData* get() {return p;}
//  CxTIrMarketData* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTIrMarketData *p;
    void deallocate() {CxIrMarketDataFree(p); p=0; }
    void copy(const CxTIrMarketData*ptr) {
        if (ptr) {
            p = CxIrMarketDataCopy((CxTIrMarketData*)ptr);
            if (!p) throw StringException ("Could not copy CxTIrMarketData");
        } else {
            p = 0;
        }
    }
};

class CxTCreditMarketDataSP {
public:
    CxTCreditMarketDataSP() : p(0) {}
    CxTCreditMarketDataSP(CxTCreditMarketData* ptr) : p(ptr) {}
    CxTCreditMarketDataSP(const CxTCreditMarketData* ptr) {copy(ptr);}
    CxTCreditMarketDataSP(const CxTCreditMarketDataSP &src) {copy(src.p);}
    const CxTCreditMarketDataSP& operator= (const CxTCreditMarketDataSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTCreditMarketDataSP() {deallocate();}
    CxTCreditMarketData* get() {return p;}
//  CxTCreditMarketData* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTCreditMarketData *p;
    void deallocate() {CxCreditMarketDataFree(p); p=0; }
    void copy(const CxTCreditMarketData*ptr) {
        if (ptr) {
            p = CxCreditMarketDataCopy((CxTCreditMarketData*)ptr);
            if (!p) throw StringException ("Could not copy CxTCreditMarketData");
        } else {
            p = 0;
        }
    }
};

class CxTCdsCarryOutputsSP {
public:
    CxTCdsCarryOutputsSP() : p(0) {}
    CxTCdsCarryOutputsSP(CxTCdsCarryOutputs* ptr) : p(ptr) {}
    CxTCdsCarryOutputsSP(const CxTCdsCarryOutputs* ptr) {copy(ptr);}
    CxTCdsCarryOutputsSP(const CxTCdsCarryOutputsSP &src) {copy(src.p);}
    const CxTCdsCarryOutputsSP& operator= (const CxTCdsCarryOutputsSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTCdsCarryOutputsSP() {deallocate();}
    CxTCdsCarryOutputs* get() {return p;}
//  CxTCdsCarryOutputs* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTCdsCarryOutputs *p;
    void deallocate() {CxCdsCarryOutputsFree(p); p=0; }
    void copy(const CxTCdsCarryOutputs*ptr) {
        if (ptr) {
            p = CxCdsCarryOutputsCopy((CxTCdsCarryOutputs*)ptr);
            if (!p) throw StringException ("Could not copy CxTCdsCarryOutputs");
        } else {
            p = 0;
        }
    }
};

class CxTContingentLegSP {
public:
    CxTContingentLegSP() : p(0) {}
    CxTContingentLegSP(CxTContingentLeg* ptr) : p(ptr) {}
    CxTContingentLegSP(const CxTContingentLeg* ptr) {copy(ptr);}
    CxTContingentLegSP(const CxTContingentLegSP &src) {copy(src.p);}
    const CxTContingentLegSP& operator= (const CxTContingentLegSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTContingentLegSP() {deallocate();}
    CxTContingentLeg* get() {return p;}
//  CxTContingentLeg* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTContingentLeg *p;
    void deallocate() {CxContingentLegFree(p); p=0; }
    void copy(const CxTContingentLeg*ptr) {
        if (ptr) {
            p = CxContingentLegCopy((CxTContingentLeg*)ptr);
            if (!p) throw StringException ("Could not copy CxTContingentLeg");
        } else {
            p = 0;
        }
    }
};

class CxTCreditCurveSP {
public:
    CxTCreditCurveSP() : p(0) {}
    CxTCreditCurveSP(CxTCreditCurve* ptr) : p(ptr) {}
    CxTCreditCurveSP(const CxTCreditCurve* ptr) {copy(ptr);}
    CxTCreditCurveSP(const CxTCreditCurveSP &src) {copy(src.p);}
    const CxTCreditCurveSP& operator= (const CxTCreditCurveSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTCreditCurveSP() {deallocate();}
    CxTCreditCurve* get() {return p;}
//  CxTCreditCurve* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTCreditCurve *p;
    void deallocate() {CxCreditCurveFree(p); p=0; }
    void copy(const CxTCreditCurve*ptr) {
        if (ptr) {
            p = CxCreditCurveCopy((CxTCreditCurve*)ptr);
            if (!p) throw StringException ("Could not copy CxTCreditCurve");
        } else {
            p = 0;
        }
    }
};

class CxTCdsTypeSP {
public:
    CxTCdsTypeSP() : p(0) {}
    CxTCdsTypeSP(CxTCdsType* ptr) : p(ptr) {}
    CxTCdsTypeSP(const CxTCdsType* ptr) {copy(ptr);}
    CxTCdsTypeSP(const CxTCdsTypeSP &src) {copy(src.p);}
    const CxTCdsTypeSP& operator= (const CxTCdsTypeSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTCdsTypeSP() {deallocate();}
    CxTCdsType* get() {return p;}
//  CxTCdsType* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTCdsType *p;
    void deallocate() {CxCdsTypeFree(p); p=0; }
    void copy(const CxTCdsType*ptr) {
        if (ptr) {
            p = CxCdsTypeCopy((CxTCdsType*)ptr);
            if (!p) throw StringException ("Could not copy CxTCdsType");
        } else {
            p = 0;
        }
    }
};

class CxTCdsConventionsSP {
public:
    CxTCdsConventionsSP() : p(0) {}
    CxTCdsConventionsSP(CxTCdsConventions* ptr) : p(ptr) {}
    CxTCdsConventionsSP(const CxTCdsConventions* ptr) {copy(ptr);}
    CxTCdsConventionsSP(const CxTCdsConventionsSP &src) {copy(src.p);}
    const CxTCdsConventionsSP& operator= (const CxTCdsConventionsSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTCdsConventionsSP() {deallocate();}
    CxTCdsConventions* get() {return p;}
//  CxTCdsConventions* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTCdsConventions *p;
    void deallocate() {CxCdsConventionsFree(p); p=0; }
    void copy(const CxTCdsConventions*ptr) {
        if (ptr) {
            p = CxCdsConventionsCopy((CxTCdsConventions*)ptr);
            if (!p) throw StringException ("Could not copy CxTCdsConventions");
        } else {
            p = 0;
        }
    }
};
#endif /* _CX_CX_HPP */
