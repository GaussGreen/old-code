/*
***************************************************************************
** HEADER FILE: crxdata.hpp
**
** Provides mini C++ classes wrapping up the C data structures.
**
** Defines data structures for crxflow library.
***************************************************************************
*/

#ifndef _CRX_CRXDATA_HPP
#define _CRX_CRXDATA_HPP

#include "crxdata.h"
#include "crxmacros.h"

#include <alib/bondcnst.h>
#include <alib/convert.h>
#include <alib/dtivlo.h>
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

class CrxTBondSpreadVolCurveSP {
public:
    CrxTBondSpreadVolCurveSP() : p(0) {}
    CrxTBondSpreadVolCurveSP(CrxTBondSpreadVolCurve* ptr) : p(ptr) {}
    CrxTBondSpreadVolCurveSP(const CrxTBondSpreadVolCurve* ptr) {copy(ptr);}
    CrxTBondSpreadVolCurveSP(const CrxTBondSpreadVolCurveSP &src) {copy(src.p);}
    const CrxTBondSpreadVolCurveSP& operator= (const CrxTBondSpreadVolCurveSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTBondSpreadVolCurveSP() {deallocate();}
    CrxTBondSpreadVolCurve* get() {return p;}
//  CrxTBondSpreadVolCurve* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTBondSpreadVolCurve *p;
    void deallocate() {CrxBondSpreadVolCurveFree(p); p=0; }
    void copy(const CrxTBondSpreadVolCurve*ptr) {
        if (ptr) {
            p = CrxBondSpreadVolCurveCopy((CrxTBondSpreadVolCurve*)ptr);
            if (!p) throw StringException ("Could not copy CrxTBondSpreadVolCurve");
        } else {
            p = 0;
        }
    }
};

class CrxTCdsOptionQOptimizationSP {
public:
    CrxTCdsOptionQOptimizationSP() : p(0) {}
    CrxTCdsOptionQOptimizationSP(CrxTCdsOptionQOptimization* ptr) : p(ptr) {}
    CrxTCdsOptionQOptimizationSP(const CrxTCdsOptionQOptimization* ptr) {copy(ptr);}
    CrxTCdsOptionQOptimizationSP(const CrxTCdsOptionQOptimizationSP &src) {copy(src.p);}
    const CrxTCdsOptionQOptimizationSP& operator= (const CrxTCdsOptionQOptimizationSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTCdsOptionQOptimizationSP() {deallocate();}
    CrxTCdsOptionQOptimization* get() {return p;}
//  CrxTCdsOptionQOptimization* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTCdsOptionQOptimization *p;
    void deallocate() {CrxCdsOptionQOptimizationFree(p); p=0; }
    void copy(const CrxTCdsOptionQOptimization*ptr) {
        if (ptr) {
            p = CrxCdsOptionQOptimizationCopy((CrxTCdsOptionQOptimization*)ptr);
            if (!p) throw StringException ("Could not copy CrxTCdsOptionQOptimization");
        } else {
            p = 0;
        }
    }
};

class CrxTCdsOptionQOptModelSP {
public:
    CrxTCdsOptionQOptModelSP() : p(0) {}
    CrxTCdsOptionQOptModelSP(CrxTCdsOptionQOptModel* ptr) : p(ptr) {}
    CrxTCdsOptionQOptModelSP(const CrxTCdsOptionQOptModel* ptr) {copy(ptr);}
    CrxTCdsOptionQOptModelSP(const CrxTCdsOptionQOptModelSP &src) {copy(src.p);}
    const CrxTCdsOptionQOptModelSP& operator= (const CrxTCdsOptionQOptModelSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTCdsOptionQOptModelSP() {deallocate();}
    CrxTCdsOptionQOptModel* get() {return p;}
//  CrxTCdsOptionQOptModel* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTCdsOptionQOptModel *p;
    void deallocate() {CrxCdsOptionQOptModelFree(p); p=0; }
    void copy(const CrxTCdsOptionQOptModel*ptr) {
        if (ptr) {
            p = CrxCdsOptionQOptModelCopy((CrxTCdsOptionQOptModel*)ptr);
            if (!p) throw StringException ("Could not copy CrxTCdsOptionQOptModel");
        } else {
            p = 0;
        }
    }
};

class CrxTBondPriceOptionCalcSP {
public:
    CrxTBondPriceOptionCalcSP() : p(0) {}
    CrxTBondPriceOptionCalcSP(CrxTBondPriceOptionCalc* ptr) : p(ptr) {}
    CrxTBondPriceOptionCalcSP(const CrxTBondPriceOptionCalc* ptr) {copy(ptr);}
    CrxTBondPriceOptionCalcSP(const CrxTBondPriceOptionCalcSP &src) {copy(src.p);}
    const CrxTBondPriceOptionCalcSP& operator= (const CrxTBondPriceOptionCalcSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTBondPriceOptionCalcSP() {deallocate();}
    CrxTBondPriceOptionCalc* get() {return p;}
//  CrxTBondPriceOptionCalc* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTBondPriceOptionCalc *p;
    void deallocate() {CrxBondPriceOptionCalcFree(p); p=0; }
    void copy(const CrxTBondPriceOptionCalc*ptr) {
        if (ptr) {
            p = CrxBondPriceOptionCalcCopy((CrxTBondPriceOptionCalc*)ptr);
            if (!p) throw StringException ("Could not copy CrxTBondPriceOptionCalc");
        } else {
            p = 0;
        }
    }
};

class CrxTBondSpreadOptionSP {
public:
    CrxTBondSpreadOptionSP() : p(0) {}
    CrxTBondSpreadOptionSP(CrxTBondSpreadOption* ptr) : p(ptr) {}
    CrxTBondSpreadOptionSP(const CrxTBondSpreadOption* ptr) {copy(ptr);}
    CrxTBondSpreadOptionSP(const CrxTBondSpreadOptionSP &src) {copy(src.p);}
    const CrxTBondSpreadOptionSP& operator= (const CrxTBondSpreadOptionSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTBondSpreadOptionSP() {deallocate();}
    CrxTBondSpreadOption* get() {return p;}
//  CrxTBondSpreadOption* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTBondSpreadOption *p;
    void deallocate() {CrxBondSpreadOptionFree(p); p=0; }
    void copy(const CrxTBondSpreadOption*ptr) {
        if (ptr) {
            p = CrxBondSpreadOptionCopy((CrxTBondSpreadOption*)ptr);
            if (!p) throw StringException ("Could not copy CrxTBondSpreadOption");
        } else {
            p = 0;
        }
    }
};

class CrxTCdsOptionSP {
public:
    CrxTCdsOptionSP() : p(0) {}
    CrxTCdsOptionSP(CrxTCdsOption* ptr) : p(ptr) {}
    CrxTCdsOptionSP(const CrxTCdsOption* ptr) {copy(ptr);}
    CrxTCdsOptionSP(const CrxTCdsOptionSP &src) {copy(src.p);}
    const CrxTCdsOptionSP& operator= (const CrxTCdsOptionSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTCdsOptionSP() {deallocate();}
    CrxTCdsOption* get() {return p;}
//  CrxTCdsOption* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTCdsOption *p;
    void deallocate() {CrxCdsOptionFree(p); p=0; }
    void copy(const CrxTCdsOption*ptr) {
        if (ptr) {
            p = CrxCdsOptionCopy((CrxTCdsOption*)ptr);
            if (!p) throw StringException ("Could not copy CrxTCdsOption");
        } else {
            p = 0;
        }
    }
};

class CrxTBondRepoCurveSP {
public:
    CrxTBondRepoCurveSP() : p(0) {}
    CrxTBondRepoCurveSP(CrxTBondRepoCurve* ptr) : p(ptr) {}
    CrxTBondRepoCurveSP(const CrxTBondRepoCurve* ptr) {copy(ptr);}
    CrxTBondRepoCurveSP(const CrxTBondRepoCurveSP &src) {copy(src.p);}
    const CrxTBondRepoCurveSP& operator= (const CrxTBondRepoCurveSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTBondRepoCurveSP() {deallocate();}
    CrxTBondRepoCurve* get() {return p;}
//  CrxTBondRepoCurve* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTBondRepoCurve *p;
    void deallocate() {CrxBondRepoCurveFree(p); p=0; }
    void copy(const CrxTBondRepoCurve*ptr) {
        if (ptr) {
            p = CrxBondRepoCurveCopy((CrxTBondRepoCurve*)ptr);
            if (!p) throw StringException ("Could not copy CrxTBondRepoCurve");
        } else {
            p = 0;
        }
    }
};

class CrxTQDistSP {
public:
    CrxTQDistSP() : p(0) {}
    CrxTQDistSP(CrxTQDist* ptr) : p(ptr) {}
    CrxTQDistSP(const CrxTQDist* ptr) {copy(ptr);}
    CrxTQDistSP(const CrxTQDistSP &src) {copy(src.p);}
    const CrxTQDistSP& operator= (const CrxTQDistSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTQDistSP() {deallocate();}
    CrxTQDist* get() {return p;}
//  CrxTQDist* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTQDist *p;
    void deallocate() {CrxQDistFree(p); p=0; }
    void copy(const CrxTQDist*ptr) {
        if (ptr) {
            p = CrxQDistCopy((CrxTQDist*)ptr);
            if (!p) throw StringException ("Could not copy CrxTQDist");
        } else {
            p = 0;
        }
    }
};

class CrxTBondPriceVolCurveSP {
public:
    CrxTBondPriceVolCurveSP() : p(0) {}
    CrxTBondPriceVolCurveSP(CrxTBondPriceVolCurve* ptr) : p(ptr) {}
    CrxTBondPriceVolCurveSP(const CrxTBondPriceVolCurve* ptr) {copy(ptr);}
    CrxTBondPriceVolCurveSP(const CrxTBondPriceVolCurveSP &src) {copy(src.p);}
    const CrxTBondPriceVolCurveSP& operator= (const CrxTBondPriceVolCurveSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTBondPriceVolCurveSP() {deallocate();}
    CrxTBondPriceVolCurve* get() {return p;}
//  CrxTBondPriceVolCurve* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTBondPriceVolCurve *p;
    void deallocate() {CrxBondPriceVolCurveFree(p); p=0; }
    void copy(const CrxTBondPriceVolCurve*ptr) {
        if (ptr) {
            p = CrxBondPriceVolCurveCopy((CrxTBondPriceVolCurve*)ptr);
            if (!p) throw StringException ("Could not copy CrxTBondPriceVolCurve");
        } else {
            p = 0;
        }
    }
};

class CrxTBondSpreadOptionCalcSP {
public:
    CrxTBondSpreadOptionCalcSP() : p(0) {}
    CrxTBondSpreadOptionCalcSP(CrxTBondSpreadOptionCalc* ptr) : p(ptr) {}
    CrxTBondSpreadOptionCalcSP(const CrxTBondSpreadOptionCalc* ptr) {copy(ptr);}
    CrxTBondSpreadOptionCalcSP(const CrxTBondSpreadOptionCalcSP &src) {copy(src.p);}
    const CrxTBondSpreadOptionCalcSP& operator= (const CrxTBondSpreadOptionCalcSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTBondSpreadOptionCalcSP() {deallocate();}
    CrxTBondSpreadOptionCalc* get() {return p;}
//  CrxTBondSpreadOptionCalc* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTBondSpreadOptionCalc *p;
    void deallocate() {CrxBondSpreadOptionCalcFree(p); p=0; }
    void copy(const CrxTBondSpreadOptionCalc*ptr) {
        if (ptr) {
            p = CrxBondSpreadOptionCalcCopy((CrxTBondSpreadOptionCalc*)ptr);
            if (!p) throw StringException ("Could not copy CrxTBondSpreadOptionCalc");
        } else {
            p = 0;
        }
    }
};

class CrxTCdsOptionCalcSP {
public:
    CrxTCdsOptionCalcSP() : p(0) {}
    CrxTCdsOptionCalcSP(CrxTCdsOptionCalc* ptr) : p(ptr) {}
    CrxTCdsOptionCalcSP(const CrxTCdsOptionCalc* ptr) {copy(ptr);}
    CrxTCdsOptionCalcSP(const CrxTCdsOptionCalcSP &src) {copy(src.p);}
    const CrxTCdsOptionCalcSP& operator= (const CrxTCdsOptionCalcSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTCdsOptionCalcSP() {deallocate();}
    CrxTCdsOptionCalc* get() {return p;}
//  CrxTCdsOptionCalc* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTCdsOptionCalc *p;
    void deallocate() {CrxCdsOptionCalcFree(p); p=0; }
    void copy(const CrxTCdsOptionCalc*ptr) {
        if (ptr) {
            p = CrxCdsOptionCalcCopy((CrxTCdsOptionCalc*)ptr);
            if (!p) throw StringException ("Could not copy CrxTCdsOptionCalc");
        } else {
            p = 0;
        }
    }
};

class CrxTBondPriceOptionSP {
public:
    CrxTBondPriceOptionSP() : p(0) {}
    CrxTBondPriceOptionSP(CrxTBondPriceOption* ptr) : p(ptr) {}
    CrxTBondPriceOptionSP(const CrxTBondPriceOption* ptr) {copy(ptr);}
    CrxTBondPriceOptionSP(const CrxTBondPriceOptionSP &src) {copy(src.p);}
    const CrxTBondPriceOptionSP& operator= (const CrxTBondPriceOptionSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTBondPriceOptionSP() {deallocate();}
    CrxTBondPriceOption* get() {return p;}
//  CrxTBondPriceOption* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTBondPriceOption *p;
    void deallocate() {CrxBondPriceOptionFree(p); p=0; }
    void copy(const CrxTBondPriceOption*ptr) {
        if (ptr) {
            p = CrxBondPriceOptionCopy((CrxTBondPriceOption*)ptr);
            if (!p) throw StringException ("Could not copy CrxTBondPriceOption");
        } else {
            p = 0;
        }
    }
};

class CrxTBondPriceSP {
public:
    CrxTBondPriceSP() : p(0) {}
    CrxTBondPriceSP(CrxTBondPrice* ptr) : p(ptr) {}
    CrxTBondPriceSP(const CrxTBondPrice* ptr) {copy(ptr);}
    CrxTBondPriceSP(const CrxTBondPriceSP &src) {copy(src.p);}
    const CrxTBondPriceSP& operator= (const CrxTBondPriceSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CrxTBondPriceSP() {deallocate();}
    CrxTBondPrice* get() {return p;}
//  CrxTBondPrice* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CrxTBondPrice *p;
    void deallocate() {CrxBondPriceFree(p); p=0; }
    void copy(const CrxTBondPrice*ptr) {
        if (ptr) {
            p = CrxBondPriceCopy((CrxTBondPrice*)ptr);
            if (!p) throw StringException ("Could not copy CrxTBondPrice");
        } else {
            p = 0;
        }
    }
};
#endif /* _CRX_CRXDATA_HPP */
