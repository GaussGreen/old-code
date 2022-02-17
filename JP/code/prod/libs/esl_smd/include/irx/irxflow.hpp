/*
***************************************************************************
** HEADER FILE: irxflow.hpp
**
** Provides mini C++ classes wrapping up the C data structures.
**
** Defines data structures and enums used in the irxflow library.
***************************************************************************
*/

#ifndef _IRX_IRXFLOW_HPP
#define _IRX_IRXFLOW_HPP

#include "irxflow.h"
#include "mktconv.h"
#include "rate.h"
#include "swap.h"

#include <irxutils/include/macros.h>
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

class IrxTSwapSP {
public:
    IrxTSwapSP() : p(0) {}
    IrxTSwapSP(IrxTSwap* ptr) : p(ptr) {}
    IrxTSwapSP(const IrxTSwap* ptr) {copy(ptr);}
    IrxTSwapSP(const IrxTSwapSP &src) {copy(src.p);}
    const IrxTSwapSP& operator= (const IrxTSwapSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~IrxTSwapSP() {deallocate();}
    IrxTSwap* get() {return p;}
//  IrxTSwap* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    IrxTSwap *p;
    void deallocate() {irxSwapFree(p); p=0; }
    void copy(const IrxTSwap*ptr) {
        if (ptr) {
            p = irxSwapCopy((IrxTSwap*)ptr);
            if (!p) throw StringException ("Could not copy IrxTSwap");
        } else {
            p = 0;
        }
    }
};

class IrxTSwapPaymentsSP {
public:
    IrxTSwapPaymentsSP() : p(0) {}
    IrxTSwapPaymentsSP(IrxTSwapPayments* ptr) : p(ptr) {}
    IrxTSwapPaymentsSP(const IrxTSwapPayments* ptr) {copy(ptr);}
    IrxTSwapPaymentsSP(const IrxTSwapPaymentsSP &src) {copy(src.p);}
    const IrxTSwapPaymentsSP& operator= (const IrxTSwapPaymentsSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~IrxTSwapPaymentsSP() {deallocate();}
    IrxTSwapPayments* get() {return p;}
//  IrxTSwapPayments* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    IrxTSwapPayments *p;
    void deallocate() {irxSwapPaymentsFree(p); p=0; }
    void copy(const IrxTSwapPayments*ptr) {
        if (ptr) {
            p = irxSwapPaymentsCopy((IrxTSwapPayments*)ptr);
            if (!p) throw StringException ("Could not copy IrxTSwapPayments");
        } else {
            p = 0;
        }
    }
};

class IrxTZeroCurveSP {
public:
    IrxTZeroCurveSP() : p(0) {}
    IrxTZeroCurveSP(IrxTZeroCurve* ptr) : p(ptr) {}
    IrxTZeroCurveSP(const IrxTZeroCurve* ptr) {copy(ptr);}
    IrxTZeroCurveSP(const IrxTZeroCurveSP &src) {copy(src.p);}
    const IrxTZeroCurveSP& operator= (const IrxTZeroCurveSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~IrxTZeroCurveSP() {deallocate();}
    IrxTZeroCurve* get() {return p;}
//  IrxTZeroCurve* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    IrxTZeroCurve *p;
    void deallocate() {irxZeroCurveFree(p); p=0; }
    void copy(const IrxTZeroCurve*ptr) {
        if (ptr) {
            p = irxZeroCurveCopy((IrxTZeroCurve*)ptr);
            if (!p) throw StringException ("Could not copy IrxTZeroCurve");
        } else {
            p = 0;
        }
    }
};

class IrxTCashFlowListSP {
public:
    IrxTCashFlowListSP() : p(0) {}
    IrxTCashFlowListSP(IrxTCashFlowList* ptr) : p(ptr) {}
    IrxTCashFlowListSP(const IrxTCashFlowList* ptr) {copy(ptr);}
    IrxTCashFlowListSP(const IrxTCashFlowListSP &src) {copy(src.p);}
    const IrxTCashFlowListSP& operator= (const IrxTCashFlowListSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~IrxTCashFlowListSP() {deallocate();}
    IrxTCashFlowList* get() {return p;}
//  IrxTCashFlowList* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    IrxTCashFlowList *p;
    void deallocate() {irxCashFlowListFree(p); p=0; }
    void copy(const IrxTCashFlowList*ptr) {
        if (ptr) {
            p = irxCashFlowListCopy((IrxTCashFlowList*)ptr);
            if (!p) throw StringException ("Could not copy IrxTCashFlowList");
        } else {
            p = 0;
        }
    }
};

class IrxTMarketConvSP {
public:
    IrxTMarketConvSP() : p(0) {}
    IrxTMarketConvSP(IrxTMarketConv* ptr) : p(ptr) {}
    IrxTMarketConvSP(const IrxTMarketConv* ptr) {copy(ptr);}
    IrxTMarketConvSP(const IrxTMarketConvSP &src) {copy(src.p);}
    const IrxTMarketConvSP& operator= (const IrxTMarketConvSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~IrxTMarketConvSP() {deallocate();}
    IrxTMarketConv* get() {return p;}
//  IrxTMarketConv* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    IrxTMarketConv *p;
    void deallocate() {irxMarketConvFree(p); p=0; }
    void copy(const IrxTMarketConv*ptr) {
        if (ptr) {
            p = irxMarketConvCopy((IrxTMarketConv*)ptr);
            if (!p) throw StringException ("Could not copy IrxTMarketConv");
        } else {
            p = 0;
        }
    }
};

class IrxTSwapZeroCurveSP {
public:
    IrxTSwapZeroCurveSP() : p(0) {}
    IrxTSwapZeroCurveSP(IrxTSwapZeroCurve* ptr) : p(ptr) {}
    IrxTSwapZeroCurveSP(const IrxTSwapZeroCurve* ptr) {copy(ptr);}
    IrxTSwapZeroCurveSP(const IrxTSwapZeroCurveSP &src) {copy(src.p);}
    const IrxTSwapZeroCurveSP& operator= (const IrxTSwapZeroCurveSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~IrxTSwapZeroCurveSP() {deallocate();}
    IrxTSwapZeroCurve* get() {return p;}
//  IrxTSwapZeroCurve* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    IrxTSwapZeroCurve *p;
    void deallocate() {irxSwapZeroCurveFree(p); p=0; }
    void copy(const IrxTSwapZeroCurve*ptr) {
        if (ptr) {
            p = irxSwapZeroCurveCopy((IrxTSwapZeroCurve*)ptr);
            if (!p) throw StringException ("Could not copy IrxTSwapZeroCurve");
        } else {
            p = 0;
        }
    }
};
#endif /* _IRX_IRXFLOW_HPP */
