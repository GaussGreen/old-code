/*
***************************************************************************
** HEADER FILE: alib.hpp
**
** Provides mini C++ classes wrapping up the C data structures.
**
** Defines ALIB structures that we re-use within CX.
***************************************************************************
*/

#ifndef _CX_ALIB_HPP
#define _CX_ALIB_HPP


#include <alib/bastypes.h>
#include <alib/bondcnst.h>
#include <alib/bonds.h>
#include <alib/cdate.h>
#include <alib/convert.h>
#include <alib/datatab.h>
#include <alib/dtivlo.h>
#include <alib/dtlist.h>
#include <alib/gtomat.h>
#include <alib/rtypeo.h>
#include <alib/svolfx.h>
#include <alib/tcurve.h>
#include <alib/voldef.h>
#include <alib/zcurve.h>
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

class TZeroCurveSP {
public:
    TZeroCurveSP() : p(0) {}
    TZeroCurveSP(TZeroCurve* ptr) : p(ptr) {}
    TZeroCurveSP(const TZeroCurve* ptr) {copy(ptr);}
    TZeroCurveSP(const TZeroCurveSP &src) {copy(src.p);}
    const TZeroCurveSP& operator= (const TZeroCurveSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TZeroCurveSP() {deallocate();}
    TZeroCurve* get() {return p;}
//  TZeroCurve* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TZeroCurve *p;
    void deallocate() {GtoZeroCurveDelete(p); p=0; }
    void copy(const TZeroCurve*ptr) {
        if (ptr) {
            p = GtoZeroCurveMakeCopy((TZeroCurve*)ptr);
            if (!p) throw StringException ("Could not copy TZeroCurve");
        } else {
            p = 0;
        }
    }
};

class TFXSVolCurveSP {
public:
    TFXSVolCurveSP() : p(0) {}
    TFXSVolCurveSP(TFXSVolCurve* ptr) : p(ptr) {}
    TFXSVolCurveSP(const TFXSVolCurve* ptr) {copy(ptr);}
    TFXSVolCurveSP(const TFXSVolCurveSP &src) {copy(src.p);}
    const TFXSVolCurveSP& operator= (const TFXSVolCurveSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TFXSVolCurveSP() {deallocate();}
    TFXSVolCurve* get() {return p;}
//  TFXSVolCurve* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TFXSVolCurve *p;
    void deallocate() {GtoFXSVolCurveDelete(p); p=0; }
    void copy(const TFXSVolCurve*ptr) {
        if (ptr) {
            p = GtoFXSVolCurveCopy((TFXSVolCurve*)ptr);
            if (!p) throw StringException ("Could not copy TFXSVolCurve");
        } else {
            p = 0;
        }
    }
};

class TDataTableSP {
public:
    TDataTableSP() : p(0) {}
    TDataTableSP(TDataTable* ptr) : p(ptr) {}
    TDataTableSP(const TDataTable* ptr) {copy(ptr);}
    TDataTableSP(const TDataTableSP &src) {copy(src.p);}
    const TDataTableSP& operator= (const TDataTableSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TDataTableSP() {deallocate();}
    TDataTable* get() {return p;}
//  TDataTable* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TDataTable *p;
    void deallocate() {GtoDataTableDelete(p); p=0; }
    void copy(const TDataTable*ptr) {
        if (ptr) {
            p = GtoDataTableCopy((TDataTable*)ptr);
            if (!p) throw StringException ("Could not copy TDataTable");
        } else {
            p = 0;
        }
    }
};

class TRateTypeSP {
public:
    TRateTypeSP() : p(0) {}
    TRateTypeSP(TRateType* ptr) : p(ptr) {}
    TRateTypeSP(const TRateType* ptr) {copy(ptr);}
    TRateTypeSP(const TRateTypeSP &src) {copy(src.p);}
    const TRateTypeSP& operator= (const TRateTypeSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TRateTypeSP() {deallocate();}
    TRateType* get() {return p;}
//  TRateType* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TRateType *p;
    void deallocate() {GtoRateTypeDelete(p); p=0; }
    void copy(const TRateType*ptr) {
        if (ptr) {
            p = GtoRateTypeMakeCopy((TRateType*)ptr);
            if (!p) throw StringException ("Could not copy TRateType");
        } else {
            p = 0;
        }
    }
};

class TDateListSP {
public:
    TDateListSP() : p(0) {}
    TDateListSP(TDateList* ptr) : p(ptr) {}
    TDateListSP(const TDateList* ptr) {copy(ptr);}
    TDateListSP(const TDateListSP &src) {copy(src.p);}
    const TDateListSP& operator= (const TDateListSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TDateListSP() {deallocate();}
    TDateList* get() {return p;}
//  TDateList* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TDateList *p;
    void deallocate() {GtoFreeDateList(p); p=0; }
    void copy(const TDateList*ptr) {
        if (ptr) {
            p = GtoCopyDateList((TDateList*)ptr);
            if (!p) throw StringException ("Could not copy TDateList");
        } else {
            p = 0;
        }
    }
};

class TMatrix2DSP {
public:
    TMatrix2DSP() : p(0) {}
    TMatrix2DSP(TMatrix2D* ptr) : p(ptr) {}
    TMatrix2DSP(const TMatrix2D* ptr) {copy(ptr);}
    TMatrix2DSP(const TMatrix2DSP &src) {copy(src.p);}
    const TMatrix2DSP& operator= (const TMatrix2DSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TMatrix2DSP() {deallocate();}
    TMatrix2D* get() {return p;}
//  TMatrix2D* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TMatrix2D *p;
    void deallocate() {GtoMatrixFree(p); p=0; }
    void copy(const TMatrix2D*ptr) {
        if (ptr) {
            p = GtoMatrixCopy((TMatrix2D*)ptr);
            if (!p) throw StringException ("Could not copy TMatrix2D");
        } else {
            p = 0;
        }
    }
};

class TDateIntervalSP {
public:
    TDateIntervalSP() : p(0) {}
    TDateIntervalSP(TDateInterval* ptr) : p(ptr) {}
    TDateIntervalSP(const TDateInterval* ptr) {copy(ptr);}
    TDateIntervalSP(const TDateIntervalSP &src) {copy(src.p);}
    const TDateIntervalSP& operator= (const TDateIntervalSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TDateIntervalSP() {deallocate();}
    TDateInterval* get() {return p;}
//  TDateInterval* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TDateInterval *p;
    void deallocate() {GtoDateIntervalDelete(p); p=0; }
    void copy(const TDateInterval*ptr) {
        if (ptr) {
            p = GtoDateIntervalMakeCopy((TDateInterval*)ptr);
            if (!p) throw StringException ("Could not copy TDateInterval");
        } else {
            p = 0;
        }
    }
};

class TBondSP {
public:
    TBondSP() : p(0) {}
    TBondSP(TBond* ptr) : p(ptr) {}
    TBondSP(const TBond* ptr) {copy(ptr);}
    TBondSP(const TBondSP &src) {copy(src.p);}
    const TBondSP& operator= (const TBondSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TBondSP() {deallocate();}
    TBond* get() {return p;}
//  TBond* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TBond *p;
    void deallocate() {GtoFreeTBond(p); p=0; }
    void copy(const TBond*ptr) {
        if (ptr) {
            p = GtoCopyTBond((TBond*)ptr);
            if (!p) throw StringException ("Could not copy TBond");
        } else {
            p = 0;
        }
    }
};

class TCurveSP {
public:
    TCurveSP() : p(0) {}
    TCurveSP(TCurve* ptr) : p(ptr) {}
    TCurveSP(const TCurve* ptr) {copy(ptr);}
    TCurveSP(const TCurveSP &src) {copy(src.p);}
    const TCurveSP& operator= (const TCurveSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TCurveSP() {deallocate();}
    TCurve* get() {return p;}
//  TCurve* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TCurve *p;
    void deallocate() {GtoFreeTCurve(p); p=0; }
    void copy(const TCurve*ptr) {
        if (ptr) {
            p = GtoCopyCurve((TCurve*)ptr);
            if (!p) throw StringException ("Could not copy TCurve");
        } else {
            p = 0;
        }
    }
};

class TCashFlowListSP {
public:
    TCashFlowListSP() : p(0) {}
    TCashFlowListSP(TCashFlowList* ptr) : p(ptr) {}
    TCashFlowListSP(const TCashFlowList* ptr) {copy(ptr);}
    TCashFlowListSP(const TCashFlowListSP &src) {copy(src.p);}
    const TCashFlowListSP& operator= (const TCashFlowListSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TCashFlowListSP() {deallocate();}
    TCashFlowList* get() {return p;}
//  TCashFlowList* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TCashFlowList *p;
    void deallocate() {GtoFreeCFL(p); p=0; }
    void copy(const TCashFlowList*ptr) {
        if (ptr) {
            p = GtoCopyCFL((TCashFlowList*)ptr);
            if (!p) throw StringException ("Could not copy TCashFlowList");
        } else {
            p = 0;
        }
    }
};

class TBondTradeDataSP {
public:
    TBondTradeDataSP() : p(0) {}
    TBondTradeDataSP(TBondTradeData* ptr) : p(ptr) {}
    TBondTradeDataSP(const TBondTradeData* ptr) {copy(ptr);}
    TBondTradeDataSP(const TBondTradeDataSP &src) {copy(src.p);}
    const TBondTradeDataSP& operator= (const TBondTradeDataSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TBondTradeDataSP() {deallocate();}
    TBondTradeData* get() {return p;}
//  TBondTradeData* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TBondTradeData *p;
    void deallocate() {GtoFreeTBondTradeData(p); p=0; }
    void copy(const TBondTradeData*ptr) {
        if (ptr) {
            p = GtoCopyTBondTradeData((TBondTradeData*)ptr);
            if (!p) throw StringException ("Could not copy TBondTradeData");
        } else {
            p = 0;
        }
    }
};

class TVolDefIRSP {
public:
    TVolDefIRSP() : p(0) {}
    TVolDefIRSP(TVolDefIR* ptr) : p(ptr) {}
    TVolDefIRSP(const TVolDefIR* ptr) {copy(ptr);}
    TVolDefIRSP(const TVolDefIRSP &src) {copy(src.p);}
    const TVolDefIRSP& operator= (const TVolDefIRSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~TVolDefIRSP() {deallocate();}
    TVolDefIR* get() {return p;}
//  TVolDefIR* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    TVolDefIR *p;
    void deallocate() {GtoVolDefIRFree(p); p=0; }
    void copy(const TVolDefIR*ptr) {
        if (ptr) {
            p = GtoVolDefIRCopy((TVolDefIR*)ptr);
            if (!p) throw StringException ("Could not copy TVolDefIR");
        } else {
            p = 0;
        }
    }
};
#endif /* _CX_ALIB_HPP */
