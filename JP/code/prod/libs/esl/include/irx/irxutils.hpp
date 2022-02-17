/*
***************************************************************************
** HEADER FILE: irxutils.hpp
**
** Provides mini C++ classes wrapping up the C data structures.
**
** Defines data structures and enums used in the irxutils library.
***************************************************************************
*/

#ifndef _IRX_IRXUTILS_HPP
#define _IRX_IRXUTILS_HPP

#include "irxutils.h"
#include "matrix2d.h"
#include "macros.h"
#include "calendar.h"
#include "dateutils.h"
#include "mdmin.h"
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

class IrxTMultiDimMinStateSP {
public:
    IrxTMultiDimMinStateSP() : p(0) {}
    IrxTMultiDimMinStateSP(IrxTMultiDimMinState* ptr) : p(ptr) {}
    IrxTMultiDimMinStateSP(const IrxTMultiDimMinState* ptr) {copy(ptr);}
    IrxTMultiDimMinStateSP(const IrxTMultiDimMinStateSP &src) {copy(src.p);}
    const IrxTMultiDimMinStateSP& operator= (const IrxTMultiDimMinStateSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~IrxTMultiDimMinStateSP() {deallocate();}
    IrxTMultiDimMinState* get() {return p;}
//  IrxTMultiDimMinState* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    IrxTMultiDimMinState *p;
    void deallocate() {irxMultiDimMinStateFree(p); p=0; }
    void copy(const IrxTMultiDimMinState*ptr) {
        if (ptr) {
            p = irxMultiDimMinStateCopy((IrxTMultiDimMinState*)ptr);
            if (!p) throw StringException ("Could not copy IrxTMultiDimMinState");
        } else {
            p = 0;
        }
    }
};

class IrxTCalendarSP {
public:
    IrxTCalendarSP() : p(0) {}
    IrxTCalendarSP(IrxTCalendar* ptr) : p(ptr) {}
    IrxTCalendarSP(const IrxTCalendar* ptr) {copy(ptr);}
    IrxTCalendarSP(const IrxTCalendarSP &src) {copy(src.p);}
    const IrxTCalendarSP& operator= (const IrxTCalendarSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~IrxTCalendarSP() {deallocate();}
    IrxTCalendar* get() {return p;}
//  IrxTCalendar* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    IrxTCalendar *p;
    void deallocate() {irxCalendarFree(p); p=0; }
    void copy(const IrxTCalendar*ptr) {
        if (ptr) {
            p = irxCalendarCopy((IrxTCalendar*)ptr);
            if (!p) throw StringException ("Could not copy IrxTCalendar");
        } else {
            p = 0;
        }
    }
};

class IrxTMatrix2DSP {
public:
    IrxTMatrix2DSP() : p(0) {}
    IrxTMatrix2DSP(IrxTMatrix2D* ptr) : p(ptr) {}
    IrxTMatrix2DSP(const IrxTMatrix2D* ptr) {copy(ptr);}
    IrxTMatrix2DSP(const IrxTMatrix2DSP &src) {copy(src.p);}
    const IrxTMatrix2DSP& operator= (const IrxTMatrix2DSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~IrxTMatrix2DSP() {deallocate();}
    IrxTMatrix2D* get() {return p;}
//  IrxTMatrix2D* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    IrxTMatrix2D *p;
    void deallocate() {irxMatrixFree(p); p=0; }
    void copy(const IrxTMatrix2D*ptr) {
        if (ptr) {
            p = irxMatrixCopy((IrxTMatrix2D*)ptr);
            if (!p) throw StringException ("Could not copy IrxTMatrix2D");
        } else {
            p = 0;
        }
    }
};

class IrxTDateListSP {
public:
    IrxTDateListSP() : p(0) {}
    IrxTDateListSP(IrxTDateList* ptr) : p(ptr) {}
    IrxTDateListSP(const IrxTDateList* ptr) {copy(ptr);}
    IrxTDateListSP(const IrxTDateListSP &src) {copy(src.p);}
    const IrxTDateListSP& operator= (const IrxTDateListSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~IrxTDateListSP() {deallocate();}
    IrxTDateList* get() {return p;}
//  IrxTDateList* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    IrxTDateList *p;
    void deallocate() {irxDateListFree(p); p=0; }
    void copy(const IrxTDateList*ptr) {
        if (ptr) {
            p = irxDateListCopy((IrxTDateList*)ptr);
            if (!p) throw StringException ("Could not copy IrxTDateList");
        } else {
            p = 0;
        }
    }
};
#endif /* _IRX_IRXUTILS_HPP */
