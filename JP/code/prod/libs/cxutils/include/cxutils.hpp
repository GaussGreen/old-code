/*
***************************************************************************
** HEADER FILE: cxutils.hpp
**
** Provides mini C++ classes wrapping up the C data structures.
**
** Defines data structures used in the cxutils library.
***************************************************************************
*/

#ifndef _CX_CXUTILS_HPP
#define _CX_CXUTILS_HPP

#include "cxutils.h"
#include "cxmacros.h"
#include "calendar.h"
#include "zerocurve.h"
#include "dateutils.h"
#include "surface.h"
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

class CxTCalendarSP {
public:
    CxTCalendarSP() : p(0) {}
    CxTCalendarSP(CxTCalendar* ptr) : p(ptr) {}
    CxTCalendarSP(const CxTCalendar* ptr) {copy(ptr);}
    CxTCalendarSP(const CxTCalendarSP &src) {copy(src.p);}
    const CxTCalendarSP& operator= (const CxTCalendarSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTCalendarSP() {deallocate();}
    CxTCalendar* get() {return p;}
//  CxTCalendar* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTCalendar *p;
    void deallocate() {CxCalendarFree(p); p=0; }
    void copy(const CxTCalendar*ptr) {
        if (ptr) {
            p = CxCalendarCopy((CxTCalendar*)ptr);
            if (!p) throw StringException ("Could not copy CxTCalendar");
        } else {
            p = 0;
        }
    }
};

class CxTMultiDimMinStateSP {
public:
    CxTMultiDimMinStateSP() : p(0) {}
    CxTMultiDimMinStateSP(CxTMultiDimMinState* ptr) : p(ptr) {}
    CxTMultiDimMinStateSP(const CxTMultiDimMinState* ptr) {copy(ptr);}
    CxTMultiDimMinStateSP(const CxTMultiDimMinStateSP &src) {copy(src.p);}
    const CxTMultiDimMinStateSP& operator= (const CxTMultiDimMinStateSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTMultiDimMinStateSP() {deallocate();}
    CxTMultiDimMinState* get() {return p;}
//  CxTMultiDimMinState* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTMultiDimMinState *p;
    void deallocate() {CxMultiDimMinStateFree(p); p=0; }
    void copy(const CxTMultiDimMinState*ptr) {
        if (ptr) {
            p = CxMultiDimMinStateCopy((CxTMultiDimMinState*)ptr);
            if (!p) throw StringException ("Could not copy CxTMultiDimMinState");
        } else {
            p = 0;
        }
    }
};

class CxTSurfaceSP {
public:
    CxTSurfaceSP() : p(0) {}
    CxTSurfaceSP(CxTSurface* ptr) : p(ptr) {}
    CxTSurfaceSP(const CxTSurface* ptr) {copy(ptr);}
    CxTSurfaceSP(const CxTSurfaceSP &src) {copy(src.p);}
    const CxTSurfaceSP& operator= (const CxTSurfaceSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTSurfaceSP() {deallocate();}
    CxTSurface* get() {return p;}
//  CxTSurface* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTSurface *p;
    void deallocate() {CxSurfaceFree(p); p=0; }
    void copy(const CxTSurface*ptr) {
        if (ptr) {
            p = CxSurfaceCopy((CxTSurface*)ptr);
            if (!p) throw StringException ("Could not copy CxTSurface");
        } else {
            p = 0;
        }
    }
};

class CxTSurfaceInterpSP {
public:
    CxTSurfaceInterpSP() : p(0) {}
    CxTSurfaceInterpSP(CxTSurfaceInterp* ptr) : p(ptr) {}
    CxTSurfaceInterpSP(const CxTSurfaceInterp* ptr) {copy(ptr);}
    CxTSurfaceInterpSP(const CxTSurfaceInterpSP &src) {copy(src.p);}
    const CxTSurfaceInterpSP& operator= (const CxTSurfaceInterpSP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~CxTSurfaceInterpSP() {deallocate();}
    CxTSurfaceInterp* get() {return p;}
//  CxTSurfaceInterp* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    CxTSurfaceInterp *p;
    void deallocate() {CxSurfaceInterpFree(p); p=0; }
    void copy(const CxTSurfaceInterp*ptr) {
        if (ptr) {
            p = CxSurfaceInterpCopy((CxTSurfaceInterp*)ptr);
            if (!p) throw StringException ("Could not copy CxTSurfaceInterp");
        } else {
            p = 0;
        }
    }
};
#endif /* _CX_CXUTILS_HPP */
