/*
***************************************************************************
** HEADER FILE: mqdata.hpp
**
** Provides mini C++ classes wrapping up the C data structures.
**
** Data manipulation functions for the MQDATA structure.
***************************************************************************
*/

#ifndef _CRX_MQDATA_HPP
#define _CRX_MQDATA_HPP

#include "mqdata.h"

#include <crxflow/include/crxmacros.h>
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

class MQDATASP {
public:
    MQDATASP() : p(0) {}
    MQDATASP(MQDATA* ptr) : p(ptr) {}
    MQDATASP(const MQDATA* ptr) {copy(ptr);}
    MQDATASP(const MQDATASP &src) {copy(src.p);}
    const MQDATASP& operator= (const MQDATASP &src) {
        if (src.p != p) {
            deallocate(); }
        copy(src.p);
        return *this; }
    ~MQDATASP() {deallocate();}
    MQDATA* get() {return p;}
//  MQDATA* operator->() {return p;}
    bool isNull() const {return p==0;}
private:
    MQDATA *p;
    void deallocate() {CrxMqdataFree(p); p=0; }
    void copy(const MQDATA*ptr) {
        if (ptr) {
            p = CrxMqdataCopy((MQDATA*)ptr);
            if (!p) throw StringException ("Could not copy MQDATA");
        } else {
            p = 0;
        }
    }
};
#endif /* _CRX_MQDATA_HPP */
