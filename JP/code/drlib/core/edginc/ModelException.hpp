//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ModelException.hpp
//
//   Description : Prototype exception
//
//   Author      : Andrew J Swain
//
//   Date        : 8 November 2000
//
//
//----------------------------------------------------------------------------

#ifndef MODELEXCEPTION_HPP
#define MODELEXCEPTION_HPP

#include <exception>
#include <string>
#include "Format.hpp"

CORE_BEGIN_NAMESPACE

using namespace std;


class ModelException_i;

class CORE_DLL ModelException : public exception
{
public:
    ModelException();
    ModelException(const char* msg);
    ModelException(const string& msg);
    ModelException(const string& routine, const string& msg);
    ModelException(const ModelException& e);
    ModelException(const exception& e);

    /** used to wrap assert - strips path from file name & turns assert into
        an exception
    */
    static ModelException fromAssert(const char* file, int lineNumber, const char* msg);

    /** Adds routine: Failed to call stack */
    // deprecated - prefer one that takes reference.
    ModelException(exception* e, const string &routine);
    ModelException(exception& e, const string &routine);

    /** Creates a ModelException from an exception and adds a message to it */
    static ModelException addTextToException(exception& e, const string &message);

    /** adds method: msg and "method: Failed" to the call stack */
    // deprecated - prefer one that takes reference
    ModelException(exception* e, const string &method, const string &msg);
    ModelException(exception& e, const string &method, const string &msg);



    ModelException & operator= (const ModelException& rhs);
    virtual ~ModelException() throw ();

    virtual const char *what() const throw ();

    /** Returns true if the ModelException is empty. This is the case if
        it has been constructed using the default constructor and
        addMsg() has not been invoked */
    bool empty() const;

    // methods that make this useful
    /** echo stack trace to stdout */
    void printStackTrace() const;

    /** echo error trace to stdout - like printStackTrace but
        without all the Failed's */
    void printErrorTrace() const;

    /** write stack trace to the error log */
    void errorLog() const;

    /** add a message to the stack trace */
    void addMsg(const string& msg);

    /** return stack trace as a dynamic C style string -
        use free() to cleanup */
    char* stackTrace() const;

    /** returns the original exception, if any */
    virtual ModelException* getCause();

protected:
    /* these functions are to try to allow an infrastructure that supports
       different types of exceptions but also ensures that they don't get
       lost when people do "catch (exception& e){ throw ModelException(e);}"
    */

    /** creates a deep copy */
    virtual ModelException* clone() const;

    /** indicates whether this exception is derived from ModelException -
        used to drive whether this is stored as a 'cause' when a new
        exception is created */
    virtual bool isDerived() const;

    /** sets the first message in the stack to the supplied string */
    void setDescription(const string& description);

    /** like clone but does not copy the cause exception over */
    static ModelException create(const ModelException& e);

private:
    //auto_ptr<ModelException_i> mine;
    ModelException_i* mine;

    /** add a message to the log */
    void addLog(const string& method, const string& message);

    /* [possibly] store cause */
    void setCause(const ModelException* e);
};

/** Helper macro QLIB_VERIFY to be used as a replacement for 'assert' */
namespace QLIBVERIFY
{
//    inline void checkAndThrow(
//        bool condition, 
//        const string& method, 
//        const string& msg)
//    {
//        if (!condition)
//            throw ModelException(method, msg);
//    }
} // end of namespace 

#define checkAndThrow(condition, method, msg) if (condition) ; else throw ModelException((method), (msg))

#define QLIB_VERIFY(cond, msg) \
    /*QLIBVERIFY::*/ checkAndThrow( (cond), Format::toString("%s(%i)", (char *)__FUNCTION__, (int)__LINE__), (msg))


CORE_END_NAMESPACE

#endif

