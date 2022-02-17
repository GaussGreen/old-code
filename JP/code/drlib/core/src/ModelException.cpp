//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ModelException.cpp
//
//   Description : Prototype exception
//
//   Author      : Andrew J Swain
//
//   Date        : 8 November 2000
//
//
//----------------------------------------------------------------------------

#include <iostream>

#include "edginc/coreConfig.hpp"
#include "edginc/ErrorHandler.hpp"
#include "edginc/Malloc.hpp"
#include "edginc/ModelException.hpp"

#if defined(_MSC_VER)
#define snprintf _snprintf
#endif

CORE_BEGIN_NAMESPACE

class CORE_DLL Error
{
public:
    Error(const string& method, const string& message) :
        method(method), message(message) {}

    Error(){}

    string method;
    string message;


};

class CORE_DLL ModelException_i
{
public:
    // fields
    ModelException* cause;   // originating exception, possibly null
    vector<string>  stack;    // the full error stack
    vector<Error>   logger;   // like the stack , but without all the Failed's

    ModelException_i(): cause(0){}

    ~ModelException_i(){
        delete cause;
    }
};

ModelException::ModelException() : mine(new ModelException_i)
{
    // empty
}

ModelException::ModelException(const char* msg) : mine(new ModelException_i)
{
    addMsg(msg);
}

ModelException::ModelException(const string& msg) : mine(new ModelException_i)
{
    addMsg(msg);
}

ModelException::ModelException(const string& routine, const string& msg):
    mine(new ModelException_i){
    addMsg(routine + ": " + msg);
    addMsg(routine + ": Failed");
    addLog(routine, msg);
}

// used to wrap assert - strips path from file name
ModelException ModelException::fromAssert(const char* file, int lineNumber, const char* msg) {
    // want to strip out any path from file name handling both fwd & back slash
    const char* stripfile = strrchr(file, '/'); // look for / first
    if (stripfile) {
        // got it - move along to file name
        stripfile++;
    }
    else {
        // look for \ now
        stripfile = strrchr(file, '\\');
        if (stripfile) {
            // got it - move along to file name
            stripfile++;
        }
        else {
            stripfile = file;
        }
    }

    //512 bytes is plenty of space for a filename and line number, 
    //but we use snprintf to be 100% safe.
    char fileAndLine[512];
    snprintf(fileAndLine, 511, "%s(%i)", stripfile, lineNumber);
    fileAndLine[511] = 0;
    return ModelException(fileAndLine, msg);
}

/* [possibly] store cause */
void ModelException::setCause(const ModelException* e){
    if (e->isDerived()){
        // if we have an exception derived from ModelException then store
        // a copy of it (would rather not copy but not much choice)
        ModelException* copy = e->clone();
        delete mine->cause;
        mine->cause = copy;
        // and clear the stack etc
        mine->stack.resize(0);
        mine->logger.resize(0);
    }
}

//ModelException::ModelException(const ModelException& e) : mine(e.mine) {
ModelException::ModelException(const ModelException& e) :
    mine(new ModelException_i) {
    *this = e;
    setCause(&e);
}

ModelException::ModelException(const exception& e) :
    mine(new ModelException_i) {
    const ModelException* existing = dynamic_cast<const ModelException*>(&e);
    if (existing) {
        *this = *existing;
        setCause(existing);
    }
    else {
        addMsg(e.what());
    }
}

ModelException::ModelException(exception& e, const string &routine) :
    mine(new ModelException_i) {
    ModelException* existing = dynamic_cast<ModelException*>(&e);
    if (existing) {
        *this = *existing;
        setCause(existing);
    }
    else {
        addMsg(e.what());
    }
    addMsg(routine+": Failed");
}

ModelException::ModelException(exception* e, const string &routine) :
    mine(new ModelException_i) {
    ModelException* existing = dynamic_cast<ModelException*>(e);
    if (existing) {
        *this = *existing;
        setCause(existing);
    }
    else {
        addMsg(e->what());
    }
    addMsg(routine+": Failed");
}

ModelException ModelException::addTextToException(exception& e, const string &message) {
    ModelException newException(e);
    newException.addMsg(message);
    return newException;
}

// adds method: msg and method: Failed to the call stack
ModelException::ModelException(exception&    e,
                               const string& method,
                               const string &msg) :
    mine(new ModelException_i) {
    ModelException* existing = dynamic_cast<ModelException*>(&e);
    if (existing) {
        *this = *existing;
        setCause(existing);
    }
    else {
        addMsg(e.what());
    }
    addMsg(method + ": " + msg);
    addMsg(method + ": Failed");
    addLog(method, msg);
}

ModelException::ModelException(exception* e,
                               const string& method,
                               const string &msg) :
    mine(new ModelException_i) {
    ModelException* existing = dynamic_cast<ModelException*>(e);
    if (existing) {
        *this = *existing;
        setCause(existing);
    }
    else {
        addMsg(e->what());
    }
    addMsg(method + ": " + msg);
    addMsg(method + ": Failed");
    addLog(method, msg);
}

/** Returns true if the ModelException is empty. This is the case if
    it has been constructed using the default constructor and
    addMsg() has not been invoked */
bool ModelException::empty() const{
    return mine->stack.empty();
}

ModelException & ModelException::operator= (const ModelException& rhs) {
    mine->stack  = rhs.mine->stack;
    mine->logger = rhs.mine->logger;
    ModelException* copy = rhs.mine->cause? rhs.mine->cause->clone(): 0;
    delete mine->cause;
    mine->cause = copy;
    return (*this);
}

ModelException::~ModelException() throw () {
    delete mine;
}

const char * ModelException::what() const throw () {
    ModelException_i* relMine = mine->cause? mine->cause->mine: mine;
    return relMine->stack.size() ? relMine->stack[0].c_str() : "";
}

/** sets the first message in the stack to the supplied string */
void ModelException::setDescription(const string& description){
    if (mine->stack.empty()){
        mine->stack.push_back(description);
    } else {
        mine->stack[0] = description;
    }
}


// echo stack trace to stdout
void ModelException::printStackTrace() const {
    if (mine->cause){
        mine->cause->printStackTrace();
    }
    for (unsigned int i = 0; i < mine->stack.size(); ++i) {
        cout << mine->stack[i] << endl;
    }
}

// echo error trace to stdout
void ModelException::printErrorTrace() const {
    if (mine->cause){
        mine->cause->printErrorTrace();
    }
    for (unsigned int i = 0; i < mine->logger.size(); ++i) {
        cout << mine->logger[i].method << ": " << mine->logger[i].message << endl;
    }
}


// write stack trace to the error log
void ModelException::errorLog() const {
    if (mine->cause){
        mine->cause->errorLog();
    }
    // maybe ErrorHandler should be a parameter ?
    for (unsigned int i = 0; i < mine->stack.size(); ++i) {
        ErrorHandler::writeMsg(mine->stack[i]);
    }
    ErrorHandler::flush();
}

// add a message to the stack trace
void ModelException::addMsg(const string& msg) {
    mine->stack.push_back(msg);
}

// add a message to the log
void ModelException::addLog(const string& method, const string& message) {
    mine->logger.push_back(Error(method, message));
}


/** return stack trace as a dynamic C style string -
    use free() to cleanup */
char* ModelException::stackTrace() const {
    unsigned int size;
    unsigned int i;
    char*        message = 0;

    char*        causeMessage = mine->cause? mine->cause->stackTrace(): 0;
    const char*  causeToUse = causeMessage? causeMessage: "";
    for (i = 0, size = 0; i < mine->stack.size(); ++i) {
        size += mine->stack[i].length();
    }

    // room for N "\n" and 1 null terminator
    size += strlen(causeToUse) + mine->stack.size() + 1;

    message = NEW_ARRAY(char, size);

    strcpy(message, causeToUse);
    FREE(causeMessage);
    for (i = 0; i < mine->stack.size(); ++i) {
        strcat(message, mine->stack[i].c_str());
        strcat(message, "\n");
    }

    return message;
}

/** creates a deep copy */
ModelException* ModelException::clone() const{
    return new ModelException(*this);
}

/** like clone but does not copy the cause exception over */
ModelException ModelException::create(const ModelException& e){
    ModelException newEx;
    *newEx.mine = *e.mine;
    newEx.mine->cause = 0;
    return newEx;
}


/** indicates whether this exception is derived from ModelException */
bool ModelException::isDerived() const{
    return false;
}

/** returns the original exception, if any */
ModelException* ModelException::getCause(){
    return mine->cause;
}

CORE_END_NAMESPACE
