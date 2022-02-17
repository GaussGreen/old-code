//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRParseException.cpp
//
//   Description : Defines exception thrown when a parse area is encountered
//                 when parsing 'flex rules'
//
//   Author      : Mark A Robson
//
//   Date        : 30 July 2002
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/FRParseException.hpp"

DRLIB_BEGIN_NAMESPACE

//// copy constructor
FRParseException::FRParseException(const FRParseException& e): 
    ModelException(create(e)), /* avoid infinite recursion  - does not copy
                                  cause */
    numErrors(e.numErrors){} 

//// constructor: string is used to build parent
FRParseException::FRParseException(const string& what): 
    ModelException(what), numErrors(1){}

//// constructor: string is used to build parent
FRParseException::FRParseException(int numErrors, const string& what):
    ModelException(what), numErrors(numErrors){}

/** constructor: string is formed from what + "Failed parsing: "+
    expression + "< pos space> ^" */
FRParseException FRParseException::make(
    const string& what, const string& expression, int pos){
    string err1("Failed parsing: "+expression);
    string err2("                "+string(pos, ' ')+"^");
    return FRParseException(what+"\n"+err1+"\n"+err2);
}

/** constructor: string is formed from "Parse error, unexpected "+
    what + "Failed parsing: "+ expression + "< pos space> ^" */
FRParseException FRParseException::make(
    char unexpected, const string& expression, int pos){
    string what("Parse error, unexpected '"+string(1, unexpected)+"'");
    return make(what, expression, pos);
}


/** creates a [deep] copy of the exception */
ModelException* FRParseException::clone() const{
    return new FRParseException(*this);
}

/** indicates whether this exception is derived from ModelException -
    used to drive whether this is stored as a 'cause' when a new
    exception is created */
bool FRParseException::isDerived() const{
    return true;
}

/** Returns the number of errors that this exception holds */
int FRParseException::getNumErrors() const{
    return numErrors;
}

const int FRParseException::MAX_NUM_ERRORS= 30;

DRLIB_END_NAMESPACE
