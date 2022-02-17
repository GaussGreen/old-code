//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : FRParseException.hpp
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

#ifndef EDR_FRPARSEEXCEPTION_HPP
#define EDR_FRPARSEEXCEPTION_HPP
#include "edginc/ModelException.hpp"

DRLIB_BEGIN_NAMESPACE

/** Defines exception thrown when a parse area is encountered
    when parsing 'flex rules' */ 
class PRODUCTS_DLL FRParseException: public ModelException{
public:
    static const int MAX_NUM_ERRORS; // = 30
    //// copy constructor
    FRParseException(const FRParseException& e);

    //// constructor: string is used to build parent
    FRParseException(int numErrors, const string& what);

    //// constructor: string is used to build parent
    FRParseException(const string& what);

    /** constructor: string is formed from what + "Failed parsing: "+
        expression + "< pos space> ^" */
    static FRParseException make(
        const string& what, const string& expression, int pos);

    /** constructor: string is formed from "Parse error, unexpected "+
        what + "Failed parsing: "+ expression + "< pos space> ^" */
    static FRParseException make(
        char unexpected, const string& expression, int pos);

    /** creates a [deep] copy of the exception */
    virtual ModelException* clone() const;
        
    /** indicates whether this exception is derived from ModelException -
        used to drive whether this is stored as a 'cause' when a new
        exception is created */
    virtual bool isDerived() const;

    /** Returns the number of errors that this exception holds */
    int getNumErrors() const;

private:
    int numErrors;
};

DRLIB_END_NAMESPACE
#endif
