//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : ResultsSet.hpp
//
//   Description : Return type of a composite instrument
//
//   Author      : Andrew J Swain
//
//   Date        : 24 May 2001
//
//
//----------------------------------------------------------------------------


#ifndef RESULTSSET_HPP
#define RESULTSSET_HPP
#include <string>
#include "edginc/Object.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/Results.hpp"

DRLIB_BEGIN_NAMESPACE
/** Return type of a composite instrument  */
class ADDINS_DLL ResultsSet: public CObject {
public:
    static CClassConstSP const TYPE;
   
    ResultsSet(Results*       composite,
               CResultsArray* components);


private:
    friend class ResultsSetHelper;
    friend class ResultsSetAddin;
    ResultsSet();
    ResultsSet(const ResultsSet &rhs);
    ResultsSet& operator=(const ResultsSet& rhs);

    ResultsSP       composite;
    CResultsArraySP components;
};

// typedef for smart pointers to ResultsSet
typedef smartConstPtr<ResultsSet> ResultsSetConstSP;
typedef smartPtr<ResultsSet> ResultsSetSP;

DRLIB_END_NAMESPACE
#endif
