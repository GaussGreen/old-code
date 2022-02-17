/**
 * @file NotApplicableException.hpp
 */

#ifndef QLIB_NoSubjectsFoundException_H
#define QLIB_NoSubjectsFoundException_H

#include "edginc/ModelException.hpp"

DRLIB_BEGIN_NAMESPACE

/**
 * ModelException thrown by SensMgr::theFirst()
 */

class RISKMGR_DLL NoSubjectsFoundException: public ModelException {

public:

    NoSubjectsFoundException(const string& message);
    NoSubjectsFoundException();
};

DRLIB_END_NAMESPACE

// Local Variables: ***
// c-basic-offset:4 ***
// End: ***

#endif
