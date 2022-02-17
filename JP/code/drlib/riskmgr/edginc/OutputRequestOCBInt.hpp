//
//

#ifndef EDG_OUTPUT_REQUEST_OCB_INT_HPP
#define EDG_OUTPUT_REQUEST_OCB_INT_HPP
#include "edginc/OutputName.hpp"
#include "edginc/Array.hpp"
#include "edginc/AtomicArray.hpp"
#include "edginc/OutputRequest.hpp"

DRLIB_BEGIN_NAMESPACE
class Results;

/** Identifies additional outputs which don't require tweaking - for example,
    indicative vol, fwd at maturity, effective strike etc. These will usually
    be added to the result set during the first pricing. */
class RISKMGR_DLL OutputRequestOCBInt: public OutputRequest{
public:
    friend class OutputRequestOCBIntHelper;
    static CClassConstSP const TYPE;

    //// creates deep copy
    virtual IObject* clone() const;

    //// ensures that the output request is a valid one
//    virtual void validatePop2Object();
    double getMTMPrice() {return mtmPrice;}
    double getQuoteFace() {return quoteFace;}
    bool getIsClean() {return isClean;}
private:
//    string               requestName;
    double               mtmPrice;
    double               quoteFace;
    bool                 isClean;
        
    // field not registered

    bool                 hasFinished;    // indicates whether request has been  $unregistered
                                         // calculated successfully
    /** default constructor for reflection */
    OutputRequestOCBInt();

    // not implemented
    OutputRequestOCBInt(const OutputRequestOCBInt &rhs);
    OutputRequestOCBInt& operator=(const OutputRequestOCBInt& rhs);
};

typedef smartPtr<OutputRequestOCBInt>  OutputRequestOCBIntSP;

DRLIB_END_NAMESPACE

#endif
