//----------------------------------------------------------------------------
//                                                                           
// Group       : CH Quantitative Research                                    
//                                                                           
// Description : Interface for objects that are to be used for pricing 
//               ICreditLossConfig objects. The precise methods are engine
//               dependent (eg a 'Effective Loss Curve' model would have 
//               different methods to a monte carlo), so there is nothing
//               in this interface
//                                                                           
// Date        : July 2006                                                   
//                                                                           
//----------------------------------------------------------------------------

#ifndef QLIB_ICREDITLOSSGEN_HPP
#define QLIB_ICREDITLOSSGEN_HPP

#include "edginc/Object.hpp"
#include "edginc/VirtualDestructorBase.hpp"
#include "edginc/DECLARE.hpp"


DRLIB_BEGIN_NAMESPACE


/** Interface for objects that are to be used for pricing ICreditLossConfig
    objects. The precise methods are engine dependent (eg a 'Effective Loss
    Curve' model would have different methods to a monte carlo), so there
    is nothing in this interface */
class MARKET_DLL ICreditLossGen : virtual public VirtualDestructorBase {

public:
    virtual ~ICreditLossGen();

protected:
    ICreditLossGen();

private:    
    ICreditLossGen(const ICreditLossGen& rhs); // don't use
    ICreditLossGen& operator=(const ICreditLossGen& rhs); // don't use
};

DECLARE_REF_COUNT(ICreditLossGen);

DRLIB_END_NAMESPACE

#endif
