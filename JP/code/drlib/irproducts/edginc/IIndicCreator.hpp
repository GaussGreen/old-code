#ifndef _IIndicCreator_HPP
#define _IIndicCreator_HPP

#include "edginc/Model.hpp"
#include "edginc/DECLARE.hpp"

DRLIB_BEGIN_NAMESPACE

// generic indicator function. The product getValue must return values within range [0;1]
class IIndicCreator: virtual public IProdCreator {
public:
    static CClassConstSP const TYPE;

protected:
    IIndicCreator(){}
private:
    static void load(CClassSP& clazz);
};
DECLARE(IIndicCreator);

DRLIB_END_NAMESPACE

#endif
