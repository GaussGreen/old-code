#ifndef IMCINTOPRODUCT_HPP
#define IMCINTOPRODUCT_HPP
DRLIB_BEGIN_NAMESPACE

#include "edginc/Object.hpp"
#include "edginc/Model.hpp"

class MonteCarlo;
class IMCProduct; // note that clients will have to include "MCProduct.hpp"

/** interface that the instrument must implement */
class MCARLO_DLL IMCIntoProduct: virtual public CModel::IModelIntoProduct {
public:
    static CClassConstSP const TYPE;

    virtual ~IMCIntoProduct() {};

    /** Creates an instance of an IMCProduct */
    virtual IMCProduct* createProduct(const MonteCarlo* model) const = 0;
};

DRLIB_END_NAMESPACE

#endif
