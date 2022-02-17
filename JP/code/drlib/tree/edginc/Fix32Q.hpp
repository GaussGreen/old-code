//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : Fix32Q.hpp
//
//----------------------------------------------------------------------------

#ifndef _FIX32Q_HPP
#define _FIX32Q_HPP

#include "edginc/SingleCcyTree.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL Fix32Q : public SingleCcyTree {
public:

    /*************** FDModel interface ****************/
protected:
    virtual void retrieveFactor();

    /********************** local methods  **************************/
protected:
    Fix32Q(const CClassConstSP &type = TYPE) : SingleCcyTree(type) {}
private:
    static IObject* defaultConstructor(void) { return new Fix32Q(); }
    static void load(CClassSP& clazz); /* declare what is being exposed */

    /***************************** variables ********************************/
public:
    static CClassConstSP const TYPE;

};

typedef smartPtr<Fix32Q> Fix32QSP;

DRLIB_END_NAMESPACE

#endif

