//----------------------------------------------------------------------------
//
//   Group       : Global QR&D
//
//   Filename    : TMX.hpp
//
//----------------------------------------------------------------------------

#ifndef _TMX_HPP
#define _TMX_HPP

#include "edginc/SingleCcyTree.hpp"

DRLIB_BEGIN_NAMESPACE

class TREE_DLL TMX : public SingleCcyTree {
public:

    /*************** FDModel interface ****************/
protected:
    virtual void retrieveFactor();

    /********************** local methods  **************************/
protected:
    TMX(const CClassConstSP &type = TYPE) : SingleCcyTree(type) {}
private:
    static IObject* defaultConstructor(void) { return new TMX(); }
    static void load(CClassSP& clazz); /* declare what is being exposed */

    /***************************** variables ********************************/
public:
    static CClassConstSP const TYPE;

};

typedef smartPtr<TMX> TMXSP;

DRLIB_END_NAMESPACE

#endif

