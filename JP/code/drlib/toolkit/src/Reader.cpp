//----------------------------------------------------------------------------
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE
/** build an object from a given node. This default implementation
    simply calls read(Node*, this) */
smartPtr<IObject> Reader::read(Node* node){
    return read(node, this);
}

// default constructor/destructor empty
Reader::Reader() {}

Reader::~Reader() {}

Reader::Node::Node() {}

Reader::Node::~Node() {}

DRLIB_END_NAMESPACE
