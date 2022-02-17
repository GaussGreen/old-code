//----------------------------------------------------------------------------
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Writer.hpp"

DRLIB_BEGIN_NAMESPACE

// default constructor/destructor empty
Writer::Writer() {}

Writer::~Writer() {}

/** default implementation of writing out an object */
void Writer::write(const string& id, const IObject* object) {
    // route back to object
    object->write(id, this);
}

DRLIB_END_NAMESPACE
