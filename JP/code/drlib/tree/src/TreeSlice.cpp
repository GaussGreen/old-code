//----------------------------------------------------------------------------
//
//   Group       : Derivatives Research
//
//   Filename    : TreeSliceBasic.cpp
//
//   Description : Basic tree slice implementation.
//
//   Date        : Sep 8, 2005
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/TreeSlice.hpp"

const char* printTreeSliceP(const NAMESPACE::TreeSlice *p) {
    static char buf[10000];
    buf[0]=0;
    if (p)
        p->printDetails(buf);
    else
        strcat(buf, "NULL");
    return buf;
}

const char* printTreeSliceRC(const NAMESPACE::TreeSlice &p) {
    static char buf[10000];
    buf[0]=0;
    p.printDetails(buf);
    return buf;
}

const char* printTreeSliceR(NAMESPACE::TreeSlice &p) {
    return printTreeSliceRC(p);
}

DRLIB_BEGIN_NAMESPACE

/********************* ExposedSliceOverride *********************/

void TreeSlice::ExposedSliceOverride::clear() { 
    if( slice ) {
        slice->exposedSlice = oldValue; 
    }
    slice = 0;
}

void TreeSlice::ExposedSliceOverride::set( const TreeSlice * source, const TreeSlice * target ) {
    // don't override same slice twice
    if( target != source->exposedSlice ) {
        clear();
        slice = source;
        oldValue = source->exposedSlice;
        source->exposedSlice = target;
    }
}


/********************* TreeSlice *********************/

TreeSlice::TreeSlice() : treeStep(-1), iter(0) 
{ 
    exposedSlice = this; 
}

TreeSlice::TreeSlice(const TreeSlice &s) : treeStep(-1), iter(0)
{ 
    exposedSlice = this;
    operator=(s);
}

void TreeSlice::printDebug(char *s) const {
    if (iter) {
        char buf[30];
        sprintf(buf, "%f",*iter);
        strcat(s,buf);
    } else {
        strcat(s,"NULL");
    }
    strcat(s,"{");
    {
        static const string::size_type maxLen=500;
        char buf[maxLen+1];
        int len = (name.size() < maxLen ? name.size() : maxLen);
        strncpy(buf, name.data(), len);
        buf[len]=0;
        strcat(s, buf);
    }
    strcat(s,"}");
}

void TreeSlice::copyTreeSlice(const TreeSlice &s) {
    name = s.name;
    treeStep = s.treeStep;
}

TreeSliceSP TreeSlice::calcSmoothStep() const {
    throw ModelException("TreeSlice::calcSmoothStep", "Not implemented");
}

double TreeSlice::getCentre() const {
    throw ModelException("TreeSlice::getCentre", "Not implemented");
}

void TreeSlice::resizeVector(
    int newSize,
    const TreeSlice & proto,
    vector< TreeSliceSP > & slices,
    const string & namePrefix,
    bool copyValues )
{
    int oldSize = slices.size();
    if( oldSize == newSize )
        return;

    slices.resize( newSize );
    for( int i = oldSize; i < newSize; ++i )
    {
        slices[ i ] = proto.clone( copyValues );
#ifdef DEBUG
        slices[ i ]->name = namePrefix + "[" + Format::toString( i ) + "]";
#endif
    }
}

void TreeSlice::printDetails(char *s) const {
    static char buf[30];
    strcat(s, typeName().c_str());
    strcat(s, " \"");
    strcat(s, name.c_str());
    strcat(s, "\" step=");
    sprintf(buf, "%d", treeStep);
    strcat(s, buf);
}


DRLIB_END_NAMESPACE
