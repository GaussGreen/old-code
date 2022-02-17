//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XMLReader.hpp
//
//   Description : Reads XML streams
//
//   Author      : Andrew J Swain
//
//   Date        : 18 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_XMLREADER_HPP
#define EDR_XMLREADER_HPP

#include "edginc/Reader.hpp"


DRLIB_BEGIN_NAMESPACE

/** How we read XML - key method is to read an object from the XML file */
class TOOLKIT_DLL XMLReader: public Reader{
public:
     ~XMLReader();

   /** build an object */
    virtual smartPtr<IObject> read();

    /** build an object from a given node */
    virtual smartPtr<IObject> read(Node* node, Reader* subNodeReader);

    /** build an object from a given node. This default implementation
        simply calls read(Node*, this). We have to call the default
        implementation ourselves becauase defining the above 'read' methods
        hides all 'read' methods (yes C++ is just dumb at times). In theory
        you can work around it using a 'using' but then not all compilers 
        support it ... */
    virtual smartPtr<IObject> read(Node* node);

    /** return the top level node in the document */
    virtual Node* root();

    /** call once before doing anything */
    static void initialize();

    /** call once when all done */
    static void terminate();

    /** Reads [first] cdata section from this node. Fails if none found */
    string cdata(Node* node);

    /** constructor */
    XMLReader(const string& xmlin, bool isFile);
    
private:
    XMLReader(XMLReader& rhs); // not wanted
    XMLReader& operator=(const XMLReader& rhs);// not wanted
    class Imp;
    // fields
    Imp*  my; // hides implementation
    // internal classes
    class MyNode;
    class Error;
    typedef refCountPtr<Error> ErrorSP;
};

DRLIB_END_NAMESPACE

#endif
