//----------------------------------------------------------------------------
//
//
//----------------------------------------------------------------------------

#ifndef EDR_READER_HPP
#define EDR_READER_HPP
#include "edginc/smartPtr.hpp"

DRLIB_BEGIN_NAMESPACE

class IObject;
class CClass;
typedef const CClass* CClassConstSP;
typedef CClass* CClassSP;

/** The Reader class defines an interface for reading objects from a stream */
class TOOLKIT_DLL Reader {
public:

    class Node;
    typedef refCountPtr<Node> NodeSP;
    typedef vector<NodeSP> NodeList;
    typedef refCountPtr<NodeList> NodeListSP;

    /** A Reader returns data as a 'document' of connected 'nodes' */
    class TOOLKIT_DLL Node {
    public:
        /** what's the name of this node ? */
        virtual const string name() = 0;

        /** what class does this node represent ? */
        virtual CClassConstSP getClass() = 0;

        /** is this node a null object ? */
        virtual bool isNull() = 0;

        /** is this node an array ? */
        virtual bool isArray() = 0;

        /** if so, how long is it ? */
        virtual int arrayLength() = 0;
        
        //// does the specified attribute exist. If so returns it
        virtual bool attributeExists(const string& name,
                                     string&       theAttribute) = 0;

        /** return named attribute. Fails if it does not exist */
        virtual string attribute(const string& name) = 0;

        /** return value */
        virtual string value() = 0;

        /** return child nodes */
        virtual NodeListSP children() = 0;

        virtual ~Node();
    protected:
        Node();
    };

    /** build an object */
    virtual smartPtr<IObject> read() = 0;

    /** build an object from a given node. This default implementation
        simply calls read(Node*, this) */
    virtual smartPtr<IObject> read(Node* node);

    /** build an object from a given node but using supplied reader for
        subNodes (this allows Readers to be wrapped since typically they
        pass themselves down to subNodes - by passing the wrapped reader
        in as the subNodeReader the use of the wrapped reader is preserved). */
    virtual smartPtr<IObject> read(Node* node, Reader* subNodeReader) = 0;

    /** return the top level node in the document */
    virtual Node* root() = 0;

    virtual ~Reader();
protected:
    Reader();
};

typedef refCountPtr<Reader> ReaderSP;

DRLIB_END_NAMESPACE
#endif
