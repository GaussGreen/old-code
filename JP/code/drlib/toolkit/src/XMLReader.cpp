//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XMLReader.cpp
//
//   Description : Reads XML streams
//
//   Author      : Andrew J Swain
//
//   Date        : 18 January 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/Object.hpp"
#include "edginc/XMLReader.hpp"
#include "edginc/Format.hpp"
#include "edginc/ModelException.hpp"
#include "edginc/Socket.hpp"
#include "edginc/Class.hpp"
#include "edginc/Null.hpp"
#include "edginc/Hashtable.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/XMLWriter.hpp"
#include "edginc/Array.hpp"

#include "xercesc/dom/DOMElement.hpp"
#include "xercesc/dom/DOMNodeList.hpp"
#include "xercesc/dom/DOMNamedNodeMap.hpp"
#include "xercesc/dom/DOMException.hpp"
#include "xercesc/util/XMLException.hpp"
#include "xercesc/util/PlatformUtils.hpp"
#include "xercesc/parsers/XercesDOMParser.hpp"
#include "xercesc/sax/ErrorHandler.hpp"
#include "xercesc/sax/SAXParseException.hpp"
#include "xercesc/framework/LocalFileInputSource.hpp"
#include "xercesc/framework/MemBufInputSource.hpp"
XERCES_CPP_NAMESPACE_USE
DRLIB_BEGIN_NAMESPACE

// defined just below
static string transcode(const XMLCh* s);
// for performance cache the XMLCh for common attributes. Not in Imp
// class as need to declare them for use in MyNode class
static XMLCh* XML_REF_TYPE_ATTRIBUTE = 0;
static XMLCh* XML_OBJECT_TYPE_ATTRIBUTE = 0;
static XMLCh* XML_ARRAY_LENGTH = 0;

// create a ModelException based on a XMLException
static ModelException xmlException(
    const string&       caller,
    const XMLException& e) {
    const XMLCh* err = e.getMessage();
    return ModelException(caller, transcode(err)); 
}

// convert a DOMString into a C++ string
static string transcode(const XMLCh* s) {
    static const string method("XMLReader::Node::transcode");
    try {
        char  *raw = XMLString::transcode(s);
        string wrapped(raw);
        delete[] raw;
        return wrapped;
    }
    catch (XMLException& e) {
        throw xmlException(method, e);
    }
}

// create an ModelException based on a DOMException
static ModelException domException(
    const string&       caller,
    const DOMException& e) {
    const XMLCh* err = e.msg;
    return ModelException(caller, transcode(err)); 
}

static CClassConstSP classLookupHack(const string& className){
    // for backward compatibility for dr wrappers until
    // 4.1 is released into prod. 
    string name(className);
    if (name == "ObjectArray"){
        name = "IObjectArray"; // Map ObjectArray to IObjectArray
    } else if (name == "CAssetWrapperArray"){
        name = "AssetWrapperArray"; // CAssetWrapperArray to AssetWrapperArray
    }
    return CClass::forName(name);
}

// error handler for parser
class XMLReader::Error : public XERCES_CPP_NAMESPACE::ErrorHandler {
public:
    Error() {}

    ~Error() {} 

    static string format(int i) {
        char buffer[50];
        sprintf(buffer, "%d", i);
        return string(buffer);
    }

    // -----------------------------------------------------------------------
    //  Implementation of the error handler interface
    // -----------------------------------------------------------------------
    void warning(const SAXParseException& toCatch) {
        // live with warnings
    }

    void error(const SAXParseException& toCatch) {
        // apparently should live with errors
        // only "fatal" errors are fatal
    }

    void fatalError(const SAXParseException& toCatch) {
        string sysId = transcode(toCatch.getSystemId());
        string msg   = transcode(toCatch.getMessage());
        string line  = format((int)toCatch.getLineNumber());
        string col   = format((int)toCatch.getColumnNumber());

        throw ModelException("xml parse error (" + msg + ") in " + sysId +
                             " at line " + line + " and column " + col);
    }

    void resetErrors() {}
};

class XMLReader::MyNode: public Node{
    DOMElement*   node;
    CClassConstSP clazz; // lazy look up
    friend class XMLReader;
public:
    ~MyNode(){} // apparently no need to delete node

    /** what's the name of this node ? */
    const string name() {
        static const string method("XMLReader::MyNode::name");
        try {
            return transcode(node->getNodeName());
        }
        catch (exception& e) {
            throw ModelException(e, method);
        } 
        catch (XMLException& e) {
            throw xmlException(method, e);
        }
    }
    
    /** what class does this node represent ? */
    CClassConstSP getClass() {
        static const string method("XMLReader::MyNode::getClass");
        if (!clazz){
            try {
                string typeKey;
                if (!attributeExists(XML_OBJECT_TYPE_ATTRIBUTE, typeKey)){
                    typeKey = 
                        attribute(CObject::OBJECT_TYPE_ATTRIBUTE); // will fail
                }
                clazz = classLookupHack(typeKey);
            }
            catch (exception& e) {
                throw ModelException(e, method);
            }
        }
        return clazz;
    }
    
    /** is this node a null object ? */
    bool isNull() {
        try {
            return CNull::TYPE->isAssignableFrom(getClass());
        }
        catch (exception& e) {
            throw ModelException(e, "XMLReader::Node::isNull");
        }
    }
    
    /** is this node an array ? */
    bool isArray() {
        return getClass()->isArray();
    }
    
    /** if so, how long is it ? */
    int arrayLength() {
        string length;
        if (!attributeExists(XML_ARRAY_LENGTH, length)){
            length = attribute(IArray::ARRAY_LENGTH); // will fail
        }
        return(atoi(length.c_str()));
    }
    
    //// does the specified attribute exist
    virtual bool attributeExists(const string& name,
                                 string&       theAttribute){
        static const string method("XMLReader::Node::attributeExists");
        XMLCh* attrib = 0;
        try {
            DOMNamedNodeMap* map = node->getAttributes();
            if (map){
                attrib = XMLString::transcode(name.c_str());
                DOMNode* lnode = map? map->getNamedItem(attrib): 0;
                delete[] attrib;
                attrib = 0;
                if (lnode){
                    theAttribute = (transcode(lnode->getNodeValue()));
                    return true;
                }
            }
            return false;
        } catch (XMLException& e) {
            delete[] attrib;
            throw xmlException(method, e);
        }
    }

    /** return named attribute */
    string attribute(const string& name) {
        string theAttribute;
        if (!attributeExists(name, theAttribute)){
            throw ModelException("XMLReader::Node::attribute",
                                 "no attribute named: " + name);
        }
        return theAttribute;
    }

    /** return value */
    string value() {
        static const string method("XMLReader::Node::value");
        try {
            DOMNode* child = node->getFirstChild();
            if (child && child->getNodeType() == DOMNode::ELEMENT_NODE){
                // shouldn't be here - should have called children()
                throw ModelException(method, "Node ("+name()+") contains an"
                                     " element node");
            }
            string s(child? transcode(child->getNodeValue()): "");
            return s;
        } catch (XMLException& e) {
            throw xmlException(method, e);
        }
    }

    MyNode(DOMElement* node): node(node), clazz(0) {}
    
    /** return child nodes */
    NodeListSP children() {
        static const string method("XMLReader::MyNode::children");
        try {
            DOMNodeList*        nl = node->getChildNodes();       
            unsigned int numChildren = nl->getLength();
            NodeListSP kids(new NodeList());
            kids->reserve(numChildren);
            for (unsigned int i = 0; i < numChildren; i++) {
                DOMNode*  child = nl->item(i);
                if (child->getNodeType() == DOMNode::ELEMENT_NODE) {
                    NodeSP node(new MyNode(static_cast<DOMElement*>(child)));
                    kids->push_back(node);
                }
            }
            return kids;
        }
        catch (exception& e) {
            throw ModelException(e, method);
        }    
    }

    /** Reads [first] cdata section from this node. Fails if none found */
    string cdata(){
        static const string method("XMLReader::MyNode::cdata");
        try{
            DOMNodeList*        nl = node->getChildNodes();       
            unsigned int numChildren = nl->getLength();
            for (unsigned int i = 0; i < numChildren; i++) {
                DOMNode*  child = nl->item(i);
                if (child->getNodeType() == DOMNode::CDATA_SECTION_NODE) {
                    string s = transcode(child->getNodeValue());
                    return s;
                }
            }
        } catch (XMLException& e) {
            throw xmlException(method, e);
        }
        throw ModelException(method," No cdata section found");
    }
private:
    //// does the specified attribute exist
    bool attributeExists(XMLCh*        attrib,
                         string&       theAttribute){
        static const string method("XMLReader::Node::attributeExists");
        try {
            DOMNamedNodeMap* map = node->getAttributes();
            if (map){
                DOMNode* lnode = map? map->getNamedItem(attrib): 0;
                if (lnode){
                    theAttribute = (transcode(lnode->getNodeValue()));
                    return true;
                }
            }
            return false;
        } catch (XMLException& e) {
            throw xmlException(method, e);
        }
    }
};

typedef refCountPtr<InputSource> InputSourceSP;
typedef refCountPtr<XercesDOMParser> ParserSP;
class XMLReader::Imp{
public:
    static bool   activated;
    ParserSP      parser;
    ErrorSP       handler;
    InputSourceSP source;
    Hashtable     cache;
    string        urlBuffer;

    /** call once before doing anything */
    static void initialize(){
        static const string method("XMLReader::initialize");
        if (!activated){
            try {
                XMLPlatformUtils::Initialize();
                XML_REF_TYPE_ATTRIBUTE =
                    XMLString::transcode(XMLWriter::REF_TYPE_ATTRIBUTE.c_str());
                XML_OBJECT_TYPE_ATTRIBUTE = XMLString::transcode(
                    CObject::OBJECT_TYPE_ATTRIBUTE.c_str());
                XML_ARRAY_LENGTH = XMLString::transcode(
                    IArray::ARRAY_LENGTH.c_str());
                activated = true;
            } catch (XMLException& e) {
                throw xmlException(method, e);
            }
        }
    }
    static void terminate(){
        static const string method("XMLReader::terminate");
        try {
            XMLPlatformUtils::Terminate();
            delete[] XML_REF_TYPE_ATTRIBUTE;
            XML_REF_TYPE_ATTRIBUTE = 0;
            delete[] XML_OBJECT_TYPE_ATTRIBUTE;
            XML_OBJECT_TYPE_ATTRIBUTE = 0;
            delete[] XML_ARRAY_LENGTH;
            XML_ARRAY_LENGTH = 0;
            activated = false;
        } catch (XMLException& e) {
            throw xmlException(method, e);
        }
    }
        

    ~Imp(){}

    Imp(const string& xmlIn, bool isFile){
        static const char* id = "xml input buffer";
        XMLCh* tmpBuf = 0;
        try{
            if (CString::equalsIgnoreCase(xmlIn, "http://", 7)) {
                // magic hack for reading from URL
                urlBuffer = Socket::urlRead(xmlIn);
                source = InputSourceSP(
                    new MemBufInputSource((const XMLByte*)urlBuffer.c_str(),
                                          urlBuffer.length(),
                                          id));
            }
            else if (isFile) {
                tmpBuf = XMLString::transcode(xmlIn.c_str());
                source = InputSourceSP(new LocalFileInputSource(tmpBuf));
            }
            else {
                source = InputSourceSP(
                    new MemBufInputSource((const XMLByte*)xmlIn.c_str(),
                                          xmlIn.length(),
                                          id));
            }
        } catch (XMLException& e){
            delete[] tmpBuf;
            throw xmlException(id, e);
        }
        delete[] tmpBuf;
    }
    /** return the top level node in the document */
    Node* root(){
        static const string method("XMLReader::root");
        try {
            if (!activated) {
                throw ModelException(method, "xml parser not activated");
            }
            
            parser = ParserSP(new XercesDOMParser);
            
            handler = ErrorSP(new Error);
            parser->setErrorHandler(handler.get());
            
            // parse the XML into a DOM document
            parser->parse(*source);
            DOMDocument* doc = parser->getDocument();
            
            // get the root of the document
            DOMElement* rootNode = doc->getDocumentElement();      
            
            if (!rootNode) {
                throw ModelException(method, "failed to parse xml");
            }
            
            return new MyNode(rootNode);
        }
        catch (XMLException& e) {
            throw xmlException(method, e);
        }   
        catch (DOMException& e) {
            throw domException(method, e);
        }   
        catch (exception& e) {
            throw ModelException(e, method);
        }
    }
};

bool XMLReader::Imp::activated = false;

/** build an object */
IObjectSP XMLReader::read(){
    static const string method("XMLReader::read");
    try {
        // get the root of the document
        NodeSP rootNode(root());      
        
        // assemble the object
        return (read(rootNode.get()));
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

smartPtr<IObject> XMLReader::read(Node* node){
    return read(node, this);
}

/** build an object from a given node */
IObjectSP XMLReader::read(Node* node, Reader* subNodeReader){
    MyNode* myNode = static_cast<MyNode*>(node);
    bool isType = false;
    string ref;
    // is there a ref attribute?
    bool isRef = myNode->attributeExists(XML_REF_TYPE_ATTRIBUTE, ref);
    if (isRef){
        // reference exists - is there a TYPE attribute
        string typeAttrib;
        isType = myNode->attributeExists(XML_OBJECT_TYPE_ATTRIBUTE, typeAttrib);
        if (isType){
            // cache this whilst we're here
            static_cast<MyNode*>(node)->clazz = classLookupHack(typeAttrib);
        }
    }
    IObjectSP obj;
    if (!isRef || isType){
        // if no ref or type attribute exists call CObject::xmlRead
        obj = IObjectSP(CObject::read(node, subNodeReader));
        // if ref attribute, store object in hash table
        if (isRef){
            my->cache.put(ref, obj);
        }
    } else {
        // else get value from hash table, return
        obj = my->cache.get(ref);
        if (!obj){
            throw ModelException("XMLInputStream::read",
                                 "Node refers to unknown REF \""+ref+"\"");
        }
    }
    return obj;
}

/** return the top level node in the document */
Reader::Node* XMLReader::root(){
    return my->root();
}

/** call once before doing anything */
void XMLReader::initialize(){
    Imp::initialize();
}

/** Reads [first] cdata section from this node. Fails if none found */
string XMLReader::cdata(Node* node){
    // it had better be our node
    MyNode& myNode = dynamic_cast<MyNode&>(*node);
    return myNode.cdata();
}

XMLReader::~XMLReader(){
    delete my;
}

/** call once when all done */
void XMLReader::terminate(){
    Imp::terminate();
}

/** constructor */
XMLReader::XMLReader(const string& xmlIn, bool isFile): 
    my(new Imp(xmlIn, isFile)){}

DRLIB_END_NAMESPACE
