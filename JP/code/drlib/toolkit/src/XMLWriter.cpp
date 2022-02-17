//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XMLWriter.cpp
//
//   Description : Writer for XML
//
//   Author      : Andrew J Swain
//
//   Date        : 17 January 2001
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include <sstream>
#include "edginc/XMLWriter.hpp"
#include "edginc/Class.hpp"
#include "edginc/Version.hpp"
#include "edginc/Null.hpp"
#include "edginc/Format.hpp"
#include "edginc/Hashtable.hpp"
#include <fstream>
#include ext_hash_map

DRLIB_BEGIN_NAMESPACE
const string XMLWriter::REF_TYPE_ATTRIBUTE  = "REF";
static const string MARKUP_AMPERSAND    = "amp;";
static const string MARKUP_LESS_THAN    = "lt;";
static const string MARKUP_GREATER_THAN = "gt;";
static const string MARKUP_QUOTE        = "quot;";
static const string MARKUP_APOS         = "apos;";

static void writeDisclaimer(ostream* stream){
    *stream << "<!-- WARNING: The format of this XML does NOT constitute a "
            << "public interface." << endl 
            << "Quantitative Research reserve the right "
            << "to change the format at will and do" << endl
            << "not support any applications that attempt to "
            << "parse, interpret or in any way" << endl
            << "manipulate this XML other than via publicly "
            << "supported QR interfaces." << endl
            <<"Anyone assuming the style, format or content "
            << "of this XML does so entirely at their own risk." 
            << " -->" << endl;
}

static void writeVersion(ostream* stream){
    writeDisclaimer(stream);
    *stream << "<!-- QLIB Version " << CVersion::DRLibVersion() << 
        " -->" << endl;
}

// replace any & with entity reference &amp;
static string purgeAmpersand(const string& s, std::string::size_type x = 0) {
    std::string::size_type i;

    // does it have any &'s ?
    i = s.find('&', x);
    if (i == std::string::npos) {
        // no &, nothing to do
        return s;
    }

    // is it markup ?
    if (!strncmp(s.c_str()+i+1, MARKUP_AMPERSAND.c_str(), 4) ||
        !strncmp(s.c_str()+i+1, MARKUP_LESS_THAN.c_str(), 3) ||
        !strncmp(s.c_str()+i+1, MARKUP_GREATER_THAN.c_str(), 3) ||
        !strncmp(s.c_str()+i+1, MARKUP_QUOTE.c_str(), 5) ||
        !strncmp(s.c_str()+i+1, MARKUP_APOS.c_str(), 5)) {
        // continue to do rest of string
        return purgeAmpersand(s, i+1);
    }

    // replace it with markup
    string p = s;

    p.insert(i+1, MARKUP_AMPERSAND);
    
    return purgeAmpersand(p, i+1);
}

// replace any < with entity reference &lt;
static string purgeLessThan(const string& s, std::string::size_type x = 0) {
    std::string::size_type i;

    // does it have any <'s ?
    i = s.find("<", x);
    if (i == std::string::npos) {
        // no <, nothing to do
        return s;
    }

    // replace it with markup
    string p = s;

    p[i] = '&';
    p.insert(i+1, MARKUP_LESS_THAN);
    
    return purgeLessThan(p, i+1);
}

// replace any > with entity reference &gt;
static string purgeGreaterThan(const string& s, std::string::size_type x = 0) {
    std::string::size_type i;

    // does it have any <'s ?
    i = s.find(">", x);
    if (i == std::string::npos) {
        // no <, nothing to do
        return s;
    }

    // replace it with markup
    string p = s;

    p[i] = '&';
    p.insert(i+1, MARKUP_GREATER_THAN);
    
    return purgeGreaterThan(p, i+1);
}

class XMLWriter::Imp{
    struct PtrHash {
        // hash func
        size_t operator()(const IObjectConstSP& obj) const {
            return (size_t)obj.get();}
        // equality func
        bool operator()(const IObjectConstSP& obj1,
                        const IObjectConstSP& obj2) const {
            return (obj1.get() == obj2.get());
        }
    };
    struct StringHash {
        size_t operator()(const string& str) const {
            return (hash_string(str.c_str()));
        }
    };
public:
    Imp(): isCDATA(false), allowBackReferences(true) {};

    typedef hash_map<IObjectConstSP, string, PtrHash, PtrHash> ObjectHashMap;
    typedef hash_map<string, int, StringHash> HandleIndexMap;
    ObjectHashMap  writtenObjects;
    HandleIndexMap handleIndex;
    bool           isCDATA;
    bool           allowBackReferences;

    /** Stores the reference name to use for the given object in refName.
        Returns true if the object has been written out before */
    bool getRefName(const string&   tagName, 
                    const IObject*  obj,
                    string&         refName) {
        /* this method is only called if refCount > 2 (so obj must already
           be associated with a SP */
        IObjectConstSP theObj(IObjectConstSP::attachToRef(obj));
        ObjectHashMap::const_iterator iter = writtenObjects.find(theObj);
        if (iter != writtenObjects.end()){
            refName = iter->second;
            return true;
        }
        // work out base handle name
        const string& clazz = theObj->getClass()->getName();
        string baseName = tagName+"."+clazz;
        int& index = handleIndex[baseName]; // create if not there
        index++; // move to next one
        baseName += "."+Format::toString(index);
        writtenObjects[theObj] = baseName; // and record
        refName = baseName;
        return false;
    }
};

XMLWriter::XMLWriter(ostream* stream) : 
    stream(stream), delStream(false), buffer(0), my(new Imp()) {
    writeVersion(stream);
}

XMLWriter::XMLWriter(string& xmlOut, bool isFile): 
    stream(0), delStream(true), buffer(isFile? 0: &xmlOut), my(new Imp())  {
    if (isFile){
        stream = new ofstream(xmlOut.c_str());
        writeVersion(stream);
    } else {
        this->buffer->erase();
    }
}

XMLWriter::XMLWriter(const string& fileName):
    stream(0), delStream(true), buffer(0), my(new Imp())  {
    stream = new ofstream(fileName.c_str());
    writeVersion(stream);
}

XMLWriter::~XMLWriter() {
    if (delStream && stream) {
        delete stream;
    }
}

//// write string data to the stream
void XMLWriter::write(const string& data, bool filterAllMarkups) {

    // Always filter ampersand - but lessthan is part of syntax so
    // filter that only for external data, not tag delineation
    string d = filterAllMarkups ? 
        purgeAmpersand(purgeLessThan(purgeGreaterThan(data))) : 
        purgeAmpersand(data);
    try {
        if (stream) {           
            *stream << d;
        }
        else {
            buffer->append(d);
        }
    }
    catch (...){
        throw ModelException("XMLWriter::write", 
                             "output stream error writing " + d);
    }
}

/** write an XML start tag <tag TYPE='objectType' attributes>
    where objectType = 
    convertToPublicRep(object)->getClass()->getName(). Returns non null
    IObjectSP if the object has not already been written out.
    Additionally, the returned object will be the public representation
    of the object if there is one. A new line is
    appended if appendNewLine is true and the object has not already
    been written out. 
    If object is null then no TYPE tag is written. */
IObjectConstSP XMLWriter::objectStart(const string&  id,
                                      const string&  attributes,
                                      const IObject* object,
                                      bool           appendNewLine){
   const static string routine = "XMLWriter::objectStart";
    try{
        IObjectConstSP pubObj;
        string data;
        if (!object){
            data = "<" + id + (attributes.empty()?
                               "": " "+attributes)+ ">";            
        } else {
            string modAttributes = attributes;
            /* why 2 in next line? 1 for the field and 1 for the SP which has
               been used to access the field. */
            if (my->allowBackReferences && object->getRefCount() > 2){
                // duplicates exist.
                string refName;
                bool alreadyWritten = my->getRefName(id, object, refName);
                string refAttr = REF_TYPE_ATTRIBUTE + "='" +
                    purgeLessThan(purgeGreaterThan(refName)) + "'";
                // Has it been written out before
                if (alreadyWritten){
                    string data = "<" + id + " " + refAttr + ">";
                    write(data, false);
                    return IObjectConstSP(); // return null
                }
                modAttributes += attributes.empty()? refAttr: " "+refAttr;
            }
            IObjectConstSP obj(IObjectConstSP::attachToRef(object)); 
            pubObj = IObjectConstSP(CObject::convertToPublicRep(obj));
            data = "<" + id + " " +  CObject::OBJECT_TYPE_ATTRIBUTE + "='" +
                purgeLessThan(purgeGreaterThan(pubObj->getClass()->getName())) +
                "'" + (modAttributes.empty()? "": " "+modAttributes)+ ">";
        }
        write(data, false);
        if (appendNewLine){
            if (stream) {
                *stream << endl;
            } else {
                write("\n");
            }
        }
        return pubObj;
    } catch (exception& e){
        throw ModelException(e, routine);
    }
}  
  
/** finish writing out an object. object can be null */
void XMLWriter::objectEnd(const string& id, const IObject* object){
    string label = "</" + id + ">\n";
    write(label, false);
    if (stream){
        stream->flush();
    }
}
    

/** called before writing out individual objects - starts process */
void XMLWriter::documentStart(const string& id) {
    // do nothing
}

/** called after writing out individual objects - ends process */
void XMLWriter::documentEnd(const string& id) {
    if (stream){
        stream->flush();
    }
}

/** start a comment */
void XMLWriter::commentStart() {
    write("<!--\n", false);
}

/** end a comment */ 
void XMLWriter::commentEnd() {
    write("-->\n", false);
}

//// simply write string data to the stream
void XMLWriter::write(const string& data) {
    // don't purge "markup" from a CDATA section
    write(data, !my->isCDATA);
}

/** write a null object to the stream */
void XMLWriter::writeNull(const string& id) {
    const static string routine = "XMLWriter::writeNullObject";
    try {
        string data = "<" + id + " " +  CObject::OBJECT_TYPE_ATTRIBUTE 
            + "='" + CNull::TYPE->getName() + "'>";
        write(data, false);
        objectEnd(id, 0);
    }
    catch (exception& e){
        throw ModelException(&e, routine);
    }       
}

/** start a CDATA section */
void XMLWriter::cdataStart() {
    write("<![CDATA[\n", false);
    my->isCDATA = true;
}

/** end a CDATA section */
void XMLWriter::cdataEnd() {
    write("]]>\n", false);
    my->isCDATA = false;
}


/** Invokes the outputWrite method on object using the stream with
    which the XMLOutputStream was created with */
void XMLWriter::outputWrite(const string&  linePrefix,
                            const string&  prefix,
                            const IObject* object) const {
    static const string routine = "XMLWriter::outputWrite";
    try{
        if (!object){
            throw ModelException(routine, "NULL Object");
        }
        ostringstream buf;
        object->outputWrite(linePrefix, prefix, buf);
        /* ffs */ const_cast<XMLWriter*>(this)->write(buf.str(), true);
    }
    catch (exception& e){
        throw ModelException(e, routine);
    }
}

/** way of disabling REF=... to allow easier comparison between XML files */
void XMLWriter::allowBackReferences(bool allow) {
    my->allowBackReferences = allow;
}

DRLIB_END_NAMESPACE
