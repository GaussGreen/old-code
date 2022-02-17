//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : XMLWriter.hpp
//
//   Description : Writer for XML
//
//   Author      : Andrew J Swain
//
//   Date        : 17 January 2001
//
//
//----------------------------------------------------------------------------

#ifndef EDR_XMLWRITER_HPP
#define EDR_XMLWRITER_HPP

#include "edginc/Writer.hpp"

DRLIB_BEGIN_NAMESPACE

class TOOLKIT_DLL XMLWriter: public Writer{
public:
    // used for identifying 'references'
    static const string REF_TYPE_ATTRIBUTE; // = "REF"

    XMLWriter(ostream* stream);
    // constructor from file name or buffer
    XMLWriter(string& xmlOut, bool isFile);
    // constructor from file name
    XMLWriter(const string& fileName);

    ~XMLWriter();

    /** called before writing out individual objects - starts process */
    virtual void documentStart(const string& id);

    /** called after writing out individual objects - ends process */
    virtual void documentEnd(const string& id);

    /** write an XML start tag <tag TYPE='objectType' attributes>
        where objectType = 
        convertToPublicRep(object)->getClass()->getName(). Returns non null
        IObjectSP if the object has not already been written out.
        Additionally, the returned object will be the public representation
        of the object if there is one. A new line is
        appended if appendNewLine is true and the object has not already
        been written out. 
        If object is null then no TYPE tag is written. */
    virtual IObjectConstSP objectStart(const string&  id,
                                       const string&  attributes,
                                       const IObject* object,
                                       bool           appendNewLine);

    /** finish writing out an object. object can be null */
    virtual void objectEnd(const string& id, const IObject* object);

    /** start a comment */
    virtual void commentStart();
 
    /** end a comment */
    virtual void commentEnd();

    /** simply write string data to the stream (all 'markups' are filtered) */
    void write(const string& data);

    /** write a null object to the stream */
    void writeNull(const string& id);

    /** xml specific - start a cdata section */
    void cdataStart();
 
    /** end cdata section */
    void cdataEnd();

    /** Invokes the outputWrite method on object using the stream with
        which the XMLWriter was created with */
    void outputWrite(const string&  linePrefix,
                     const string&  prefix,
                     const IObject* object) const;

    /** way of disabling REF=... to allow easier comparison between XML files */
    void allowBackReferences(bool allow);
    
private:
    class Imp;
    friend class Imp; // hide hash table headers
    //XMLWriter(XMLWriter&); // not wanted
    XMLWriter(const XMLWriter &rhs);
    XMLWriter& operator=(const XMLWriter& rhs);
    void write(const string&, bool);
    // fields
    ostream* stream;
    bool     delStream;
    string*  buffer;
    auto_ptr<Imp> my; // hides implementation
    
};

DRLIB_END_NAMESPACE

#endif
