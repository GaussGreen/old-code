//----------------------------------------------------------------------------
//
//
//----------------------------------------------------------------------------

#ifndef EDR_WRITER_HPP
#define EDR_WRITER_HPP

#include "edginc/Object.hpp"

DRLIB_BEGIN_NAMESPACE

/** The Writer class defines an interface for writing objects to a stream */
class TOOLKIT_DLL Writer {
public:
    /** called before writing out individual objects - starts process */
    virtual void documentStart(const string& id) = 0;

    /** called after writing out individual objects - ends process */
    virtual void documentEnd(const string& id) = 0;

    /** default implementation of writing out an object */
    virtual void write(const string& id, const IObject* object);

    /** write an XML start tag <tag TYPE='objectType' attributes>
        where objectType = 
        convertToPublicRep(object)->getClass()->getName(). Returns non null
        IObjectSP if the object has not already been written out.
        Additionally, the returned object will be the public representation
        of the object if there is one. A new line is
        appended if appendNewLine is true and the object has not already
        been written out. */
   virtual IObjectConstSP objectStart(const string&  id,
                                      const string&  attributes,
                                      const IObject* object,
                                      bool           appendNewLine) = 0;

    /** finish writing out an object */
    virtual void objectEnd(const string& id, const IObject* object) = 0;

    /** start a comment */
    virtual void commentStart() = 0;
 
    /** end a comment */
    virtual void commentEnd() = 0;
    
    /** write free-form data */
    virtual void write(const string& data) = 0;

    /** write out a null object */
    virtual void writeNull(const string& id) = 0;

    /** Invokes the outputWrite method on object using the stream with
        which the XMLOutputStream was created with */
    virtual void outputWrite(const string&  linePrefix,
                             const string&  prefix,
                             const IObject* object) const = 0;

    virtual ~Writer();
protected:
    Writer();
};

typedef refCountPtr<Writer> WriterSP;

DRLIB_END_NAMESPACE
#endif
