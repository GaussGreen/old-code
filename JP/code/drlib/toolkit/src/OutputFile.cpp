//----------------------------------------------------------------------------
//
//   Group       : EDG Derivatives Research
//
//   Filename    : OutputFile.cpp
//
//   Description : Responsible for creation of .out files
//
//   Author      : Mark A Robson
//
//   Date        : 8 Feb 2001
//
//----------------------------------------------------------------------------


#include "edginc/config.hpp"
#include "edginc/OutputFile.hpp"
#include "edginc/XMLReader.hpp"

DRLIB_BEGIN_NAMESPACE

const string OutputFile::OUTPUT_FILE_TAG = "OutputFile";
const string OutputFile::OUTPUT_TAG = "Output";
const string OutputFile::SUMMARY_TAG = "Summary";
const string OutputFile::RESULTS_TAG = "Results";
const string OutputFile::EXCEPTION_TAG = "Exception";
const string OutputFile::PERMUTATION = "Permutation ";
const string OutputFile::LINE_PREFIX = "Output ";

/** Open output file, ready for writing results to it */
OutputFile::OutputFile(const string& filename):
    writer(filename), permNum(0), filename(filename)
{
    try{
        writer.objectStart(OUTPUT_FILE_TAG, "", 0, true);
    } catch (exception& e){
        throw ModelException(e, "OutputFile::OutputFile");
    }
}

string OutputFile::getFilename(void) const { return filename; }

void OutputFile::write(ModelException& e) {
    static const string method = "OutputFile::write";
    try {
        char* msg = e.stackTrace();
        try {
            char buffer[10];
            sprintf(buffer, "%d:\n", ++permNum);
            writer.write(PERMUTATION+buffer);
            writer.write(msg);
        }
        catch (exception& e) { 
            free(msg);
            throw ModelException(e, method);
        }
        free(msg);
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}

/** Writes object to the stream in the standard output format
    together with a string array which represents any exceptions
    thrown during the regression test run. The string array is
    written out within a <Exception> .. </Exception> block. If the
    string array is empty this block is skipped. Multiple objects
    can be written by repeated calls to this function. 
    The supplied object can be null provided exceptions is not empty */
void OutputFile::write(const IObject* object, const CStringArray& exceptions){
    static const string routine = "OutputFile::write";
    try{
        if (!object && exceptions.empty()){
            throw ModelException(routine, "NULL Object and no exceptions");
        }
        permNum++;
        writer.objectStart(OUTPUT_TAG, "", 0, true);
        writer.objectStart(SUMMARY_TAG, "", 0, true);
        char buffer[10];
        sprintf(buffer, "%d:\n", permNum);
        writer.write(PERMUTATION+buffer);
        if (!exceptions.empty()){
            writer.objectStart(EXCEPTION_TAG, "", 0, true);
            writer.cdataStart();
            for (int i = 0; i < exceptions.size(); i++){
                writer.write(exceptions[i]);
            }
            writer.cdataEnd();
            writer.objectEnd(EXCEPTION_TAG, 0);
        }
        if (object){
            writer.outputWrite(LINE_PREFIX, "", object);
        }
        writer.objectEnd(SUMMARY_TAG, 0);
        if (object){
            object->write(RESULTS_TAG, &writer); // todo: problematic
        }
        writer.objectEnd(OUTPUT_TAG, 0);
    }
    catch (exception& e){
        throw ModelException(e, routine);
    }    
}

/** Read an object back in from an XML serialised OutputFile */
IObjectSP OutputFile::xmlReadResults(const string& filename, int perm) {
    static const string method = "OutputFile::xmlReadResults";
    IObjectSP results;
    try{
        // get the root of the document
        XMLReader  xml(filename, true);
        Reader::NodeSP root(xml.root());

        // <Results> section is inside <Output> inside the root
        Reader::NodeListSP nl(root->children());
        int permSoFar = 0;

        // walk through root children
        for (unsigned int i = 0; i < nl->size() && !results; i++) {
            Reader::NodeSP child((*nl)[i]);

            string fieldname(child->name());

            // if we hit the <Output> section
            if (fieldname == OUTPUT_TAG) {
                Reader::NodeListSP nl2(child->children());
                permSoFar++;  // increase count of permutations read in

                if (permSoFar == perm) {
                    // walk through until we find <Results>
                    for (unsigned int j = 0; j < nl2->size() && !results; j++) {
                        Reader::NodeSP  ochild((*nl2)[j]);
                        if (ochild->name() == RESULTS_TAG) {
                            // read it in
                            results = xml.read(ochild.get());
                        }
                    }
                    // was it there at all ?
                    if (!results) {
                        throw ModelException(method, 
                                             "couldn't find "+RESULTS_TAG+ 
                                             " tag inside "+OUTPUT_TAG + 
                                             " tag");
                    }
                }
            }
        }
        // was it there at all ?
        if (!results) {
            throw ModelException(method, "couldn't find " + OUTPUT_TAG + 
                                 " tag inside document root");
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }  
    return results;
}

static Reader::NodeSP findChild(const Reader::NodeSP& node, 
                                const string&         toFind) {
    static const string method = "OutputFile::findChild";

    try {
        Reader::NodeSP child;

        // walk through children
        Reader::NodeListSP nl(node->children());
        for (unsigned int i = 0; i < nl->size(); i++) {
            child = (*nl)[i];
            string fieldname = child->name();
            if (fieldname == toFind) {
                return child;
            }
        }

        string inside = node->name();
        throw ModelException(method, "couldn't find " + toFind + 
                             " tag inside " + inside + " tag");
    }
    catch (exception& e) {
        throw ModelException(e, method);
    }
}


/** Read errors back in from an XML serialised OutputFile */
// this is pretty dreadful code right now
string OutputFile::xmlReadErrors(const string& filename, int perm){
    static const string method = "OutputFile::xmlReadErrors";
    string errors;
    try {
        // get the root of the document
        XMLReader  xml(filename, true);
        Reader::NodeSP root(xml.root());

        // <Exception> section is inside <Summary> inside <Output> 
        // inside the root

        Reader::NodeListSP nl(root->children());
        int permSoFar = 0;

        // walk through root children
        bool   found = false;
        for (unsigned int i = 0; i < nl->size() && !found; i++) {
            Reader::NodeSP child((*nl)[i]);

            string fieldname(child->name());

            // if we hit the <Output> section
            if (fieldname == OUTPUT_TAG) {
                found = true;
                Reader::NodeListSP nl2(child->children());
                permSoFar++;  // increase count of permutations read in
                
                if (permSoFar == perm) {
                    Reader::NodeSP summary(findChild(child, SUMMARY_TAG));
                    Reader::NodeSP error(findChild(summary, EXCEPTION_TAG));
                    errors = xml.cdata(error.get());
                }
            }
        }
        // was it there at all ?
        if (!found) {
            throw ModelException(method, "couldn't find " + OUTPUT_TAG + 
                                 " tag inside document root");
        }
    }
    catch (exception& e){
        throw ModelException(e, method);
    }  
    return errors;
}


OutputFile::~OutputFile(){
    writer.objectEnd(OUTPUT_FILE_TAG, 0);
}

/** Creates an output file name from an input file name */
string OutputFile::createOutputFileName(const string& inputFileName){
    // see if it ends in . something
    std::string::size_type idx = inputFileName.rfind(".");

    if (idx == std::string::npos) {
        return inputFileName + ".out";
    }
    else {
        return inputFileName.substr(0, idx) + ".out";
    }
}



DRLIB_END_NAMESPACE

