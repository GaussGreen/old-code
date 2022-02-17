//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : MCWriter.cpp
//
//   Description : Writer of MC-generated data for RM applications
//
//   Author      : Anatoly Morosov
//
//   Date        : October 2005
//
//
//
//----------------------------------------------------------------------------

#include "edginc/config.hpp"
#include "edginc/XMLWriter.hpp"
#include "edginc/MCWriter.hpp"


#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>
#include <cstring>
#include <iomanip>

// Note that Cygwin defines "unix" as well
#if defined(unix)
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>
#endif

using namespace std;

DRLIB_BEGIN_NAMESPACE

namespace RM_Assets {

string GetErrorMsg()
{
#ifdef  WIN32
    DWORD code = ::GetLastError();
    const char* buff;
    ::FormatMessage(FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_ALLOCATE_BUFFER, 0, code, 0, (LPSTR) &buff, 0, 0);
    return buff;
#else
    ostringstream ss;
    ss << "errno = " << errno;
    return ss.str();
#endif
}

MCWriterImpl::MCWriterImpl(const MCWriterDataHolder &_assets,
                   const std::vector<TDate> &_dates) :
    m_assets(_assets), m_offsetmap(), m_full_size(0)
{
    vector<size_t> asizes(_dates.size() + 1, 0);
    for(size_t i=0; i<m_assets.size(); ++i)
    {
        for(size_t d=0; d<_dates.size(); ++d)
            asizes[d] += m_assets[i]->get_size_for_date(_dates[d]);
    }

    for(size_t d=0; d<_dates.size(); ++d)
        m_offsetmap[_dates[d]] = (m_full_size += asizes[d]);
}

void MCWriterImpl::writeHeader(std::ostream& out_data, 
        std::ostream& out_offsets)
{
    for(size_t i=0; i<m_assets.size(); ++i)
        m_assets[i]->write_object_header(out_data);

    for(map<TDate,long>::const_iterator it = m_offsetmap.begin();
        it != m_offsetmap.end(); 
        ++it)
        out_offsets<<it->first<<" "<<it->second<<endl;
}

size_t MCWriterImpl::writeInitialState(vector<double> &out_buf)
{
    size_t start_size = out_buf.size();
    for(size_t i=0; i<m_assets.size(); ++i)
    {
        m_assets[i]->set_initial_record(m_offsetmap.begin()->first, m_assets);
        m_assets[i]->write_one_record(out_buf);
    }
    return out_buf.size() - start_size;
}

size_t MCWriterImpl::writeTheState(vector<double> &out_buf)
{
    size_t start_size = out_buf.size();
    for(size_t i=0; i<m_assets.size(); ++i)
    {
//         size_t start_size = out_buf.size();
        m_assets[i]->write_one_record(out_buf);
//         cout << "Asset " << i << " += " << (out_buf.size() - start_size) << " records" << endl;
    }
    return out_buf.size() - start_size;
}

void MCWriterImpl::writeInitialDebugInfo(std::ostream& out, bool sane)
{
    for(size_t i=0; i<m_assets.size(); ++i)
    {
        m_assets[i]->set_initial_record(m_offsetmap.begin()->first, m_assets);
        m_assets[i]->write_debug_record(out, sane);
    }
}

void MCWriterImpl::writeStateDebugInfo(std::ostream& out, bool sane)
{
    for(size_t i=0; i<m_assets.size(); ++i)
        m_assets[i]->write_debug_record(out, sane);
}

/// collect diffusion snapshots accross all assets.
StateRecords MCWriterImpl::getAssetStateRecords() const {
    StateRecords result;
    for(size_t i=0; i<m_assets.size(); ++i)
        result.push_back(m_assets[i]->create_record());
    return result;
}


size_t MCWriterImpl::getTotalStoragePerPath() const
{
    return m_full_size;
}

size_t MCWriterImpl::getLongestSingleDateStorage() const 
{
    // the longest one-date memory is required to store the very first
    // date, so the offsets between 1st and second dates will give it.

    if (m_offsetmap.size()<2)
        return m_full_size;
    else
    {
        map<TDate, long>::const_iterator it = m_offsetmap.begin();
        // assuming the offset of the first date is ALWAYS 0;
        return (++it)->second; 
    }
}


MCWriterDecorator::MCWriterDecorator(IMCWriterSP someObj, unsigned long m) :
    obj(someObj), mode(m)
{
}

MCWriterImpl * MCWriterDecorator::getImpl()
{
    if (obj.get())
        return obj->getImpl();
    else 
        return NULL;
}

void MCWriterDecorator::notifyStartHeader()
{
    if (obj.get())
        obj->notifyStartHeader();
}

void MCWriterDecorator::notifyStartPath(int pathIdx)
{
    if (obj.get())
        obj->notifyStartPath(pathIdx);
}

void MCWriterDecorator::notifyEndPath(int pathIdx)
{
    if (obj.get())
        obj->notifyEndPath(pathIdx);
}

void MCWriterDecorator::notifyEndDate(int pathIdx, int iStep,  const DateTime& iDate)
{
    if (obj.get())
        obj->notifyEndDate(pathIdx, iStep, iDate);
}

void MCWriterDecorator::finalize()
{
    if (mode & SliceMode)
    {
        // close the files if open
        for (size_t i=0; i<fd.size(); ++i)
        {
#ifdef WIN32
            CloseHandle(fd[i]);
#else
            close(fd[i]);
#endif
        }
    }

    if (obj.get())
        obj->finalize();
}

/*****************************************************************************
    MCBinaryWriter
    writes every path in a separate file
    Every path is concatenation of blocks for each diffusion timepoint;
******************************************************************************/
    
MCBinaryWriter::MCBinaryWriter(string basename, IMCWriterSP next, unsigned long m) :
    MCWriterDecorator(next, m),
    baseName(basename)
{
}

void MCBinaryWriter::notifyStartHeader(void)
{
    static const string method = "void MCBinaryWriter::notifyStartHeader(void)";
    string name1(baseName + "_header.dat");
    ofstream h1(name1.c_str());
    if (!h1)
        throw ModelException(method, "Could not create " + name1);
    h1.exceptions (std::ios_base::badbit | std::ios_base::failbit); // detect further IO errors
    
    string name2(baseName + "_offsets.dat");
    ofstream h2(name2.c_str());
    if (!h2)
        throw ModelException(method, "Could not create " + name2);
    h2.exceptions (std::ios_base::badbit | std::ios_base::failbit); // detect further IO errors
    
    getImpl()->writeHeader(h1, h2);
}

void MCBinaryWriter::notifyStartPath(int pathIdx) // new path started
{
    buf.clear();
    buf.reserve(getImpl()->getTotalStoragePerPath());
    notifyEndDate( pathIdx, -1, DateTime(/* don't care */) );

    MCWriterDecorator::notifyStartPath(pathIdx);
}

void MCBinaryWriter::notifyEndPath(int pathIdx) // end of the path
{
    if (mode & PathMode)
    {
        ostringstream name;
        name << baseName << pathIdx << ".dat";
        string s = name.str();
        
        FILE* bfile = fopen(s.c_str(), "wb");
        
        if (bfile == NULL)
        {
            int err = errno;
            string errMsg("ERROR: fopen(" + s + ") failed with error: " + strerror(err));
            throw ModelException(errMsg);
        }
        
        // This can't be always true due to maxCurveMat trim ?
        if (buf.size() != getImpl()->getTotalStoragePerPath()) {
            clog << "Error: pathIdx= " << pathIdx << " buf.size= " << buf.size() << " total_storage= " <<  getImpl()->getTotalStoragePerPath() << endl;
        }
        
        ASSERT(buf.size() == getImpl()->getTotalStoragePerPath());
        
        size_t ret = fwrite(& buf[0], sizeof(double), buf.size(), bfile);
        if (ret < buf.size() ||  ferror(bfile))
        {
            int err = errno;
            string errMsg(string("ERROR: fwrite() to " + s + " failed with error: ") + strerror(err));
            throw ModelException(errMsg);
        }
        
#if defined __USE_BSD || defined __USE_XOPEN
        if (fsync(fileno(bfile)) != 0)  // fclose() will not commit data to disk immediately; but we want to make sure there is no nasty surprises.
        {
            int err = errno;
            string errMsg(string("ERROR: fsync(") + s + ") failed with error: " + strerror(err));
            throw ModelException(errMsg);
        }
#endif

        if (fclose(bfile) != 0)
        {
            int err = errno;
            string errMsg(string("ERROR: fclose(") + s + ") failed with error: " + strerror(err));
            throw ModelException(errMsg);
        }
        
        bfile = NULL;
    }

    MCWriterDecorator::notifyEndPath(pathIdx);
}

void MCBinaryWriter::notifyEndDate(int pathIdx, int step, const DateTime& iDateTime) // called w
{
    size_t loc = buf.size();    // start of data for this date, needed for slice mode output
    if (step < 0)  // time 0 data
        getImpl()->writeInitialState(buf);
    else
        getImpl()->writeTheState(buf);

    if (mode & SliceMode)
    {
        size_t iStep = step + 1;

        // extend the fd vector if necessary
        if (fd.size() <= iStep)
        {
            size_t i = fd.size();
            fd.resize(iStep + 1);
            for ( ; i<fd.size(); ++i)
                fd[i] = INVALID_HANDLE_VALUE;
        }

        // open the file if not yet open
        if (fd[iStep] == INVALID_HANDLE_VALUE)
        {
            ostringstream ss;
            ss << baseName << "-slice" << iStep << ".dat";
            string filename = ss.str();
            fd[iStep] = 
#ifdef  WIN32
                CreateFile (filename.c_str(),
                            GENERIC_WRITE,
                            0,
                            0,
                            CREATE_ALWAYS,
                            FILE_ATTRIBUTE_NORMAL,
                            0);
            if (fd[iStep] == INVALID_HANDLE_VALUE)
                throw ModelException(string("cannot create file " + filename + " - " + GetErrorMsg()));
#else
            fd[iStep] = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
            if (fd[iStep] < 0)
            {
	    	    ostringstream ss;
		        ss << errno;
                    throw ModelException(string("cannot create file " + filename + " errno: " + ss.str()));
            }
#endif
        }

        // write data to file
#ifdef  WIN32
        DWORD bytes;
        BOOL success = WriteFile(fd[iStep], &buf[loc], sizeof(double) * (buf.size() - loc), &bytes, 0);
        if (success == 0 || sizeof(double) * (buf.size() - loc) != bytes)
            throw ModelException(string("MCBinaryWriter::notifyEndDate(): writing binary data failed - ") + GetErrorMsg());
#else
        size_t bytes = write(fd[iStep], &buf[loc], sizeof(double) * (buf.size() - loc));
        if (sizeof(double) * (buf.size() - loc) != bytes)
            throw ModelException(string("MCBinaryWriter::notifyEndDate(): writing binary data failed"));
#endif
    }

    if (step < 0)
        return;     // this function is called from notifyStartPath(), which itself is already chained, so no chaining this call

    MCWriterDecorator::notifyEndDate(pathIdx, step, iDateTime);
}

MCDebugWriter::MCDebugWriter(string basename, IMCWriterSP next, unsigned long m, bool s) :
    MCWriterDecorator(next, m),
    baseName(basename),
    sane(s)
{
}
    
void MCDebugWriter::notifyStartPath(int pathIdx)
{
    if (mode & PathMode)
    {
        ostringstream name;
        name << baseName << pathIdx << (sane ? ".sck" : ".txt");
        string str = name.str();

        dout.open(str.c_str());
        dout.setf(ios::fixed);
        dout.precision(12);

        notifyEndDate( pathIdx, -1, DateTime(/* don't care */) );
    }
    
    MCWriterDecorator::notifyStartPath(pathIdx);
}

void MCDebugWriter::notifyEndPath(int pathIdx) // end of the path
{
    if (mode & PathMode)
    {
        dout.close(); //TODO: check I/O
    }
    MCWriterDecorator::notifyEndPath(pathIdx);
}

void MCDebugWriter::notifyEndDate(int pathIdx, int step, const DateTime& iDateTime) // called w
{
    if (mode & PathMode)
    {
        if (step < 0)
            getImpl()->writeInitialDebugInfo(dout, sane);
        else
            getImpl()->writeStateDebugInfo(dout, sane); // TODO: check I/O
    }

    if (mode & SliceMode)
    {
        size_t iStep = step + 1;
        ostringstream  ss;
        ss.setf(ios::fixed);
        ss.precision(12);

        ss << "Path: " << pathIdx << endl;
        if (step < 0)
            getImpl()->writeInitialDebugInfo(ss, sane);
        else
            getImpl()->writeStateDebugInfo(ss, sane);

        string data = ss.str();

        // extend the fd vector if necessary
        if (fd.size() <= iStep)
        {
            unsigned int i = fd.size();
            fd.resize(iStep + 1);
            for ( ; i<fd.size(); ++i)
                fd[i] = INVALID_HANDLE_VALUE;
        }

        // open the file if not yet open
        if (fd[iStep] == INVALID_HANDLE_VALUE)
        {
            ostringstream ss;
            ss << baseName << "-slice" << iStep << (sane ? ".sck" : ".txt");
            string filename = ss.str();
            fd[iStep] = 
#ifdef  WIN32
                CreateFile (filename.c_str(),
                            GENERIC_WRITE,
                            0,
                            0,
                            CREATE_ALWAYS,
                            FILE_ATTRIBUTE_NORMAL,
                            0);
            if (! fd[iStep])
                throw ModelException(string("cannot create file " + filename));
#else
            fd[iStep] = open(filename.c_str(), O_WRONLY | O_CREAT | O_TRUNC, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
            if (fd[iStep] < 0)
                throw ModelException(string("cannot create file " + filename));
#endif
        }

        // write data to file
#ifdef  WIN32
        DWORD bytes;
        BOOL success = WriteFile(fd[iStep], data.c_str(), sizeof(char) * data.length(), &bytes, 0);
        if (success == 0 || sizeof(char) * data.length() != bytes)
            throw ModelException(string("writing debug data failed"));
#else
        size_t bytes = write(fd[iStep], data.c_str(), sizeof(char) * data.length());
        if (sizeof(char) * data.length() != bytes)
            throw ModelException(string("writing debug data failed"));
#endif
    }

    if (step < 0)
        return;     // this function is called from notifyStartPath(), which itself is already chained, so no chaining this call

    MCWriterDecorator::notifyEndDate(pathIdx, step, iDateTime);
}

USING_DRLIB_NAMESPACE 

//////////////////// write XML output ///////////////

MCXMLWriter::MCXMLWriter(string basename, IMCWriterSP next) :
    MCWriterDecorator(next),
    xmlOut(basename),
    xml(CStringArraySP(new CStringArray))
{
}
    
void MCXMLWriter::notifyStartHeader(void)
{
    ostringstream h1; // ofstream h1("QSPIheader.dat");
    ostringstream h2; // ofstream h2("QSPIoffsets.dat");

    getImpl()->writeHeader(h1, h2);

    xml->push_back(h1.str());
    xml->push_back(h2.str());

    MCWriterDecorator::notifyStartHeader();

}
void MCXMLWriter::notifyStartPath(int pathIdx)
{
    dout.str(""); // reset dout
    
    MCWriterDecorator::notifyStartPath(pathIdx);
}

void MCXMLWriter::notifyEndPath(int pathIdx) // end of the path
{
    xml->push_back(dout.str());
    MCWriterDecorator::notifyEndPath(pathIdx);
}

void MCXMLWriter::notifyEndDate(int pathIdx, int iStep, const DateTime& iDateTime) // called w
{

    getImpl()->writeStateDebugInfo(dout); // TODO: check I/O
    MCWriterDecorator::notifyEndDate(pathIdx, iStep, iDateTime);
}

MCXMLWriter::~MCXMLWriter()
{
    xml->write("SimpathIDiffusion", & xmlOut); // write dout to the XML stream
}
//////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////////////////
// forward declaration doesn't work with smartPtr<> :( Have to put definition here.

class MCARLO_DLL MCRegTestWriterImpl : public CObject
{
    MCWriterImpl *w;
    map<int, map< int, StateRecords> > fullDiffusion;
    static CClassConstSP const TYPE;
public:
        // Loop over all 
    MCRegTestWriterImpl(MCWriterImpl *_w) : CObject(TYPE), w(_w) {}
    virtual ~MCRegTestWriterImpl() {}
    virtual void outputWrite(const string& linePrefix, const string& prefix, ostream& stream) const
    {
        for(map<int, map< int, StateRecords> >::const_iterator itPath = fullDiffusion.begin();
            itPath != fullDiffusion.end();
            ++itPath)
            for(map< int, StateRecords>::const_iterator itTime = itPath->second.begin();
                itTime != itPath->second.end();
                ++itTime)
                for(StateRecords::const_iterator itAsset= itTime->second.begin();
                    itAsset != itTime->second.end();
                    ++itAsset)
                {
                    IDiffusionRecordSP rec = *itAsset;
                    string pref = prefix + "path=" + Format::toString(itPath->first) + "_step=" + Format::toString(itTime->first);
                    rec->outputWrite(linePrefix, pref, stream);
                }
    }
    
    virtual void notifyStartHeader(void) {}
    virtual void notifyStartPath(int pathIdx) {} // new path started
    virtual void notifyEndPath(int pathIdx) {} // end of the path
    virtual void notifyEndDate(int pathIdx, int iStep, const DateTime& iDateTime)
    {
        fullDiffusion[pathIdx][iStep] = w->getAssetStateRecords();
    }
    
};
DECLARE(MCRegTestWriterImpl);

CClassConstSP const MCRegTestWriterImpl::TYPE = CClass::registerInterfaceLoadMethod(
    "MCRegTestWriterImpl", typeid(MCRegTestWriterImpl), 0);

///////////////////////////////////////////////////////////////////////////////////////////////


MCRegTestWriter::MCRegTestWriter(IMCWriterSP next) :
    MCWriterDecorator(next),
    impl(new MCRegTestWriterImpl(getImpl()))
{
}
    

void MCRegTestWriter::notifyStartHeader(void) {
    impl->notifyStartHeader();
    MCWriterDecorator::notifyStartHeader();
}
        
void MCRegTestWriter::notifyStartPath(int pathIdx) {
    impl->notifyStartPath(pathIdx);
    MCWriterDecorator::notifyStartPath(pathIdx);
}

void MCRegTestWriter::notifyEndPath(int pathIdx) {
    impl->notifyEndPath(pathIdx);
    MCWriterDecorator::notifyEndPath(pathIdx);
    
}
void MCRegTestWriter::notifyEndDate(int pathIdx, int iStep, const DateTime& iDateTime) {
    impl->notifyEndDate(pathIdx, iStep, iDateTime);
    MCWriterDecorator::notifyEndDate(pathIdx, iStep, iDateTime);
}

IObjectSP   MCRegTestWriter::getResults() const {
    return impl;
}

MCRegTestWriter::~MCRegTestWriter()
{
}

} // end of RMAsset namespace



DRLIB_END_NAMESPACE
