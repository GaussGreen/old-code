//	CFastStream
//
//	This is a minimum standard replacement for std::stringstream.
//	I do not claim that it has all the functionality of the latter.
//  It is however much faster for how we are using it.
//
//  Author :			David Cuin
//
//////////////////////////////////////////////////////////////////////

#ifndef _FASTSTREAM_H
#define _FASTSTREAM_H

class CFastStream;

typedef CFastStream						xmlstreamer;	// You can replace this with a class derived from std::stringstream, so long as that class implements a (char*) operator

class CFastStream
{
public:
	CFastStream(void);
	CFastStream(const std::string& sz);
	virtual								~CFastStream(void);	
	CFastStream&						operator<<(const std::string& s);
	CFastStream&						operator<<(long n);
	operator char*(void) const;

	// Throwaway functions that perhaps should be in XmlStreamer but that would
	// mean moving XmlStreamer into a more global include space which causes
	// other problems.
	static void							FromXml(estring* psz, bool bLeaveAmpersand);
	static void							ToXml(estring* psz, bool bLeaveAmpersand);
	static HRESULT						WriteScalar(const std::string& szPrefix, const CComVariant& vIn, xmlstreamer& ssOut);
	static const std::string			s_szEncoding;
	
protected:	
	char*								m_sz;
	size_t								m_nWatermark;
	size_t								m_nSize;
};


#endif