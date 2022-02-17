//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : IDiffusionRecord
//
//   Description : Provides information about one step of difussion (i.e. diffusion results for specific (path, asset, timepoint)
//
//----------------------------------------------------------------------------

#ifndef IDIFFUSIONRECORD_HPP
#define IDIFFUSIONRECORD_HPP

#include "edginc/DECLARE.hpp"
#include "edginc/Object.hpp"
#include "edginc/Atomic.hpp"
#include "edginc/AssetData.hpp"
#include "edginc/ConvenienceMacros.hpp"
#include <string>
#include <vector>
#include <ostream>

DRLIB_BEGIN_NAMESPACE

using std::vector;
using std::string;
using std::ostream;

using RM_Assets::TDate;


// namespace RM_Assets {



class MCARLO_DLL IDiffusionRecord : public CObject
{
	
public:
    IDiffusionRecord(CClass const *   cl = TYPE) :  CObject(cl) {}
	virtual void outputWrite(const string& linePrefix, const string& prefix, ostream& stream) const
    {
        QLIB_VERIFY(0, "Has to be implemented");
    }
	QUICK_CLASS_REFLECTION(IDiffusionRecord);
};

DECLARE(IDiffusionRecord);
DECLARE(IDiffusionRecordArray);
DECLARE(IDiffusionRecordArrayArray);


typedef vector<IDiffusionRecordSP> StateRecords;



class MCARLO_DLL DiffusionRecordBase : public IDiffusionRecord
{

	/*TDate*/  CIntSP        date;
	
	CStringSP       type;
	CStringSP       name;
protected:
    DiffusionRecordBase(CClass const *   cl = TYPE) : IDiffusionRecord(cl) {}
public:
	
	DiffusionRecordBase(const TDate& date, const string& type, const string& name, CClass const *   cl = TYPE);
	
	string header() const;
	
	virtual void outputWrite(const string& linePrefix, const string& prefix, ostream& stream) const;
	QUICK_CLASS_REFLECTION(DiffusionRecordBase);
};


class MCARLO_DLL DiffusionRecordCurrency : public DiffusionRecordBase
{
    
    CDoubleSP          spot;
	CDoubleSP          discount;
    CDoubleArray       lnFutureDF;
    CIntArray          offsets;
protected:	
    DiffusionRecordCurrency(CClass const *   cl = TYPE) : DiffusionRecordBase(cl) {}
public:
	
	DiffusionRecordCurrency(
		const TDate& date,
		const string& type,
		const string& name,
		double  spot,
		double  discount,
		const vector<double>&  lnFutureDF,
		const vector<int>&  offsets,
		size_t  npoints,
		CClass const *   cl = TYPE
		);
	virtual void outputWrite(const string& linePrefix, const string& prefix, ostream& stream) const 
	{
		string start = linePrefix + prefix + header();
		stream << start << "_spot: " << spot->toString() << std::endl;
		stream << start << "_discount: " << discount->toString() << std::endl;
		for(int i=0; i != lnFutureDF.size(); ++i)
			stream << start << "_lnFutureDF["<< Format::toString(offsets[i])<< "]: " <<  lnFutureDF[i] << std::endl;
	}
	QUICK_CLASS_REFLECTION(DiffusionRecordCurrency);
	
};

class MCARLO_DLL DiffusionRecordBasis : public DiffusionRecordBase
{
    
    CDoubleSP          spot,  discount;
    CDoubleArray  lnFutureDF;
    CIntArray    offsets;
protected:    
    DiffusionRecordBasis(CClass const *   cl = TYPE) : DiffusionRecordBase(cl) {}
		
public:
	DiffusionRecordBasis(
		const TDate& date,
		const string& type,
		const string& name,
		double  spot,
		double  discount,
		const vector<double>&  lnFutureDF,
		const vector<int>&  offsets,
		size_t  npoints,
		CClass const *   cl = TYPE);
	
	virtual void outputWrite(const string& linePrefix, const string& prefix, ostream& stream) const 
	{
	    string start = linePrefix + prefix + header();
	    stream << start << "_spot: " << spot->toString() << std::endl;
		stream << start << "_discount: " << discount->toString() << std::endl;
		for(int i=0; i != lnFutureDF.size(); ++i)
			stream << start << "_lnFutureDF["<< Format::toString(offsets[i])<< "]: " <<  lnFutureDF[i] << std::endl; 
		
	}
	QUICK_CLASS_REFLECTION(DiffusionRecordBasis);
	
};

class MCARLO_DLL DiffusionRecordEquity : public DiffusionRecordBase
{
    CDoubleSP spot;
protected:    
    DiffusionRecordEquity(CClass const *   cl = TYPE) : DiffusionRecordBase(cl) {}
public:
	DiffusionRecordEquity(
		const TDate& date,
		const string& type,
		const string& name,
		double spot,
		CClass const *   cl = TYPE);
	
	virtual void outputWrite(const string& linePrefix, const string& prefix, ostream& stream) const 
	{
	    string start = linePrefix + prefix + header();
		stream << start << "_spot: " << spot->toString() << std::endl;
	}
	QUICK_CLASS_REFLECTION(DiffusionRecordEquity);

};

class MCARLO_DLL DiffusionRecordCredit : public DiffusionRecordBase
{
    CDoubleSP sdf;
    CDoubleArray lnFutureSpread;
    CIntArray   offsets;
protected:    
	DiffusionRecordCredit(CClass const *   cl = TYPE) : DiffusionRecordBase(cl) {}
public:
	DiffusionRecordCredit(
		const TDate& date,
		const string& type,
		const string& name,
		double sdf,
		const vector<double>&  lnFutureSpread,
		const vector<int>&  offsets,
		size_t  npoints,
		CClass const *   cl = TYPE);
	
	virtual void outputWrite(const string& linePrefix, const string& prefix, ostream& stream) const 
	{
		string start = linePrefix + prefix + header();
		stream << start << "_sdf: " << sdf->toString() << std::endl;
		for(int i=0; i != lnFutureSpread.size(); ++i)
			stream << start << "_lnFutureSpread["<< Format::toString(offsets[i])<< "]: " <<  lnFutureSpread[i] << std::endl;
	}

	QUICK_CLASS_REFLECTION(DiffusionRecordCredit);
	
};

class MCARLO_DLL DiffusionRecordEnergy  : public DiffusionRecordBase
{
    CDoubleArray lnExpDF;
    CIntArray   offsets;
protected:
    DiffusionRecordEnergy(CClass const *   cl = TYPE) : DiffusionRecordBase(cl) {}
public:
	DiffusionRecordEnergy(
		const TDate& date,
		const string& type,
		const string& name,
		const vector<double>&  lnExpDF,
		const vector<int>&  offsets,
		size_t  npoints,
        CClass const *   cl = TYPE);
	
	virtual void outputWrite(const string& linePrefix, const string& prefix, ostream& stream) const 
	{
		string start = linePrefix + prefix + header();
		for(int i=0; i != lnExpDF.size(); ++i)
			stream << start << "_lnExpDF["<< Format::toString(offsets[i])<< "]: " <<  lnExpDF[i] << std::endl;
	}
	
	QUICK_CLASS_REFLECTION(DiffusionRecordEnergy);
	
};


// } // end of namespace RM_Aseets

DRLIB_END_NAMESPACE
        
#endif //MCWRITER_HPP
