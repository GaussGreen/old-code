#include "edginc/config.hpp"
#include "edginc/IDiffusionRecord.hpp"

DRLIB_BEGIN_NAMESPACE

void IDiffusionRecord::load(CClassSP& clazz)
{
	QUICK_CLASS_LOAD(IDiffusionRecord); // FIXME QUICK_IFACE_LOAD ? 
    SUPERCLASS(CObject);
}

REGISTER_MYCLASS(IDiffusionRecord);

DEFINE_TEMPLATE_TYPE_WITH_NAME("IDiffusionRecordArray", IDiffusionRecordArray);
INSTANTIATE_TEMPLATE(class MCARLO_DLL array<IDiffusionRecordSP _COMMA_ IDiffusionRecord>);

//////////////////////////////////////////
// define TYPEs for 2 and 3D arrays
DEFINE_TEMPLATE_TYPE_WITH_NAME("IDiffusionRecordArrayArray", IDiffusionRecordArrayArray);
DEFINE_TEMPLATE_TYPE_WITH_NAME("IDiffusionRecordArrayArrayArray", IDiffusionRecordArrayArrayArray);
//////////////////////////////////////////////////////////

DiffusionRecordBase::DiffusionRecordBase(const TDate& date, const string& type, const string& name, CClass const *   cl ) :
	IDiffusionRecord(cl),
	date(CInt::create(date)),
	type(CString::create(type)),
	name(CString::create(name))
{
}

string DiffusionRecordBase::header() const
{
	return string("_Date_") + date->toString() + "_" + type->toString() + "_" + name->toString();
}

void DiffusionRecordBase::load(CClassSP& clazz)
{
	QUICK_CLASS_LOAD(DiffusionRecordBase);
	SUPERCLASS(IDiffusionRecord);
	FIELD(date, "Tdate value of date");
	FIELD(type, "Asset's type");
	FIELD(name, "Asset's name");
}

REGISTER_MYCLASS(DiffusionRecordBase);

void DiffusionRecordBase::outputWrite(const string& linePrefix, const string& prefix, ostream& stream) const
{
	stream << linePrefix << prefix << header() << endl;
}


/////////////////////////////////////////////////////////////

DiffusionRecordCurrency::DiffusionRecordCurrency(
		const TDate& date,
		const string& type,
		const string& name,
		double  spot,
		double  discount,
		const vector<double>&  lnFutureDF,
		const vector<int>&  offsets,
		size_t  npoints,
		CClass const *   cl
		) :
		DiffusionRecordBase(date, type, name, cl),
		spot(CDouble::create(spot)),
		discount(CDouble::create(discount)),
        lnFutureDF(lnFutureDF.begin(), lnFutureDF.begin()+npoints),
        offsets(offsets.begin(), offsets.begin()+npoints)
        {
        }

void DiffusionRecordCurrency::load(CClassSP& clazz) {
	QUICK_CLASS_LOAD(DiffusionRecordCurrency);
	SUPERCLASS(DiffusionRecordBase);
	FIELD(spot, "FX spot value");
	FIELD(discount, "1/MM");
	FIELD(lnFutureDF, "expected factors");
	FIELD(offsets, "date offsets for expected factors");
}

REGISTER_MYCLASS(DiffusionRecordCurrency);

/////////////////////////////////////////////////////////////////////////
DiffusionRecordBasis::DiffusionRecordBasis(
	const TDate& date,
	const string& type,
	const string& name,
	double  spot,
	double  discount,
	const vector<double>&  lnFutureDF,
	const vector<int>&  offsets,
	size_t  npoints,
	CClass const *   cl ) :
	DiffusionRecordBase(date, type, name, cl),
	spot(CDouble::create(spot)),
	discount(CDouble::create(discount)),
	lnFutureDF(lnFutureDF.begin(), lnFutureDF.begin()+npoints),
	offsets(offsets.begin(), offsets.begin()+npoints)
{
}

void DiffusionRecordBasis::load(CClassSP& clazz) {
	QUICK_CLASS_LOAD(DiffusionRecordBasis);
	SUPERCLASS(DiffusionRecordBase);
	FIELD(spot, "FX spot value");
	FIELD(discount, "1/MM");
	FIELD(lnFutureDF, "expected factors");
	FIELD(offsets, "date offsets for expected factors");
}

REGISTER_MYCLASS(DiffusionRecordBasis);

/////////////////////////////////////////////////////////////////////////
DiffusionRecordEquity::DiffusionRecordEquity(
		const TDate& date,
		const string& type,
		const string& name,
		double spot,
		CClass const *  cl) :
	DiffusionRecordBase(date, type, name, cl),
	spot(CDouble::create(spot))
{
}

void DiffusionRecordEquity::load(CClassSP& clazz) {
	QUICK_CLASS_LOAD(DiffusionRecordEquity);
	SUPERCLASS(DiffusionRecordBase);
	FIELD(spot, "EQ spot value");
}

REGISTER_MYCLASS(DiffusionRecordEquity);


/////////////////////////////////////////////////////////////////////////
DiffusionRecordCredit::DiffusionRecordCredit(
		const TDate& date,
		const string& type,
		const string& name,
		double sdf,
		const vector<double>&  lnFutureSpread,
		const vector<int>&  offsets,
		size_t  npoints,
		CClass const *   cl) :
        DiffusionRecordBase(date, type, name, cl),
		sdf(CDouble::create(sdf)),
        lnFutureSpread(lnFutureSpread.begin(), lnFutureSpread.begin()+npoints),
        offsets(offsets.begin(), offsets.begin()+npoints)
{
}
void DiffusionRecordCredit::load(CClassSP& clazz) {
	QUICK_CLASS_LOAD(DiffusionRecordCredit);
	SUPERCLASS(DiffusionRecordBase);
	FIELD(sdf, "survival disc factor");
	FIELD(lnFutureSpread, "ln of the future spread");
	FIELD(offsets, "future date offsets");
}

REGISTER_MYCLASS(DiffusionRecordCredit);

/////////////////////////////////////////////////////////////////////////
DiffusionRecordEnergy::DiffusionRecordEnergy(
		const TDate& date,
		const string& type,
		const string& name,
		const vector<double>&  lnExpDF,
		const vector<int>&  offsets,
		size_t  npoints,
        CClass const *   cl ) :
        DiffusionRecordBase(date, type, name, cl),
        lnExpDF(lnExpDF.begin() + npoints, lnExpDF.end()),
        offsets(offsets.begin() + npoints, offsets.end())
{
}

void DiffusionRecordEnergy::load(CClassSP& clazz) {
	QUICK_CLASS_LOAD(DiffusionRecordEnergy);
	SUPERCLASS(DiffusionRecordBase);
	FIELD(lnExpDF, "lnExpDF");
	FIELD(offsets, "future date offsets");
}

REGISTER_MYCLASS(DiffusionRecordEnergy);


DRLIB_END_NAMESPACE
