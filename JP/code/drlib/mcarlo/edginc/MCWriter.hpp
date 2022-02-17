//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : MCWriter.hpp
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

#ifndef MCWRITER_HPP
#define MCWRITER_HPP

#include "edginc/AssetData.hpp"
#include "edginc/FastTDateMap.hpp"
#include "edginc/smartPtr.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/XMLWriter.hpp"
#include "edginc/AtomicArray.hpp"

#include "edginc/IDiffusionRecord.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/FORWARD_DECLARE.hpp"

#include "edginc/SimpleEquity.hpp"
#include "edginc/FXAsset.hpp"
#include "edginc/BasisIndexCurve.hpp"
#include "edginc/UntweakableBasisIndexCurve.hpp"
#include "edginc/ICDSParSpreads.hpp"
#include "edginc/EnergyFuturesCurve.hpp"

#include <map>
#include <sstream>
#include <string>
#include <fstream>
#include <iostream>

#ifdef WIN32
#include <windows.h>
#endif

DRLIB_BEGIN_NAMESPACE
/// AssetInfoWritable hierarchy
namespace RM_Assets{

class MCARLO_DLL IAssetWritable : public IAsset
{
public:

    virtual ~IAssetWritable() {}

    virtual void                write_object_header(std::ostream& out_stream) = 0;
    virtual void                set_initial_record(TDate today, const std::vector<IAssetWritable*>& m_assets) = 0;    // this seems simpler
    virtual void                write_one_record( vector<double> &buf ) = 0;
    virtual size_t              get_size_for_date(TDate d) = 0;

    virtual void                write_debug_record( std::ostream& out_stream, bool sane = false ) = 0;
    virtual IDiffusionRecordSP  create_record() = 0; // create one record of diffusion results at a specific (path, time, asset)
};

DECLARE_REF_COUNT(IAssetWritable);


class MCARLO_DLL IAssetInfoWritable : public IAssetWritable, public AssetData
{
public:

    IAssetInfoWritable(const AssetData& data) : AssetData(data), valid(false) {}

    virtual ~IAssetInfoWritable() {}

    // top level base interface partially implemented
    const AssetName&    get_name()   const  { return m_name; }

    // writing interface partially implemented
    void write_asset_header(std::ostream& out_stream)
    {
        out_stream<<m_name<<" "<<m_asset_maturity<<" ";
    }
    bool    valid; // as we reuse assets, keep track if assets is relevant for the given timepoint
};

DECLARE_REF_COUNT(IAssetInfoWritable);

class MCARLO_DLL AssetInfoWritable_Currency : public IAssetInfoWritable, public AssetData_Currency
{
public:
    AssetInfoWritable_Currency(const AssetData& base_asset, const AssetData_Currency& ccy_asset, const YieldCurve* y, const FXAsset* f) :
      IAssetInfoWritable(base_asset), AssetData_Currency(ccy_asset), yc(y), fx(f) {}

    const AssetType    get_type() const { return get_typeA(); }
    static std::string typeAsString()   { return get_typeS(); }

    virtual void       set_initial_record(TDate today, const std::vector<IAssetWritable*>& m_assets) // compute once, write many times
    {
        set_cur_date(today);
        valid = true;
        m_ssSpotFX = fx ? fx->getSpot() : 1.;
        m_ssDiscount = 1.;
        for (size_t i=0, npoints = get_npoints(today); i < npoints; ++i)
            m_ssLnFutureDF[i] = log(yc->pv(DateTime(today + m_offsets[i], 0)));

    }

    // writing interface implemented
    void   write_one_record( vector<double> &buf )
    {
        if (get_cur_date() > get_maturity())
            return;

        buf.push_back(m_ssSpotFX);
        buf.push_back(m_ssDiscount);
        for (size_t i=0, npoints = get_npoints(get_cur_date()); i < npoints; ++i)
            buf.push_back(m_ssLnFutureDF[i]);
    }
    void   write_object_header(std::ostream& out_stream)
    {
        out_stream<<typeAsString()<<" ";
        write_asset_header(out_stream);
        out_stream<<m_curve_end<<" "<<m_offsets.size();
        for(size_t i=0; i<m_offsets.size(); ++i)
            out_stream<<" "<<m_offsets[i];
        out_stream<<std::endl;
    }
    size_t get_size_for_date(TDate d)
    {
        if (d > get_maturity())
            return 0;
        return get_npoints(d) + 2;
    }
    void write_debug_record( std::ostream& out, bool sane )
    {
        if (get_cur_date() > get_maturity())
            return;
        AssetData::dump(out, typeAsString());
        AssetData_Currency::dump(out, get_name(), get_npoints(get_cur_date()), sane);
    }
    
    IDiffusionRecordSP create_record() {
        // stores(Date_now, Asset, SpotFX, PathDF, LnExpDF, offsets, npoints (#valid points)
        return IDiffusionRecordSP( new DiffusionRecordCurrency(
                get_cur_date(), 
                typeAsString(), get_name(), 
                m_ssSpotFX, 
                m_ssDiscount, 
                m_ssLnFutureDF,
				vector<int>(m_offsets.begin(), m_offsets.end()), // FIXME: performance waste
				get_npoints(get_cur_date())));
    }
private:
    const YieldCurve* yc;
    const FXAsset* fx;
};

DECLARE_REF_COUNT(AssetInfoWritable_Currency);

class MCARLO_DLL AssetInfoWritable_Basis : public IAssetInfoWritable, public AssetData_Basis
{
public:
    AssetInfoWritable_Basis(const AssetData& base_asset, const AssetData_Basis& ccy_asset, const IBasisIndexCurve* b) :
      IAssetInfoWritable(base_asset), AssetData_Basis(ccy_asset), basis(b)
      {
          const UntweakableBasisIndexCurve* bs = dynamic_cast<const UntweakableBasisIndexCurve*>(basis);
          m_isAdditive = bs ? bs->isAdditiveBasis() : false;
      }

      const AssetType    get_type() const { return get_typeA(); }
      static std::string typeAsString()   { return get_typeS(); }

      virtual void        set_initial_record(TDate today, const std::vector<IAssetWritable*>& m_assets)
      {
          set_cur_date(today);
          valid = true;
          for (size_t i=0, npoints = get_npoints(today); i < npoints; ++i)
              m_ssFutureSpread[i] = basis->parSpread(DateTime(today + m_offsets[i], 0));
      }

      // writing interface implemented
      void   write_one_record( vector<double> &buf )
      {
          if (get_cur_date() > get_maturity())
              return;

          for (size_t i=0, npoints = get_npoints(get_cur_date()); i < npoints; ++i)
              buf.push_back(m_ssFutureSpread[i]);
      }
      void   write_object_header(std::ostream& out_stream)
      {
          out_stream<<typeAsString()<<" ";
          write_asset_header(out_stream);
          const YieldCurveWrapper& wrapper = basis->getRefCurve();
          out_stream << wrapper.get()->getCcy() << " " << (m_isAdditive ? 1 : 0) << " ";
          
          out_stream<<m_curve_end<<" "<<m_offsets.size();
          for(size_t i=0; i<m_offsets.size(); ++i)
              out_stream<<" "<<m_offsets[i];
          out_stream<<std::endl;
      }
      size_t get_size_for_date(TDate d)
      {
          if (d > get_maturity())
              return 0;
          return get_npoints(d);
      }
      void write_debug_record( std::ostream& out, bool sane )
      {
          if (get_cur_date() > get_maturity())
              return;
          AssetData::dump(out, typeAsString());
          AssetData_Basis::dump(out, get_name(), get_npoints(get_cur_date()), sane);
      }

    IDiffusionRecordSP create_record() {
        // stores(Date_now, Asset, SpotFX, PathDF, LnExpDF, offsets, npoints (#valid points)
        return IDiffusionRecordSP( new DiffusionRecordBasis(
                get_cur_date(), 
                typeAsString(), get_name(), 
                0, 
                0, 
                m_ssFutureSpread,
				vector<int>(m_offsets.begin(), m_offsets.end()), // FIXME: performance waste
				get_npoints(get_cur_date())));
    }
private:
    const IBasisIndexCurve* basis;
};

typedef refCountPtr<AssetInfoWritable_Basis> AssetInfoWritable_BasisSP;

class MCARLO_DLL AssetInfoWritable_Equity : public IAssetInfoWritable, public AssetData_Equity
{
public:
    AssetInfoWritable_Equity(const AssetData& base_asset, const AssetData_Equity& eq_asset, const SimpleEquity* e) :
      IAssetInfoWritable(base_asset), AssetData_Equity(eq_asset), eq(e)
      {
          yc = eq->getYCIsoCode();
      }

    const AssetType    get_type() const { return get_typeA(); }
    static std::string typeAsString()   { return get_typeS(); }

    virtual void       set_initial_record(TDate today, const std::vector<IAssetWritable*>& m_assets)
    {
        set_cur_date(today);
        valid = true;
        m_ssSpotEQ = eq->getSpot();
    }
    void   write_one_record( vector<double> &buf )
    {
        if (get_cur_date() > get_maturity())
            return;
        buf.push_back(m_ssSpotEQ);
    }
    void   write_object_header(std::ostream& out_stream)
    {
        out_stream  <<  typeAsString()  <<  " ";
        write_asset_header(out_stream);
        out_stream  <<  yc  << " ";

        // dividend info: borrow curve and dividend list
        const BorrowCurve* bc = eq->getEquity()->getBorrow().get();
        const DoubleArray&   rates = bc->getRates();
        const DateTimeArray& dates = bc->getDates();
        ASSERT(dates.size() == rates.size());

        const DividendList* div = eq->getDivList().get(); // not included for now

        out_stream << dates.size();
        out_stream << std::endl;
        for (int i=0; i<dates.size(); ++i)
            out_stream << dates[i].getDate() << " " << rates[i] << std::endl;

    }
    size_t get_size_for_date(TDate d)
    {
        if (d > get_maturity())
            return 0;
        return 1;
    }
    void write_debug_record( std::ostream& out, bool sane )
    {
        if (get_cur_date() > get_maturity())
            return;
        AssetData::dump(out, typeAsString());
        AssetData_Equity::dump(out, get_name(), sane);
    }

    IDiffusionRecordSP create_record() {
    // stores(Date_now, Asset, SpotFX, PathDF, LnExpDF, offsets, npoints (#valid points)
        return IDiffusionRecordSP( new DiffusionRecordEquity(
                    get_cur_date(), 
                    typeAsString(), get_name(), 
                    m_ssSpotEQ));
    }

    const SimpleEquity* eq; // for writing initial data
};

DECLARE_REF_COUNT(AssetInfoWritable_Equity);

class MCARLO_DLL AssetInfoWritable_Credit : public IAssetInfoWritable, public AssetData_Credit
{
public:
    AssetInfoWritable_Credit(const AssetData& base_asset, const AssetData_Credit& cr_asset, const ICDSParSpreads* a) :
          IAssetInfoWritable(base_asset), AssetData_Credit(cr_asset), cds(a) {}

    const AssetType    get_type() const { return get_typeA(); }
    static std::string typeAsString()   { return get_typeS(); }

    // writing interface implemented
    virtual void       set_initial_record(TDate today, const std::vector<IAssetWritable*>& m_assets)
    {
        set_cur_date(today);
        valid = true;
        m_ssNDP = 1.;
        DefaultRatesSP def = cds->defaultRates();
        DateTime valueDate(today, 0);
        for (size_t i=0, npoints = get_npoints(today); i < npoints; ++i)
            m_ssLnExpSDF[i] = log ( def->calcDefaultPV(valueDate, DateTime(today + m_offsets[i], 0)) );
    }
    void   write_one_record( vector<double> &buf )
    {
        if (!valid || get_cur_date() > get_maturity())
            return;
        buf.push_back(m_ssNDP);
        for (size_t i=0, npoints = get_npoints(get_cur_date()); i < npoints; ++i)
            buf.push_back(m_ssLnExpSDF[i]);
    }
    void   write_object_header(std::ostream& out_stream)
    {
        out_stream<<typeAsString()<<" ";
        write_asset_header(out_stream);
        out_stream<<m_curve_end<<" "<<m_offsets.size();
        for(size_t i=0; i<m_offsets.size(); ++i)
            out_stream<<" "<<m_offsets[i];
        out_stream<<std::endl;
    }
    size_t get_size_for_date(TDate d)
    {
        if (d > get_maturity())
            return 0;
        return get_npoints(d) + 1;
    }
    void write_debug_record( std::ostream& out, bool sane )
    {
        if (!valid || get_cur_date() > get_maturity())
            return;
        AssetData::dump(out, typeAsString());
        AssetData_Credit::dump(out, get_name(), get_npoints(get_cur_date()), sane);
    }
    IDiffusionRecordSP create_record() {
        // stores(Date_now, Asset, SpotFX, PathDF, LnExpDF, offsets, npoints (#valid points)
        return IDiffusionRecordSP( new DiffusionRecordCredit(
                get_cur_date(), 
                typeAsString(), get_name(), 
                m_ssNDP, 
                m_ssLnExpSDF,
				vector<int>(m_offsets.begin(), m_offsets.end()), // FIXME: performance waste
				get_npoints(get_cur_date())));
    }

    const ICDSParSpreads* cds;
};

DECLARE_REF_COUNT(AssetInfoWritable_Credit);


class MCARLO_DLL AssetInfoWritable_Energy : public IAssetInfoWritable, public AssetData_Energy
{
public:
    AssetInfoWritable_Energy(const AssetData& base_asset, const AssetData_Energy& energy_asset, const EnergyFuturesCurve* en) :
          IAssetInfoWritable(base_asset), AssetData_Energy(energy_asset), enrg(en) {}

    const AssetType    get_type() const { return get_typeA(); }
    static std::string typeAsString()   { return get_typeS(); }

    virtual void       set_initial_record(TDate today, const std::vector<IAssetWritable*>& m_assets)
    {
        set_cur_date(today);
        valid = true;
        EnergyUnderlyerConstSP underlyer = enrg->getEnergyUnderlyer();
        const EnergyFuturesCurve* tier1 = 0;
        if (! underlyer->getParentName().empty())
        {
            // this is really ugly and inefficient
            // TO DO: find better ways to do this
            for (std::vector<IAssetWritable*>::const_iterator i = m_assets.begin(); i != m_assets.end(); ++i)
            {
                const AssetInfoWritable_Energy* a = dynamic_cast<const AssetInfoWritable_Energy*>(*i);
                if (a && a->get_name() == underlyer->getParentName())
                {
                    tier1 = a->enrg;
                    break;
                }
            }
            ASSERT(tier1 != 0);
        }
        for (size_t i = get_startingpoint(today); i < m_ssLnFuturePrices.size(); ++i)
        {
            DateTime futdat(m_dates[i], 0);
            m_ssLnFuturePrices[i] = log( enrg->interpolate(futdat) + (tier1 ? tier1->interpolate(futdat) : 0) );
        }
    }

    // writing interface implemented
    void   write_one_record( vector<double> &buf )
    {
        if (get_cur_date() > get_maturity())
            return;
        for (size_t i = get_startingpoint(get_cur_date());
            i < m_ssLnFuturePrices.size(); ++i)
            buf.push_back(m_ssLnFuturePrices[i]);
    }
    void   write_object_header(std::ostream& out_stream)
    {
        out_stream<<typeAsString()<<" ";
        write_asset_header(out_stream);
        out_stream<<m_dates.size();
        for(size_t i=0; i<m_dates.size(); ++i)
            out_stream<<" "<<m_dates[i];
        out_stream<<std::endl;
    }
    size_t get_size_for_date(TDate d)
    {
        if (d > get_maturity())
            return 0;
        return m_ssLnFuturePrices.size() - get_startingpoint(d);
    }
    void write_debug_record( std::ostream& out, bool sane )
    {
        if (get_cur_date() > get_maturity())
            return;
        AssetData::dump(out, typeAsString());
        AssetData_Energy::dump(out, get_name(), get_startingpoint(get_cur_date()), sane);
    }
    IDiffusionRecordSP create_record() {
    // stores(Date_now, Asset, SpotFX, PathDF, LnExpDF, offsets, npoints (#valid points)
        return IDiffusionRecordSP( new DiffusionRecordEnergy(
            get_cur_date(), 
            typeAsString(), 
            get_name(), 
            m_ssLnFuturePrices, 
            vector<int>(m_dates.begin(), m_dates.end()), // FIXME: performance problem  
            get_startingpoint(get_cur_date()))); // pass in the first index
    }
    private:
    const EnergyFuturesCurve* enrg;

};

DECLARE_REF_COUNT(AssetInfoWritable_Energy);


} // end of namespace



// MCWriter
namespace RM_Assets{

typedef std::vector<IAssetWritable*> MCWriterDataHolder;



// basic class that know how to extract states and print them
// Currently it is the only available implementation suitable for MCWriter class initialization
class MCARLO_DLL MCWriterImpl
{
public:
    MCWriterImpl(const MCWriterDataHolder &_assets,
             const std::vector<TDate> &_dates);

    virtual ~MCWriterImpl() // does not assume ownership
    {}

    void writeHeader(std::ostream& out_data,
                     std::ostream& out_offsets);

    size_t writeInitialState(vector<double>& out_buf);
    size_t writeTheState(vector<double>& out_buf);       // move to a more specific class
    void   writeInitialDebugInfo(std::ostream& out, bool sane = false);
    void   writeStateDebugInfo(std::ostream& out, bool sane = false); // move to a more specific class

    size_t getTotalStoragePerPath() const;
    size_t getLongestSingleDateStorage() const;
    
    
    StateRecords  getAssetStateRecords() const; // return diffusion results for all assets

private:

    MCWriterImpl();

    MCWriterDataHolder  m_assets;
    std::map<TDate, long>         m_offsetmap;
    size_t                        m_full_size;

};

DECLARE_REF_COUNT(MCWriterImpl);


/** Root of the  MCWriter hirerachy
    The hirerachy has three main components:
    - the IMCWriter class that
      1. defines interface to the events we support
      2. defines an abstract  method to get an instance of implementation that provides methods needed to handle the events (i.e. MCWriterImpl class)
    - The MCWriter class that has sole responsibility to return the current implementation strategy
    - The MCWriterDecorator and classes derived from it.
      1. every class derived from Decorator provides a constructor with a smart pointer to another IMCWriter instance
      2. Every overriden method should call at the end the same-named method in the MCWriterDecorator base class so other classes in the chain have chance to react.

    In this way we can dynamically create a chain of writers: we want binary, text and XML dumps of diffusion process. The filename is usually passed in the constructor.
    TODO:
    - check if it's better to switch to the even based approach
    - redefine MCWriterImpl interface to expose class neutral methods to look into the SVs (currently there are two methods: for binary and text outputs)

*/



/** IMCWriter responsibilities:
   1. be the root class for MCWriter hierarchy
   2. provide list of notification events and the getImpl() method
*/

class MCARLO_DLL IMCWriter {
public:
    virtual MCWriterImpl*     getImpl(void) = 0;      // get state reader
    virtual void notifyStartHeader(void) = 0;
    virtual void notifyStartPath(int pathIdx) = 0; // new path started
    virtual void notifyEndPath(int pathIdx) = 0;
    virtual void notifyEndDate(int pathIdx, int iStep, const DateTime& iDateTime) = 0;
    virtual void finalize() = 0;
    virtual ~IMCWriter() {}

    enum
    {
        PathMode    = 1,
        SliceMode   = 1 << 1
    };
};

DECLARE_REF_COUNT(IMCWriter);

#ifdef WIN32
    typedef HANDLE  FileHandle;
#else
    typedef int     FileHandle;
    #define INVALID_HANDLE_VALUE    (-1)
#endif

/** MCWriter responsibilities
    1. Provide getImpl() method that returns the needed way to look into diffusion state
    2. Provide default reaction to notify() events (default is no action)
    3. Supposed to be the inner most class in the chain of decoration
 */

class MCARLO_DLL MCWriter : public IMCWriter
{
    MCWriterImpl impl;
public:
    MCWriter(const MCWriterDataHolder &_assets,
             const std::vector<TDate> &_dates) :
        impl(_assets, _dates)
    {
    }
    virtual MCWriterImpl*     getImpl(void) {return & impl;}      // get state reader
    virtual void notifyStartHeader(void) {}
    virtual void notifyStartPath(int pathIdx) {} // new path started
    virtual void notifyEndPath(int pathIdx) {}
    virtual void notifyEndDate(int pathIdx, int iStep,  const DateTime& iDateTime) {}
    virtual void finalize() {}
};

DECLARE_REF_COUNT(MCWriter);


// Decorator class, to derive all decorators from If derived class
// overrides one of the notify methods, it should call the same-named
// method from this (MCWriterDecorator) class to pass it along the
// chain.

/** MCWriterDecorator responsibilities:
    - helps to implement Decorator pattern
    - keeps pointer to the next IMCWriter object
    - forwards requests to the next object in chain
*/

class MCARLO_DLL MCWriterDecorator : public IMCWriter
{
    IMCWriterSP         obj;

protected:
    unsigned long       mode;
    vector<FileHandle>  fd; // only used for slice mode writer
                            // raw OS fd instead of FILE or ostream is used to avoid the limit imposed in those libraries
public:
    MCWriterDecorator(IMCWriterSP someobj, unsigned long mode=PathMode);
    virtual MCWriterImpl*  getImpl(void);
    virtual void notifyStartHeader(void);
    virtual void notifyStartPath(int pathIdx); // new path started
    virtual void notifyEndPath(int pathIdx);
    virtual void notifyEndDate(int pathIdx, int iStep, const DateTime& iDateTime);
    virtual void finalize();
};

// Binary Writer -- writes everything to a file as a binary stream
/** MCBinaryWriter responsibilities:
    - dump diffusion in binary format to files: headers and a stream of double[]. The numbers, currently, are in the host byte order
*/
class MCARLO_DLL MCBinaryWriter : public MCWriterDecorator {

    vector<double> buf; // binary buffer
    string         baseName;

public:
    MCBinaryWriter(string basename, IMCWriterSP next, unsigned long mode=PathMode);
    virtual void notifyStartHeader(void);       // writes header with binary file annotation
    virtual void notifyStartPath(int pathIdx);  // new path started
    virtual void notifyEndPath(int pathIdx);    // end of the path
    virtual void notifyEndDate(int pathIdx, int iStep,  const DateTime& iDateTime); // called w
    //virtual void finalize();
};

DECLARE_REF_COUNT(MCBinaryWriter);

// DebugWriter -- writes everything to a file as a readable text file
/** MCDebugWriter responsibilities:
    - dump diffusion in a human readable verbose format
*/
class MCARLO_DLL MCDebugWriter : public MCWriterDecorator {

    std::ofstream dout;
    string  baseName;
    bool    sane;

public:
    MCDebugWriter(string basename, IMCWriterSP next, unsigned long mode=PathMode, bool sane = false);

    virtual void notifyStartPath(int pathIdx); // new path started
    virtual void notifyEndPath(int pathIdx); // end of the path
    virtual void notifyEndDate(int pathIdx, int iStep, const DateTime& iDateTime); // called w
};

DECLARE_REF_COUNT(MCDebugWriter);

/** MCXMLWriter responsibilities :
    - dump diffusion results in an XML format
*/
class MCARLO_DLL MCXMLWriter : public MCWriterDecorator {

    CStringArraySP xml;
    XMLWriter xmlOut;
    std::ostringstream dout;

public:
    MCXMLWriter(string basename, IMCWriterSP next);

    virtual void notifyStartHeader(void);
    virtual void notifyStartPath(int pathIdx); // new path started
    virtual void notifyEndPath(int pathIdx); // end of the path
    virtual void notifyEndDate(int pathIdx, int iStep, const DateTime& iDateTime); // called w

    IObjectSP   getResults() const { return xml;}

    virtual ~MCXMLWriter();
};

DECLARE_REF_COUNT(MCXMLWriter);

//////////////////////////////////////////////////////////////

FORWARD_DECLARE(MCRegTestWriterImpl);

/** Regression  test compatible writer
*/


class MCARLO_DLL MCRegTestWriter : public MCWriterDecorator {

    MCRegTestWriterImplSP impl;
    
public:
    MCRegTestWriter(IMCWriterSP next);
    virtual ~MCRegTestWriter();
    virtual void notifyStartHeader(void);
    virtual void notifyStartPath(int pathIdx); // new path started
    virtual void notifyEndPath(int pathIdx); // end of the path
    virtual void notifyEndDate(int pathIdx, int iStep, const DateTime& iDateTime); // called w

    IObjectSP   getResults() const;
};

DECLARE_REF_COUNT(MCRegTestWriter);

} // end of namespace RM_Aseets



DRLIB_END_NAMESPACE
        
#endif //MCWRITER_HPP
