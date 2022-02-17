//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : MCReader.hpp
//
//   Description : Reader of MC-generated data for RM applications
//                 (not to be #included inside QLib, but to be used by clients)
//
//   Author      : Anatoly Morosov
//
//   Date        : October 2005
//
//
//
//----------------------------------------------------------------------------
#ifndef MCREADER_HPP
#define MCREADER_HPP

#ifdef DRLIB_BEGIN_NAMESPACE
#    define CORE_CONFIG_HH      // to bypass coreConfig.hpp
#else
#include "edginc/config.hpp"
#endif

#include "edginc/AssetData.hpp"
#include "edginc/FastTDateMap.hpp"
//#include "edginc/smartPtr.hpp"
#include <map>
#include <sstream>
#include <string>
#include <algorithm>

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

using std::map;

DRLIB_BEGIN_NAMESPACE

/// AssetInfoReadable hierarchy
namespace RM_Assets{


class MomentsBase
{
public:
    virtual ~MomentsBase() {}
    virtual void calcAbsMoments(int num_sample) = 0;
    virtual void calcCentralMoments() = 0;
    virtual void calcStatistics()     = 0;
    virtual void dump(ostream& os)    = 0;
};

template<typename T>
class Moments : public MomentsBase, public vector<boost::shared_ptr<T> >
{
public:
    Moments(int maxOrderOfMoments) : vector<boost::shared_ptr<T> >(maxOrderOfMoments) {}
    void calcAbsMoments(int num_sample);
    void calcCentralMoments();
    void calcStatistics();
    void dump(ostream& os);
};


typedef double (PTRFUN3d_to_d)(double*& , double*& , double ) ;


inline double LinearInterpolation (double*& p1, double*& p2, double f)
{
    return (1.0-f)* *p1++ + f* *p2++;
}

inline double NoInterpolation (double*& p1, double*& p2, double f)
{
    return *p1++;
}

inline double NoInterpolation2 (double*& p1, double*& p2, double f)
{
    ++ p2;
    return *p1++;
}


class MCARLO_DLL IAssetReadable : public IAsset
{
public:

    virtual ~IAssetReadable() {}

    virtual IAssetReadable*     build_new_object(std::istream& in_stream) = 0;

    virtual void                read_for_date(    TDate d1, double* &buf1,
                                                  TDate d2, double* &buf2,
                                                  TDate d) = 0;

    virtual void                read(TDate d, const double*& buf) = 0;

    virtual void                dump(ostream& os, bool sane=false) = 0;

    // for statistical computations
    virtual void                accumulate(MomentsBase*) = 0;
    virtual void                reset() = 0;
};

typedef boost::shared_ptr<IAssetReadable> IAssetReadableSP;

class MCARLO_DLL IAssetInfoReadable : public IAssetReadable, public AssetData
{
public:

    IAssetInfoReadable(std::istream& in_stream) 
    {        
        in_stream>>m_name>>m_asset_maturity;
    }

    IAssetInfoReadable(){}
    
    virtual ~IAssetInfoReadable() {}

    // top level base interface partially implemented
    const AssetName&    get_name()   const  { return m_name; }

    // reading interface partially implemented
    void  read_for_date( TDate d1, double* &buf1,
                         TDate d2, double* &buf2,
                         TDate d)
    {
        set_cur_date(d2);
        
        if (d1 > get_maturity()) return; // do nothing past maturity
        if (d2 > get_maturity())         // do not extrapolate past maturity
            read_interpolated(  buf1, 
                                buf2,    // should not be touched in this function
                                0, 
                                NoInterpolation);

        else if (d1 == d2)               // interpolation not necessary
            read_interpolated(  buf1, 
                                buf2, 
                                0,
                                NoInterpolation2); 
        else                             // interpolate
            read_interpolated(  buf1, 
                                buf2, 
                                d, d1, d2, 
                                LinearInterpolation);

        set_cur_date(d);
    }


    virtual void read_interpolated( double* &buf1,
                                    double* &buf2,
                                    double f,
                                    PTRFUN3d_to_d *func) = 0;

    virtual void read_interpolated( double* &buf1,
                                    double* &buf2,
                                    TDate d, TDate d1, TDate d2,
                                    PTRFUN3d_to_d *func) = 0;
};

typedef boost::shared_ptr<IAssetInfoReadable> IAssetInfoReadableSP; 

#define ACCUMULATE(moments, field)                      \
    {                                                   \
        double term = field;                            \
        for (size_t i_=0; i_<(moments).size(); ++i_)    \
        {                                               \
            (moments)[i_]->field += term;               \
            term *= field;                              \
        }                                               \
    }


class MCARLO_DLL AssetInfo_Currency : public IAssetInfoReadable, public AssetData_Currency
{
public:
    AssetInfo_Currency(std::istream& in_stream) : IAssetInfoReadable(in_stream)
    {
        long npts = 0;
        in_stream>>m_curve_end>>npts;
        m_offsets.assign(npts,0);
        for(int i=0; i<npts; ++i)
            in_stream>>m_offsets[i];
        m_ssLnFutureDF.assign(npts,0.0);
    }
    AssetInfo_Currency(){}
    IAssetReadable* build_new_object(std::istream& in_stream)
    {
        return new AssetInfo_Currency(in_stream);
    }

    const AssetType    get_type() const { return get_typeA(); }
    static std::string typeAsString()   { return get_typeS(); } 

    // reading interface
    void read_interpolated( double* &buf1,
                            double* &buf2,
                            double f,
                            PTRFUN3d_to_d *func)
    {
        m_ssSpotFX = func(buf1, buf2, f);
        m_ssDiscount  = func(buf1, buf2, f);
        size_t npoints = get_npoints(get_cur_date());
        m_ssLnFutureDF.resize(npoints);     // throw away any invalid entry at end
        for (size_t i=0; i < npoints; ++i)
            m_ssLnFutureDF[i] = func(buf1, buf2, f);
    }

    void read_interpolated( double* &buf1,
                            double* &buf2,
                            TDate d, TDate d1, TDate d2,
                            PTRFUN3d_to_d *func)
    {
        double f = (d-d1)/(double)(d2-d1);
        m_ssSpotFX = func(buf1, buf2, f);
        m_ssDiscount  = func(buf1, buf2, f);
        size_t npoints = get_npoints(get_cur_date());
        m_ssLnFutureDF.resize(npoints);     // throw away any invalid entry at end
        for (size_t i=0; i < npoints; ++i)
            m_ssLnFutureDF[i] = func(buf1, buf2, f);
        buf1 += get_npoints(d1) - npoints;  // this is ugly but correct
    }

    // read without interpolation, useful for slice-wise data access
    void read(TDate d, const double*& buf)
    {
        set_cur_date(d);
        if (m_asset_maturity < d)
            return;     // nothing to read
        m_ssSpotFX      = *buf ++;
        m_ssDiscount    = *buf ++;
        size_t npoints = get_npoints(get_cur_date());
        m_ssLnFutureDF.resize(npoints);     // throw away any invalid entry at end
        for (size_t i=0; i < npoints; ++i)
            m_ssLnFutureDF[i] = *buf ++;
    }

    virtual void dump(ostream& os, bool sane=false)
    {
        if (m_asset_maturity < m_ssCurDate)
            return;
        AssetData::dump(os, get_typeS());
        AssetData_Currency::dump(os, m_name, m_ssLnFutureDF.size(), sane);
    }

    void resize()
    {
        m_ssLnFutureDF.resize(get_npoints(m_ssCurDate));
    }

    virtual void accumulate(MomentsBase* b)
    {
        Moments<AssetInfo_Currency>* m = dynamic_cast<Moments<AssetInfo_Currency>*>(b);
        if (m == 0)
            return;
        ACCUMULATE(*m, m_ssSpotFX)
        ACCUMULATE(*m, m_ssDiscount)
        for (size_t i=0; i < m_ssLnFutureDF.size(); ++i)
            ACCUMULATE(*m, m_ssLnFutureDF[i])
    }
    virtual void                reset()
    {
        m_ssSpotFX      = 0;
        m_ssDiscount    = 0;
        for (size_t i=0; i < m_ssLnFutureDF.size(); ++i)
            m_ssLnFutureDF[i] = 0;
    }
};

typedef boost::shared_ptr<AssetInfo_Currency> AssetInfo_CurrencySP;

class MCARLO_DLL AssetInfo_Basis : public IAssetInfoReadable, public AssetData_Basis
{
public:
    AssetInfo_Basis(std::istream& in_stream) : IAssetInfoReadable(in_stream)
    {
        int isAdditive;
        in_stream   >> m_underlyer >> isAdditive;
        m_isAdditive = isAdditive ? true : false;

        long npts = 0;
        in_stream   >> m_curve_end>>npts;
        m_offsets.assign(npts,0);
        for(int i=0; i<npts; ++i)
            in_stream >> m_offsets[i];
        m_ssFutureSpread.assign(npts,0.0);
    }
    AssetInfo_Basis(){}
    IAssetReadable* build_new_object(std::istream& in_stream)
    {
        return new AssetInfo_Basis(in_stream);
    }

    const AssetType    get_type() const { return get_typeA(); }
    static std::string typeAsString()   { return get_typeS(); } 

    // reading interface
    void read_interpolated( double* &buf1,
        double* &buf2,
        double f,
        PTRFUN3d_to_d *func)
    {
        size_t npoints = get_npoints(get_cur_date());
        m_ssFutureSpread.resize(npoints);     // throw away any invalid entry at end
        for (size_t i=0; i < npoints; ++i)
            m_ssFutureSpread[i] = func(buf1, buf2, f);
    }

    void read_interpolated( double* &buf1,
        double* &buf2,
        TDate d, TDate d1, TDate d2,
        PTRFUN3d_to_d *func)
    {
        double f = (d-d1)/(double)(d2-d1);
        size_t npoints = get_npoints(get_cur_date());
        m_ssFutureSpread.resize(npoints);     // throw away any invalid entry at end
        for (size_t i=0; i < npoints; ++i)
            m_ssFutureSpread[i] = func(buf1, buf2, f);
        buf1 += get_npoints(d1) - npoints;  // this is ugly but correct
    }

    // read without interpolation, useful for slice-wise data access
    void read(TDate d, const double*& buf)
    {
        set_cur_date(d);
        if (m_asset_maturity < d)
            return;     // nothing to read
        size_t npoints = get_npoints(get_cur_date());
        m_ssFutureSpread.resize(npoints);     // throw away any invalid entry at end
        for (size_t i=0; i < npoints; ++i)
            m_ssFutureSpread[i] = *buf ++;
    }

    virtual void dump(ostream& os, bool sane=false)
    {
        if (m_asset_maturity < m_ssCurDate)
            return;
        AssetData::dump(os, get_typeS());
        AssetData_Basis::dump(os, m_name, m_ssFutureSpread.size(), sane);
    }

    void resize()
    {
        m_ssFutureSpread.resize(get_npoints(m_ssCurDate));
    }

    virtual void accumulate(MomentsBase* b)
    {
        Moments<AssetInfo_Basis>* m = dynamic_cast<Moments<AssetInfo_Basis>*>(b);
        if (m == 0)
            return;
        for (size_t i=0; i < m_ssFutureSpread.size(); ++i)
            ACCUMULATE(*m, m_ssFutureSpread[i])
    }
    virtual void                reset()
    {
        for (size_t i=0; i < m_ssFutureSpread.size(); ++i)
            m_ssFutureSpread[i] = 0;
    }

    string              m_underlyer;
};

typedef boost::shared_ptr<AssetInfo_Basis> AssetInfo_BasisSP;

class MCARLO_DLL AssetInfo_Equity : public IAssetInfoReadable, public AssetData_Equity
{
public:
    AssetInfo_Equity(std::istream& in_stream) : IAssetInfoReadable(in_stream)
    {
        int num;
        in_stream   >> yc;
        in_stream   >> num;
        borrowCurve.resize(num);
        for (int i=0; i<num; ++i)
            in_stream >> borrowCurve[i].first >> borrowCurve[i].second;
    }
    AssetInfo_Equity(){}
    IAssetReadable* build_new_object(std::istream& in_stream)
    {
        return new AssetInfo_Equity(in_stream);
    }
    
    const AssetType    get_type() const { return get_typeA(); }
    static std::string typeAsString()   { return get_typeS(); } 

    void read_interpolated( double* &buf1,
                            double* &buf2,
                            double f,
                            PTRFUN3d_to_d *func)
    {
        m_ssSpotEQ = func(buf1, buf2, f);
    }
  
    void read_interpolated( double* &buf1,
                            double* &buf2,
                            TDate d, TDate d1, TDate d2,
                            PTRFUN3d_to_d *func)
    {
        double f = (d-d1)/(double)(d2-d1);
        m_ssSpotEQ = func(buf1, buf2, f);
    }

    // read without interpolation, useful for slice-wise data access
    void read(TDate d, const double*& buf)
    {
        set_cur_date(d);
        if (m_asset_maturity < d)
            return;     // nothing to read
        m_ssSpotEQ      = *buf ++;
    }

    virtual void dump(ostream& os, bool sane=false)
    {
        if (m_asset_maturity < m_ssCurDate)
            return;
        AssetData::dump(os, get_typeS());
        AssetData_Equity::dump(os, m_name, sane);
    }

    void resize()
    {}

    virtual void accumulate(MomentsBase* b)
    {
        Moments<AssetInfo_Equity>* m = dynamic_cast<Moments<AssetInfo_Equity>*>(b);
        if (m == 0)
            return;
        ACCUMULATE(*m, m_ssSpotEQ)
    }
    virtual void                reset()
    {
        m_ssSpotEQ = 0;
    }
};

typedef boost::shared_ptr<AssetInfo_Equity> AssetInfo_EquitySP;

class MCARLO_DLL AssetInfo_Credit : public IAssetInfoReadable, public AssetData_Credit
{
public:
    AssetInfo_Credit(std::istream& in_stream) : IAssetInfoReadable(in_stream)    
    {
        long npts = 0;
        in_stream>>m_curve_end>>npts;
        m_offsets.assign(npts,0);
        for(int i=0; i<npts; ++i)
            in_stream>>m_offsets[i];
        m_ssLnExpSDF.assign(npts,0.0);
    }
    AssetInfo_Credit(){}
    IAssetReadable* build_new_object(std::istream& in_stream)
    {
        return new AssetInfo_Credit(in_stream);
    }
    
    const AssetType    get_type() const { return get_typeA(); }
    static std::string typeAsString()   { return get_typeS(); } 
   
    void read_interpolated( double* &buf1,
                            double* &buf2,
                            double f,
                            PTRFUN3d_to_d *func)
    {
        m_ssNDP = func(buf1, buf2, f);
        size_t npoints = get_npoints(get_cur_date());
        m_ssLnExpSDF.resize(npoints);     // throw away any invalid entry at end
        for (size_t i=0; i < npoints; ++i)
            m_ssLnExpSDF[i] = func(buf1, buf2, f);
    }

    void read_interpolated( double* &buf1,
                            double* &buf2,
                            TDate d, TDate d1, TDate d2,
                            PTRFUN3d_to_d *func)
    {
        double f = (d-d1)/(double)(d2-d1);
        m_ssNDP = func(buf1, buf2, f); 
        size_t npoints = get_npoints(get_cur_date());
        m_ssLnExpSDF.resize(npoints);     // throw away any invalid entry at end
        for (size_t i=0; i < npoints; ++i)
            m_ssLnExpSDF[i] = func(buf1, buf2, f);
        buf1 += get_npoints(d1) - npoints;  // this is ugly but correct since buf1 possibly needs to advance more
    }

    // read without interpolation, useful for slice-wise data access
    void read(TDate d, const double*& buf)
    {
        set_cur_date(d);
        if (m_asset_maturity < d)
            return;     // nothing to read
        m_ssNDP      = *buf ++;
        size_t npoints = get_npoints(get_cur_date());
        m_ssLnExpSDF.resize(npoints);     // throw away any invalid entry at end
        for (size_t i=0; i < npoints; ++i)
            m_ssLnExpSDF[i] = *buf ++;
    }

    virtual void dump(ostream& os, bool sane=false)
    {
        if (m_asset_maturity < m_ssCurDate)
            return;
        AssetData::dump(os, get_typeS());
        AssetData_Credit::dump(os, m_name, m_ssLnExpSDF.size(), sane);
    }

    void resize()
    {
        m_ssLnExpSDF.resize(get_npoints(m_ssCurDate));
    }

    virtual void accumulate(MomentsBase* b)
    {
        Moments<AssetInfo_Credit>* m = dynamic_cast<Moments<AssetInfo_Credit>*>(b);
        if (m == 0)
            return;
        ACCUMULATE(*m, m_ssNDP)
        for (size_t i=0; i < m_ssLnExpSDF.size(); ++i)
            ACCUMULATE(*m, m_ssLnExpSDF[i])
    }
    virtual void                reset()
    {
        m_ssNDP = 0;
        for (size_t i=0; i < m_ssLnExpSDF.size(); ++i)
            m_ssLnExpSDF[i] = 0;
    }
};

typedef boost::shared_ptr<AssetInfo_Credit> AssetInfo_CreditSP;

class MCARLO_DLL AssetInfo_Energy : public IAssetInfoReadable, public AssetData_Energy
{
public:
    AssetInfo_Energy(std::istream& in_stream) : IAssetInfoReadable(in_stream)
    {
        long npts = 0;
        in_stream>>npts;
        m_dates.assign(npts, 0);
        for(int i=0; i<npts; ++i)
            in_stream>>m_dates[i];
        m_ssLnFuturePrices.assign(npts, 0.0);
    }
    AssetInfo_Energy(){}
    IAssetReadable* build_new_object(std::istream& in_stream)
    {
        return new AssetInfo_Energy(in_stream);
    }
    
    const AssetType    get_type() const { return get_typeA(); }
    static std::string typeAsString()   { return get_typeS(); } 
    
    void read_interpolated( double* &buf1,
                            double* &buf2,
                            double f,
                            PTRFUN3d_to_d *func)
    {
        for (size_t i = get_startingpoint(get_cur_date()); 
                    i < m_ssLnFuturePrices.size(); ++i)
            m_ssLnFuturePrices[i] = func(buf1, buf2, f);
    }

    void read_interpolated( double* &buf1,
                            double* &buf2,
                            TDate d, TDate d1, TDate d2,
                            PTRFUN3d_to_d *func)
    {
        double f = (d-d1)/(double)(d2-d1);
        for (size_t i = get_startingpoint(get_cur_date()); i < m_ssLnFuturePrices.size(); ++i)
            m_ssLnFuturePrices[i] = func(buf1, buf2, f);
            // what to do here?
            //  for energy, there is no need to adjust buf since the maturity date instead of the tenor
            //  are fixed and exactly the same number futures are simulated at any date
    }

    void read(TDate d, const double*& buf)
    {
        set_cur_date(d);
        if (m_asset_maturity < d)
            return;     // nothing to read
        for (size_t i = get_startingpoint(get_cur_date()); i < m_ssLnFuturePrices.size(); ++i)
            m_ssLnFuturePrices[i] = *buf ++;
    }

    virtual void dump(ostream& os, bool sane=false)
    {
        if (m_asset_maturity < m_ssCurDate)
            return;
        AssetData::dump(os, get_typeS());
        AssetData_Energy::dump(os, m_name, get_startingpoint(get_cur_date()), sane);
    }
    virtual void accumulate(MomentsBase* b)
    {
        Moments<AssetInfo_Energy>* m = dynamic_cast<Moments<AssetInfo_Energy>*>(b);
        if (m == 0)
            return;
        for (size_t i = get_startingpoint(get_cur_date()); i < m_ssLnFuturePrices.size(); ++i)
            ACCUMULATE(*m, m_ssLnFuturePrices[i])
    }

    void resize()
    {
        m_ssLnFuturePrices.resize(m_dates.size());
    }

    virtual void                reset()
    {
        for (size_t i = get_startingpoint(get_cur_date()); i < m_ssLnFuturePrices.size(); ++i)
            m_ssLnFuturePrices[i] = 0;
    }
};

typedef boost::shared_ptr<AssetInfo_Energy> AssetInfo_EnergySP;

} // end of namespace 


// MCReader
namespace RM_Assets{


class MCARLO_DLL MCReaderBuilder
{

public:
    MCReaderBuilder()
    {
        m_typemap[AssetInfo_Currency::typeAsString()] = new AssetInfo_Currency;
        m_typemap[AssetInfo_Equity  ::typeAsString()] = new AssetInfo_Equity;
        m_typemap[AssetInfo_Credit  ::typeAsString()] = new AssetInfo_Credit;
        m_typemap[AssetInfo_Energy  ::typeAsString()] = new AssetInfo_Energy;
        m_typemap[AssetInfo_Basis   ::typeAsString()] = new AssetInfo_Basis;
    }

    ~MCReaderBuilder()
    {
        for(std::map<std::string, IAssetReadable* >::iterator i=m_typemap.begin(); i!=m_typemap.end(); ++i)
            delete i->second;
    }

    IAssetReadable* buildOneAssetInfo(std::istream& in_stream);

    IAssetReadable* buildOneAssetInfo(const std::string& in_record)
    {
        std::istringstream stringStream(in_record);
        return buildOneAssetInfo(stringStream);
    }

private:
    std::map<std::string, IAssetReadable* > m_typemap;

};


class MCARLO_DLL MCReader
{

public:
    MCReader() : m_mcreader(), m_assets(), m_offsetmap()
    {}

    ~MCReader(); // owns and needs to clean IAssetReadables
    

    void readAssetList(std::istream & in_data, std::istream & in_offsets);

    // path == -1 
    //      indicates reading path-mode data files
    //      and thus the argument path is ignored
    //      will interpolate if date does not fall on a timeline date
    // path >= 0
    //      indicates reading slice-mode data files
    //      and thus the argument path is used to compute offset
    //      if date does not fall on a timeline date,
    //          it will be set to the next timeline date
    //          but no interpolation is performed
    // returns
    //      -1 if nothing is read because of empty data or date past last date
    //      0  otherwise
    int  readTheState(TDate& d, double* in_buf, int path=-1);

    typedef map<AssetName, IAssetReadableSP> AssetSubMap;
    const IAssetReadable* lookupAsset(AssetType type, const std::string& name) const;
    const AssetSubMap*    lookupAsset(AssetType type) const;
    const std::vector<IAssetReadableSP>&    getAssets() const { return m_assets; }

    void dump(ostream& os, bool sane=false)
    {
        for(size_t i=0; i<m_assets.size(); ++i)
            m_assets[i]->dump(os, sane);
    }

private:
    MCReaderBuilder m_mcreader;

    std::vector<IAssetReadableSP>    m_assets;
    std::map<AssetType, AssetSubMap> m_assets_map;
    FastTDateMap<long>      m_offsetmap;
};

} // end of namespace
DRLIB_END_NAMESPACE

#endif //MCREADER_HPP
