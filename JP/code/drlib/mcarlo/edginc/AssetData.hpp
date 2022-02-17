//----------------------------------------------------------------------------
//
//   Group       : QR&D
//
//   Filename    : AssetData.hpp
//
//   Description : Information for Writer/Reader of MC-generated data 
//                  for RM applications
//
//   Author      : Anatoly Morosov
//
//   Date        : October 2005
//
//
//
//----------------------------------------------------------------------------
#ifndef ASSETDATA_HPP
#define ASSETDATA_HPP

#ifdef _MSC_VER
#pragma warning( disable : 4786 )
#endif

#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include <iostream>
using std::ostream;
using std::vector;
using std::pair;
using std::string;

#define SANE(cond) (cond ? " (GOOD) " : " (BAD) " )

DRLIB_BEGIN_NAMESPACE

namespace RM_Assets{


enum AssetType
{
    UNKNOWN = 0,
    CURRENCY,
    EQUITY,
    CREDIT,
    ENERGY,
    BASIS
};


// currently... but may be extended later
typedef std::string AssetName;
typedef long TDate;

/// Common abstract root of AssetInfo classes
class MCARLO_DLL IAsset
{
public:
    virtual ~IAsset() {}
    virtual const AssetType     get_type() const = 0;
    virtual const AssetName&    get_name() const = 0;
};



class MCARLO_DLL AssetData
{
public:
    AssetData(const AssetName& name, TDate maturity) : 
      m_name(name), m_asset_maturity(maturity), m_ssCurDate(0) {}

    AssetData() : 
      m_name(), m_asset_maturity(0), m_ssCurDate(0) {}

    const TDate         get_maturity()  const { return m_asset_maturity; }
    const TDate         get_cur_date()  const { return m_ssCurDate; }
    void                set_cur_date(TDate d) { m_ssCurDate = d; }

    void dump(ostream& out, const string& type)
    {
        out << "[" << m_name << "] Date now =: "<< m_ssCurDate <<std::endl;
        out << "[" << m_name << "] Asset: "<< type  <<" " << m_name << std::endl;
    }

    AssetName           m_name;
    TDate               m_asset_maturity;
    TDate               m_ssCurDate;
};

class MCARLO_DLL AssetData_Currency
{
public:

    static AssetType   get_typeA() { return CURRENCY; }
    static std::string get_typeS() { return "Currency";} 

    AssetData_Currency() : m_curve_end(0), m_ssSpotFX(0), m_ssDiscount(0) {}
    AssetData_Currency(TDate curve_end, const std::vector<long>  &offsets) : 
      m_curve_end(curve_end), m_offsets(offsets),  
      m_ssSpotFX(0), m_ssDiscount(0), m_ssLnFutureDF(offsets.size(), 0.0) {}

    // local functions
    size_t get_npoints(TDate cur_date)
    {
        return std::upper_bound(m_offsets.begin(), 
                                m_offsets.end(), 
                                m_curve_end - cur_date) 
                        - m_offsets.begin();
    }

    void dump(ostream& out, const string& name, size_t npoints, bool sane)
    {
        out<<"[" << name << "] SpotFX: " << m_ssSpotFX   << (sane ? SANE(0<=m_ssSpotFX && m_ssSpotFX<1e6) : "") << std::endl;
        out<<"[" << name << "] PathDF: " << m_ssDiscount << (sane ? SANE(0<m_ssDiscount && m_ssDiscount<=1) : "") << std::endl;
        for (size_t i=0; i < npoints; ++i)
            out<<"[" << name << "] LnExpDF(now, now+" << m_offsets[i] << "): " << m_ssLnFutureDF[i] << (sane ? SANE(-1e+3 < m_ssLnFutureDF[i] && m_ssLnFutureDF[i] <= 0 && (i==0 || m_ssLnFutureDF[i]<=m_ssLnFutureDF[i-1])) : "") << std::endl;
        out << "[" << name << "] ------- * -------" << std::endl;
    }

    // fixed data
    TDate               m_curve_end;
    std::vector<long>   m_offsets;

    // snapshot data (updatable)
    double              m_ssSpotFX;
    double              m_ssDiscount;
    std::vector<double> m_ssLnFutureDF;

};

class MCARLO_DLL AssetData_Basis
{
public:

    static AssetType   get_typeA() { return BASIS; }
    static std::string get_typeS() { return "Basis";} 

    AssetData_Basis() : m_curve_end(0), m_isAdditive(true) {}
    AssetData_Basis(TDate curve_end, const std::vector<long>  &offsets) : 
        m_curve_end(curve_end), m_offsets(offsets), m_isAdditive(true), m_ssFutureSpread(offsets.size(), 0.0) {}

        // local functions
        size_t get_npoints(TDate cur_date)
        {
            return std::upper_bound(m_offsets.begin(), 
                m_offsets.end(), 
                m_curve_end - cur_date) 
                - m_offsets.begin();
        }

        void dump(ostream& out, const string& name, size_t npoints, bool sane)
        {
            for (size_t i=0; i < npoints; ++i)
                out << "[" << name << "] FutureSpread(now, now+" << m_offsets[i] << "): " 
                    << m_ssFutureSpread[i] << (sane ? SANE((m_isAdditive && -.5 < m_ssFutureSpread[i] && m_ssFutureSpread[i] < .5) || ( 0. < m_ssFutureSpread[i] && m_ssFutureSpread[i] < 5. && ! m_isAdditive)) : "")
                    << std::endl;
            out << "[" << name <<"] ------- * -------" << std::endl;
        }

        // fixed data
        TDate               m_curve_end;
        std::vector<long>   m_offsets;
        bool                m_isAdditive;

        // snapshot data (updatable)
        std::vector<double> m_ssFutureSpread;

};

class MCARLO_DLL AssetData_Equity 
{
public:

    static AssetType   get_typeA() { return EQUITY; }
    static std::string get_typeS() { return "Equity";} 

#if 0
<<<<<<< .working
    AssetData_Equity(const string& ccy = "") :  yc(ccy), m_ssSpotEQ(0) {}
=======
#else
    AssetData_Equity() : m_ssSpotEQ(0)    // this will be accessible from outside qlib, so don't use any QLIB data type
    {}
#endif
//>>>>>>> .merge-right.r66336

    size_t get_npoints(TDate cur_date){ return 0; }

    void dump(ostream& out, const string& name, bool sane)
    {
        out << "[" << name <<"] SpotEQ: " << m_ssSpotEQ << (sane ? SANE(0<=m_ssSpotEQ && m_ssSpotEQ<=1e6) : "") << std::endl;
        out << "[" << name <<"] ------- * -------" << std::endl;
    }

    // static data (not updated in diffusion)
    string                       yc;            // pay ccy, used only in reader
    vector<pair<TDate, double> > borrowCurve;   // dividend info, used only in reader

    // snapshot data (updatable)
    double  m_ssSpotEQ;
};

class MCARLO_DLL AssetData_Credit 
{
public:
    static AssetType   get_typeA() { return CREDIT; }
    static std::string get_typeS() { return "Credit";} 

    AssetData_Credit() : m_curve_end(0), m_ssNDP(0) {}
    AssetData_Credit(TDate curve_end, const std::vector<long>  &offsets) : 
        m_curve_end(curve_end), m_offsets(offsets),  
        m_ssNDP(0), m_ssLnExpSDF(offsets.size(), 0.0) {}

    size_t get_npoints(TDate cur_date)
    {
        return std::upper_bound(m_offsets.begin(), 
                                m_offsets.end(), 
                                m_curve_end-cur_date) 
                        - m_offsets.begin();
    }

    void dump(ostream& out, const string& name, size_t npoints, bool sane)
    {
        out << "[" << name <<"] NDP(t0, now): " << m_ssNDP << (sane ? SANE(0<=m_ssNDP && m_ssNDP<=1) : "") <<std::endl;
        for (size_t i=0; i < npoints; ++i)
            out << "[" << name <<"] LnExpSDF(now, now+" << m_offsets[i] << "): "
            <<m_ssLnExpSDF[i] << (sane ? SANE(-1e3 < m_ssLnExpSDF[i] && m_ssLnExpSDF[i] <=0 && (i==0 || m_ssLnExpSDF[i]<=m_ssLnExpSDF[i-1])) : "") <<std::endl;
        out << "[" << name <<"] ------- * -------" << std::endl;
    }

    // fixed data
    TDate               m_curve_end;
    std::vector<long>   m_offsets;

    // snapshot data (updatable)
    double              m_ssNDP;
    std::vector<double> m_ssLnExpSDF;
};

class MCARLO_DLL AssetData_Energy 
{
public:
    static AssetType   get_typeA() { return ENERGY; }
    static std::string get_typeS() { return "Energy";} 

    AssetData_Energy() {}
    AssetData_Energy(const std::vector<TDate>  &dates) : 
    m_dates(dates), m_ssLnFuturePrices(dates.size(),0.0) {}

    size_t get_startingpoint(TDate cur_date)
    {
        return std::lower_bound(m_dates.begin(), 
                                m_dates.end(), 
                                cur_date) 
                - m_dates.begin();
    }

    void dump(ostream& out, const string& name, size_t starting_point, bool sane)
    {
        for (size_t i = starting_point; i < m_ssLnFuturePrices.size(); ++i)
            out << "[" << name <<"] LnFuturesPrice[" << m_dates[i] << "] : " << m_ssLnFuturePrices[i] << std::endl;
        out << "[" << name <<"] ------- * -------" << std::endl;
    }

    // fixed data
    std::vector<TDate>   m_dates;
   
    // snapshot data (updatable)
    std::vector<double> m_ssLnFuturePrices;

};


} // end of namespace RM_Assets

DRLIB_END_NAMESPACE

#endif //ASSETDATA_HPP
