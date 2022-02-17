//
// C++ Interface: IQMCHelperDateTimeCache
//
// Description:
//
//
// Author: Vladimir A Grebinskiy <vladimir.a.grebinskiy@jpmorgan.com>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef IQMCHelperDateTimeCache_HPP
#define IQMCHelperDateTimeCache_HPP

//#include "edginc/config.hpp"
#include "edginc/DateTime.hpp"
#include "edginc/DECLARE.hpp"
#include "edginc/TemplateIdx.hpp"

DRLIB_BEGIN_NAMESPACE

class IQMCHelperDateTimeCache
{
    public:
        static const int npos = SpotIdx::npos; // returned if getIdx() fails

        virtual void add(const DateTimeArray& dates) = 0;
        virtual void add(DateTimeArrayConstSP dates) { add(*dates);}
        virtual DateTimeArray       getDates()  = 0; //< returns current merge of all dates
        virtual DateTime            getDate(int idx) = 0; //< returns reference(!) to the getDates()[idx]
        virtual int                 getIdx(const DateTime& date) = 0; //< logically it is offset of date inside getDates() or npos
        virtual void                cleanup() {} ///< hook for memory critical situatiuons; visible state is not changed
        virtual void                debug() {} ///< aux hook
        virtual void                trim(const DateTime& maxDate) {} /// more memory tricks

};

DECLARE_REF_COUNT(IQMCHelperDateTimeCache);

DRLIB_END_NAMESPACE

#endif // QMCHelperDateTimeCache_HPP
