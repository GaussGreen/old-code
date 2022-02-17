#ifndef FASTTDATEMAP_HPP
#define FASTTDATEMAP_HPP

#include "edginc/AssetData.hpp"
#include <vector>

DRLIB_BEGIN_NAMESPACE

namespace RM_Assets{


// A very fast lookup, but using additional (Dlast - Dfirst)*sizeof(short) bytes
template <typename T>
class FastTDateMap
{
public:

    FastTDateMap() {}
    ~FastTDateMap() {}

    void populate(  const std::vector<TDate> & dates,
                    const std::vector<T>     & values)
    {
        if (!dates.empty() && dates.size() == values.size() 
            && dates.back() >= dates.front())
        {
            m_dates   = dates;
            m_values  = values;
            m_dtstart = dates.at(0);
            m_idx.assign(dates.back() - dates.front() + 1, -1);

            for(size_t i=0, j=0; i<m_idx.size(); ++i)
            {
                if ((long)i + m_dtstart == dates.at(j))
                    ++j;

                m_idx.at(i) = j - 1;
            }

        }

    }

    void get_prev(TDate d,                        // (I)
                  TDate &d_prev, T &val_prev)     // (O)
    {
        size_t idx  = m_idx.at(d-m_dtstart);
        d_prev      = m_dates.at(idx);
        val_prev    = m_values.at(idx);
    }
    void get_next(TDate d,                        // (I)
                  TDate &d_next, T &val_next)     // (O)
    {
        size_t idx  = m_idx.at(d-m_dtstart);

        if(m_dates.at(idx) < d) // no scrolling if the date is exact
            ++idx;

        d_next      = m_dates.at(idx);
        val_next    = m_values.at(idx);
    }

    void get_bracket(TDate d,                        // (I)
                     TDate &d_prev, T &val_prev,     // (O)
                     TDate &d_next, T &val_next)     // (O)
    {
        size_t idx  = m_idx.at(d-m_dtstart);
        // note the offsets data file actually stores the end of data for this date
        // which is also the start of next date
        d_prev      = m_dates.at(idx);
        val_prev    = idx==0 ? 0 : m_values.at(idx-1);

        if(m_dates.at(idx) < d) // no scrolling if the date is exact
            ++idx;

        d_next      = m_dates.at(idx);
        val_next    = idx==0 ? 0 : m_values.at(idx-1);
    }

    int get_size(TDate& d)
    {
        std::vector<TDate>::iterator found = std::lower_bound(m_dates.begin(), m_dates.end(), d);
        if (found == m_dates.end())
            return -1;  // beyond end
        if (*found > d)
            d = *found;
        int idx = found - m_dates.begin();
        return m_values[idx] - (idx ? m_values[idx-1] : 0);
    }

    const std::vector<TDate>& get_dates() const { return m_dates; }

    bool is_valid(TDate d)
    {
        return !m_dates.empty() && d >= m_dtstart && d <= m_dates.back();
    }

private:
    TDate m_dtstart;
    std::vector<size_t> m_idx;
    std::vector<TDate> m_dates;
    std::vector<T>     m_values;

};


} // end of namespace RM_Assets
DRLIB_END_NAMESPACE


#endif //FASTTDATEMAP_HPP
