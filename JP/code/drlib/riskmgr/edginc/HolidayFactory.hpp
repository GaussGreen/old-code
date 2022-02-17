#ifndef HOLIDAYFACTORY_HPP
#define HOLIDAYFACTORY_HPP

#include <map>
#include <string>

#include <edginc/Holiday.hpp>

DRLIB_BEGIN_NAMESPACE

class RISKMGR_DLL HolidayFactory
{
public:

    static HolidayFactory& instance();

    // registers the holiday with the internal cache
    void registerHoliday(const std::string& name, const HolidaySP& calendar);

    void clear();

    // creates a holiday based on the name
    HolidaySP create(const std::string& name) const;

private:
    // Singleton semantics
    HolidayFactory();
    HolidayFactory(const HolidayFactory& rhs);
    HolidayFactory& operator=(const HolidayFactory& rhs);

    HolidaySP retrieve(const std::string& name) const;

    static const std::string delims_;

    std::map<std::string, HolidaySP> cache_;
};

DRLIB_END_NAMESPACE

#endif
