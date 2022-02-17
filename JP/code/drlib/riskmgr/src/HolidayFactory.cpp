#include "edginc/config.hpp"

#include <edginc/HolidayFactory.hpp>

#include <edginc/Addin.hpp>
#include <edginc/Atomic.hpp>
#include <edginc/ModelException.hpp>

#include <boost/lexical_cast.hpp>
#include <fstream>
#include <set>

DRLIB_BEGIN_NAMESPACE

// helper funcions dec
DateTime dateTimeFromChaseString(const std::string& date);
std::vector<std::string> tokenise(std::string line, const std::string& delims);
template<class T> void unique_sort(T& container);

HolidayFactory::HolidayFactory() {}
HolidayFactory::HolidayFactory(const HolidayFactory& rhs) {}

HolidayFactory&
HolidayFactory::operator=(const HolidayFactory& rhs)
{
    return *this;
}

const std::string HolidayFactory::delims_ = "+";

HolidayFactory&
HolidayFactory::instance()
{
    static HolidayFactory instance_;

    return instance_;
}

void
HolidayFactory::registerHoliday(const std::string& name, const HolidaySP& calendar)
{
    if (calendar.get())
    {
        std::map<std::string, HolidaySP>::const_iterator it = cache_.find(name);
        if (it == cache_.end())
        {
            cache_.insert(std::make_pair(name, calendar));
        }
        else
        {
            throw ModelException("HolidayFactory::registerHoliday", "Holiday calendar'" + name + "' already exists!");
        }
    }
}

void
HolidayFactory::clear()
{
    cache_.clear();
}

HolidaySP
HolidayFactory::create(const std::string& name) const
{
    HolidaySP calendar;

    // handle composite calendars (e.g. 'NYC+LON')
    std::vector<std::string> calNames = tokenise(name, HolidayFactory::delims_);
    unique_sort(calNames);  // sorting makes NYC+TOK+LON == TOK+LON+NYC == LON+NYC+TOK, unique makes LON+LON+LON == LON
    std::string compositeName;
    for (std::size_t i = 0; i < calNames.size(); ++i)
    {
        compositeName += calNames[i];
    }

    calendar = retrieve(compositeName);

    if (!calendar.get())
    {
        // attempt to construct the composite calendar
        if (calNames.size() > 1)
        {
            DateTimeArray hols;
            bool weekends = false;

            for (std::size_t i = 0; i < calNames.size(); ++i)
            {
                calendar = retrieve(calNames[i]);
                if (calendar.get())
                {
                    bool we;
                    std::copy(calendar->toALIB(we)->begin(), calendar->toALIB(we)->end(), hols.back_inserter());
                    weekends |= we;
                }
                else
                {
                    throw ModelException("HolidayFactory::create", "Holiday calendar '" + calNames[i] + "' does not exist!");
                }
            }

            // We need to ignore time for holidays, which unique_sort is unable to do
            // so we call this function to set all dates to have the same times
            // and do a unique_sort on that!
            DateTime::setTimeOfDay(hols, DateTime::START_OF_DAY_TIME);

            unique_sort(hols);
            calendar = HolidaySP(new Holiday(compositeName, hols, weekends));
            HolidayFactory::instance().registerHoliday(compositeName, calendar);
        }
        else
        {
            throw ModelException("HolidayFactory::create", "Holiday calendar '" + name + "' does not exist!");
        }
    }

    return calendar;
}

HolidaySP
HolidayFactory::retrieve(const std::string& name) const
{
    HolidaySP calendar;

    std::map<std::string, HolidaySP>::const_iterator it = cache_.find(name);
    if (it != cache_.end())
    {
        calendar = (*it).second;
    }

    return calendar;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// Addin functions...

////////////////////////////////////////////////////////
// Create a calendar
struct CalendarAddin : public CObject
{
    static CClassConstSP const TYPE;

    std::string name;

    IObjectSP create()
    {
        HolidaySP calendar;

        try
        {
            calendar = HolidayFactory::instance().create(name);
        }
        catch (std::exception& e)
        {
            throw ModelException(e, "CalendarAddin::create");
        }

        return calendar;
    }

    /** for reflection */
    CalendarAddin() : CObject(TYPE) {}

    static void load(CClassSP& clazz)
    {
        REGISTER(CalendarAddin, clazz);
        SUPERCLASS(CObject);
        EMPTY_SHELL_METHOD(defaultConstructor);
        FIELD(name, "Calendar Name");
        Addin::registerObjectMethod("Calendar",
                                    Addin::UTILITIES,
                                    "Create a holiday calendar object from a string e.g. LON",
                                    false,
                                    Addin::returnHandle,
                                    &CalendarAddin::create);
    }

    static IObject* defaultConstructor()
    {
        return new CalendarAddin;
    }
};

CClassConstSP const CalendarAddin::TYPE = CClass::registerClassLoadMethod("CalendarAddin", typeid(CalendarAddin), load);

////////////////////////////////////////////////////////
// Load calendars from a file
struct LoadCalendarsAddin : public CObject
{
    static CClassConstSP const TYPE;
    std::string filename;

    std::string loadCalendars();

    /** for reflection */
    LoadCalendarsAddin() : CObject(TYPE) {}

    static void load(CClassSP& clazz);

    static IObject* defaultConstructor();
};

CClassConstSP const LoadCalendarsAddin::TYPE = CClass::registerClassLoadMethod("LoadCalendarsAddin",
                                                                               typeid(LoadCalendarsAddin),
                                                                               load);

void
LoadCalendarsAddin::load(CClassSP& clazz)
{
    REGISTER(LoadCalendarsAddin, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(filename, "Filename");
    Addin::registerStringMethod("LoadCalendars",
                                Addin::UTILITIES,
                                "Loads holidays from a file containing lines "
                                    "of the form: mm/dd/yyyy,calendarName "
                                    "e.g. 12/25/2006,LON",
                                &LoadCalendarsAddin::loadCalendars);
}

IObject*
LoadCalendarsAddin::defaultConstructor()
{
    return new LoadCalendarsAddin;
}

std::string
LoadCalendarsAddin::loadCalendars()
{
    static bool initialised = false;
    static const std::string funcName("LoadCalendarsAddin::loadCalendars");

    if (!initialised)
    {
        std::ifstream ifs(filename.c_str());
        if (!ifs)
        {
            throw ModelException(funcName, "Could not open input file '" + filename + "'");
        }

        try
        {
            // file is of the form : mm/dd/yyyy,calCode
            //    e.g. 12/25/2006,LON
            std::string line;
            std::vector<std::string> tokens;

            std::map<std::string, DateTimeArray> holidays;
            while (std::getline(ifs, line))
            {
                if (!line.empty() && line[0] != '#')
                {
                    tokens = tokenise(line, ", \t#");

                    if (tokens.size() < 2)
                    {
                        throw ModelException(funcName, "Each line must be of the form : 'mm/dd/yyyy,calendarName'");
                    }

                    holidays[tokens[1]].push_back(dateTimeFromChaseString(tokens[0]));
                }
            }

            std::map<std::string, DateTimeArray>::iterator it  = holidays.begin();
            std::map<std::string, DateTimeArray>::iterator end = holidays.end();
            for (; it != end; ++it)
            {
                if (!(*it).second.empty())
                {
                    unique_sort((*it).second);
                    HolidaySP calendar = HolidaySP(new Holiday((*it).first, (*it).second, true));
                    HolidayFactory::instance().registerHoliday((*it).first, calendar);
                }
            }

            ifs.close();
            initialised = true;
        }
        catch (const ModelException& e)
        {
            HolidayFactory::instance().clear();
            throw e;
        }
        catch (std::exception& e)
        {
            HolidayFactory::instance().clear();
            throw ModelException(e, funcName);
        }
    }

    return "Done!";
}

////////////////////////////////////////////////////////
// Create and register a calendar from scratch
struct CreateCalendarAddin : CObject
{
    static CClassConstSP const TYPE;

    std::string name;
    DateTimeArraySP holidays;
    bool useWeekends;
    
    IObjectSP createCalendar();

    /** for reflection */
    CreateCalendarAddin() : CObject(TYPE) {}

    static void load(CClassSP& clazz);

    static IObject* defaultConstructor();
};

CClassConstSP const CreateCalendarAddin::TYPE = CClass::registerClassLoadMethod("CreateCalendarAddin",
                                                                                typeid(CreateCalendarAddin),
                                                                                load);

void
CreateCalendarAddin::load(CClassSP& clazz)
{
    clazz->setPublic(); // make visible to EAS/spreadsheet
    REGISTER(CreateCalendarAddin, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(name, "Calendar name");
    FIELD(holidays, "holidays");
    FIELD(useWeekends, "Whether weekends are holidays");

    Addin::registerObjectMethod("CreateCalendar",
                                Addin::UTILITIES,
                                "Creates and registers a Calendar object",
                                false,
                                Addin::returnHandle,
                                &CreateCalendarAddin::createCalendar);
}

IObject*
CreateCalendarAddin::defaultConstructor()
{
    return new CreateCalendarAddin;
}

IObjectSP
CreateCalendarAddin::createCalendar()
{
    static const std::string funcName("CreateCalendarAddin::createCalendar");

    static std::set<std::string> createdCalendars;

    HolidaySP cal;

    if (name.find("+") != std::string::npos)
    {
        throw ModelException(funcName, "Cannot use character '+' as part of the calendar name");
    }

    if (createdCalendars.find(name) == createdCalendars.end())
    {
        cal = HolidaySP(new Holiday(name, *holidays, useWeekends));
        HolidayFactory::instance().registerHoliday(name, cal);

        createdCalendars.insert(name);
    }
    else
    {
        cal = HolidayFactory::instance().create(name);
    }

    return cal;
}

////////////////////////////////////////////////////////
// Register a pre-existing calendar with the HolidayFactory
struct RegisterCalendarAddin : CObject
{
    static CClassConstSP const TYPE;

    HolidaySP calendar;
    
    IObjectSP registerCalendar();

    /** for reflection */
    RegisterCalendarAddin() : CObject(TYPE) {}

    static void load(CClassSP& clazz);

    static IObject* defaultConstructor();
};

CClassConstSP const RegisterCalendarAddin::TYPE = CClass::registerClassLoadMethod("RegisterCalendarAddin",
                                                                                  typeid(RegisterCalendarAddin),
                                                                                  load);

void
RegisterCalendarAddin::load(CClassSP& clazz)
{
    clazz->setPublic();
    REGISTER(RegisterCalendarAddin, clazz);
    SUPERCLASS(CObject);
    EMPTY_SHELL_METHOD(defaultConstructor);
    FIELD(calendar, "Calendar");

    Addin::registerObjectMethod("RegisterCalendar",
                                Addin::UTILITIES,
                                "Registers a pre-created Holiday object with the holiday factory",
                                false,
                                Addin::returnHandle,
                                &RegisterCalendarAddin::registerCalendar);
}

IObject*
RegisterCalendarAddin::defaultConstructor()
{
    return new RegisterCalendarAddin;
}

IObjectSP
RegisterCalendarAddin::registerCalendar()
{
    static const std::string funcName("RegisterCalendarAddin::registerCalendar");

    static std::set<std::string> registeredCalendars;

    if (calendar.get())
    {
        if (calendar->getName().find("+") != std::string::npos)
        {
            throw ModelException(funcName, "Cannot register calendar which has '+' as part of the calendar name");
        }

        if (registeredCalendars.find(calendar->getName()) == registeredCalendars.end())
        {
            HolidayFactory::instance().registerHoliday(calendar->getName(), calendar);

            registeredCalendars.insert(calendar->getName());
        }
    }

    return calendar;
}

// helper functions
DateTime
dateTimeFromChaseString(const std::string& date)
{
    // Convert from Chase US style date "mm/dd/yyyy" to QLib DateTime
    std::vector<std::string> dt = tokenise(date, "/");
    if (dt.size() != 3)
    {
        throw ModelException("Date must be of the form: mm/dd/yyyy");
    }

    return DateTime::MonthDayYear(boost::lexical_cast<int>(dt[1]),
                                  boost::lexical_cast<int>(dt[0]),
                                  boost::lexical_cast<int>(dt[2])).toDateTime();
            }

std::vector<std::string>
tokenise(std::string line, const std::string& delims)
{
    std::vector<std::string> tokens;
    std::size_t pos;

    while (!line.empty())
    {
        pos = line.find_first_of(delims);

        if (pos == std::string::npos)
        {
            tokens.push_back(line);
            line.clear();
        }
        else
        {
            tokens.push_back(line.substr(0, pos));
            line = line.substr(++pos);
        }
    }

    return tokens;
}

template<class T>
void
unique_sort(T& container)
{
    std::sort(container.begin(), container.end());
    container.erase(std::unique(container.begin(), container.end()), container.end());
}

DRLIB_END_NAMESPACE
