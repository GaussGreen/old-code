/////////////////////////////////////////////////////////////////////////////
// Name:        tests/xll/xll.cpp
// Purpose:     test Excel add-in
// Author:      Vadim Zeitlin
// Created:     2006-03-21
// RCS-ID:      $Id: xll.cpp,v 1.20 2006/05/05 13:31:00 zeitlin Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/XL/addin.h"
#include "ito33/XL/sheetobject.h"

#include "ito33/numeric/interpolation.h"

#include <map>
#include <vector>

using namespace ito33;

class TempCurveManager;

// ----------------------------------------------------------------------------
// the main add-in class
// ----------------------------------------------------------------------------

class TestAddIn : public XL::AddIn
{
public:
  TestAddIn() : XL::AddIn("Test Excel Add-in", "Temperature")
  {
    m_tempCurves = NULL;
  }

  virtual ~TestAddIn();

  TempCurveManager& GetTempCurveManager();

  virtual void RegisterAllFunctions();

private:
  TempCurveManager *m_tempCurves;
};

XL_IMPLEMENT_ADDIN(TestAddIn)

// ============================================================================
// temperature curves
// ============================================================================

class TempCurve
{
public:
  TempCurve()
  {
    XLL_TRACE("Creating TempCurve(%p)", this);
  }

  double GetValueAt(double d) const;

  ~TempCurve()
  {
    XLL_TRACE("Destroying TempCurve(%p)", this);
  }

  std::vector<double> m_dates,
                      m_temps;
};

class TempCurveManager : public XL::SheetObjectManager<TempCurve> { };

double TempCurve::GetValueAt(double d) const
{
  double t;
  numeric::LinearInterpolate(&m_dates[0], &m_temps[0], m_dates.size(), &d, &t, 1);
  return t;
}

// ============================================================================
// exported functions implementation
// ============================================================================

XLLEXPORT XLOPER *WhatIs(XLOPER *xlo_)
{
  const XL::Oper * const xlo = static_cast<XL::Oper *>(xlo_);

  XLL_TRACE("Input parameter is \"%s\"", xlo->Dump().c_str());
  return XL::Oper::Return(xlo->Dump());
}

XLLEXPORT XLOPER *WhatIsDeref(const XL::Oper *xlo)
{
  XLL_TRACE("Input parameter is \"%s\"", xlo->Dump().c_str());
  return XL::Oper::Return(xlo->Dump());
}

XLLEXPORT double Celsius(double f)
{
  return ((f - 32)*5)/9.;
}

XLLEXPORT bool IsHotter(double c, double f)
{
  return c > Celsius(f);
}

XLLEXPORT XLOPER *Hottest(const XL::Oper *xlo)
{
  try
  {
    XLL_TRACE("Input parameter is \"%s\"", xlo->Dump().c_str());

    XL::Oper xloArray;
    if ( !xlo->AsArray(&xloArray) )
      return XL::Oper::ErrValue();

    unsigned w, h;
    const XL::Oper *data = xloArray.GetAsArray(w, h);

    XLL_TRACE("Input range is %lux%lu", w, h);

    double hottest = 0.;
    for ( unsigned n = 0; n < w*h; n++, data++ )
    {
      double t;
      if ( data->As(&t) && t > h )
        hottest = t;
    }

    return XL::Oper::Return(hottest);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT const char *DegreesScale(const char *abbrev)
{
  const char ch = abbrev[0];
  if ( ch != '\0' && abbrev[1] == '\0' )
  {
    if ( ch == 'f' || ch == 'F' )
      return "Fahrenheit";
    if ( ch == 'c' || ch == 'C' )
      return "Celsius";
  }

  return "?";
}

XLLEXPORT XLOPER *CurrentTempC()
{
  const XL::Oper xloFormula("RTD(\"ITO33RTDTest.Server\", \"\", \"spot\", 20)");

  XL::Oper *xloRC = new XL::Oper;
  XL::Call(xlfEvaluate, XL::ByRef(*xloRC), xloFormula);
  return xloRC->Detach();
}

// this is a bit simplistic but will do for the test purposes
static const double tempInMonth[Date::MONTHS_IN_YEAR] =
{
  -5, -10, 5, 10, 15, 20, 25, 30,  15, 10, 5, 0,
};

XLLEXPORT double TempAtDate(XL::Oper *date)
{
  Date dateVal;
  if ( !date || !date->As(&dateVal) )
    return 0.;

  return tempInMonth[dateVal.GetMonth() - Date::Jan];
}

XLLEXPORT XLOPER *TempCurveCreate(XL::Oper *curve)
{
  try
  {
    TempCurveManager& tempCurveManager = GetTheTestAddIn().GetTempCurveManager();

    tempCurveManager.PruneInactive();

    TempCurve& tcurve = tempCurveManager.GetForThisCell();

    if ( !curve->AsTable(tcurve.m_dates, tcurve.m_temps) )
      return XL::Oper::ErrValue();

    return tempCurveManager.HandleOf(tcurve);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *TempCurveAt(XL::SheetObjectHandle hcurve, double date)
{
  TempCurveManager& tempCurveManager = GetTheTestAddIn().GetTempCurveManager();
  TempCurve * const curve = tempCurveManager.FromHandle(hcurve);
  if ( !curve )
    return XL::Oper::ErrValue();

  return XL::Oper::Return(curve->GetValueAt(date));
}

XLLEXPORT int TempCurveLen(XL::SheetObjectHandle hcurve)
{
  return (int)(GetTheTestAddIn().GetTempCurveManager()[hcurve].m_dates.size());
}

XLLEXPORT XLOPER *TempInCell()
{
  try
  {
    XL::Oper xloRef;
    XL::Call(xlfCaller, XL::ByRef(xloRef));

    unsigned x, y;
    if ( !xloRef.GetAsSRef(x, y) )
      return XL::Oper::ErrValue();

    return XL::Oper::Return(String::Printf("Cell (%d, %d) is hot", x, y));
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *MyT(const XL::Oper *xlo)
{
  try
  {
    XL::Oper xloArray;
    if ( !xlo->AsArray(&xloArray) )
      return XL::Oper::ErrValue();

    unsigned w, h;
    const XL::Oper *src = xloArray.GetAsArray(w, h);

    const unsigned size = w*h;
    XL::Oper *data = new XL::Oper[size];

    XL::Oper *dst = data;
    for ( unsigned x = 0; x < w; ++x )
    {
      const XL::Oper *col = src + x;
      for ( unsigned y = 0; y < h; ++y, ++dst, col += w )
      {
        *dst = *col;
      }
    }

    return XL::Oper::ReturnArray(h, w, data);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}

XLLEXPORT XLOPER *MyA(int size)
{
  try
  {
    XL::Oper *data = new XL::Oper[size];
    for ( int n = 0; n < size; ++n )
    {
      data[n] = XL::Oper(n + 1);
    }

    return XL::Oper::ReturnArray(1, size, data);
  }
  catch ( ... )
  {
    return XL::ReturnException();
  }
}
// ============================================================================
// TestAddIn implementation
// ============================================================================

TempCurveManager& TestAddIn::GetTempCurveManager()
{
  if ( !m_tempCurves )
    m_tempCurves = new TempCurveManager;

  return *m_tempCurves;
}

void TestAddIn::RegisterAllFunctions()
{
  XL_REGISTER_0(CurrentTempC,
                "Returns current room temperature in degrees Celsius");

  XL_REGISTER_0(TempInCell,
                "Returns temperature in the calling cell");

  XL_REGISTER_1(WhatIs,
                "Returns the string showing the parameter value",
                "param", "Anything at all");

  XL_REGISTER_1(WhatIsDeref,
                "Returns the string showing the dereferenced parameter value",
                "param", "Anything at all");

  XL_REGISTER_1(TempAtDate,
                "Returns temperature at the given date in degrees Celsius",
                "date", "The date for which to return the temperature");

  XL_REGISTER_1(Celsius,
                "Returns the temperature expressed in degrees Celsius",
                "fahrenheit",
                "Temperature in degrees Fahrenheit");

  XL_REGISTER_2(IsHotter,
                "Returns true if celsius is hotter than fahrenheit",
                "celsius",
                "Temperature in degrees Celsius",
                "fahrenheit",
                "Temperature in degrees Fahrenheit");

  XL_REGISTER_1(Hottest,
                "Returns the hottest temperature in the range",
                "temps",
                "Range of temperature values");

  XL_REGISTER_1(DegreesScale,
                "Returns \"Celsius\" or \"Fahrenheit\"",
                "abbrev",
                "Abbreviated scale name, either \"C\" or \"F\"");

  XL_REGISTER_1(TempCurveCreate,
                "Creates a new temperature curve object and returns its handle",
                "curve",
                "Range of cells containing dates and temperatures");
  XL_REGISTER_1(TempCurveLen,
                "Returns the number of points in the given curve",
                "curve",
                "Handle of the curve object");
  XL_REGISTER_2(TempCurveAt,
                "Returns N-th point of the temperature curve",
                "curve",
                "Handle of the curve object",
                "date",
                "Date for which to return the temperature");

  XL_REGISTER_1(MyT, "Transposes a range", "range", "Range to be transposed");
  XL_REGISTER_1(MyA, "Returns an array of given size", "N", "Array size");
}

TestAddIn::~TestAddIn()
{
  delete m_tempCurves;
}

