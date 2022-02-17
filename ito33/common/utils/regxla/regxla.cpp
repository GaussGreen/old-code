/////////////////////////////////////////////////////////////////////////////
// Name:        regxla/regxla.cpp
// Purpose:     register an Excel add-in (XLA file)
// Author:      Vadim Zeitlin
// Created:     2004-12-30
// RCS-ID:      $Id: regxla.cpp,v 1.7 2006/07/20 23:09:02 zeitlin Exp $
// Copyright:   (c) 2004 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/*
    This tool installs or uninstalls an Excel plug-in.

    Syntax: regxla [/s] [/u] <filename.xla>
    /s      silent, don't give success messages
    /u      uninstall

    The specified plugin is [un]installed for all Excel versions starting from
    Excel 2002 (a.k.a. Excel XP, a.k.a. Excel 10) installed on the system.

    TODO:
      - allow selecting Excel version
 */

// ============================================================================
// declarations
// ============================================================================

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/string.h"
#include "ito33/exception.h"
#include "ito33/vector.h"

#include "ito33/win32/regkey.h"

#ifdef _MSC_VER
  // these warning have to be disabled and not just temporarily disabled
  // because they will be given at the end of the compilation of the current
  // source -- and there is absolutely nothing we can do about them

  // 'foo': copy constructor could not be generated
  #pragma warning(disable:4511)

  // 'foo': assignment operator could not be generated
  #pragma warning(disable:4512)

  // conditional expression is constant
  #pragma warning(disable:4127)

  // local variable 'foo' may be used without having been initialized
  #pragma warning(disable:4701)
#endif

#include <boost/program_options.hpp>

using namespace ito33;
using ito33::Win32::RegKey;

// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------

static const char *EXCEL_KEY = "Excel";
static const char *EXCEL_OPTIONS_SUBKEY = "Options";
static const char *EXCEL_ADDIN_MANAGER_SUBKEY = "Add-in Manager";

// ----------------------------------------------------------------------------
// globals
// ----------------------------------------------------------------------------

static bool gs_silent = false;
static bool gs_uninstall = false;

static std::string gs_filenameXLA;

// ----------------------------------------------------------------------------
// ExcelVersionIterator: iterate over all installed Excel versions
// ----------------------------------------------------------------------------

class ExcelVersionIterator
{
public:
  // initialize the iterator, use IsOfficeInstalled() to test if any Office
  // versions are available
  ExcelVersionIterator();

  // call ProcessExcelVersion() with the root registry key of this version as
  // parameter
  void ProcessAll();

  // return true if any Office seems to be present
  bool IsOfficeInstalled() const { return m_keyOffice != NULL; }

  // return true if we have actually found an installed Excel version in the
  // registry
  bool FoundExcel() const { return m_foundExcel; }

  // return true if everything succeeded and false if at least one call of
  // ProcessExcelVersion() has failed (should be only used if FoundExcel())
  bool AllSucceeded() const { return m_allOk; }

  // close any open registry keys
  virtual ~ExcelVersionIterator();

protected:
  // this function should perform the requested action for the given Excel
  // version and return true on success, false on failure
  //
  // name is the user-recognizable name of this Excel version (2002, 2003, ...)
  virtual bool ProcessExcelVersion(RegKey& keyVer, const char *name) = 0;

private:
  // the top level Microsoft Office key, created in the ctor if Office is
  // installed, is NULL otherwise
  RegKey *m_keyOffice;

  // flag which is set to true if we have found any Excel versions
  bool m_foundExcel;

  // flag set to false if at least one call to ProcessExcelVersion() failed
  bool m_allOk;


  NO_COPY_CLASS(ExcelVersionIterator);
};

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// message output functions
// ----------------------------------------------------------------------------

static void
Output(DWORD mbFlags, const char *title, const char *format, va_list argptr)
{
  ::MessageBoxA(NULL, String::PrintfV(format, argptr).c_str(), title, mbFlags);
}

static void ErrorMsg(const char *format, ...)
{
  va_list argptr;
  va_start(argptr, format);

  Output(MB_ICONERROR, "Excel Add-In Registration Error", format, argptr);

  va_end(argptr);
}

static void Warning(const char *format, ...)
{
  va_list argptr;
  va_start(argptr, format);

  Output(MB_ICONWARNING, "Excel Add-In Registration Warning", format, argptr);

  va_end(argptr);
}

static void Message(const char *format, ...)
{
  if ( gs_silent )
    return;

  va_list argptr;
  va_start(argptr, format);

  Output(MB_ICONINFORMATION, "Excel Add-In Registration", format, argptr);

  va_end(argptr);
}

static void Debug(const char *format, ...)
{
  va_list argptr;
  va_start(argptr, format);

  ::OutputDebugStringA((String::PrintfV(format, argptr) + "\r\n").c_str());

  va_end(argptr);
}

// ----------------------------------------------------------------------------
// cmd line parsing
// ----------------------------------------------------------------------------

bool ParseCommandLine(const char *lpCmdLine)
{
  namespace popts = boost::program_options;

  popts::options_description options;
  options.add_options()
    ( "silent,s",     "silent, don't give success messages" )
    ( "uninstall,u",  "uninstall Excel add-in instead of installing it" )
    ( "xla",          "" );

  popts::positional_options_description params;
  params.add("xla", 1);

  std::string msg;
  std::vector<std::string> args = popts::split_winmain(lpCmdLine);
  popts::variables_map variables;
  try
  {
    popts::store(
        popts::command_line_parser(args)
        .options(options)
        .positional(params)
        .style(
          popts::command_line_style::unix_style |
          popts::command_line_style::allow_slash_for_short)
        .run(), variables);
    popts::notify(variables);
  }
  catch ( std::exception& e )
  {
    msg = "Wrong command line syntax: ";
    msg += e.what();
  }

  if ( msg.empty() && !variables.count("xla") )
  {
    msg = "Missing XLA file name.";
  }

  if ( !msg.empty() )
  {
    ErrorMsg("%s\n"
          "\n"
          "Usage: regxla [/s] [/u] <filename.xla>\n"
          "/s\tsilent, don't give success message\n"
          "/u\tuninstall the add-in instead of installing it", msg.c_str());

    return false;
  }

  gs_silent = variables.count("silent") != 0;
  gs_uninstall = variables.count("uninstall") != 0;

  gs_filenameXLA = variables["xla"].as<std::string>();

  return true;
}

// ----------------------------------------------------------------------------
// ExcelVersionIterator implementation
// ----------------------------------------------------------------------------

ExcelVersionIterator::ExcelVersionIterator()
{
  m_foundExcel = false;
  m_allOk = true;
  m_keyOffice = NULL;

  // check if we have Microsoft Office at all
  static const char *OFFICE_KEY = "Software\\Microsoft\\Office";
  if ( RegKey::Exists(HKEY_CURRENT_USER, OFFICE_KEY) )
  {
    m_keyOffice = new RegKey(HKEY_CURRENT_USER, OFFICE_KEY);
  }
  //else: leave m_keyOffice unopened
}

void ExcelVersionIterator::ProcessAll()
{
  if ( !m_keyOffice )
    return;

  // we could iterate over all subkeys of keyOffice and just install the plugin
  // for all versions but as we currently support only Excel 2002/2003/2007, it
  // is easier to check just for them -- but this will need to be updated for
  // any new Excel version
  static const char *versions[] = { "10.0", "11.0", "12.0" };
  static const char *names[] = { "2002", "2003", "2007" };
  for ( size_t n = 0; n < SIZEOF(versions); ++n )
  {
    if ( RegKey::Exists(*m_keyOffice, versions[n]) )
    {
      m_foundExcel = true;

      RegKey regVer(*m_keyOffice, versions[n]);

      m_allOk &= ProcessExcelVersion(regVer, names[n]);
    }
  }
}

ExcelVersionIterator::~ExcelVersionIterator()
{
  delete m_keyOffice;
}

// ----------------------------------------------------------------------------
// miscellaneous helpers
// ----------------------------------------------------------------------------

// quote the string if it contains any embedded spaces or other non alnum
// characters
static std::string GetPossiblyQuotedPath(const std::string& unquoted)
{
  std::string quoted;
  if ( unquoted.find(' ') == std::string::npos )
  {
    quoted = unquoted;
  }
  else // do quote it
  {
    quoted = '"' + unquoted + '"';
  }

  return quoted;
}

// ----------------------------------------------------------------------------
// installation
// ----------------------------------------------------------------------------

static std::string GetExcelOpenValue(int n)
{
  // Excel numbers the "open" commands in a typically logical way: first one is
  // "open", second one is "open1", then "open2" and so on
  std::string open = "open";
  if ( n )
    open += String::Printf("%d", n);

  return open;
}

// install the XLA for a single Excel version
//
// the keyExcel is supposed to be HKCU\Software\Microsoft\Office\<version>\Excel
static bool InstallXLAForVersion(RegKey& keyExcel, const char *name)
{
  // if the path contains any quotes it should be quoted as otherwise Excel
  // won't be able to open it
  std::string filenameXLA(GetPossiblyQuotedPath(gs_filenameXLA));

  RegKey keyOptions(keyExcel, EXCEL_OPTIONS_SUBKEY);

  // find the first "free" open value checking that we're not already installed
  // in the meanwhile
  std::string open;
  for ( int n = 0; keyOptions.HasValue(open = GetExcelOpenValue(n)); ++n )
  {
    if ( keyOptions.GetValue(open) == filenameXLA )
    {
      Warning("Add-in \"%s\" is already installed for Excel %s.",
              gs_filenameXLA.c_str(),
              name);

      // don't install it again
      return true;
    }
  }

  // add the new key
  keyOptions.Set(open.c_str(), filenameXLA);


  // also add it to the list of plug-ins which should be shown in XL add-in
  // manager dialog box
  RegKey keyAddinManager(keyExcel, EXCEL_ADDIN_MANAGER_SUBKEY);
  if ( !keyAddinManager.HasValue(filenameXLA) )
    keyAddinManager.Set(filenameXLA.c_str(), "");

  return true;
}

static bool InstallXLA()
{
  class Installer : public ExcelVersionIterator
  {
  public:
    Installer() { m_installed = false; }

    // return true if we have managed to install for at least some Excel
    // version
    bool WasInstalled() const { return m_installed; }

  protected:
    virtual bool ProcessExcelVersion(RegKey& keyVer, const char *name)
    {
      RegKey keyExcel(keyVer, EXCEL_KEY);
      if ( !InstallXLAForVersion(keyExcel, name) )
        return false;

      Message("Successfully installed \"%s\" add-in for Microsoft Excel %s.",
              gs_filenameXLA.c_str(),
              name);

      m_installed = true;

      return true;
    }

  private:
    bool m_installed;
  };

  // check if we have Office at all
  Installer installer;
  if ( !installer.IsOfficeInstalled() )
  {
    ErrorMsg("No Microsoft Office versions found on this machine, please "
          "install Microsoft Office before installing the add-ins.");

    return false;
  }

  installer.ProcessAll();

  if ( !installer.FoundExcel() )
  {
    ErrorMsg("Microsoft Excel version is not supported, please "
             "install Microsoft Excel 2002, 2003 or 2007 before "
             "installing the add-ins.");
    return false;
  }

  if ( !installer.WasInstalled() )
  {
    // error message already given
    return false;
  }

  return true;
}

// ----------------------------------------------------------------------------
// uninstallation
// ----------------------------------------------------------------------------

static bool UninstallXLAForVersion(RegKey& keyExcel, const char *name)
{
  RegKey keyOptions(keyExcel, EXCEL_OPTIONS_SUBKEY);

  bool found = false;

  const std::string filenameXLA = GetPossiblyQuotedPath(gs_filenameXLA);

  // check the values of all "open<x>" values
  std::string open;
  for ( int n = 0; keyOptions.HasValue(open = GetExcelOpenValue(n)); )
  {
    if ( keyOptions.GetValue(open) == filenameXLA )
    {
      // check that it's the first time we find this add-in
      if ( found )
      {
        Warning("Multiple installations of Excel add-in \"%s\" found "
                "for Excel %s",
                gs_filenameXLA.c_str(),
                name);
      }
      else
      {
        found = true;
      }

      // remove this one
      keyOptions.DeleteValue(open);

      // and rename all the remaining keys so that they remain sequential
      std::string open2;
      for ( int n2 = n + 1;
            keyOptions.HasValue(open2 = GetExcelOpenValue(n2));
            ++n2 )
      {
        keyOptions.Set(GetExcelOpenValue(n2 - 1).c_str(),
                            keyOptions.GetValue(open2));
        keyOptions.DeleteValue(open2);
      }

      // remain on the same value as it has now changed and should be examined
      // again
    }
    else // go to the next value
    {
      n++;
    }
  }

  // also remove this add-in from the list of all add-ins (active or not)
  RegKey keyAddinManager(keyExcel, EXCEL_ADDIN_MANAGER_SUBKEY);
  if ( keyAddinManager.HasValue(filenameXLA) )
  {
    keyAddinManager.DeleteValue(filenameXLA);
  }
  else // not fatal but still log this
  {
    Debug("Add-in \"%s\" not found under %s\\%s",
          filenameXLA.c_str(), EXCEL_KEY, EXCEL_ADDIN_MANAGER_SUBKEY);
  }

  return found;
}

static bool UninstallXLA()
{
  class Uninstaller : public ExcelVersionIterator
  {
  public:
    Uninstaller() { }

  protected:
    virtual bool ProcessExcelVersion(RegKey& keyVer, const char *name)
    {
      RegKey keyExcel(keyVer, EXCEL_KEY);
      if ( !UninstallXLAForVersion(keyExcel, name) )
        return false;

      Message("Successfully uninstalled \"%s\" add-in for Microsoft Excel %s.",
              gs_filenameXLA.c_str(),
              name);

      return true;
    }
  };

  // check if we have Office at all
  Uninstaller uninstaller;
  if ( uninstaller.IsOfficeInstalled() )
  {
    uninstaller.ProcessAll();

    if ( uninstaller.FoundExcel() )
      return uninstaller.AllSucceeded();
  }

  Warning("Excel add-in \"%s\" hasn't been found, not uninstalled.",
          gs_filenameXLA.c_str());

  return true;
}

// ----------------------------------------------------------------------------
// entry point
// ----------------------------------------------------------------------------

int WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR lpCmdLine, int)
{
  try
  {
    if ( !ParseCommandLine(lpCmdLine) )
      return 1;

    if ( !(gs_uninstall ? UninstallXLA() : InstallXLA()) )
      return 2;

    return 0;
  }
  catch ( ito33::Exception& e )
  {
    ErrorMsg("Unexpected exception has occured: %s", e.GetFullMessage().c_str());
  }
  catch ( std::exception& e )
  {
    ErrorMsg("Unexpected exception has occured: %s", e.what());
  }
  catch ( ... )
  {
    ErrorMsg("Unexpected exception has occured.");
  }

  return -1;
}

