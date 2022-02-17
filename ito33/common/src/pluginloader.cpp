///////////////////////////////////////////////////////////////////////////////
// Name:        src/common/pluginloader.cpp
// Purpose:     Plugins interface
// Author:      Vaclav Slavik
// Created:     2005-01-06
// RCS-ID:      $Id: pluginloader.cpp,v 1.13 2006/05/27 20:13:16 zhang Exp $
// Copyright:   (c) 2004 Trilemma LLP
///////////////////////////////////////////////////////////////////////////////

// ============================================================================
// declarations
// ============================================================================

#define RLOG_COMPONENT "pluginload"

// ----------------------------------------------------------------------------
// headers
// ----------------------------------------------------------------------------

#include "ito33/common.h"
#include "ito33/useexception.h"
#include "ito33/getenv.h"
#include "ito33/log.h"
#include "ito33/pluginloader.h"

#include <memory>

#ifdef _MSC_VER
  #define putenv _putenv
#endif

using namespace std;
using namespace ito33;

// ----------------------------------------------------------------------------
// macros
// ----------------------------------------------------------------------------

ITO33_DEFINE_LOG_CATEGORY(LogPlugin, "plugin");

#define TRACE_PLUGIN    ITO33_TRACE_CATEGORY(LogPlugin)

// ----------------------------------------------------------------------------
// constants
// ----------------------------------------------------------------------------

static const char *PATH_ENV_VAR =
#ifdef _WIN32
      "PATH"
#elif defined(__MACOSX__)
      "DYLD_LIBRARY_PATH"
#else // other Unix
      "LD_LIBRARY_PATH"
#endif
    ;

static const char PATH_SEPARATOR =
#ifdef _WIN32
          ';'
#else // Unix
          ':'
#endif
      ;

// ============================================================================
// implementation
// ============================================================================

// ----------------------------------------------------------------------------
// PluginLoader
// ----------------------------------------------------------------------------

/* static */
string
PluginLoaderBase::GetPluginFileName(const char *prefix, const char *name)
{
  // name of the plugin has plg_<prefix><name>.dll format but we don't add the
  // extension here, it is done by DynamicLibrary
  return String::Printf("plg_%s%s", prefix, name);
}

/* static */
void PluginLoaderBase::AddSearchPath(const char *pathNew)
{
  CHECK_VOID( pathNew && *pathNew != '\0', "invalid path" );

  std::string path(pathNew);
  std::string pathOld = GetEnv(PATH_ENV_VAR);
  if ( !pathOld.empty() )
  {
    path += PATH_SEPARATOR;
    path += pathOld;
  }

  #if ( defined(_MSC_VER) && _MSC_VER >= 1400 )
  {
    _putenv_s(PATH_ENV_VAR, path.c_str());
  }
  #else
  {
    // NB: putenv() makes the string pointed to by its argument part of the
    //     environment, it does _not_ make a copy. This means that we must
    //     ensure that the string will be kept in memory for as long as the
    //     value is assigned to PATH_ENV_VAR (i.e. until the next call to
    //     AddSearchPath() or putenv() changing it), hence the use of static
    //     variable.
    static std::string s_envContent;
    s_envContent = String::Printf("%s=%s", PATH_ENV_VAR, path.c_str());
    putenv((char*)s_envContent.c_str());
  }
  #endif
}

Plugin *PluginLoaderBase::DoLoad(const char *prefix,
                                 const char *name,
                                 const char *category,
                                 int apiver,
                                 DllHandle_t *handle)
{
  extern const ito33::Error ITO33_DLL_ERROR;

  try
  {
    TRACE_PLUGIN("trying to load plugin %s using path \"%s\"",
                 name, GetEnv(PATH_ENV_VAR).c_str());

    const string dllname = GetPluginFileName(prefix, name);

    DynamicLibrary dll(dllname.c_str());

    const string path = dll.GetFullPath();

    TRACE_PLUGIN("plugin %s loaded from \"%s\"", dllname.c_str(), path.c_str());

    const string symbol(String::Printf(PLUGIN_FACTORY_FUNC_NAME_STR(%s),
                        category));

    TRACE_PLUGIN("getting symbol %s from \"%s\"", symbol.c_str(), path.c_str());

    FactoryFuncType factory;

    try
    {
      // NB: there is no C++ cast (not even reinterpret_cast<>!) which can cast
      //     from "void *" returned by GetSymbol() to function pointer
      factory = (FactoryFuncType)dll.GetSymbol(symbol.c_str());
    }
    catch ( ito33::Exception& )
    {
      throw EXCEPTION_MSG
            (
              ITO33_DLL_ERROR,
              String::Printf
              (
                TRANS("DLL \"%s\" is not a valid ITO 33 plugin."),
                path.c_str()
              )
            );
    }

    TRACE_PLUGIN("calling factory %s in \"%s\"", symbol.c_str(), path.c_str());

    std::string error;
    Plugin *obj = (*factory)(&error, apiver, NULL /*unused*/);

    if ( !obj )
    {
      throw EXCEPTION_MSG
            (
              ITO33_DLL_ERROR,
              String::Printf
              (
                TRANS("Error loading plugin %s from \"%s\": %s."),
                name,
                path.c_str(),
                error.c_str()
              )
            );
    }

    TRACE_PLUGIN("successfully loaded %s from %s", category, path.c_str());

    DllHandle_t h = dll.Detach();

    if ( handle )
      *handle = h;
    //else: unload the DLL on program shutdown, not sooner - we can't unload
    //      it from Plugin dtor, which is the only place where we know it's
    //      no longer needed, because execution would return to now-unloaded
    //      code in Plugin dtor. We can't put it to global variable with list
    //      of DLLs to unload on program shutdown either because
    //      datastore::plugin::Reader is used as global object also destroyed
    //      at program shutdown and we have no way of knowing which of the
    //      two global objects would be cleared first (it would crash if the
    //      list was cleared first). The OS unloads DLLs when the process
    //      exits anyway, so leave it to the OS:

    return obj;
  }
  catch ( ito33::Exception& e )
  {
    if ( e.GetErrorCode() == ITO33_DLL_ERROR )
      throw;

    throw EXCEPTION_MSG
          (
            ITO33_DLL_ERROR,
            String::Printf(TRANS("Failed to load plugin \"%s\"."), name)
          );
  }
}
