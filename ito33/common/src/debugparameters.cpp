/////////////////////////////////////////////////////////////////////////////
// Name:        common/src/debugparameters.cpp
// Purpose:     Struct holding some debug parameters
// Created:     2005/02/01
// RCS-ID:      $Id: debugparameters.cpp,v 1.5 2006/08/19 23:18:54 wang Exp $
// Copyright:   (c) 2005 -  Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

#include "ito33/useexception.h"
#include "ito33/debugparameters.h"
#include "ito33/getenv.h"

extern const ito33::Error ITO33_UNDEFINED_ENV_VAR;

namespace ito33
{

std::string DebugParameters::GetDebugOutputFile(const char* defaultName) const
{
  std::string outputfile;

  // do we need the default directory?
  if ( filename.empty() ||
       filename.find_first_of("/\\") == std::string::npos )
  {
    outputfile = GetDebugOutputDir()
               + ( filename.empty() ? defaultName : filename.c_str() );
  }
  else // we have a full filename (maybe not absolute, but not basename neither)
  {
    outputfile = filename;
  }

  return outputfile;
}

/* static */
std::string DebugParameters::GetDebugOutputDir()
{
  std::string dir(GetEnv("ITO33InstDir"));

  if ( dir.empty() )
    throw EXCEPTION_MSG
          (
            ITO33_UNDEFINED_ENV_VAR,
            TRANS("Evironment variable 'ITO33InstDir' undefined."
                  "It is required for outputting debug files.")
          );

  static const char FILENAME_SEPARATOR =
                    #ifdef _WIN32
                      '\\'
                    #else
                      '/'
                    #endif
                    ;

  if ( *dir.rbegin() != FILENAME_SEPARATOR )
    dir += FILENAME_SEPARATOR;

  dir += std::string("output") + FILENAME_SEPARATOR;

  return dir;
}


/* static */
std::string DebugParameters::GetDebugDir()
{
  std::string dir(GetEnv("ITO33DebugDir"));

  if ( !dir.empty() )
    return dir +
           #ifdef _WIN32
             '\\';
           #else
              '/';
           #endif
  
  return GetDebugOutputDir();
}

} // namespace ito33
