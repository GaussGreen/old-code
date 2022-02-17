/*************************************************************************** 
*  Name:        ddlmaker/src/tools.cpp                              
*  Purpose:	Part of the ddlmaker utility                                 
*  Author:      Vadim Zeitlin                                            
*  Created:     2005-06-07                                                   
*  RCS-ID:      $Id: tools.cpp,v 1.2 2005/07/06 09:03:02 cosmin Exp $    
*  Copyright:   (c) 2005 Trilemma LLP                                        
***************************************************************************/
#include "tools.h"

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>

#include <string>

// output message to stdout
void Output(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  std::string format(fmt);
  format += '\n';
  vfprintf(stdout, format.c_str(), ap);
  fflush(stdout);
  va_end(ap);
}

// output the message to stderr and terminates the application
void Error(const char *fmt, ...)
{
  va_list ap;
  va_start(ap, fmt);
  std::string format(fmt);
  format += '\n';
  vfprintf(stderr, format.c_str(), ap);
  va_end(ap);

  exit(-1);
}
