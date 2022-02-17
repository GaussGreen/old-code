/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/bll.h
// Purpose:     wrapper for boost::lambda headers
// Author:      Vadim Zeitlin
// Created:     2005-10-23
// RCS-ID:      $Id: bll.h,v 1.2 2006/06/28 12:34:36 vaclav Exp $
// Copyright:   (c) 2005 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/bll.h
    @brief Wrapper around boost::lambda library (BLL) headers.

    Including &lt;boost/lambda/lambda.hpp&gt; directly results in many
    warnings under VC++ so this header disables them. It also brings the
    placeholders and bind() into the including namespace: although contrary to
    the usual practice, this has to be done to use boost::lambda in practice.
 */

#ifndef _ITO33_BLL_H_
#define _ITO33_BLL_H_

#ifdef _MSC_VER
  #pragma warning(push)
  // 'class' : assignment operator could not be generated
  #pragma warning(disable:4512)
  // 'class' : destructor could not be generated
  #pragma warning(disable:4513)
#endif

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#ifdef _MSC_VER
  #pragma warning(pop)
#endif

// NB: We can't simply have "using boost::lambda::_1" in global namespace,
//     because it would conflict with boost::bind as soon as you included
//     both boost/bind.hpp and ito33/bll.h. But if we put the placeholders into 
//     a namespace instead and use "using namespace", the conflict wouldn't
//     happen until some code actually uses them -- and that can be resolved
//     by writing e.g. "using boost::lambda::_1;" above the conflicting code.
namespace Private
{

namespace bll_placeholders
{

using boost::lambda::_1;
using boost::lambda::_2;
using boost::lambda::_3;
using boost::lambda::bind;

} // namespace bll_placeholders

} // namespace Private

using namespace Private::bll_placeholders;

#endif // _ITO33_BLL_H_
