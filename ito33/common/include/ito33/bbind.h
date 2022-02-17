/////////////////////////////////////////////////////////////////////////////
// Name:        ito33/bbind.h
// Purpose:     wrapper for boost::bind headers
// Author:      Vaclav Slavik
// Created:     2006-03-13
// RCS-ID:      $Id: bbind.h,v 1.1 2006/06/28 12:34:36 vaclav Exp $
// Copyright:   (c) 2006 Trilemma LLP
/////////////////////////////////////////////////////////////////////////////

/**
    @file ito33/bbind.h
    @brief Wrapper around boost::bind headers.

    Including &lt;boost/bind.hpp&gt; directly results in in _1, _2, ..., _9
    symbols being defined in anonymous namespace. This conflicts with
    ito33/bll.h, which contains "using boost::lambda::_1" etc. statements.

    This header includes &lt;boost/bind.hpp&lt;, but disables these symbols
    and defines them in the boost:: namespace instead. It also adds similar
    using statements to these in ito33/bll.h ("using boost::_1;"), which makes
    it possible to resolve conflicts between these two headers with "using".
 */

#ifndef _ITO33_BBIND_H_
#define _ITO33_BBIND_H_

// include boost/bind.hpp, but disable declarations of _1, _2 etc. placeholders
// in global namespace: we'll do them ourselves below
#define BOOST_BIND_NO_PLACEHOLDERS
#include <boost/bind.hpp>

// NB: See ito33/bll.h for explanation of why we can't put the _1, _2 etc.
//     placeholders directly into the global namespace.
//
//     We put them into anonymous namespace (just like boost does) in order
//     to avoid linker errors about multiple definitions of the symbols. We
//     use boost_bind_placeholders namespace for reasons explained in bll.h.
namespace
{

namespace boost_bind_placeholders
{

#ifdef BOOST_BIND_PLACEHOLDERS_HPP_INCLUDED
  #error "<boost/bind/placeholders.hpp> already included, include ito33/bbind.h sooner"
#endif

// this code is copied from <boost/bind/placeholders.hpp> (and slightly
// simplified for our use by removing code for N>3 and compilers we don't use):

#if defined(BOOST_MSVC)

static boost::arg<1> _1;
static boost::arg<2> _2;
static boost::arg<3> _3;

#else

boost::arg<1> _1;
boost::arg<2> _2;
boost::arg<3> _3;

#endif

using boost::bind;

} // namespace boost_bind_placeholders

} // anonymous namespace

// this makes the placeholders available in global namespace:
using namespace boost_bind_placeholders;

// and this makes them available as boost::_1 etc. so that conflicts with
// ito33/bll.h can be resolved by typing "using boost::_1;":
namespace boost
{
  using namespace boost_bind_placeholders;
}

#endif // _ITO33_BOOSTBIND_H_
