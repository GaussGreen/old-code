/************************************************************************
 * Module:	DRL - drl_getopt
 * Function:	A drl_getopt() routine for portability
 * 		(usually available in UNIX only)
 * Author:	GNU (modified by C. Daher)
 * Revision:	$Header$
 ************************************************************************/
#include "drlstd.h"		/* platform compatibility */

#include "drlgetop.h"		/* */

/*#if defined(_INCLUDE_DRL_GETOPT)*/

#define	STDC_HEADERS

/* Getopt for GNU.
   NOTE: drl_getopt is now part of the C library, so if you don't know what
   "Keep this file name-space clean" means, talk to roland@gnu.ai.mit.edu
   before changing it!

   Copyright (C) 1987, 88, 89, 90, 91, 1992 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.  */


/* AIX requires this to be the first thing in the file.  */
#ifdef __GNUC__
#define alloca __builtin_alloca
#else /* not __GNUC__ */
#if defined (HAVE_ALLOCA_H) || (defined(sparc) && (defined(sun) || (!defined(USG) && !defined(SVR4) && !defined(__svr4__))))
#include <alloca.h>
#else
#ifdef _AIX
#pragma alloca
#else
char *alloca ();
#endif
#endif /* alloca.h */
#endif /* not __GNUC__ */




/* This tells Alpha OSF/1 not to define a drl_getopt prototype in <stdio.h>.  */
#ifndef _NO_PROTO
#  define _NO_PROTO
#endif

#include <stdio.h>

#if defined(USG) || defined(STDC_HEADERS) || defined(__GNU_LIBRARY__)
#include <string.h>
#endif

/* This needs to come after some library #include
   to get __GNU_LIBRARY__ defined.  */
#ifdef	__GNU_LIBRARY__
#undef	alloca
/* Don't include stdlib.h for non-GNU C libraries because some of them
   contain conflicting prototypes for drl_getopt.  */
#include <stdlib.h>
#else	/* Not GNU C library.  */

# ifdef	UNIX
#  define	__alloca	alloca
# else
#  define	__alloca	malloc
#endif
#endif	/* GNU C library.  */

/*$$$ They can't write something portable ! */
#define	alloca	malloc

#if !__STDC__
#define const
#endif

/* If GETOPT_COMPAT is defined, `+' as well as `--' can introduce a
   long-named option.  Because this is not POSIX.2 compliant, it is
   being phased out.  */
#define GETOPT_COMPAT

/* This version of `drl_getopt' appears to the caller like standard Unix `drl_getopt'
   but it behaves differently for the user, since it allows the user
   to intersperse the options with the other arguments.

   As `drl_getopt' works, it permutes the elements of ARGV so that,
   when it is done, all the options precede everything else.  Thus
   all application programs are extended to handle flexible argument order.

   Setting the environment variable POSIXLY_CORRECT disables permutation.
   Then the behavior is completely standard.

   GNU application programs can use a third alternative mode in which
   they can distinguish the relative order of options and other arguments.  */

/*----------------------------------------------------------------------
 * drl_getopt.h: from the GNU library
 */


/* Declarations for drl_getopt.
   Copyright (C) 1989, 1990, 1991, 1992 Free Software Foundation, Inc.

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.  */

#ifndef _GETOPT_H
#define _GETOPT_H 1

#ifdef	__cplusplus
extern "C" {
#endif

/* For communication from `drl_getopt' to the caller.
   When `drl_getopt' finds an option that takes an argument,
   the argument value is returned here.
   Also, when `ordering' is RETURN_IN_ORDER,
   each non-option ARGV-element is returned here.  */

extern char *drl_optarg;

/* Index in ARGV of the next element to be scanned.
   This is used for communication to and from the caller
   and for communication between successive calls to `drl_getopt'.

   On entry to `drl_getopt', zero means this is the first call; initialize.

   When `drl_getopt' returns EOF, this is the index of the first of the
   non-option elements that the caller should itself scan.

   Otherwise, `drl_optind' communicates from one call to the next
   how much of ARGV has been scanned so far.  */

extern int drl_optind;

/* Callers store zero here to inhibit the error message `drl_getopt' prints
   for unrecognized options.  */

extern int drl_opterr;

/* Describe the long-named options requested by the application.
   The LONG_OPTIONS argument to drl_getopt_long or drl_getopt_long_only is a vector
   of `struct option' terminated by an element containing a name which is
   zero.

   The field `has_arg' is:
   no_argument		(or 0) if the option does not take an argument,
   required_argument	(or 1) if the option requires an argument,
   optional_argument 	(or 2) if the option takes an optional argument.

   If the field `flag' is not NULL, it points to a variable that is set
   to the value given in the field `val' when the option is found, but
   left unchanged if the option is not found.

   To have a long-named option do something other than set an `int' to
   a compiled-in constant, such as set a value from `drl_optarg', set the
   option's `flag' field to zero and its `val' field to a nonzero
   value (the equivalent single-letter option character, if there is
   one).  For long options that have a zero `flag' field, `drl_getopt'
   returns the contents of the `val' field.  */

struct option
{
#if	__STDC__
  const char *name;
#else
  char *name;
#endif
  /* has_arg can't be an enum because some compilers complain about
     type mismatches in all the code that assumes it is an int.  */
  int has_arg;
  int *flag;
  int val;
};

/* Names for the values of the `has_arg' field of `struct option'.  */

enum _argtype
{
  no_argument,
  required_argument,
  optional_argument
};

#if __STDC__ || defined(PROTO)
#if defined(__GNU_LIBRARY__)
/* Many other libraries have conflicting prototypes for drl_getopt, with
   differences in the consts, in stdlib.h.  To avoid compilation
   errors, only prototype drl_getopt for the GNU C library.  */
extern int drl_getopt (int argc, char *const *argv, const char *shortopts);
#endif /* not __GNU_LIBRARY__ */
extern int drl_getopt_long (int argc, char *const *argv, const char *shortopts,
		        const struct option *longopts, int *longind);
extern int drl_getopt_long_only (int argc, char *const *argv,
			     const char *shortopts,
		             const struct option *longopts, int *longind);

/* Internal only.  Users should not call this directly.  */
extern int _drl_getopt_internal (int argc, char *const *argv,
			     const char *shortopts,
		             const struct option *longopts, int *longind,
			     int long_only);
#else /* not __STDC__ */
extern int drl_getopt ();
extern int drl_getopt_long ();
extern int drl_getopt_long_only ();

extern int _drl_getopt_internal ();
#endif /* not __STDC__ */

#ifdef	__cplusplus
}
#endif

#endif /* _GETOPT_H */

/*----------------------------------------------------------------------
 * Back to drl_getopt.c
 */

/* For communication from `drl_getopt' to the caller.
   When `drl_getopt' finds an option that takes an argument,
   the argument value is returned here.
   Also, when `ordering' is RETURN_IN_ORDER,
   each non-option ARGV-element is returned here.  */

char *drl_optarg = 0;

/* Index in ARGV of the next element to be scanned.
   This is used for communication to and from the caller
   and for communication between successive calls to `drl_getopt'.

   On entry to `drl_getopt', zero means this is the first call; initialize.

   When `drl_getopt' returns EOF, this is the index of the first of the
   non-option elements that the caller should itself scan.

   Otherwise, `drl_optind' communicates from one call to the next
   how much of ARGV has been scanned so far.  */

int drl_optind = 0;

/* The next char to be scanned in the option-element
   in which the last option character we returned was found.
   This allows us to pick up the scan where we left off.

   If this is zero, or a null string, it means resume the scan
   by advancing to the next ARGV-element.  */

static char *nextchar;

/* Callers store zero here to inhibit the error message
   for unrecognized options.  */

int drl_opterr = 1;

/* Describe how to deal with options that follow non-option ARGV-elements.

   If the caller did not specify anything,
   the default is REQUIRE_ORDER if the environment variable
   POSIXLY_CORRECT is defined, PERMUTE otherwise.

   REQUIRE_ORDER means don't recognize them as options;
   stop option processing when the first non-option is seen.
   This is what Unix does.
   This mode of operation is selected by either setting the environment
   variable POSIXLY_CORRECT, or using `+' as the first character
   of the list of option characters.

   PERMUTE is the default.  We permute the contents of ARGV as we scan,
   so that eventually all the non-options are at the end.  This allows options
   to be given in any order, even with programs that were not written to
   expect this.

   RETURN_IN_ORDER is an option available to programs that were written
   to expect options and other ARGV-elements in any order and that care about
   the ordering of the two.  We describe each non-option ARGV-element
   as if it were the argument of an option with character code 1.
   Using `-' as the first character of the list of option characters
   selects this mode of operation.

   The special argument `--' forces an end of option-scanning regardless
   of the value of `ordering'.  In the case of RETURN_IN_ORDER, only
   `--' can cause `drl_getopt' to return EOF with `drl_optind' != ARGC.  */

static enum
{
  REQUIRE_ORDER, PERMUTE, RETURN_IN_ORDER
} ordering;

#ifdef	__GNU_LIBRARY__
#include <string.h>
#define	my_index	strchr
#define	my_bcopy(src, dst, n)	memcpy ((dst), (src), (n))
#else

/* Avoid depending on library functions or files
   whose names are inconsistent.  */

char *getenv ();

static char *
my_index (string, chr)
     char *string;
     int chr;
{
  while (*string)
    {
      if (*string == chr)
	return string;
      string++;
    }
  return 0;
}

static void
my_bcopy (from, to, size)
     char *from, *to;
     int size;
{
  int i;
  for (i = 0; i < size; i++)
    to[i] = from[i];
}
#endif				/* GNU C library.  */

/* Handle permutation of arguments.  */

/* Describe the part of ARGV that contains non-options that have
   been skipped.  `first_nonopt' is the index in ARGV of the first of them;
   `last_nonopt' is the index after the last of them.  */

static int first_nonopt;
static int last_nonopt;

/* Exchange two adjacent subsequences of ARGV.
   One subsequence is elements [first_nonopt,last_nonopt)
   which contains all the non-options that have been skipped so far.
   The other is elements [last_nonopt,drl_optind), which contains all
   the options processed since those non-options were skipped.

   `first_nonopt' and `last_nonopt' are relocated so that they describe
   the new indices of the non-options in ARGV after they are moved.  */

static void
exchange (argv)
     char **argv;
{
  int nonopts_size = (last_nonopt - first_nonopt) * sizeof (char *);
#ifdef _CRAY
  char *temp[last_nonopt - first_nonopt];
#else
  char **temp = (char **) __alloca (nonopts_size);
#endif

  /* Interchange the two blocks of data in ARGV.  */

  my_bcopy ((char *) &argv[first_nonopt], (char *) temp, nonopts_size);
  my_bcopy ((char *) &argv[last_nonopt], (char *) &argv[first_nonopt],
	    (drl_optind - last_nonopt) * sizeof (char *));
  my_bcopy ((char *) temp,
	    (char *) &argv[first_nonopt + drl_optind - last_nonopt],
	    nonopts_size);

  /* Update records for the slots the non-options now occupy.  */

  first_nonopt += (drl_optind - last_nonopt);
  last_nonopt = drl_optind;
}

/* Scan elements of ARGV (whose length is ARGC) for option characters
   given in OPTSTRING.

   If an element of ARGV starts with '-', and is not exactly "-" or "--",
   then it is an option element.  The characters of this element
   (aside from the initial '-') are option characters.  If `drl_getopt'
   is called repeatedly, it returns successively each of the option characters
   from each of the option elements.

   If `drl_getopt' finds another option character, it returns that character,
   updating `drl_optind' and `nextchar' so that the next call to `drl_getopt' can
   resume the scan with the following option character or ARGV-element.

   If there are no more option characters, `drl_getopt' returns `EOF'.
   Then `drl_optind' is the index in ARGV of the first ARGV-element
   that is not an option.  (The ARGV-elements have been permuted
   so that those that are not options now come last.)

   OPTSTRING is a string containing the legitimate option characters.
   If an option character is seen that is not listed in OPTSTRING,
   return '?' after printing an error message.  If you set `drl_opterr' to
   zero, the error message is suppressed but we still return '?'.

   If a char in OPTSTRING is followed by a colon, that means it wants an arg,
   so the following text in the same ARGV-element, or the text of the following
   ARGV-element, is returned in `drl_optarg'.  Two colons mean an option that
   wants an optional arg; if there is text in the current ARGV-element,
   it is returned in `drl_optarg', otherwise `drl_optarg' is set to zero.

   If OPTSTRING starts with `-' or `+', it requests different methods of
   handling the non-option ARGV-elements.
   See the comments about RETURN_IN_ORDER and REQUIRE_ORDER, above.

   Long-named options begin with `--' instead of `-'.
   Their names may be abbreviated as long as the abbreviation is unique
   or is an exact match for some defined option.  If they have an
   argument, it follows the option name in the same ARGV-element, separated
   from the option name by a `=', or else the in next ARGV-element.
   When `drl_getopt' finds a long-named option, it returns 0 if that option's
   `flag' field is nonzero, the value of the option's `val' field
   if the `flag' field is zero.

   The elements of ARGV aren't really const, because we permute them.
   But we pretend they're const in the prototype to be compatible
   with other systems.

   LONGOPTS is a vector of `struct option' terminated by an
   element containing a name which is zero.

   LONGIND returns the index in LONGOPT of the long-named option found.
   It is only valid when a long-named option has been found by the most
   recent call.

   If LONG_ONLY is nonzero, '-' as well as '--' can introduce
   long-named options.  */

int
_drl_getopt_internal (argc, argv, optstring, longopts, longind, long_only)
     int argc;
     char *const *argv;
     const char *optstring;
     const struct option *longopts;
     int *longind;
     int long_only;
{
  int option_index;

  drl_optarg = 0;

  /* Initialize the internal data when the first call is made.
     Start processing options with ARGV-element 1 (since ARGV-element 0
     is the program name); the sequence of previously skipped
     non-option ARGV-elements is empty.  */

  if (drl_optind == 0)
    {
      first_nonopt = last_nonopt = drl_optind = 1;

      nextchar = NULL;

      /* Determine how to handle the ordering of options and nonoptions.  */

      if (optstring[0] == '-')
	{
	  ordering = RETURN_IN_ORDER;
	  ++optstring;
	}
      else if (optstring[0] == '+')
	{
	  ordering = REQUIRE_ORDER;
	  ++optstring;
	}
      else if (getenv ("POSIXLY_CORRECT") != NULL)
	ordering = REQUIRE_ORDER;
      else
	ordering = PERMUTE;
    }

  if (nextchar == NULL || *nextchar == '\0')
    {
      if (ordering == PERMUTE)
	{
	  /* If we have just processed some options following some non-options,
	     exchange them so that the options come first.  */

	  if (first_nonopt != last_nonopt && last_nonopt != drl_optind)
	    exchange ((char **) argv);
	  else if (last_nonopt != drl_optind)
	    first_nonopt = drl_optind;

	  /* Now skip any additional non-options
	     and extend the range of non-options previously skipped.  */

	  while (drl_optind < argc
		 && (argv[drl_optind][0] != '-' || argv[drl_optind][1] == '\0')
#ifdef GETOPT_COMPAT
		 && (longopts == NULL
		     || argv[drl_optind][0] != '+' || argv[drl_optind][1] == '\0')
#endif				/* GETOPT_COMPAT */
		 )
	    drl_optind++;
	  last_nonopt = drl_optind;
	}

      /* Special ARGV-element `--' means premature end of options.
	 Skip it like a null option,
	 then exchange with previous non-options as if it were an option,
	 then skip everything else like a non-option.  */

      if (drl_optind != argc && !strcmp (argv[drl_optind], "--"))
	{
	  drl_optind++;

	  if (first_nonopt != last_nonopt && last_nonopt != drl_optind)
	    exchange ((char **) argv);
	  else if (first_nonopt == last_nonopt)
	    first_nonopt = drl_optind;
	  last_nonopt = argc;

	  drl_optind = argc;
	}

      /* If we have done all the ARGV-elements, stop the scan
	 and back over any non-options that we skipped and permuted.  */

      if (drl_optind == argc)
	{
	  /* Set the next-arg-index to point at the non-options
	     that we previously skipped, so the caller will digest them.  */
	  if (first_nonopt != last_nonopt)
	    drl_optind = first_nonopt;
	  return EOF;
	}

      /* If we have come to a non-option and did not permute it,
	 either stop the scan or describe it to the caller and pass it by.  */

      if ((argv[drl_optind][0] != '-' || argv[drl_optind][1] == '\0')
#ifdef GETOPT_COMPAT
	  && (longopts == NULL
	      || argv[drl_optind][0] != '+' || argv[drl_optind][1] == '\0')
#endif				/* GETOPT_COMPAT */
	  )
	{
	  if (ordering == REQUIRE_ORDER)
	    return EOF;
	  drl_optarg = argv[drl_optind++];
	  return 1;
	}

      /* We have found another option-ARGV-element.
	 Start decoding its characters.  */

      nextchar = (argv[drl_optind] + 1
		  + (longopts != NULL && argv[drl_optind][1] == '-'));
    }

  if (longopts != NULL
      && ((argv[drl_optind][0] == '-'
	   && (argv[drl_optind][1] == '-' || long_only))
#ifdef GETOPT_COMPAT
	  || argv[drl_optind][0] == '+'
#endif				/* GETOPT_COMPAT */
	  ))
    {
      const struct option *p;
      char *s = nextchar;
      int exact = 0;
      int ambig = 0;
      const struct option *pfound = NULL;
      int indfound = 0;

      while (*s && *s != '=')
	s++;

      /* Test all options for either exact match or abbreviated matches.  */
      for (p = longopts, option_index = 0; p->name;
	   p++, option_index++)
	if (!strncmp (p->name, nextchar, s - nextchar))
	  {
	    if (s - nextchar == strlen (p->name))
	      {
		/* Exact match found.  */
		pfound = p;
		indfound = option_index;
		exact = 1;
		break;
	      }
	    else if (pfound == NULL)
	      {
		/* First nonexact match found.  */
		pfound = p;
		indfound = option_index;
	      }
	    else
	      /* Second nonexact match found.  */
	      ambig = 1;
	  }

      if (ambig && !exact)
	{
	  if (drl_opterr)
	    fprintf (stderr, "%s: option `%s' is ambiguous\n",
		     argv[0], argv[drl_optind]);
	  nextchar += strlen (nextchar);
	  drl_optind++;
	  return '?';
	}

      if (pfound != NULL)
	{
	  option_index = indfound;
	  drl_optind++;
	  if (*s)
	    {
	      /* Don't test has_arg with >, because some C compilers don't
		 allow it to be used on enums.  */
	      if (pfound->has_arg)
		drl_optarg = s + 1;
	      else
		{
		  if (drl_opterr)
		    {
		      if (argv[drl_optind - 1][1] == '-')
			/* --option */
			fprintf (stderr,
				 "%s: option `--%s' doesn't allow an argument\n",
				 argv[0], pfound->name);
		      else
			/* +option or -option */
			fprintf (stderr,
			     "%s: option `%c%s' doesn't allow an argument\n",
			     argv[0], argv[drl_optind - 1][0], pfound->name);
		    }
		  nextchar += strlen (nextchar);
		  return '?';
		}
	    }
	  else if (pfound->has_arg == 1)
	    {
	      if (drl_optind < argc)
		drl_optarg = argv[drl_optind++];
	      else
		{
		  if (drl_opterr)
		    fprintf (stderr, "%s: option `%s' requires an argument\n",
			     argv[0], argv[drl_optind - 1]);
		  nextchar += strlen (nextchar);
		  return '?';
		}
	    }
	  nextchar += strlen (nextchar);
	  if (longind != NULL)
	    *longind = option_index;
	  if (pfound->flag)
	    {
	      *(pfound->flag) = pfound->val;
	      return 0;
	    }
	  return pfound->val;
	}
      /* Can't find it as a long option.  If this is not drl_getopt_long_only,
	 or the option starts with '--' or is not a valid short
	 option, then it's an error.
	 Otherwise interpret it as a short option.  */
      if (!long_only || argv[drl_optind][1] == '-'
#ifdef GETOPT_COMPAT
	  || argv[drl_optind][0] == '+'
#endif				/* GETOPT_COMPAT */
	  || my_index ((char*)optstring, *nextchar) == NULL)
	{
	  if (drl_opterr)
	    {
	      if (argv[drl_optind][1] == '-')
		/* --option */
		fprintf (stderr, "%s: unrecognized option `--%s'\n",
			 argv[0], nextchar);
	      else
		/* +option or -option */
		fprintf (stderr, "%s: unrecognized option `%c%s'\n",
			 argv[0], argv[drl_optind][0], nextchar);
	    }
	  nextchar = (char *) "";
	  drl_optind++;
	  return '?';
	}
    }

  /* Look at and handle the next option-character.  */

  {
    char c = *nextchar++;
    char *temp = my_index ((char*)optstring, c);

    /* Increment `drl_optind' when we start to process its last character.  */
    if (*nextchar == '\0')
      ++drl_optind;

    if (temp == NULL || c == ':')
      {
	if (drl_opterr)
	  {
	    if (c < 040 || c >= 0177)
	      fprintf (stderr, "%s: unrecognized option, character code 0%o\n",
		       argv[0], c);
	    else
	      fprintf (stderr, "%s: unrecognized option `-%c'\n", argv[0], c);
	  }
	return '?';
      }
    if (temp[1] == ':')
      {
	if (temp[2] == ':')
	  {
	    /* This is an option that accepts an argument optionally.  */
	    if (*nextchar != '\0')
	      {
		drl_optarg = nextchar;
		drl_optind++;
	      }
	    else
	      drl_optarg = 0;
	    nextchar = NULL;
	  }
	else
	  {
	    /* This is an option that requires an argument.  */
	    if (*nextchar != '\0')
	      {
		drl_optarg = nextchar;
		/* If we end this ARGV-element by taking the rest as an arg,
		   we must advance to the next element now.  */
		drl_optind++;
	      }
	    else if (drl_optind == argc)
	      {
		if (drl_opterr)
		  fprintf (stderr, "%s: option `-%c' requires an argument\n",
			   argv[0], c);
		c = '?';
	      }
	    else
	      /* We already incremented `drl_optind' once;
		 increment it again when taking next ARGV-elt as argument.  */
	      drl_optarg = argv[drl_optind++];
	    nextchar = NULL;
	  }
      }
    return c;
  }
}

int
drl_getopt (argc, argv, optstring)
     int argc;
     char *const *argv;
     const char *optstring;
{
  return _drl_getopt_internal (argc, argv, optstring,
			   (const struct option *) 0,
			   (int *) 0,
			   0);
}

int
drl_getopt_long (argc, argv, options, long_options, opt_index)
     int argc;
     char *const *argv;
     const char *options;
     const struct option *long_options;
     int *opt_index;
{
  return _drl_getopt_internal (argc, argv, options, long_options, opt_index, 0);
}

/* Like drl_getopt_long, but '-' as well as '--' can indicate a long option.
   If an option that starts with '-' (not '--') doesn't match a long option,
   but does match a short option, it is parsed as a short option
   instead. */

int 
drl_getopt_long_only (argc, argv, options, long_options, opt_index)
     int argc;
     char *const *argv;
     const char *options;
     const struct option *long_options;
     int *opt_index;
{
  return _drl_getopt_internal (argc, argv, options, long_options, opt_index, 1);
}


#ifdef TEST

/* Compile with -DTEST to make an executable for use in testing
   the above definition of `drl_getopt'.  */

int
main (argc, argv)
     int argc;
     char **argv;
{
  int c;
  int digit_drl_optind = 0;

  while (1)
    {
      int this_option_drl_optind = drl_optind ? drl_optind : 1;

      c = drl_getopt (argc, argv, "abc:d:0123456789");
      if (c == EOF)
	break;

      switch (c)
	{
	case '0':
	case '1':
	case '2':
	case '3':
	case '4':
	case '5':
	case '6':
	case '7':
	case '8':
	case '9':
	  if (digit_drl_optind != 0 && digit_drl_optind != this_option_drl_optind)
	    printf ("digits occur in two different argv-elements.\n");
	  digit_drl_optind = this_option_drl_optind;
	  printf ("option %c\n", c);
	  break;

	case 'a':
	  printf ("option a\n");
	  break;

	case 'b':
	  printf ("option b\n");
	  break;

	case 'c':
	  printf ("option c with value `%s'\n", drl_optarg);
	  break;

	case '?':
	  break;

	default:
	  printf ("?? drl_getopt returned character code 0%o ??\n", c);
	}
    }

  if (drl_optind < argc)
    {
      printf ("non-option ARGV-elements: ");
      while (drl_optind < argc)
	printf ("%s ", argv[drl_optind++]);
      printf ("\n");
    }

  exit (0);
}

#endif /* TEST */
/*#endif	defined(_INCLUDE_DRL_GETOPT) */
