/* -------------------------------------------------------------------------

   File: CCdate.c
   Path: /home/julia/projets/dev/Common/libdate/SCCS/s.CCdate.c
   Description: bib de date
   Created: 96/10/06
   Author: Jacques WERNERT
   Modified: 00/03/20 15:24:02
   Last maintained by: Jacques WERNERT
   Revision: 1.23

   -------------------------------------------------------------------------

   Note:

   ------------------------------------------------------------------------- */

#include "CCcommon.h"
SCCS_ID (CCdate_c_SccsId, "@(#)CCdate.c	1.23, modified 00/03/20");

#define _CCdate_c
#include "CCdate.h"

#define BADDATE	-1L
#define MINDATE	4749L
#define MAXDATE	10852487L 
#define GREGOR	2299161L
#define IGREG	(15 + 31L * (10 + 12L * 1582))

/*
**      Function name: DAT_is_week_end
**
**      Description: est-ce que le timestamp GMT correspond a un week-end ?
**      Input: timestamp GMT
**      Output: 1: Oui / 0: Non / -1: Erreur
*/
#ifdef __STDC__
int DAT_is_week_end (time_t timestamp)
#else
int DAT_is_week_end (timestamp)
     time_t timestamp;
#endif
{
  
  struct tm *tm;
  
#ifdef Solaris
  
  struct tm w_tm;
  
  tm = localtime_r (&timestamp, &w_tm);
  
#else
  
  tm = localtime (&timestamp);
  
#endif
  
  if (tm) {
    return ((tm->tm_wday == 0) || (tm->tm_wday == 6));
  } else {
    return -1;
  }
  
} /* DAT_is_week_end */

/*
**      Function name: DAT_gmt_to_syb_string
**
**      Description: convertit un timestamp GMT en une chaine utilisable 
**                   dans SYBASE via convert
**      Input: timestamp GMT
**      Output: chaine STATIQUE
**      MT-Level: MT-unsafe
*/
#ifdef __STDC__
char *DAT_gmt_to_syb_string (time_t timestamp)
#else
char *DAT_gmt_to_syb_string (timestamp)
     time_t timestamp;
#endif
{
  static char string[DAT_FSTRING_SIZE];
  struct tm *tm;
  
#ifdef Solaris
  
  struct tm w_tm;
  
  tm = localtime_r (&timestamp, &w_tm);
  
#else
  
  tm = localtime (&timestamp);
  
#endif
  
  strftime (string, DAT_FSTRING_SIZE, DAT_FSTRING_FORMAT, tm);
  
  return (string);
  
} /* DAT_gmt_to_syb_string */

/*
**      Function name: DAT_gmt_to_syb_string_r
**
**      Description: convertit un timestamp GMT en une chaine utilisable 
**                   dans SYBASE via convert
**      Input: timestamp GMT
**      Output: chaine allouee
**      MT-Level: MT-safe
*/
#ifdef __STDC__
char *DAT_gmt_to_syb_string_r (time_t timestamp)
#else
char *DAT_gmt_to_syb_string_r (timestamp)
     time_t timestamp;
#endif
{

  char *string = (char *)malloc (sizeof (char) * DAT_FSTRING_SIZE);
  struct tm *tm;
  
#ifdef Solaris
  
  struct tm w_tm;
  
  tm = localtime_r (&timestamp, &w_tm);
  
#else
  
  tm = localtime (&timestamp);
  
#endif
  
  strftime (string, DAT_FSTRING_SIZE, DAT_FSTRING_FORMAT, tm);
  
  return (string);
  
} /* DAT_gmt_to_syb_string_r */

/*
**      Function name: DAT_gmt_to_XLstring
**
**      Description: converti un timestamp GMT en une chaine XL
**      Input: timestamp GMT
**      Output: chaine allouee
**      MT-Level: MT-unsafe
*/
#ifdef __STDC__
char *DAT_gmt_to_XLstring (time_t timestamp)
#else
char *DAT_gmt_to_XLstring (timestamp)
     time_t timestamp;
#endif
{
  char *res, *work = DAT_gmt_to_string (timestamp);
  
  res = (char *)malloc(sizeof (char) * strlen (work) + 2);
  if (! res) return (NULL);
  
  strcpy (res + 1, work);
  res[0] = ' ';
  
  return (res);
  
} /* DAT_gmt_to_XLstring */

/*
**      Function name: DAT_gmt_to_string
**
**      Description: convertie un timestamp GMT en une chaine
**      Input: timestamp GMT
**      Output: chaine STATIQUE
**      MT-Level: MT-unsafe
*/
#ifdef __STDC__
char *DAT_gmt_to_string (time_t timestamp)
#else
char *DAT_gmt_to_string (timestamp)
     time_t timestamp;
#endif
{
  char *heurelocale;
  int slen;
  
  heurelocale = ctime (&timestamp);
  if (! heurelocale) return (NULL);
  
  /*
   * suppression du \n
   */
  ((slen=strlen(heurelocale)) > 0) ?
    (heurelocale[slen-1]='\0') : (heurelocale[0]='\0');
  
  return (heurelocale);
  
} /* DAT_gmt_to_string */

#ifdef Solaris

/*
**      Function name: DAT_gmt_to_string_r
**
**      Description: convertie un timestamp GMT en une chaine
**      Input: timestamp GMT
**      Output: chaine allouee
**      MT-Level: MT-safe
*/
#ifdef __STDC__
char *DAT_gmt_to_string_r (time_t timestamp)
#else
char *DAT_gmt_to_string_r (timestamp)
     time_t timestamp;
#endif
{
  
  int slen;
  char *heurelocale;
  char *w_heurelocale = (char *)malloc (sizeof (char) * DAT_STRING_SIZE);
  struct tm *tm;
  
  if (! w_heurelocale) return (NULL);

  heurelocale = ctime_r (&timestamp, w_heurelocale, DAT_STRING_SIZE);
  if (! heurelocale) return (NULL);
  
  /*
   * suppression du \n
   */
  ((slen=strlen(heurelocale)) > 0) ?
    (heurelocale[slen-1]='\0') : (heurelocale[0]='\0');
  
  return (heurelocale);
  
} /* DAT_gmt_to_string_r */

#endif /* Solaris */

/*
 **      Function name: DAT_gmt_to_mystring
 **
 **      Description: convertit un timestamp GMT en une chaine
 **      Input: timestamp GMT
 **      Output: buf
 **      MT-Level: MT-safe
 */
#ifdef __STDC__
char *DAT_gmt_to_mystring (time_t timestamp, char *buf, int bufsize)
#else
char *DAT_gmt_to_mystring (timestamp, buf, bufsize)
     time_t timestamp;
     char *buf;
     int bufsize;
#endif /* __STDC__ */
{
  char *heurelocale;
  int slen;

#ifdef Solaris

  if (bufsize < DAT_STRING_SIZE) return (NULL);
  heurelocale = ctime_r (&timestamp, buf, bufsize);

#elif defined SunOS

  heurelocale = ctime (&timestamp);
  if (heurelocale) {
    if (bufsize >= strlen (heurelocale) + 1) {
      CC_STRCPY (buf, heurelocale, bufsize);
    }
  }

#elif defined WIN32

  heurelocale = ctime (&timestamp);
  if (heurelocale) {
    if (bufsize >= (int)strlen (heurelocale) + 1) {
      CC_STRCPY (buf, heurelocale, bufsize);
    }
  }
  
#endif

  if (! heurelocale) return (NULL);
  
  /*
   * suppression du \n
   */
  ((slen=strlen(buf)) > 0) ?
    (buf[slen-1]='\0') : (buf[0]='\0');
  
  return (buf);
  
} /* DAT_gmt_to_mystring */

/*
**      Function name: DAT_gmt_to_secs_since_midnight
**
**      Description: converti un timestamp GMT en nb de secondes
**                   depuis 0h00 en heure locale
**      Input: timestamp GMT
**      Output: Nb secs depuis 0h00 (heure locale)
**      MT-Level: MT-safe
*/
time_t DAT_gmt_to_secs_since_midnight (timestamp)
     time_t timestamp;
{
  struct tm *tm;
  
#ifdef Solaris
  
  struct tm w_tm;
  
  tm = localtime_r (&timestamp, &w_tm);
  
#else
  
  tm = localtime (&timestamp);
  
#endif
  
  if (! tm) return (DAT_error);
  
  return (3600*tm->tm_hour+60*tm->tm_min+tm->tm_sec);
  
} /* DAT_gmt_to_secs_since_midnight */

/*
**      Function name: DAT_secs_since_midnight
**
**      Description: Converti un nombre de secondes depuis 0h00
**                   en heure locale en un timestamp GMT
**      Input: Nb secs depuis 0h00 (heure locale)
**      Output: timestamp GMT ou DAT_error
**      MT-Level: MT-unsafe (a porter)
*/
time_t DAT_secs_since_midnight_to_gmt (secs, today)
     time_t secs, today;
{
  struct tm *tm;
  
  if (today == 0)
    today = DAT_now;
  
  if (! (tm = localtime (&today)))
    return (DAT_error);
  
  tm->tm_sec = secs % 60;
  secs /= 60;
  tm->tm_min = secs % 60;
  secs /= 60;
  tm->tm_hour = secs % 60;
  
#ifdef SunOS
  today = timelocal (tm);
#endif SunOS
  
#ifdef Solaris
  /*
   * Positionne le champ isdst / cf man mktime
   * if tm_isdst  is  zero,  the  original  values  are
   * assumed  to be in the main timezone and are converted to the
   * alternate timezone if the main timezone is  not  valid
   */
  tm->tm_isdst = -1;

  today = mktime (tm);
#endif Solaris
  
  return (today);
  
} /* DAT_secs_since_midnight_to_gmt */

/*
**      Function name: DAT_ssdate_to_gmt
**
**      Description: Convertit un timestamp local en timestamp GMT
**      Input: timestamp local (pb possible d'arrondi sur 1 seconde)
**      Output: timestamp GMT ou DAT_error
**      MT-Level: MT-unsafe
*/
time_t DAT_ssdate_to_gmt (ssdate)
     double ssdate;
{
  
  /*
   * Les dates tableurs sont exprimees en fraction de jours
   * -> passage en secondes
   */
  ssdate *= 24 * 3600;
  
  /*
   * Tout est un probleme de referentiel
   */
  ssdate -= 25569.0 * 24 * 3600;
  
  return (DAT_local_to_gmt ((time_t)ssdate));
  
} /* DAT_ssdate_to_gmt */

/*
**      Function name: DAT_gmt_to_ssdate
**
**      Description: Converti un timestamp local en timestamp GMT
**      Input: timestamp local
**      Output: timestamp GMT
**      MT-Level: MT-unsafe
*/
double DAT_gmt_to_ssdate (timestamp)
     time_t timestamp;
{
  
  double ssdate;
  
  ssdate = DAT_gmt_to_local (timestamp);
  
  /*
   * Tout est un probleme de referentiel
   */
  ssdate += 25569.0 * 24 * 3600;
  
  /*
   * Les dates tableurs sont exprimees en fraction de jours
   */
  ssdate /= 24 * 3600;
  
  return (ssdate);
  
} /* DAT_gmt_to_ssdate */

/*
**      Function name: DAT_local_to_gmt
**
**      Description: Convertie un timestamp local en timestamp GMT
**      Input: timestamp local
**      Output: timestamp GMT ou DAT_error
**      MT-Level: MT-unsafe
*/
time_t DAT_local_to_gmt (local)
     time_t local;
{
  
  struct tm *tm;
  time_t gmt = DAT_error;
  
  /*
   * Eclate le timestamp local sans appliquer de DST ni de decalage
   */
  if ( !(tm = gmtime (&local))) 
    return ((time_t)DAT_error);
  
  /*
   * Transformation de la structure en timestamp GMT
   */
#ifdef SunOS
  gmt = timelocal (tm);
#endif SunOS
  
#ifdef Solaris
  /*
   * Positionne le champ isdst / cf man mktime
   * if tm_isdst  is  zero,  the  original  values  are
   * assumed  to be in the main timezone and are converted to the
   * alternate timezone if the main timezone is  not  valid
   */
  tm->tm_isdst = -1;
  
  gmt = mktime (tm);
#endif /* Solaris */
  
#ifdef WIN32
  tm->tm_isdst = 0;
  gmt = mktime (tm);
#endif /* WIN32 */
  
  return (gmt);
  
} /* DAT_local_to_gmt */

/*
**      Function name: DAT_gmt_to_local
**
**      Description: Convertit un timestamp local en timestamp GMT
**      Input: timestamp local
**      Output: timestamp GMT ou DAT_error
**      MT-Level: MT-unsafe
*/
time_t DAT_gmt_to_local (gmt)
     time_t gmt;
{
  
  struct tm *tm;
  time_t local = DAT_error;
  
  /*
   * Conversion du timestamp en structure gmt
   * On considere la structure de sortie comme locale,
   * reconversion en timestamp
   */
  if (! (tm = gmtime (&gmt))) 
    return (DAT_error);
  
#ifdef SunOS
  local = timelocal (tm);
#endif SunOS
  
#ifdef Solaris
  /*
   * Positionne le champ isdst / cf man mktime
   * if tm_isdst  is  zero,  the  original  values  are
   * assumed  to be in the main timezone and are converted to the
   * alternate timezone if the main timezone is  not  valid
   */
  tm->tm_isdst = -1;

  local = mktime (tm);
#endif /* Solaris */
  
#ifdef WIN32
  tm->tm_isdst = 0;
  local = mktime (tm);
#endif /* WIN32 */
  
  if (local == -1) {
    
    /*
     * cas particulier ou la structure gmt indique une date non valide en local
     *  ex: 01/01/70 00h00 (timestamp 0) donc probleme a la reconversion
     * -> rappelle en decalant l'heure gmt si le timestamp est valide en local
     */
    if (! (tm = localtime (&gmt)))
      return (DAT_error);
    
    return (DAT_gmt_to_local (gmt + 3600) - 3600);
    
  }
  
  return (gmt + gmt - local);
  
} /* DAT_gmt_to_local */

/*
**      Function name: DAT_gmt_to_struct
**
**      Description: Converti un time_t en un ensemble {ymdhms}
**      Input: timestamp GMT * {annee, mois, jour, heure, minute, seconde}
**                           ou null
**      Output: 1: Ok / 0: Erreur
**      MT-Level: MT-safe
*/
#ifdef __STDC__
int DAT_gmt_to_struct (time_t gmt, 
		       int *year, int *month, int *day, 
		       int *hour, int *minute, int *second)
#else
     int DAT_gmt_to_struct (gmt, year, month, day, hour, minute, second)
     time_t gmt;
     int *year, *month, *day, *hour, *minute, *second;
#endif
{
  
  struct tm *tm;
  
#ifdef Solaris
  
  struct tm w_tm;
  
  tm = localtime_r (&gmt, &w_tm);
  
#else
  
  tm = localtime (&gmt);
  
#endif
  
  if (! tm) return (0);
  
  if (year) *year = tm->tm_year + 1900;
  if (month) *month = tm->tm_mon + 1;
  if (day) *day = tm->tm_mday;
  if (hour) *hour = tm->tm_hour;
  if (minute) *minute = tm->tm_min;
  if (second) *second = tm->tm_sec;
  
  return (1);
  
} /* DAT_gmt_to_struct */

/*
**      Function name: DAT_struct_to_gmt
**
**      Description: Convertit un ensemble {ymdhms} en time_t
**      Input: annee, mois, jour, heure, minute, seconde
**      Output: timestamp GMT ou DAT_error
**      MT-Level: MT-unsafe
*/
#ifdef __STDC__
time_t DAT_struct_to_gmt (int year, int month, int day, 
			  int hour, int minute, int second)
#else
     time_t DAT_struct_to_gmt (year, month, day, hour, minute, second)
     int year, month, day, hour, minute, second;
#endif
{
  
  static struct tm tm;
  time_t gmt;
  
  tm.tm_year = year - 1900;
  tm.tm_mon = month - 1;
  tm.tm_mday = day;
  tm.tm_hour = hour;
  tm.tm_min = minute;
  tm.tm_sec = second;
  
  /*
   * Transformation de la structure en timestamp GMT
   */
#ifdef SunOS
  gmt = timelocal (tm);
#endif /* SunOS */
  
#ifdef Solaris
  /*
   * Positionne le champ isdst / cf man mktime
   * if tm_isdst  is  zero,  the  original  values  are
   * assumed  to be in the main timezone and are converted to the
   * alternate timezone if the main timezone is  not  valid
   */
  tm.tm_isdst = -1;

  gmt = mktime (&tm);
#endif /* Solaris */
  
#ifdef WIN32
  tm.tm_isdst = 0;
  gmt = mktime (&tm);
#endif /* WIN32 */
  
  return (gmt);
  
} /* DAT_struct_to_gmt */

#ifdef Solaris

/*
**      Function name: DAT_gmt_to_fstring_r
**
**      Description: Convertie un timestamp en chaine
**      Input: le timestamp et son format (NULL = format par defaut)
**      Output: chaine ALLOUEE convertie ou NULL
**      MT-Level: MT-safe
*/
#ifdef __STDC__
char *DAT_gmt_to_fstring_r (time_t gmt, char *format)
#else
char *DAT_gmt_to_fstring_r (gmt, format)
     time_t gmt;
     char *format;
#endif
{
  
  char *buf;
  struct tm tm;
  
  if (! (buf = (char *)malloc (sizeof (char) * DAT_STRING_SIZE)))
    return (NULL);
  
  if (! localtime_r (&gmt, &tm)) return (NULL);
  
  if (strftime (buf, DAT_STRING_SIZE, format, &tm) == 0)
    return (NULL);
  else
    return (buf);
  
} /* DAT_gmt_to_fstring_r */

#endif /* Solaris */

/*
**      Function name: DAT_gmt_to_fstring
**
**      Description: Convertit un timestamp en chaine
**      Input: le timestamp et son format (NULL = format par defaut)
**      Output: chaine STATIQUE convertie ou NULL
**      MT-Level: MT-unsafe
*/
#ifdef __STDC__
char *DAT_gmt_to_fstring (time_t gmt, char *format)
#else
char *DAT_gmt_to_fstring (gmt, format)
     time_t gmt;
     char *format;
#endif
{
  
  static char buf[DAT_STRING_SIZE];
  struct tm *tm;
  
  tm = localtime (&gmt);
  
  if (! tm) return (NULL);
  
  if (strftime (buf, DAT_STRING_SIZE, format, tm) == 0)
    return (NULL);
  else
    return (buf);
  
} /* DAT_gmt_to_fstring */

#if Solaris
/*
**      Function name: DAT_fstring_to_gmt_r
**
**      Description: Convertie une chaine en timestamp
**      Input: la chaine et son format (NULL = format par defaut)
**      Output: timestamp GMT ou DAT_error
**      MT-Level: MT-safe
*/
#ifdef __STDC__
time_t DAT_fstring_to_gmt_r (char *source, char *format)
#else
time_t DAT_fstring_to_gmt_r (source, format)
     char *source, *format;
#endif
{

  struct tm tm;
  time_t gmt;
  static pthread_mutex_t DAT_fstring_to_gmt_r_lock = PTHREAD_MUTEX_INITIALIZER;
  
  if (! strptime (source, format, &tm))
    return (DAT_error);
  
  /*
   * debut de section crittique
   */
  mutex_lock (&DAT_fstring_to_gmt_r_lock);
  
  /*
   * Transformation de la structure en timestamp GMT
   */
#ifdef SunOS
  gmt = timelocal (&tm);
#endif SunOS
  
#ifdef Solaris
  /*
   * Positionne le champ isdst / cf man mktime
   * if tm_isdst  is  zero,  the  original  values  are
   * assumed  to be in the main timezone and are converted to the
   * alternate timezone if the main timezone is  not  valid
   */
  tm.tm_isdst = -1;

  gmt = mktime (&tm);
#endif /* Solaris */
  
#ifdef WIN32
  tm.tm_isdst = 0;
  gmt = mktime (&tm);
#endif /* WIN32 */
  
  /*
   * debut de section crittique
   */
  mutex_unlock (&DAT_fstring_to_gmt_r_lock);
  
  return (gmt);
  
} /* DAT_fstring_to_gmt_r */

#endif /* Solaris */

/*
**      Function name: DAT_fstring_to_gmt
**
**      Description: Convertit une chaine en timestamp
**      Input: la chaine et son format (NULL = format par defaut)
**      Output: timestamp GMT ou DAT_error
**      MT-Level: MT-unsafe
*/
#ifdef __STDC__
time_t DAT_fstring_to_gmt (char *source, char *format)
#else
time_t DAT_fstring_to_gmt (source, format)
     char *source, *format;
#endif
{
  
  static struct tm tm;
  time_t gmt;
  
#ifdef WIN32
  
  /*
   * TODO
   */
  return (DAT_error);
  
#else
  
  if (! strptime (source, format, &tm))
    return (DAT_error);
  
#endif /* WIN32 */
  
  /*
   * Transformation de la structure en timestamp GMT
   */
#ifdef SunOS
  gmt = timelocal (&tm);
#endif SunOS
  
#ifdef Solaris
  /*
   * Positionne le champ isdst / cf man mktime
   * if tm_isdst  is  zero,  the  original  values  are
   * assumed  to be in the main timezone and are converted to the
   * alternate timezone if the main timezone is  not  valid
   */

  tm.tm_isdst = -1;

  gmt = mktime (&tm);
#endif /* Solaris */
  
#ifdef WIN32
  tm.tm_isdst = 0;
  gmt = mktime (&tm);
#endif /* WIN32 */
  
  return (gmt);
  
} /* DAT_fstring_to_gmt */

/*
**      Function name: DAT_gmt_to_secs_until_midnight 
**
**      Description: Calcul le nombre de secondes restantes jusqu'a minuit 
**      Input: neant
**      Output: un timestamp representant le nombre de secondes restantes jusqu'a minuit
**      MT-Level: MT-safe
*/
#ifdef __STDC__
time_t DAT_gmt_to_secs_until_midnight ()
#else
time_t DAT_gmt_to_secs_until_midnight ()
#endif
{
  int y, m, d, hh, mm, ss;
  time_t now = DAT_now;
  
  DAT_gmt_to_struct (now, &y, &m, &d, &hh, &mm, &ss);
  
  if ((hh + mm + ss) == 0) return 0;
  
  return (DAT_struct_to_gmt (y, m, d + 1, 0, 0, 0) - now);
  
} /* DAT_gmt_to_secs_until_midnight */

/*
**      Function name: DAT_adgdate_to_gmt
**
**      Description: Convertit en une date Adagio (un entier au format AAAAMMJJ) en un timestamp GMT
**      Input: date Adagio 
**      Output: timestamp GMT ou DAT_error
**      MT-Level: MT-unsafe
*/
#ifdef __STDC__
time_t DAT_adgdate_to_gmt (int adgdate)
#else
time_t DAT_adgdate_to_gmt (adgdate)
        int adgdate;
#endif
{
	int y, m, d;
	div_t temp;

	temp = div (adgdate, 10000);
	y = temp.quot;	

	temp = div (temp.rem, 100);
	m = temp.quot;

	temp = div (temp.rem, 100); 
	d = temp.rem;

	return DAT_struct_to_gmt (y, m, d, 0, 0, 0);
} /* DAT_adgdate_to_gmt */

/*
**      Function name: DAT_gmt_to_adgdate
**
**      Description: Convertit un timestamp GMT en une date Adagio (un entier au format AAAAMMJJ) 
**      Input: timestamp GMT
**      Output: date Adagio (un entier au format AAAAMMJJ) ou -1 en cas d'erreur
**      MT-Level: MT-safe
*/
#ifdef __STDC__
int DAT_gmt_to_adgdate (time_t gmt)
#else
int DAT_gmt_to_adgdate (gmt)
	time_t gmt;
#endif
{
	int y, m, d, hh, mm, ss;
	int adgdate = -1;

	if(DAT_gmt_to_struct (gmt, &y, &m, &d, &hh, &mm, &ss) == 1)
	{
		adgdate = d + m * 100 + y * 10000;
	}

	return adgdate;
} /* DAT_gmt_to_adgdate */

long DAT_DMYtoJulian (int d, int m, int y, double* Julian)
{
	long Ja, Jm, Jy;
	double Jul;

	if(Julian == NULL)
	{
		return (-1);
	}

	if(y < 0)
	{
		y++;
	}

	if(m > 2)
	{
		Jy = y;
		Jm = m + 1;
	}
	else
	{
		Jy = y - 1;
		Jm = m + 13;
	}

	Jul = (floor (365.25 * (double)Jy) + floor (30.6001 * (double)Jm) + d + 1720995L);

	if(d + 31L * (m + 12L * y) >= IGREG )
	{
		Ja = (long)(0.01 * Jy);
		Jul += (2 - Ja + (long)(0.25 * Ja));
	}

	*Julian = Jul;

	return (0);
}

long DAT_JulianToDMY (double Julian, int* d, int* m, int* y)
{
	long Ja, JAlpha, Jb, Jc, Jd, Je;

	if((Julian != BADDATE) && (Julian >= MINDATE) && (Julian <= MAXDATE))
    	{
       		if( Julian >= GREGOR)
       		{
          		JAlpha = (long)(((double)(Julian - 1867216L) - 0.25) / 36524.25);
          		Ja = (long)(Julian + 1 + JAlpha - (long)(0.25 * JAlpha));
       		}
       		else
       		{
          		Ja = (long)Julian;
       		}
 
       		Jb = Ja + 1524;
       		Jc = (long)(6680.0 + ((double)(Jb - 2439870L) - 122.1) / 365.25);
       		Jd = (long)(365 * Jc + (0.25 * Jc));
       		Je = (long)((Jb - Jd) / 30.6001);
       		*d = (int)(Jb - Jd - (int)(30.6001 * Je));
       		*m = (int)Je - 1;

		if(*m > 12)
		{
          		*m -= 12;
		}
 
       		*y = (int)(Jc - 4715);
 
       		if(*m > 2)
		{
          		--(*y);
		}
 
       		if( *y <= 0)
		{
          		--(*y);
    		}
	}
    	else
    	{
		return (-1);
	}

	return (0);
}

long DAT_ssdate_to_struct (double ssdate, int* y, int* m, int* d)
{
	return (DAT_JulianToDMY( XLDateToJulian(ssdate ), d, m, y));	
}

double DAT_struct_to_ssdate (int y, int m, int d)
{
	double res;

	DAT_DMYtoJulian (d, m, y, &res);	

	return JulianToXLDate( res );
}

double XLDateToJulian(double d )
{
	return d + 2415019.0;
}

int isValidXLDate(double d)
{
	double  julDate = XLDateToJulian(d); 
	if((julDate != BADDATE) && (julDate >= MINDATE) && (julDate <= MAXDATE)) return 1; 
	return 0 ;
}
double JulianToXLDate( double d )
{
	return d - 2415019.0;
}


#ifdef TEST
/*
 * Programme de test
 */
int main ()
{
  time_t timestamp, local, timestamp2;
  struct tm *tm;
  int heure, i;
  double ssdate;
  double julian;
  int d, m, y;
  static time_t tab_stamp [] = { 0, 504954000, 512526000, 512611200, 512524801,
				 512532000, 517996800, 161420399, 887065100, 
				 889570800, 892159200, 631227930, 888710400,
				 512524801, 512528401, 811897200, 811900800,
				 811904400, 811908000, 888797400, 888883810,
				 899325000, DAT_error }; 
  
  timestamp = time(NULL);
  printf ("time_t GMT = %ld, #secs/0h00 = %li, time_t GMT = %ld\n",
	  timestamp, DAT_gmt_to_secs_since_midnight (timestamp),
	  DAT_secs_since_midnight_to_gmt (DAT_gmt_to_secs_since_midnight (timestamp)
					  , DAT_now));
  tm = localtime (&timestamp);
  printf ("-> %s %li %02li:%02li:%02li\n", asctime (tm), 
	  tm->tm_hour*3600+60*tm->tm_min+tm->tm_sec, tm->tm_hour, tm->tm_min,
	  tm->tm_sec);
  
#ifdef SunOS
  local = timestamp + tm->tm_gmtoff;
#elif Solaris
  tzset ();
  if (tm->tm_isdst)
    local = timestamp - altzone;
  else
    local = timestamp - timezone;
#endif Solaris
  heure = local % 86400;
  printf ("Methode JD, time_t GMT = %li, local = %li, heure = %i\n", 
	  timestamp, local, heure);
  
  printf ("\nEssais de conversions en dates tableurs ...\n");
  
  for (i = 0; tab_stamp[i] != DAT_error; i++) {
    ssdate = DAT_gmt_to_ssdate (tab_stamp[i]);
    timestamp2 = DAT_ssdate_to_gmt (ssdate);
    printf ("time_t = %ld, date tableur = %f, date Adg = %d, date reconvertie = %ld, diff = %d [%s]\n",
            tab_stamp[i], ssdate, DAT_gmt_to_adgdate (tab_stamp[i]), timestamp2, tab_stamp[i] - timestamp2, 
	    DAT_gmt_to_string (tab_stamp[i]));
  }

  fprintf (stderr, "Conversion ARM: %lf\n", DAT_gmt_to_ssdate (DAT_fstring_to_gmt ((char*)(const char*)"22.09.2000", "%d.%m.%Y")));

  DAT_DMYtoJulian (1, 1, 1900, &julian);
  fprintf (stderr, "01/01/1900 ==> %lf\n", julian);

  julian = 36605.0;

  DAT_ssdate_to_struct (julian, &y, &m, &d);
  fprintf (stderr, "%lf ==> %d/%d/%d\n", julian, d, m, y);
  fprintf (stderr, "%d/%d/%d ==> %lf\n", d, m, y, DAT_struct_to_ssdate (y, m, d)); 
  
  return (1);
}
/*
 * Fin du programme de test
 */
#endif TEST

/* EOF CCdate.c */
