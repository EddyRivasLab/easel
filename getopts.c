/* getopts.c
 * Implements a somewhat more powerful command line getopt interface
 * than the standard UNIX/POSIX call.
 * 
 * SVN $Id$
 * SRE, Sat Jan  1 08:50:21 2005 [Panticosa, Spain]
 */

#include <stdlib.h> 
#include <stdio.h> 
#include <string.h>
#include <ctype.h>


/* Function:  esl_getopts_Create()
 * Incept:    SRE, Tue Jan 11 11:24:16 2005 [St. Louis]
 *
 * Purpose:   
 * 
 *            Stores ptrs to <argc>, <argv>, <opt>, <usage>. Does not
 *            change these, nor their contents.
 *            
 *            Initializes.
 *            
 *            Sets default values for all options.
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
ESL_GETOPTS *
esl_getopts_Create(struct opt_s *opt, char *usage)
{
  ESL_GETOPTS *g;
  char *optname;
  int i,n;

  if ((g = malloc(sizeof(ESL_GETOPTS))) == NULL) goto FAILURE;
  g->opt       = opt;
  g->argc      = 0;
  g->argv      = NULL;
  g->usage     = usage;
  g->optind    = 1;
  g->islong    = NULL;
  g->val       = NULL;
  g->setby     = NULL;
  g->optstring = NULL;

  /* Figure out the number of options.
   */
  g->nopts = 0;
  while (g->opt[g->nopts].name != NULL)
    g->nopts++;
  
  /* Set flags on whether option is long or not,
   * for later convenience, by counting leading -'s.
   */
  if ((g->islong = malloc(sizeof(int) * g->nopts)) == NULL) goto FAILURE;
  for (i = 0; i < g->nopts; i++)
    {
      if ((n = strspn(g->opt[i].name, "-")) == 2)
	g->islong[i] = TRUE;
      else if (n==1)
	g->islong[i] = FALSE;
      else
	{
	  esl_getopts_Destroy(g);	  
	  ESL_ERROR_NULL(ESL_ECORRUPT, "internal error: malformed option");
	}
    }

  /* Set default values for all options.
   */
  if ((g->val   = malloc(sizeof(char *) * g->nopts)) == NULL) goto FAILURE;
  if ((g->setby = malloc(sizeof(int)    * g->nopts)) == NULL) goto FAILURE;
  for (i = 0; i < g->nopts; i++) 
    {
      g->val[i]   = g->opt[i].default;
      g->setby[i] = eslARG_SETBY_DEFAULT;
    }

  /* What the hell, verify the type/range of the defaults, even though it's
   * an application error (not user error) if they're invalid. 
   */
  for (i = 0; i < g->nopts; i++) 
    {
      switch (g->opt[i].type) {
	
      case eslARG_INT:
	if (! is_integer(g->val[i])) {
	  esl_getopts_Destroy(g);	  
	  esl_error(ESL_ECORRUPT, __FILE__, __LINE__, 
		    "option %s: expected integer valued default, got %s", 
		    g->opt[i].name, g->opt[i].default);
	  return NULL;
	}
	if (! verify_integer_range(g->val[i], g->opt[i].range)) {
	  esl_getopts_Destroy(g);	  
	  esl_error(ESL_ECORRUPT, __FILE__, __LINE__, 
		    "option %s: expected integer default in range %s, got %s", 
		    g->opt[i].name, g->opt[i].range, g->opt[i].default);
	  return NULL;
	}
	break;

      case eslARG_REAL:
	if /* SRE STOPPED HERE */
	  

	verify_integer(g->opt[i].name, usage, g->val[i], g->opt[i].range);
      else if (g->opt[i].type == eslARG_REAL)
	verify_real(g->opt[i].name, usage, g->val[i], g->opt[i].range);
      else if (g->opt[i].type == eslARG_CHAR)
	verify_char(g->opt[i].name, usage, g->val[i], g->opt[i].range);
    }

  /* Normal return.
   */
  return g;

  /* Abnormal return, on any alloc failure.
   */
 FAILURE:
  esl_getopts_Destroy(g);
  ESL_ERROR_NULL(ESL_EMEM, "allocation failed");
}

void
esl_getopts_Destroy(ESL_GETOPTS *g)
{
  if (g != NULL)
    {
      if (g->val   != NULL) free(g->val);
      if (g->setby != NULL) free(g->setby);
      free(g);
    }
}


int
esl_opt_ProcessCmdline(ESL_GETOPTS *g, int argc, char **argv)
{
  int   opti;
  int   togi;			/* index of a toggle-tied option */
  char *togname;		/* name of a toggle-tied option */
  char *optname;
  char *optarg;
  char  buf[64];	/* must be large enough to hold largest option + \0 */
  char  s;

  g->argc      = argc;
  g->argv      = argv;
  g->optind    = 1;
  g->optstring = NULL;

  while (esl_getopts(g, &opti, &optname, &optarg) == ESL_OK)
    {
      /* Have we already set this option? */
      if (g->setby[opti] == eslARG_SETBY_CMDLINE)
	opterror(g, "Option %s appears more than once on command line", optname);

      /* Set the option. 
       */
      g->setby[opti] = eslARG_SETBY_CMDLINE;
      if (g->opt[opti].type == eslARG_NONE)
	g->val[opti] = (char *) TRUE;
      else
	g->val[opti] = optarg;

      /* Unset toggle-tied options.
       */
      s = g->opt[opti].toggle_opts;
      while ((togi = process_optlist(g, &s, &togname)) != -1)
	{
	  if (g->setby[togi] == eslARG_SETBY_CMDLINE)
	    opterror(g, "Options %s and %s conflict, toggling each other", togname, optname);
	  
	  g->setby[togi] = eslARG_SETBY_CMDLINE;
	  if (g->opt[togi].type == eslARG_NONE)
	    g->val[togi] = (char *) FALSE;
	  else
	    g->val[opti] = NULL;
	}
    }
  return ESL_OK;
}


int
esl_opt_VerifyConfig(ESL_GETOPTS *g)
{
  int   i,j;
  char *s;
  char *reqopt;

  /* For all options that are set, 
   * overify that all their required_opts are set.
   */
  for (i = 0; i < g->nopts; i++)
    {
      if (g->setby[i] != eslARG_SETBY_DEFAULT)
	{
	  s = g->opt[i].required_opts;
	  while ((j = process_optlist(g, &s, &reqopt)) != -1)      
	    {
	      if (g->setby[j] == eslARG_SETBY_DEFAULT)
		opterror(g, "%s requires (or has no effect without) %s", 

	}


}

/* Return index <opti> of the next option in <s> up to 
 * next comma, or -1 if we're out of data (s == NULL).
 * Reset <s> to follow the ',' or to NULL if we've
 * just done the last option in the list.
 */
int 
process_optlist(ESL_GETOPTS *g, char **ret_s, char **ret_togname)
{
  char *s;
  int   i;
  int   n;
  
  if ((s = *ret_s) == NULL) return -1;
  n = strcspn(s, ",");
  if (s[1] == '-')		/* if 2nd char is a -, long option */
    {
      for (i = 0; i < g->nopts; i++)
	if (g->opt[i].longname != NULL && strncmp(g->opt[i].longname, s, n) == 0) break;
      if (i == g->nopts) return -1;
      *ret_togname = g->opt[i].longname;
    }
  else 
    {
      for (i = 0; i < g->nopts; i++)
	if (g->opt[i].stdname != NULL && strncmp(g->opt[i].stdname, s, n) == 0)  break;
      if (i == g->nopts) return -1;
      *ret_togname = g->opt[i].shortname;
    }

  if (s[n] == ',') *ret_s = s+n+1; 
  else             *ret_s = NULL;
  return i;
}
  


int
get_optidx(ESL_GETOPTS *g, char *option)
{
  int i;

  if (option[1] == '-')		/* if 2nd char is -, this is a long option */
    {
      for (i = 0; i < g->nopts; i++)
	if (g->opt[i].longname != NULL && strcmp(g->opt[i].longname, option) == 0) break;
    }
  else 				/* else we're a short option */
    {
      for (i = 0; i < g->nopts; i++)
	if (g->opt[i].stdname != NULL && strcmp(g->opt[i].stdname,  option) == 0)  break;
    }
  if (i == g->nopts) return -1;
  else               return i;
}

      
int
esl_getopts(ESL_GETOPTS *g, int *ret_opti, char **ret_optname, char **ret_optarg)
{
  char *argptr;
  int   opti;
  int   nmatch;

  *ret_optname = NULL; 
  *ret_optarg  = NULL; 

  /* Check to see if we've run out of options.
   * A '-' by itself is an argument (e.g. "read from stdin"), not an option.
   */
  if (g->optind >= g->argc || g->argv[g->optind][0] != '-' || strcmp(g->argv[g->optind], "-") == 0)
    return ESL_EOD; 		/* normal end-of-data (end of options) return  */

  /* Check to see if we're being told that this is the end
   * of the options with the special "--" flag.
   */
  if (strcmp(g->argv[g->optind], "--") == 0)
    { 
      g->optind++;
      return ESL_EOD; 		/* also a normal end-of-data return */
    }

  /* We have an option: an argv element that starts with -, but is
   * not "-" or "--".
   *
   * It's a long option if (1) we're not in the middle of a list of
   * one-char options and (2) it starts with "--". 
   */
  if (g->optstring == NULL && g->argv[optind][1] == '-')
    process_longopt(g, &opti, ret_optname, ret_optarg);
  else 
    process_stdopt(g, &opti, ret_optname, ret_optarg);

  /* Type and range checking.
   * eslARG_STRING types are unchecked, free form.
   */
  if      (g->opt[opti].type == eslARG_INT) 
    verify_integer(*ret_optname, usage, *ret_optarg, g->opt[opti].range);
  else if (g->opt[opti].type == eslARG_REAL)
    verify_real(*ret_optname, usage, *ret_optarg, g->opt[opti].range);
  else if (g->opt[opti].type == eslARG_CHAR)
    verify_char(*ret_optname, usage, *ret_optarg, g->opt[opti].range);

  /* Normal return.
   */
  *ret_opti = opti;

  return ESL_OK;
}

/* process_longopt():
 * optind is sitting on a long option, w/ syntax of one of these forms:
 *       --foo        
 *       --foo arg
 *       --foo=arg
 * (GNU getopt long option syntax.)
 * 
 * Allow abbreviations of long options, 
 * 
 * Process it, returning option number, ptr to the optname (the full copy of it
 * in opt[] structure, not the possibly abbreviated one in argv[]),
 * and optarg.
 *
 * Exits w/ mesg if:
 *   1. Option doesn't exist in opt[].
 *   2. Option abbreviation is ambiguous, matching opt[] nonuniquely.
 *   3. Option takes an argument, but no argument found.
 *   4. Option does not take an argument, but one was provided by =arg syntax.
 * 
 * g->optind is advanced to next argv element (+1, +2, +1, respectively,
 * for --foo, --foo arg, --foo=arg).
 */
static void
process_longopt(ESL_GETOPTS *g, int *ret_opti, char **ret_optname, char **ret_optarg)
{
  int   i;		/* counter for options                               */
  int   opti;		/* option number found                               */
  char *argptr;		/* ptr to arg in --foo=arg syntax                    */
  int   nmatch;		/* number of options in opt[] that match             */
  int   arglen;		/* length of argv elem's option part (up to = or \0) */

  /* Deal with options of syntax "--foo=arg" w/o modifying argv.
   */
  if ((argptr = strchr(g->argv[optind], '=')) != NULL)
    { arglen = argptr - g->argv[optind]; argptr++; } /* bump argptr off the = to the arg itself */
  else
    { arglen = strlen(g->argv[optind]); } /* and argptr == NULL from above. */

  /* Figure out which option this is.
   * One trick here is to allow abbreviations, and identify
   * ambiguities while we're doing it. (GNU getopt allows abbrevs.)
   */
  nmatch = 0;
  for (i = 0; i < g->nopts; i++)
    if (strncmp(g->opt[i].longname, g->argv[g->optind], arglen) == 0)
      {
	nmatch++;
	opti = i;
	if (arglen == strlen(g->opt[i].longname)) break; /* exact match; can stop now */
      }
  if (nmatch > 1)
    { 
      fprintf(stderr, "Abbreviated option \"%s\" is ambiguous.\n%s", g->argv[g->optind], g->usage);
      exit(1);
    }
  if (nmatch == 0)
    { 
      fprintf(stderr, "No such option \"%s\".\n%s", g->argv[g->optind], g->usage);
      exit(1);
    }
  *ret_opti    = opti;
  *ret_optname = g->opt[opti].longname;
  g->optind++;	/* optind was on the option --foo; advance the counter to next argv element */

  /* Find the argument, if there is supposed to be one.
   */
  if (g->opt[opti].type != eslARG_NONE) 
    {
      if (argptr != NULL)	/* if --foo=arg syntax, then we already found it */
	{
	  *ret_optarg = argptr;
	}
      else if (g->optind >= g->argc)
	{
	  fprintf(stderr, "Option %s requires an argument\n%s", g->opt[opti].longname, usage);
	  exit(1);
	}
      else			/* "--foo 666" style, with a space */
	{
	  *ret_optarg = g->argv[g->optind];	/* assign optind as the arg, advance counter */
	  g->optind++;
	}
    }
  else  /* if there's not supposed to be an arg, but there is, die */
    {
      if (argptr != NULL) 
	{
	  fprintf(stderr, "Option %s does not take an argument\n%s", g->opt[opti].longname, usage);
	  exit(1);
	}
      *ret_optarg = NULL;
    }

  return;
}

/* process_stdopt()
 * 
 * Standard one-char options can be concatenated.
 * If we're processing one-char options, we may either be in the middle of 
 * such a string, or starting a new one.
 * 
 * Only the last optchar in a optstring may take an argument. The argument
 * is either the remainder of the argv element (if any) or if not, the
 * next argv element.
 * 
 * Examples of syntax:
 *       -a
 *       -W arg
 *       -Warg
 *       -abc
 *       -abcW arg
 *       -abcWarg
 *       
 * Process next optchar, returning option number, ptr to the optname (the full copy of it
 * in opt[] structure, not the single character in argv[]), and optarg.
 * 
 * g->optind is advanced to the next argv element; either 0, +1, or +2,
 * depending on whether we're still processing an optstring from a prev
 * optind, starting a new optstring, or reading a "-W arg" form,
 * respectively.
 * 
 * Exits w/ mesg if:
 *   1. The option doesn't exist.
 *   2. The option takes an option, but none was found.
 */
void
process_stdopt(ESL_GETOPTS *g, int *ret_opti, char **ret_optname, char **ret_optarg)
{
  /* Do we need to start a new optstring in a new argv element?
   * (as opposed to still being in an optstring from a prev parse)
   */
  if (g->optstring == NULL) 
    {      /* assign optind as an optstring; advance counter  */
      g->optstring = g->argv[g->optind]+1;
      g->optind++;		
    }

  /* figure out what option this optchar is
   */
  for (opti = 0; opti < g->nopts; opti++)
    if (g->opt[opti].stdname != NULL && *(g->optstring) == g->opt[opti].stdname[1]) break;
  if (opti == g->nopts)
    {
      fprintf(stderr, "No such option \"-%c\".\n%s", *(g->optstring), g->usage);
      exit(1);
    }
  *ret_opti    = opti;
  *ret_optname = g->opt[opti].stdname;


  /* find the argument, if there's supposed to be one */
  if (g->opt[opti].type != eslARG_NONE) 
    {
      if (*(g->optstring+1) != '\0')   /* attached argument */
	{
	  *ret_optarg = g->optstring+1;
	}
      else if (g->optind+1 < g->argc) /* unattached argument; assign optind, advance counter  */
	{
	  *ret_optarg = g->argv[g->optind];
	  g->optind++;	      
	}
      else 
	{
	  fprintf(stderr, "Option %s requires an argument\n%s", g->opt[opti].stdname, usage);
	  exit(1);
	}
      g->optstring = NULL;   /* An optchar that takes an arg must terminate an optstring. */
    }
  else  /* if there's not supposed to be an argument, then check if we're still in an optstring */
    {
      *ret_optarg = NULL;
      if (*(g->optstring+1) != '\0')   /* yup, we're still in an optstring */
	g->optstring++; 
      else
	g->optstring = NULL;           /* nope, that's it; move to next field in args */
    }
  return;
}


/* verify_integer():
 * 
 * Verify that a string <arg> can be completely converted to an
 * integer by atoi(), and that the result is within the range
 * defined by <range> (if <range> is non-NULL). Returns 
 * on successful verification; 
 *
 * 
 * to stderr, and exit the program.
 *
 * Returns ESL_OK on success.
 * Throws ESL_EINVAL if range string is bogus.
 */
static int
verify_integer(char *option, char *usage, char *arg, char *range)
{
  int   n;
  int   upper, lower;		/* upper, lower bounds */
  char *up, *lp;		
  int   geq, leq;	        /* use >=, <= instead of >, < */
  
  if (! is_integer(arg))
    {
      fprintf(stderr, "\"%s\" requires integer argument; \"%s\" is not an integer.\n%s\n", 
	      option, arg, usage);
      exit(1);
    }
  n = atoi(arg);
  
  /* Range checking, if <range> is non-NULL
   */
  if (range != NULL) {
    if (parse_rangestring(range, 'n', &lp, &geq, &up, &leq) != ESL_OK) return ESL_EINVAL;
    if (lp != NULL)
      {
	lower = atoi(lp);
	if ((geq && ! (n >= lower)) || (!geq && !(n > lower)))
	  {
	  fprintf(stderr, "\"%s\" requires integer arg in range %s;\nArg \"%s\" is invalid.\n%s\n", 
		  option, range, arg, usage);
	  exit(1);
	}
      }
    if (up != NULL) 
      {
	upper = atoi(up);
	if ((leq && ! (n <= upper)) || (!leq && !(n < upper)))
	  {
	    fprintf(stderr, "\"%s\" requires integer arg in range %s;\nArg \"%s\" is invalid.\n%s\n", 
		    option, range, arg, usage);
	    exit(1);
	  }
      }
  }
  return ESL_OK;
}



/* verify_real():
 * 
 * Verify that a string <arg> can be completely converted to a
 * floating point real by atof(), and that the result is within the range
 * defined by <range> (if <range> is non-NULL).
 *
 * Otherwise, use <option> and <usage> to format an error message
 * to stderr, and exit the program.
 *
 * Returns ESL_OK on success.
 * Throws ESL_EINVAL if range string is bogus.
 */
static int
verify_real(char *option, char *usage, char *arg, char *range)
{
  double x;
  double upper, lower;		/* upper, lower bounds */
  char  *up, *lp;		
  int    geq, leq;	        /* use >=, <= instead of >, < */
  
  if (! is_real(arg))
    {
      fprintf(stderr, "\"%s\" requires real-valued argument; \"%s\" is not a real number.\n%s\n", 
	      option, arg, usage);
      exit(1);
    }
  x = atof(arg);
  
  /* Range checking, if <range> is non-NULL
   */
  if (range != NULL) {
    if (parse_rangestring(range, 'x', &lp, &geq, &up, &leq) != ESL_OK) return ESL_EINVAL;
    if (lp != NULL)
      {
	lower = atof(lp);
	if ((geq && ! (x >= lower)) || (!geq && !(x > lower)))
	  {
	    fprintf(stderr, "\"%s\" takes real-valued arg in range %s;\nArg \"%s\" is invalid.\n%s\n", 
		    option, range, arg, usage);
	  exit(1);
	}
      }
    if (up != NULL) 
      {
	upper = atoi(up);
	if ((leq && ! (x <= upper)) || (!leq && !(x < upper)))
	  {
	    fprintf(stderr, "\"%s\" takes real-valued arg in range %s;\nArg \"%s\" is invalid.\n%s\n", 
		    option, range, arg, usage);
	    exit(1);
	  }
      }
  }
  return ESL_OK;
}


/* verify_char():
 * 
 * Verify that a string <arg> consists solely of a single
 * char, and that it is within the range defined by ASCII <range>
 * (if <range> is non-NULL). 
 *
 * Currently, range is limited to ASCII chars that can be
 * expressed as single chars in the "range" string. Could improve
 * by allowing integer ASCII codes.
 *
 * Otherwise, use <option> and <usage> to format an error message
 * to stderr, and exit the program.
 *
 * Returns ESL_OK on success.
 * Throws ESL_EINVAL if range string is bogus.
 */
static int
verify_char(char *option, char *usage, char *arg, char *range)
{
  char   c;
  char  *upper, *lower;		
  int    geq, leq;	        /* use >=, <= instead of >, < */
  
  if (strlen(arg) > 1)
    {
      fprintf(stderr, "\"%s\" requires char argument; \"%s\" is not a single character.\n%s\n", 
	      option, arg, usage);
      exit(1);
    }
  c = *arg;
  
  /* Range checking, if <range> is non-NULL
   */
  if (range != NULL) {
    if (parse_rangestring(range, 'c', &lower, &geq, &upper, &leq) != ESL_OK) return ESL_EINVAL;
    if (lower != NULL)
      {
	if ((geq && ! (c >= *lower)) || (!geq && !(c > *lower)))
	  {
	    fprintf(stderr, "\"%s\" takes char arg in range %s;\nArg \"%s\" is invalid.\n%s\n", 
		    option, range, arg, usage);
	  exit(1);
	}
      }
    if (upper != NULL) 
      {
	if ((leq && ! (c <= *upper)) || (!leq && !(c < *upper)))
	  {
	    fprintf(stderr, "\"%s\" takes char arg in range %s;\nArg \"%s\" is invalid.\n%s\n", 
		    option, range, arg, usage);
	    exit(1);
	  }
      }
  }
  return ESL_OK;
}

/* parse_rangestring():
 * 
 * Given a range definition string in one of the following forms:
 *     "a<=c<=b",  "c>=a",  "c<=b"
 * where "a" is a lower bound, "b" is an upper bound, and the =
 * signs are optional; and c is a one-character marker for the 
 * argument value ('n' for integers, 'f' for floating-point values).
 * 
 * Sets pointers to upper and lower bound strings, for parsing by
 *  atoi() or atof() as appropriate.
 * Sets geq, leq flags to TRUE if bounds are supposed to be inclusive.
 */
static int
parse_rangestring(char *range, char c, char **ret_lowerp, int *ret_geq, char **ret_upperp, int *ret_leq)
{
  char *ptr;

  *ret_geq    = *ret_leq    = FALSE;	/* 'til proven otherwise */
  *ret_lowerp = *ret_upperp = NULL;     /* 'til proven otherwise */

  if ((ptr = strchr(range, c)) == NULL) ESL_ERROR(ESL_EINVAL, "invalid range string");  
  if (ptr == range)	/* we're "n>=a" or "n<=b" form, where n came first */  
    {
      if (ptr[1] == '>') /* "n>=a" form; lower bound */
	{
	  if (ptr[2] == '=') { *ret_geq = TRUE; *ret_lowerp = ptr+3; } 
	  else *ret_lowerp = ptr+2;
	}
      else if (ptr[1] == '<') /* "n<=a" form; upper bound */
	{
	  if (ptr[2] == '=') { *ret_leq = TRUE; *ret_upperp = ptr+3; }
	  else *ret_upperp = ptr+2;
	}
      else ESL_ERROR(ESL_EINVAL, "invalid range string");
    }
  else	/* we're in a<=n<=b form; upper bound after n, lower before */
    {
      if (*(ptr+1) != '<') ESL_ERROR(ESL_EINVAL, "invalid range string");
      if (*(ptr+2) == '=') { *ret_leq = TRUE; *ret_upperp = ptr+3; } else *ret_upperp = ptr+2;

      ptr--;
      if (*ptr == '=') { *ret_geq = TRUE; ptr--; }
      if (*ptr != '<') ESL_ERROR(ESL_EINVAL, "invalid range string");
      *ret_lowerp = range;	/* start of string */
    }
  return ESL_OK;
}



/* Function: is_integer()
 * 
 * Returns TRUE if s points to something that atoi() will parse
 * completely and convert to an integer.
 */
static int
is_integer(char *s)
{
  int hex = 0;

  if (s == NULL) return 0;
  while (isspace((int) (*s))) s++;      /* skip whitespace */
  if (*s == '-' || *s == '+') s++;      /* skip leading sign */
				        /* skip leading conversion signals */
  if ((strncmp(s, "0x", 2) == 0 && (int) strlen(s) > 2) ||
      (strncmp(s, "0X", 2) == 0 && (int) strlen(s) > 2))
    {
      s += 2;
      hex = 1;
    }
  else if (*s == '0' && (int) strlen(s) > 1)
    s++;
				/* examine remainder for garbage chars */
  if (!hex)
    while (*s != '\0')
      {
	if (!isdigit((int) (*s))) return 0;
	s++;
      }
  else
    while (*s != '\0')
      {
	if (!isxdigit((int) (*s))) return 0;
	s++;
      }
  return 1;
}


/* Function: is_real()
 * 
 * Purpose:  Returns TRUE if s is a string representation
 *           of a valid floating point number, convertable
 8           by atof().
 */
static int
is_real(char *s)
{
  int gotdecimal = 0;
  int gotexp     = 0;
  int gotreal    = 0;

  if (s == NULL) return 0;

  while (isspace((int) (*s))) s++; /* skip leading whitespace */
  if (*s == '-' || *s == '+') s++; /* skip leading sign */

  /* Examine remainder for garbage. Allowed one '.' and
   * one 'e' or 'E'; if both '.' and e/E occur, '.'
   * must be first.
   */
  while (*s != '\0')
    {
      if (isdigit((int) (*s))) 	gotreal++;
      else if (*s == '.')
	{
	  if (gotdecimal) return 0; /* can't have two */
	  if (gotexp) return 0;     /* e/E preceded . */
	  else gotdecimal++;
	}
      else if (*s == 'e' || *s == 'E')
	{
	  if (gotexp) return 0;	/* can't have two */
	  else gotexp++;
	}
      else if (isspace((int) (*s)))
	break;
      s++;
    }

  while (isspace((int) (*s))) s++;         /* skip trailing whitespace */
  if (*s == '\0' && gotreal) return 1;
  else return 0;
}


/* Function: esl_getopt()
 * 
 * Purpose:  Portable command line option parsing with abbreviated
 *           option switches. Replaces UNIX getopt(). Using UNIX getopt()
 *           hinders portability to non-UNIX platforms, and getopt()
 *           is also limited to single letter options.
 *
 *           Getopt() implements a superset of UNIX getopt().
 *           All of getopt()'s single-character switch behavior
 *           is emulated, and "--" by itself terminates the options.
 *           Additionally, Getopt() provides extended switches
 *           like "--youroptionhere", and Getopt() type checks
 *           arguments.  
 * 
 *           Extended options must start with "--", as in "--option1".
 *           Normal options must start with "-", as in "-o".
 *           Normal options may be concatenated, as in "-a -b" == "-ab".
 *           
 *           See bottom of this .c file after #fdef GETOPT_TESTDRIVER
 *           for an example of calling Getopt().
 *           
 * Args:     argc  - from main(). number of elems in argv.
 *           argv  - from main(). argv[0] is the name of the command.
 *           opt   - array of opt_s structures, defining option switches
 *           nopts - number of switches in opt
 *           usage - a (possibly long) string to print if usage error.
 *           ret_optind - RETURN: the index in argv[] of the next 
 *                        valid command-line token.
 *           ret_optname- RETURN: ptr to the name of option switch 
 *                        seen, or NULL if no option was seen.
 *           ret_optarg - RETURN: ptr to the optional argument, if any;
 *                        NULL if option takes no argument.
 *                        
 * Return:   1 if a valid option was parsed.
 *           0 if no option was found, and command-line parsing is complete.
 *           Die()'s here if an error is detected.
 */
int
Getopt(int argc, char **argv, struct opt_s *opt, int nopts, char *usage,
       int *ret_optind, char **ret_optname, char **ret_optarg)
{
  int i;
  int arglen;
  int nmatch;
  static int optind   = 1;        /* init to 1 on first call  */
  static char *optptr = NULL;     /* ptr to next valid switch */
  int opti = 0;			  /* init only to silence gcc uninit warnings */



  /* We have a real option. Find which one it is.
   * We handle single letter switches "-o" separately
   * from full switches "--option", based on the "-" vs. "--"
   * prefix -- single letter switches can be concatenated
   * as long as they don't have arguments.
   */
				/* full option */
  if (optptr == NULL && strncmp(argv[optind], "--", 2) == 0)
    {
      /* Use optptr to parse argument in options of form "--foo=666"
       */
      if ((optptr = strchr(argv[optind], '=')) != NULL)
	{ *optptr = '\0'; optptr++; }

      arglen = strlen(argv[optind]);
      nmatch = 0;
      for (i = 0; i < nopts; i++)
	if (opt[i].single == FALSE && 
	    strncmp(opt[i].name, argv[optind], arglen) == 0)
	  { 
	    nmatch++;
	    opti = i;
	    if (arglen == strlen(opt[i].name)) break; /* exact match, stop now */
	  }
      if (nmatch > 1 && arglen != strlen(opt[i].name)) 
	Die("Option \"%s\" is ambiguous; please be more specific.\n%s",
	    argv[optind], usage);
      if (nmatch == 0)
	Die("No such option \"%s\".\n%s", argv[optind], usage);

      *ret_optname = opt[opti].name;

      /* Set the argument, if there is one
       */
      if (opt[opti].argtype != sqdARG_NONE) 
	{
	  if (optptr != NULL)
	    {			/* --foo=666 style */
	      *ret_optarg = optptr;
	      optptr = NULL;
	      optind++;
	    }
	  else if (optind+1 >= argc)
	    Die("Option %s requires an argument\n%s", opt[opti].name, usage);
	  else			/* "--foo 666" style */
	    {
	      *ret_optarg = argv[optind+1];
	      optind+=2;
	    }
	}
      else  /* sqdARG_NONE */
	{
	  if (optptr != NULL) 
	    Die("Option %s does not take an argument\n%s", opt[opti].name, usage);
	  *ret_optarg = NULL;
	  optind++;
	}
    }
  else				/* else, a single letter option "-o" */
    {
				/* find the option */
      if (optptr == NULL) 
	optptr = argv[optind]+1;
      for (opti = -1, i = 0; i < nopts; i++)
	if (opt[i].single == TRUE && *optptr == opt[i].name[1])
	  { opti = i; break; }
      if (opti == -1)
	Die("No such option \"%c\".\n%s", *optptr, usage);
      *ret_optname = opt[opti].name;

				/* set the argument, if there is one */
      if (opt[opti].argtype != sqdARG_NONE) 
	{
	  if (*(optptr+1) != '\0')   /* attached argument */
	    {
	      *ret_optarg = optptr+1;
	      optind++;
	    }
	  else if (optind+1 < argc) /* unattached argument */
	    {
	      *ret_optarg = argv[optind+1];
	      optind+=2;	      
	    }
	  else Die("Option %s requires an argument\n%s", opt[opti].name, usage);

	  optptr = NULL;	/* can't concatenate after an argument */
	}
      else  /* sqdARG_NONE */
	{
	  *ret_optarg = NULL;
	  if (*(optptr+1) != '\0')   /* concatenation */
	    optptr++; 
	  else
	    {
	      optind++;                /* move to next field */
	      optptr = NULL;
	    }
	}

    }

  /* Type check the argument, if there is one
   */
  if (opt[opti].argtype != sqdARG_NONE) 
    {
      if (opt[opti].argtype == sqdARG_INT && ! IsInt(*ret_optarg))
	Die("Option %s requires an integer argument\n%s",
	    opt[opti].name, usage);
      else if (opt[opti].argtype == sqdARG_FLOAT && ! IsReal(*ret_optarg))
	Die("Option %s requires a numerical argument\n%s",
	    opt[opti].name, usage);
      else if (opt[opti].argtype == sqdARG_CHAR && strlen(*ret_optarg) != 1)
	Die("Option %s requires a single-character argument\n%s",
	    opt[opti].name, usage);
      /* sqdARG_STRING is always ok, no type check necessary */
    }

  *ret_optind = optind;
  return 1;
}



#ifdef GETOPT_TESTDRIVER 
/* cc -DGETOPT_TESTDRIVER -L ~/lib/squid.linux/ getopt.c -lsquid
 */
struct opt_s OPTIONS[] = {
  { "--test1", sqdARG_INT    },
  { "--test2", sqdARG_FLOAT  },
  { "--test3", sqdARG_STRING },
  { "--test4", sqdARG_CHAR   },
  { "-a",      sqdARG_NONE   },
  { "-b",      sqdARG_INT    },
};
#define NOPTIONS (sizeof(OPTIONS) / sizeof(struct opt_s))
    
int
main(int argc, char **argv)
{
  int   optind;
  char *optarg;
  char *optname;

  while (Getopt(argc, argv, OPTIONS, NOPTIONS, "Usage/help here",
		&optind, &optname, &optarg))
    {
      printf("Option:   index: %d name: %s argument: %s\n",
	     optind, optname, optarg);
    }
  while (optind < argc)
    {
      printf("Argument: index: %d name: %s\n", optind, argv[optind]);
      optind++;
    }


}


#endif /*GETOPT_TESTDRIVER*/

/*****************************************************************  
 * @LICENSE@
 *****************************************************************/
