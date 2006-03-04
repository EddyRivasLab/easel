/* esl_getopts.c
 * Implements a somewhat more powerful command line getopt interface
 * than the standard UNIX/POSIX call.
 * 
 * SVN $Id$
 * SRE, Sat Jan  1 08:50:21 2005 [Panticosa, Spain]
 * xref STL8/p152; STL9/p5.
 */
#include <esl_config.h>

#include <stdlib.h> 
#include <stdio.h> 
#include <string.h>
#include <ctype.h>

#include <easel.h>
#include <esl_getopts.h>

/* Forward declarations of private functions.
 */
static int set_option(ESL_GETOPTS *g, int opti, char *optarg, 
		      int setby, int do_alloc);
static int get_optidx_exactly(ESL_GETOPTS *g, char *optname, int *ret_opti);
static int get_optidx_abbrev(ESL_GETOPTS *g, char *optname, int n, 
			     int *ret_opti);
static int esl_getopts(ESL_GETOPTS *g, int *ret_opti, char **ret_optarg);
static int process_longopt(ESL_GETOPTS *g, int *ret_opti, char **ret_optarg);
static int process_stdopt(ESL_GETOPTS *g, int *ret_opti, char **ret_optarg);
static int verify_type_and_range(ESL_GETOPTS *g, int i, char *val, int setby);
static int is_integer(char *s);
static int is_real(char *s);
static int verify_integer_range(char *arg, char *range);
static int verify_real_range(char *arg, char *range);
static int verify_char_range(char *arg, char *range);
static int parse_rangestring(char *range, char c, char **ret_lowerp, 
			     int *ret_geq, char **ret_upperp, int *ret_leq);
static int process_optlist(ESL_GETOPTS *g, char **ret_s, int *ret_opti);



/*****************************************************************
 * 1. The ESL_GETOPTS object
 *****************************************************************/ 

/* Function:  esl_getopts_Create()
 * Incept:    SRE, Tue Jan 11 11:24:16 2005 [St. Louis]
 *
 * Purpose:   Creates an <ESL_GETOPTS> object, given the
 *            array of valid options <opt> (NULL-element-terminated)
 *            and a (possibly long, multiline) help/usage string
 *            in <usage>. Sets default values for all config 
 *            options (as defined in <opt>).
 *
 * Returns:   ptr to the new <ESL_GETOPTS> object.
 *
 * Throws:    NULL on failure, including allocation failures or
 *            an invalid <ESL_OPTIONS> structure.
 */
ESL_GETOPTS *
esl_getopts_Create(ESL_OPTIONS *opt, char *usage)
{
  ESL_GETOPTS *g;
  int i;

  if ((g = malloc(sizeof(ESL_GETOPTS))) == NULL) 
    ESL_ERROR_NULL(eslEMEM, "malloc failed");

  g->opt       = opt;
  g->argc      = 0;
  g->argv      = NULL;
  g->usage     = usage;
  g->optind    = 1;
  g->argi      = 1;		/* number cmdline arguments 1..n */
  g->nfiles    = 0;
  g->val       = NULL;
  g->setby     = NULL;
  g->valloc    = NULL;
  g->optstring = NULL;

  /* Figure out the number of options.
   */
  g->nopts = 0;
  while (g->opt[g->nopts].name != NULL)
    g->nopts++;
  
  /* Set default values for all options.
   * Note the valloc[] setting: we only need to dup strings
   * into allocated space if the value is volatile memory, and
   * that only happens in config files; not in defaults, cmdline,
   * or environment.
   */
  if ((g->val    = malloc(sizeof(char *) * g->nopts)) == NULL) 
    { esl_getopts_Destroy(g); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  if ((g->setby  = malloc(sizeof(int)    * g->nopts)) == NULL) 
    { esl_getopts_Destroy(g); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }
  if ((g->valloc = malloc(sizeof(int)    * g->nopts)) == NULL) 
    { esl_getopts_Destroy(g); ESL_ERROR_NULL(eslEMEM, "malloc failed"); }

  for (i = 0; i < g->nopts; i++) 
    {
      g->val[i]    = g->opt[i].defval;
      g->setby[i]  = eslARG_SETBY_DEFAULT;
      g->valloc[i] = 0;	
    }

  /* Verify type/range of the defaults, even though it's
   * an application error (not user error) if they're invalid. 
   */
  for (i = 0; i < g->nopts; i++) 
    if (verify_type_and_range(g, i, g->val[i], eslARG_SETBY_DEFAULT) != eslOK)
      { esl_getopts_Destroy(g); return NULL; }

  return g;
}

/* Function:  esl_getopts_Destroy()
 * Incept:    SRE, Thu Jan 13 08:55:10 2005 [St. Louis]
 *
 * Purpose:   Free's a created <ESL_GETOPTS> object.
 *
 * Returns:   void.
 */
void
esl_getopts_Destroy(ESL_GETOPTS *g)
{
  int i;

  if (g != NULL)
    {
      if (g->val   != NULL) 
	{
	  /* A few of our vals may have been allocated.
	   */
	  for (i = 0; i < g->nopts; i++)
	    if (g->valloc[i] > 0)
	      free(g->val[i]);
	  free(g->val);
	}
      if (g->setby  != NULL) free(g->setby);
      if (g->valloc != NULL) free(g->valloc);
      free(g);
    }
}


/* Function:  esl_getopts_Dump()
 * Incept:    SRE, Tue Jan 18 09:11:39 2005 [St. Louis]
 *
 * Purpose:   Dump the state of <g> to an output stream
 *            <ofp>, often stdout or stderr.
 */
void
esl_getopts_Dump(FILE *ofp, ESL_GETOPTS *g)
{
  int i;

  fprintf(ofp, "%12s %12s %9s\n", "Option", "Setting", "Set by");
  fprintf(ofp, "------------ ------------ ---------\n");

  for (i = 0; i < g->nopts; i++)
    {
      fprintf(ofp, "%-12s ", g->opt[i].name);

      fprintf(ofp, "%-12s ", g->val[i]);
      
      if      (g->setby[i] == eslARG_SETBY_DEFAULT) fprintf(ofp, "(default) ");
      else if (g->setby[i] == eslARG_SETBY_CMDLINE) fprintf(ofp, "cmdline   ");
      else if (g->setby[i] == eslARG_SETBY_ENV)     fprintf(ofp, "environ   ");
      else if (g->setby[i] >= eslARG_SETBY_CFGFILE) fprintf(ofp, "cfgfile   ");

      fprintf(ofp, "\n");
    }
  return;
}
  

/*****************************************************************
 * 2. Setting and testing a configuration
 *****************************************************************/ 

/* Function:  esl_opt_ProcessConfigfile()
 * Incept:    SRE, Thu Jan 13 10:25:43 2005 [St. Louis]
 *
 * Purpose:   Given an open configuration file <fp> (and
 *            its name <filename>, for error reporting),
 *            parse it and set options in <g> accordingly.
 *            Anything following a <\#> in the file is a
 *            comment. Blank (or all-comment) lines are
 *            ignored. Data lines contain one option and
 *            its optional argument: for example <--foo arg>
 *            or <-a>. All option arguments are type and
 *            range checked, as specified in <g>.
 *            
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEFORMAT> on a parse or format error in the file.
 *            <eslEINVAL> if an option argument fails a type or range
 *            check, or if an option is set twice by the same config
 *            file.
 */
int
esl_opt_ProcessConfigfile(ESL_GETOPTS *g, char *filename, FILE *fp)
{
  char *buf = NULL;
  int   n   = 0;
  char *s;
  char *optname;
  char *optarg;
  char *comment;
  int   line;
  int   opti;
  int   status;

  line = 0;
  while (esl_fgets(&buf, &n, fp) == eslOK)
    {
      line++;
      optname = NULL;
      optarg  = NULL;

      /* First token is the option, e.g. "--foo"
       */
      s = buf;
      esl_strtok(&s, " \t\n", &optname, NULL);
      if (optname   == NULL) continue; /* blank line */
      if (*optname  == '#')  continue; /* comment line */
      if (*optname  != '-') {
	esl_error(eslEFORMAT, __FILE__, __LINE__,  
		  "Parse failed at line %d of cfg file %s (saw %s, not an option)\n",
		  line, filename, optname);
	return eslEFORMAT;
      }
      
      /* Second token, if present, is the arg
       */
      esl_strtok(&s, " \t\n", &optarg, NULL);
      
      /* Anything else on the line had better be a comment
       */
      esl_strtok(&s, " \t\n", &comment, NULL);
      if (comment != NULL && *comment != '#') {
	esl_error(eslEFORMAT, __FILE__, __LINE__,  
		  "Parse failed at line %d of cfg file %s (saw %s, not a comment)\n",
		  line, filename, comment);
	return eslEFORMAT;
      }
	
      /* Now we've got an optname and an optional optarg;
       * figure out what option this is.
       */
      if (get_optidx_exactly(g, optname, &opti) != eslOK) {
	esl_error(eslEFORMAT, __FILE__, __LINE__,  
		  "%s is not a recognized option (config file %s, line %d)\n",
		  optname, filename, line);
	return eslEFORMAT;
      }

      /* Set that option.
       * Pass TRUE to set_option's do_alloc flag, because our buffer
       * is volatile memory that's going away soon - set_option needs
       * to strdup the arg, not just point to it.
       */
      status = set_option(g, opti, optarg, 
			  eslARG_SETBY_CFGFILE+g->nfiles,
			  TRUE);
      if (status != eslOK) return status;
    }

  if (buf != NULL) free(buf);
  g->nfiles++;
  return eslOK;
}




/* Function:  esl_opt_ProcessEnvironment()
 * Incept:    SRE, Thu Jan 13 10:17:58 2005 [St. Louis]
 *
 * Purpose:   For any option defined in <g> that can be modified
 *            by an environment variable, check the environment
 *            and set that option accordingly. The value provided
 *            by the environment is type and range checked.
 *            When an option is turned on that has other options 
 *            toggle-tied to it, those options are turned off.
 *            An option's state may only be changed once by the
 *            environment (even indirectly thru toggle-tying);
 *            else an error is generated.
 *            
 * Returns:   <eslOK> on success, and <g> is loaded with all
 *            options specified in the environment.
 *
 * Throws:    <eslEINVAL> on any failure, including type/range
 *            check failures.
 */
int
esl_opt_ProcessEnvironment(ESL_GETOPTS *g)
{
  int   i;
  char *optarg;
  int   status;

  for (i = 0; i < g->nopts; i++)
    if (g->opt[i].envvar != NULL &&
	(optarg = getenv(g->opt[i].envvar)) != NULL)
      {
	status = set_option(g, i, optarg, eslARG_SETBY_ENV, FALSE);
	if (status != eslOK) return status;
      }
  return eslOK;
}



/* Function:  esl_opt_ProcessCmdline()
 * Incept:    SRE, Wed Jan 12 10:12:43 2005 [St. Louis]
 *
 * Purpose:   Process a command line (<argc> and <argv>), parsing out
 *            and setting application options in <g>. Option arguments
 *            are type and range checked before they are set, if type
 *            and range information was set when <g> was created.
 *            When an option is set, if it has any other options
 *            "toggle-tied" to it, those options are also turned off.
 *            
 *            Any given option can only change state (on/off) once
 *            per command line; trying to set the same option more than
 *            once generates an error.
 *            
 *            On successful return, <g> contains settings of all
 *            command line options and their option arguments, for
 *            subsequent retrieval by <esl_opt_Get*Option()>
 *            functions.  <g> also contains an <optind> state variable
 *            pointing to the next <argv[]> element that is not an
 *            option; <esl_opt_GetArgument()> uses this to retrieves
 *            command line arguments in order of appearance.
 *            
 *            
 *            The parser starts with <argv[1]> and reads <argv[]> elements
 *            in order until it reaches an element that is not an option; 
 *            at this point, all subsequent <argv[]> elements are 
 *            interpreted as arguments to the application.
 *            
 *            Any <argv[]> element encountered in the command line that
 *            starts with <- > is an option, except <- > or <-- > by
 *            themselves. <- > by itself is interpreted as a command
 *            line argument (usually meaning ``read from stdin instead
 *            of a filename''). <-- > by itself is interpreted as
 *            ``end of options''; all subsequent <argv[]> elements are
 *            interpreted as command-line arguments even if they
 *            begin with <- >. 
 *
 * Returns:   <eslOK> on success. <g> is loaded with
 *            all option settings specified on the cmdline.
 *
 * Throws:    <eslEINVAL> on any cmdline parsing problem,
 *            including option argument type/range check failures.
 */
int
esl_opt_ProcessCmdline(ESL_GETOPTS *g, int argc, char **argv)
{
  int   opti;
  char *optarg;
  int   status;

  g->argc      = argc;
  g->argv      = argv;
  g->optind    = 1;		/* start at argv[1]             */
  g->optstring = NULL;		/* not in a -abc optstring yet  */

  /* Walk through each option in the command line using esl_getopts(),
   * which advances g->optind as the index of the next argv element we need
   * to look at.
   */
  while (esl_getopts(g, &opti, &optarg) == eslOK)
    {
      status = set_option(g, opti, optarg, eslARG_SETBY_CMDLINE, FALSE);
      if (status != eslOK) return status;
    }
  return eslOK;
}



/* Function:  esl_opt_VerifyConfig()
 * Incept:    SRE, Wed Jan 12 10:21:46 2005 [St. Louis]
 *
 * Purpose:   Given a <g> that we think is fully configured now --
 *            from config file(s), environment, and command line --
 *            verify that the configuration is self-consistent:
 *            for every option that is set, make sure that any
 *            required options are also set, and that no
 *            incompatible options are set. ``Set'' means
 *            the configured value is non-NULL (including booleans),
 *            and ``not set'' means the value is NULL. (That is,
 *            we don't go by <setby>, which refers to who
 *            determined the state of an option, even if
 *            it is turned off.)
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if a required option is not set, or
 *            if an incompatible option is set.
 */
int
esl_opt_VerifyConfig(ESL_GETOPTS *g)
{
  int   i,reqi,incompati;
  char *s;
  int   status;

  /* For all options that are set (not in default configuration,
   * and turned on with non-NULL vals), 
   * verify that all their required_opts are set.
   */
  for (i = 0; i < g->nopts; i++)
    {
      if (g->setby[i] && g->val[i] != NULL)
	{
	  s = g->opt[i].required_opts;
	  while ((status = process_optlist(g, &s, &reqi)) == eslOK)
	    if (g->val[reqi] == NULL)
	      {
		esl_error(eslEINVAL, __FILE__, __LINE__,
			  "Option %s requires (or has no effect without) option(s) %s\n\n%s", 
			  g->opt[i].name, g->opt[i].required_opts, g->usage);
		return eslEINVAL;
	      }
	  if (status != eslEOD) return status;	/* non-normal end of optlist; throw status up */
	}
    }

  /* For all options that are set (turned on with non-NULL vals),
   * verify that no incompatible options are set.
   */
  for (i = 0; i < g->nopts; i++)
    {
      if (g->setby[i] && g->val[i] != NULL)
	{
	  s = g->opt[i].incompat_opts;
	  while ((status = process_optlist(g, &s, &incompati)) == eslOK)
	    if (g->val[incompati] != NULL)
	      {
		esl_error(eslEINVAL, __FILE__, __LINE__,
			  "Option %s is incompatible with option(s) %s\n\n%s", 
			  g->opt[i].name, g->opt[i].incompat_opts, g->usage);
		return eslEINVAL;
	      }
	  /* non-normal end of optlist; throw status up */
	  if (status != eslEOD) return status;	
	}
    }

  return eslOK;
}



/*****************************************************************
 * 3. Retrieving option settings and command line args
 *****************************************************************/ 

/* Function:  esl_opt_GetBooleanOption()
 * Incept:    SRE, Wed Jan 12 13:46:09 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured TRUE/FALSE value for option <optname>
 *            from <g>, leaving it in <ret_state>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetBooleanOption(ESL_GETOPTS *g, char *optname, int *ret_state)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == eslENOTFOUND)
    ESL_ERROR(eslENOTFOUND, "no such option");

  if (g->val[opti] == NULL) *ret_state = FALSE;
  else                      *ret_state = TRUE;
  return eslOK;
}

/* Function:  esl_opt_GetIntegerOption()
 * Incept:    SRE, Wed Jan 12 11:37:28 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured value for option <optname>
 *            from <g>, leaving it in <ret_n>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetIntegerOption(ESL_GETOPTS *g, char *optname, int *ret_n)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == eslENOTFOUND)
    ESL_ERROR(eslENOTFOUND, "no such option");
  *ret_n = atoi(g->val[opti]);
  return eslOK;
}
		
/* Function:  esl_opt_GetFloatOption()
 * Incept:    SRE, Wed Jan 12 13:46:27 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured value for option <optname>
 *            from <g>, leaving it in <ret_x>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetFloatOption(ESL_GETOPTS *g, char *optname, float *ret_x)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == eslENOTFOUND)
    ESL_ERROR(eslENOTFOUND, "no such option");
  *ret_x = atof(g->val[opti]);
  return eslOK;
}

/* Function:  esl_opt_GetDoubleOption()
 * Incept:    SRE, Wed Jan 12 13:46:27 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured value for option <optname>
 *            from <g>, leaving it in <ret_x>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetDoubleOption(ESL_GETOPTS *g, char *optname, double *ret_x)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == eslENOTFOUND)
    ESL_ERROR(eslENOTFOUND, "no such option");
  *ret_x = atof(g->val[opti]);
  return eslOK;
}

/* Function:  esl_opt_GetCharOption()
 * Incept:    SRE, Wed Jan 12 13:47:36 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured value for option <optname>
 *            from <g>, leaving it in <ret_c>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetCharOption(ESL_GETOPTS *g, char *optname, char *ret_c)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == eslENOTFOUND)
    ESL_ERROR(eslENOTFOUND, "no such option");
  *ret_c = *g->val[opti];
  return eslOK;
}

/* Function:  esl_opt_GetStringOption()
 * Incept:    SRE, Wed Jan 12 13:47:36 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured value for option <optname>
 *            from <g>, leaving it in <ret_s>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetStringOption(ESL_GETOPTS *g, char *optname, char **ret_s)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == eslENOTFOUND)
    ESL_ERROR(eslENOTFOUND, "no such option");
  *ret_s = g->val[opti];
  return eslOK;
}


/* Function:  esl_opt_GetCmdlineArg()
 * Incept:    SRE, Thu Jan 13 09:21:34 2005 [St. Louis]
 *
 * Purpose:   Returns ptr to the next <argv[]> element in <g> that 
 *            is a command-line argument (as opposed to an
 *            option or an option's argument). Type check it
 *            with <type> (pass <eslARG_NONE> or <eslARG_STRING> to
 *            skip type checking), and range check it with
 *            <range> (pass NULL to skip range checking).
 *
 * Returns:   ptr to next argument.
 *
 * Throws:    NULL if we run out of arguments, or an arg
 *            fails a type/range check. On failure, prints
 *            an error message complete with application help/usage 
 *            info. 
 */
char *
esl_opt_GetCmdlineArg(ESL_GETOPTS *g, int type, char *range)
{
  char *arg;
  int   status;
  
  if (g->optind >= g->argc) 
    {
      esl_error(eslEOD, __FILE__, __LINE__,
		"Not enough command line arguments.\n\n%s", g->usage);
      return NULL;
    }
  arg = g->argv[g->optind];

  /* Type and range checking.
   */
  switch (type) 
    {
    case eslARG_NONE:	/* wouldn't make sense here, but treat as unchecked. */
    case eslARG_STRING:	/* unchecked. */
      status = eslOK;
      break;

    case eslARG_INT: 
      if (! is_integer(arg))
	{
	  esl_error(eslEINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d should be an integer; got %s\n\n%s",
		    g->argi, arg, g->usage);
	  return NULL;
	}
      status = verify_integer_range(arg, range);
      if (status == eslEINVAL)
	{
	  esl_error(eslEINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d should be integer in range %s; got %s\n\n%s", 
		    g->argi, range, arg, g->usage);
	  return NULL;
	}
      break;

    case eslARG_REAL:
      if (! is_real(arg))
	{
	  esl_error(eslEINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d should be a real-valued number; got %s\n\n%s",
		    g->argi, arg, g->usage);
	  return NULL;
	}
      status = verify_real_range(arg, range);
      if (status == eslEINVAL)
	{
	  esl_error(eslEINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d takes real number in range %s; got %s\n\n%s", 
		    g->argi, range, arg, g->usage);
	  return NULL;
	}
      break;

    case eslARG_CHAR:
      if (strlen(arg) > 1)
	{
	  esl_error(eslEINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d should be a single char; got %s\n\n%s",
		    g->argi, arg, g->usage);
	  return NULL;
	}
      status = verify_char_range(arg, range);
      if (status == eslEINVAL)
	{
	  esl_error(eslEINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d takes char in range %s; got %s\n\n%s", 
		    g->argi, range, arg, g->usage);
	  return NULL;
	}
      break;

    default: ESL_ERROR_NULL(eslEINCONCEIVABLE, "no such type");
    }

  /* We have some more possible generic errors to catch...
   */
  if (status == eslESYNTAX)
    {
      esl_error(eslESYNTAX, __FILE__, __LINE__, "range string %s for arg %d is corrupt",
		range, g->argi); 
      return NULL;
    }
  else if (status != eslOK)
    ESL_ERROR_NULL(eslEINCONCEIVABLE, "unexpected error code");

  /* Normal return. Bump the argi and optind counters.
   */
  g->optind++;
  g->argi++;
  return arg;
}

/*****************************************************************
 * 4. Formatting option help
 *****************************************************************/ 

/* Function:  esl_opt_DisplayHelp()
 * Incept:    SRE, Sun Feb 26 12:36:07 2006 [St. Louis]
 *
 * Purpose:   For each option in <go>, print one line of brief
 *            documentation for it, consisting of the option name
 *            (and argument, if any) and the help string. If space
 *            allows, default values for the options (if any) are
 *            shown in brackets. If space still allows, range restrictions 
 *            for the options (if any) are shown in parentheses.
 *
 *            If <docgroup> is non-zero, help lines are only printed
 *            for options with the matching <go->opt[i].docgrouptag>.
 *            This allows the caller to group option documentation
 *            into multiple sections with different headers. To
 *            print all options in one call, pass 0 for <docgroup>.
 *            
 *            <indent> specifies how many spaces to prefix each line with.
 *            
 *            <textwidth> specifies the maximum text width for the
 *            line.  This would typically be 80 characters. Lines are
 *            not allowed to exceed this length. If a line does exceed
 *            this length, range restriction display is silently
 *            dropped (for all options). If any line still exceeds
 *            <textwidth>, the default value display is silently dropped,
 *            for all options. If any line still exceeds <textwidth>, even 
 *            though it now consists almost solely of the option name and 
 *            its help string, an <eslEINVAL> error is thrown. The
 *            caller should either shorten the help string(s) or 
 *            increase the <textwidth>.
 *
 * Returns:   <eslOK> on success.
 *
 * Throws:    <eslEINVAL> if one or more help lines are too long for
 *            the specified <textwidth>.
 */
int
esl_opt_DisplayHelp(FILE *ofp, ESL_GETOPTS *go, int docgroup, int indent,
		    int textwidth)
{
  int optwidth     = 0;		/* maximum width for "-foo <n>" options */
  int helpwidth[3] = {0,0,0};	/* 0=everything; 1=with defaults, no range; 2=help string only */
  int show_defaults;
  int show_ranges;
  int i, n;

  /* Figure out the field widths we need in the output.
   */
  for (i = 0; i < go->nopts; i++)
    if (! docgroup || docgroup == go->opt[i].docgrouptag)
      {
	n = strlen(go->opt[i].name);                /* "--foo"  */
	if (go->opt[i].type != eslARG_NONE) n += 4; /* + " <n>" */ 
	if (n > optwidth) optwidth = n;

	n = 2;                                 /* init with " : " */
	if (go->opt[i].help != NULL) 
	  n = strlen(go->opt[i].help) + 1;     /* include " " in width */
	if (n > helpwidth[2]) helpwidth[2] = n;

	if (go->opt[i].defval != NULL)
	  n += strlen(go->opt[i].defval) + 4;  /* include "  []" in width */
	if (n > helpwidth[1]) helpwidth[1] = n;

	if (go->opt[i].range != NULL)
	  n += strlen(go->opt[i].range) + 4;   /* include "  ()" in width */
	if (n > helpwidth[0]) helpwidth[0] = n;
      }

  /* Figure out what we have room for.
   */
  if (indent + optwidth + helpwidth[0] <= textwidth)
    {
      show_defaults = TRUE;
      show_ranges   = TRUE;
    }
  else if (indent + optwidth + helpwidth[1] <= textwidth)
    {
      show_defaults = TRUE;
      show_ranges   = FALSE;
    }
  else if (indent + optwidth + helpwidth[2] <= textwidth)
    {
      show_defaults = FALSE;
      show_ranges   = FALSE;
    }
  else
    ESL_ERROR(eslEINVAL, "Help line too long");


  /* Format and print the options in this docgroup.
   */
  for (i = 0; i < go->nopts; i++)
    if (! docgroup || docgroup == go->opt[i].docgrouptag)
      {
	fprintf(ofp, "%*s", indent, "");
	n = 0;
	fprintf(ofp, "%s",  go->opt[i].name);
	n += strlen(go->opt[i].name);

	switch (go->opt[i].type) {
	case eslARG_NONE:    break;
	case eslARG_INT:     fprintf(ofp, " <n>"); n += 4; break;
	case eslARG_REAL:    fprintf(ofp, " <x>"); n += 4; break;
	case eslARG_CHAR:    fprintf(ofp, " <c>"); n += 4; break;
	case eslARG_STRING:  fprintf(ofp, " <s>"); n += 4; break;
	}

	fprintf(ofp, "%*s", optwidth-n, "");
	fprintf(ofp, " :");

	if (go->opt[i].help != NULL)
	  fprintf(ofp, " %s", go->opt[i].help);
	
	if (show_defaults && go->opt[i].defval != NULL) 
	  if (go->opt[i].type != eslARG_CHAR || *(go->opt[i].defval) != '\0')
	    fprintf(ofp, "  [%s]", go->opt[i].defval);

	if (show_ranges && go->opt[i].range != NULL)
	  fprintf(ofp, "  (%s)", go->opt[i].range);

	fprintf(ofp, "\n");
      }

  /* Fini.
   */
  return eslOK;
}



/*------------------ end of the public API -----------------------*/





/*****************************************************************
 * Miscellaneous private functions 
 *****************************************************************/ 

/* set_option()
 * 
 * Turn option <opti> ON (if it's a boolean) or set its option
 * argument to <optarg>. Record that it was set by <setby> (e.g. 
 * <eslARG_SETBY_CMDLINE>). 
 * 
 * <do_alloc> is a TRUE/FALSE flag. If <arg> is a pointer to a string
 * that isn't going to go away (e.g. into argv[], or into the
 * environment) we can get away with just pointing our option's val
 * at <arg>. But if <arg> is pointing to something volatile (e.g. 
 * the line buffer as we're reading a file) then we need to
 * strdup the arg -- and remember that we did that, so we free()
 * it when we destroy the getopts object.
 */
int
set_option(ESL_GETOPTS *g, int opti, char *optarg, int setby, int do_alloc)
{
  int   arglen;
  char *where;
  char *s;
  int   togi;
  int   status;

  if       (setby == eslARG_SETBY_DEFAULT) where = "as default";
  else if  (setby == eslARG_SETBY_CMDLINE) where = "on cmdline";
  else if  (setby == eslARG_SETBY_ENV)     where = "in env";
  else if  (setby >= eslARG_SETBY_CFGFILE) where = "in cfgfile";

  /* Have we already set this option? */
  if (g->setby[opti] == setby)
    {
      esl_error(eslEINVAL, __FILE__, __LINE__, 
		"Option %s has already been set %s.\n\n%s", 
		g->opt[opti].name, where, g->usage);
      return eslEINVAL;
    }

  /* Type and range check the option argument.
   */
  if (verify_type_and_range(g, opti, optarg, setby) != eslOK)
    return eslEINVAL;
  
  /* Set the option, being careful about when val 
   * is alloc'ed vs. not.
   */
  g->setby[opti] = setby;
  if (g->opt[opti].type == eslARG_NONE)	/* booleans: any non-NULL is TRUE... */
    g->val[opti] = (char *) TRUE;       /* so 0x1 will do fine. */
  else
    {
      /* If do_alloc is FALSE or the optarg is NULL, then:
       *    - free any previous alloc; 
       *    - set the pointer.
       */
      if (!do_alloc || optarg == NULL) 
	{
	  if (g->valloc[opti] > 0) { free(g->val[opti]); g->valloc[opti] = 0; }
	  g->val[opti] = optarg;
	}
      /* else do_alloc is TRUE, so:
       *    - make sure we have enough room, either realloc'ing or malloc'ing
       *    - copy the arg.
       */
      else
	{
	  arglen = strlen(optarg);
	  if (g->valloc[opti] < arglen+1) 
	    {
	      if (g->valloc[opti] == 0)
		g->val[opti] = malloc(sizeof(char) * (arglen+1));
	      else
		g->val[opti] = realloc(g->val[opti], sizeof(char) * (arglen+1));
	      if (g->val[opti] == NULL) ESL_ERROR(eslEMEM, "allocation failed");
	      g->valloc[opti] = arglen+1;
	    }
	  strcpy(g->val[opti], optarg);
	}
    }

  /* Unset all options toggle-tied to this one.
   */
  s = g->opt[opti].toggle_opts;
  while ((status = process_optlist(g, &s, &togi)) == eslOK)
    {
      if (g->setby[togi] == setby)
	{
	  esl_error(eslEINVAL, __FILE__, __LINE__,
		    "Options %s and %s conflict, toggling each other.\n\n%s", 
		    g->opt[togi].name, g->opt[opti].name, g->usage);
	  return eslEINVAL;
	}
	  
      g->setby[togi] = setby; /* indirectly, but still */
      if (g->valloc[togi] > 0) 	/* careful about val's that were alloc'ed */
	{ free(g->val[togi]); g->valloc[togi] = 0; }
      g->val[togi] = NULL;    /* ok for false booleans too */
    }
  if (status != eslEOD) return status; /* not a normal end of optlist */
  return eslOK;
}

/* get_optidx_exactly():
 * 
 * Find option named <optname> in <g>; set <ret_opti> to be
 * the index of the option, and return eslOK. <optname>
 * must exactly match one of the options in <g>.
 * 
 * If the option is not found, return eslENOTFOUND.
 */
static int
get_optidx_exactly(ESL_GETOPTS *g, char *optname, int *ret_opti)
{
  int i;

  for (i = 0; i < g->nopts; i++)
    if (strcmp(optname, g->opt[i].name) == 0) { *ret_opti = i; return eslOK; }
  return eslENOTFOUND;
}

/* get_optidx_abbrev():
 * 
 * Find option named <optname> in <g>; set <ret_opti> to be the index
 * of the option, and return eslOK. Allow <optname> to be an
 * abbreviation of one of the option names in <g>, so long as it is
 * unambiguous. If <n> is >0, the <optname> has an attached argument
 * (--foo=arg) and <n> is the # of characters before the = character
 * that we should match to find the option (5, in this example).
 * 
 * If the option is not found, return <eslENOTFOUND>.
 * If <optname> ambiguously matches two or more options in <g>,
 * return <eslEAMBIGUOUS>.
 */
static int
get_optidx_abbrev(ESL_GETOPTS *g, char *optname, int n, int *ret_opti)
{
  int nmatch = 0;
  int i;

  if (n == 0) 			/* unless we're told otherwise: */
    n = strlen(optname);	/* all of optname abbrev must match against the real name */

  for (i = 0; i < g->nopts; i++)
    if (strncmp(g->opt[i].name, optname, n) == 0)
      {
	nmatch++;
	*ret_opti = i;
	if (n == strlen(g->opt[i].name)) break; /* an exact match; can stop now */
      }
  if (nmatch > 1)  return eslEAMBIGUOUS;
  if (nmatch == 0) return eslENOTFOUND;
  return eslOK;
}
/*----------- end of private functions for retrieving option indices -------------*/



/*****************************************************************
 * Private functions for processing options out of a command line
 *****************************************************************/ 

/* esl_getopts():
 * 
 * Get the next option in argv[], and its argument (if any),
 * and pass this information back via <ret_opti> (index of
 * next option) and <ret_optarg).
 * 
 * Return <eslOK> on success, <eslEOD> if we're out of
 * options. 
 * 
 * Throws <eslEINVAL> if something's wrong with the options.
 */
static int
esl_getopts(ESL_GETOPTS *g, int *ret_opti, char **ret_optarg)
{
  int   opti;

  *ret_optarg  = NULL; 

  /* Check to see if we've run out of options.
   * A '-' by itself is an argument (e.g. "read from stdin"), not an option.
   */
  if (g->optstring == NULL &&
      (g->optind >= g->argc || g->argv[g->optind][0] != '-' || strcmp(g->argv[g->optind], "-") == 0))
    return eslEOD; 		/* normal end-of-data (end of options) return  */

  /* Check to see if we're being told that this is the end
   * of the options with the special "--" flag.
   */
  if (g->optstring == NULL &&
      strcmp(g->argv[g->optind], "--") == 0)
    { 
      g->optind++;
      return eslEOD; 		/* also a normal end-of-data return */
    }

  /* We have an option: an argv element that starts with -, but is
   * not "-" or "--".
   * 
   * We know the strncmp() test is ok for 2 chars, because if the option was
   * 1 char, we would've already caught it above (either it's a bare "-"
   * or it's a single non-option char, and in either case it's not an option
   * and we returned eslEOD.
   * 
   * Watch out for the case where we're in the middle of a concatenated optstring
   * of single-letter options, a la -abc
   */
  if (g->optstring == NULL && strncmp(g->argv[g->optind], "--", 2) == 0)
    process_longopt(g, &opti, ret_optarg);
  else 
    process_stdopt(g, &opti, ret_optarg);

  /* Normal return.
   */
  *ret_opti = opti;
  return eslOK;
}

/* process_longopt():
 *
 * optind is sitting on a long option, w/ syntax of one of these forms:
 *       --foo        
 *       --foo arg
 *       --foo=arg
 * (GNU getopt long option syntax.)
 * 
 * Allow unambiguous abbreviations of long options when matching;
 * e.g. --foo is ok for matching a long option --foobar.
 * 
 * Returns eslOK on success, returning the option number through
 * <ret_opti>, and a ptr to its argument through <ret_optarg> (or NULL
 * if this option takes no argument.) Internally, g->optind is
 * advanced to next argv element (+1, +2, +1, respectively, for --foo,
 * --foo arg, --foo=arg).
 *
 * Throws eslEINVAL and issues a useful error mesg if:
 *   1. Option can't be found in opt[].
 *   2. Option abbreviation is ambiguous, matching opt[] nonuniquely.
 *   3. Option takes an argument, but no argument found.
 *   4. Option does not take an argument, but one was provided by =arg syntax.
 * 
 */
static int
process_longopt(ESL_GETOPTS *g, int *ret_opti, char **ret_optarg)
{
  int   opti;		/* option number found                               */
  char *argptr;		/* ptr to arg in --foo=arg syntax                    */
  int   n;		/* length of argv elem's option part (up to = or \0) */
  int   status;

  /* Deal with options of syntax "--foo=arg" w/o modifying argv.
   */
  if ((argptr = strchr(g->argv[g->optind], '=')) != NULL)
    { n = argptr - g->argv[g->optind]; argptr++; } /* bump argptr off the = to the arg itself */
  else
    { n = strlen(g->argv[g->optind]); } /* and argptr == NULL from above. */

  /* Figure out which option this is.
   * The trick here is to allow abbreviations, and identify
   * ambiguities while we're doing it. (GNU getopt allows abbrevs.)
   */
  status = get_optidx_abbrev(g, g->argv[g->optind], n, &opti);
  if (status == eslEAMBIGUOUS)
    {
      esl_error(eslEINVAL, __FILE__, __LINE__,
		"Abbreviated option \"%s\" is ambiguous.\n\n%s", g->argv[g->optind], g->usage);
      return eslEINVAL;
    }
  if (status == eslENOTFOUND)
    { 
      esl_error(eslEINVAL, __FILE__, __LINE__,
		"No such option \"%s\".\n\n%s", g->argv[g->optind], g->usage);
      return eslEINVAL;
    }

  *ret_opti    = opti;
  g->optind++;	/* optind was on the option --foo; advance the counter to next argv element */

  /* Find the argument, if there is supposed to be one.
   */
  if (g->opt[opti].type != eslARG_NONE) 
    {
      if (argptr != NULL)	/* if --foo=arg syntax, then we already found it */
	*ret_optarg = argptr;
      else if (g->optind >= g->argc)
	{
	  esl_error(eslEINVAL, __FILE__, __LINE__,
		    "Option %s requires an argument\n\n%s", g->opt[opti].name, g->usage);
	  return eslEINVAL;
	}
      else			/* "--foo 666" style, with a space */
	*ret_optarg = g->argv[g->optind++];	/* assign optind as the arg, advance counter */
    }
  else  /* if there's not supposed to be an arg, but there is, then die */
    {
      if (argptr != NULL) 
	{
	  esl_error(eslEINVAL, __FILE__, __LINE__,
		    "Option %s does not take an argument\n\n%s", g->opt[opti].name, g->usage);
	  return eslEINVAL;
	}
      *ret_optarg = NULL;
    }

  return eslOK;
}

/* process_stdopt():
 * 
 * Either we're in the middle of working on an optstring (and optind
 * is sitting on the next argv element, which may be an argument of
 * the last char in the optstring), or optind is sitting on a "-"
 * option and we should start working on a new optstring. That is,
 * we're dealing with standard one-char options, which may be
 * concatenated into an optstring.
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
 * Process next optchar; return eslOK on success, returning option
 * number through <ret_opti> and a pointer to any argument through
 * <ret_optarg>. Internally, optind is advanced to the next argv element;
 * either 0, +1, or +2, depending on whether we're still processing an
 * optstring from a prev optind, starting a new optstring, or reading
 * a "-W arg" form, respectively.
 * 
 * Throws <eslEINVAL> and issues helpful error mesg if:
 *   1. The option doesn't exist.
 *   2. The option takes an option, but none was found.
 */
static int
process_stdopt(ESL_GETOPTS *g, int *ret_opti, char **ret_optarg)
{
  int opti;

  /* Do we need to start a new optstring in a new argv element?
   * (as opposed to still being in an optstring from a prev parse)
   */
  if (g->optstring == NULL)     
    g->optstring = g->argv[g->optind++]+1; /* init optstring on first option char; advance optind */

  /* Figure out what option this optchar is
   */
  for (opti = 0; opti < g->nopts; opti++)
    if (*(g->optstring) == g->opt[opti].name[1]) break;	/* this'll also fail appropriately for long opts. */
  if (opti == g->nopts)
    {
      esl_error(eslEINVAL, __FILE__, __LINE__,
		"No such option \"-%c\".\n\n%s", *(g->optstring), g->usage);
      return eslEINVAL;
    }
  *ret_opti    = opti;

  /* Find the argument, if there's supposed to be one */
  if (g->opt[opti].type != eslARG_NONE) 
    {
      if (*(g->optstring+1) != '\0')   /* attached argument case, a la -Warg */
	*ret_optarg = g->optstring+1;
      else if (g->optind < g->argc)  /* unattached argument; assign optind, advance counter  */
	*ret_optarg = g->argv[g->optind++];
      else 
	{
	  esl_error(eslEINVAL, __FILE__, __LINE__,
		    "Option %s requires an argument\n\n%s", g->opt[opti].name, g->usage);
	  return eslEINVAL;
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
  return eslOK;
}
/*----------- end of private functions for processing command line options -------------*/




/*****************************************************************
 * Private functions for type and range checking.
 *****************************************************************/

/* verify_type_and_range():
 *
 * Implementation of type and range checking for options.
 *
 * Given a value <val> (as a string) for option <i> in the option
 * object <g>, verify that <val> satisfies the appropriate type and
 * range.  If successful, return <eslOK>. 
 * 
 * The <setby> flag is used to help format useful error messages,
 * by saying who was responsible for a bad <val>.
 *
 * Returns: <eslOK> on success.
 *
 * Throws:  <eslEINVAL>:         <val> is not the right type.
 *          <eslERANGE>:         <val> is out of allowed range.
 *          <eslESYNTAX>:        a range string format was bogus.
 *          <eslEINCONCEIVABLE>: "can't happen" internal errors.
 */
static int
verify_type_and_range(ESL_GETOPTS *g, int i, char *val, int setby)
{
  char *where;
  int   status;

  if       (setby == eslARG_SETBY_DEFAULT) where = "as default";
  else if  (setby == eslARG_SETBY_CMDLINE) where = "on cmdline";
  else if  (setby == eslARG_SETBY_ENV)     where = "in env";
  else if  (setby >= eslARG_SETBY_CFGFILE) where = "in cfgfile";

  switch (g->opt[i].type) {

  case eslARG_NONE:	
    /* treat as unchecked, because val may be "on", 0x1, "true", etc.:
     * any non-NULL ptr means on, and NULL means off.
     */
    break;

  case eslARG_INT:
    if (! is_integer(val))
      {
	esl_error(eslEINVAL, __FILE__, __LINE__, 
		  "option %s takes integer arg; got %s %s\n\n%s", 
		  g->opt[i].name, val, where, g->usage);
	return eslEINVAL;
      }

    status = verify_integer_range(val, g->opt[i].range);
    if (status == eslERANGE)
      {
	esl_error(eslERANGE, __FILE__, __LINE__, 
		  "option %s takes integer arg in range %s; got %s %s\n\n%s", 
		  g->opt[i].name, g->opt[i].range, val, where, g->usage);
	return eslEINVAL;
      }
    else if (status == eslESYNTAX) /* eslESYNTAX, or anything else */
      {
	esl_error(eslESYNTAX, __FILE__, __LINE__, "range string %s for option %s is corrupt",
		  g->opt[i].range, g->opt[i].name); 
	return eslESYNTAX;
      }
    else if (status != eslOK) ESL_ERROR(eslEINCONCEIVABLE, "unexpected error code");
    break;

  case eslARG_REAL:
    if (! is_real(val))
      {
	esl_error(eslEINVAL, __FILE__, __LINE__, 
		  "option %s takes real-valued arg; got %s %s\n\n%s", 
		  g->opt[i].name, val, where, g->usage);
	return eslEINVAL;
      }

    status = verify_real_range(val, g->opt[i].range);
    if (status == eslERANGE)
      {
	esl_error(eslERANGE, __FILE__, __LINE__, 
		  "option %s takes real-valued arg in range %s; got %s %s\n\n%s", 
		  g->opt[i].name, g->opt[i].range, val, where, g->usage);
	return eslERANGE;
      }
    else if (status == eslESYNTAX)
      {
	esl_error(eslESYNTAX, __FILE__, __LINE__, "range string %s for option %s is corrupt",
		  g->opt[i].range, g->opt[i].name); 
	return eslESYNTAX;
      }
    else if (status != eslOK) ESL_ERROR(eslEINCONCEIVABLE, "unexpected error code");
    break;

  case eslARG_CHAR:
    if (strlen(g->val[i]) > 1)
      {
	esl_error(eslEINVAL, __FILE__, __LINE__, 
		  "option %s takes char arg; got %s %s\n\n%s", 
		  g->opt[i].name, val, where, g->usage);
	return eslEINVAL;
      }
    status = verify_char_range(val, g->opt[i].range);
    if (status == eslERANGE)
      {
	esl_error(eslERANGE, __FILE__, __LINE__, 
		  "option %s takes char arg in range %s; got %s %s\n\n%s", 
		  g->opt[i].name, g->opt[i].range, val, where, g->usage);
	return eslERANGE;
      }
    else if (status == eslESYNTAX)
      {
	esl_error(eslESYNTAX, __FILE__, __LINE__, "range string %s for option %s is corrupt",
		  g->opt[i].range, g->opt[i].name); 
	return eslESYNTAX;
      }
    else if (status != eslOK) ESL_ERROR(eslEINCONCEIVABLE, "unexpected error code");
    break;

  case eslARG_STRING: /* unchecked type. */
    if (g->opt[i].range != NULL)
      {
	esl_error(eslEINVAL, __FILE__, __LINE__, 
		  "option %s takes a string arg that can't be range checked",  g->opt[i].name);
	return eslEINVAL;
      }
    break;			
    
  default: ESL_ERROR(eslEINVAL, "no such argument type");
  }

  return eslOK;
}

/* Function: is_integer()
 * 
 * Returns TRUE if <s> points to something that atoi() will parse
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


/* is_real()
 * 
 * Returns TRUE if <s> is a string representation
 * of a valid floating point number, convertable
 * by atof().
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


/* verify_integer_range():
 * 
 * Returns <eslOK> if the string <arg>, when converted 
 * to an integer by atoi(), gives a value that lies within
 * the given <range>, if <range> is non-NULL. (If
 * <range> is NULL, there is no constraint on the range
 * of this <arg>, so return TRUE.) Else, <arg> does
 * not lie in the <range>; return <eslERANGE>. If
 * <range> is misformatted, return <eslESYNTAX>, so caller
 * can print a reasonable error message.
 * 
 * Range must be in one of three formats, matched
 * by these regexps (though regexps aren't used by the
 * parser):
 *        n>=?(\d+)           lower bound 
 *        n<=?(\d+)           upper bound
 *        (\d+)<=?n<=?(\d+)   lower and upper bound
 * Optional = signs indicate whether a bound is 
 * inclusive or not. The "n" character indicates the
 * given integer value.       
 * 
 * Returns:  <eslOK>:      <arg> is within allowed <range>.
 *           <eslERANGE>:  <arg> is not within allowed <range>.
 *           <eslESYNTAX>: something wrong with <range> string.
 */
static int
verify_integer_range(char *arg, char *range)
{
  int   n;
  int   upper, lower;		/* upper, lower bounds */
  char *up, *lp;		
  int   geq, leq;	        /* use >=, <= instead of >, < */
  
  if (range == NULL) return eslOK;
  n = atoi(arg);

  if (parse_rangestring(range, 'n', &lp, &geq, &up, &leq) != eslOK) 
    return eslESYNTAX;

  if (lp != NULL)
    {
      lower = atoi(lp);
      if ((geq && ! (n >= lower)) || (!geq && !(n > lower)))
	return eslERANGE;
    }

  if (up != NULL) 
    {
      upper = atoi(up);
      if ((leq && ! (n <= upper)) || (!leq && !(n < upper)))
	return eslERANGE;
    }
  return eslOK;
}



/* verify_real_range():
 * 
 * Verify that a string <arg>, when converted to a
 * double-precision real by atof(), gives a value that lies
 * within the range defined by <range>. If <range> is NULL,
 * there is no range constraint, and any <arg> is valid.
 *
 * Returns:  <eslOK>:      <arg> is within allowed <range>.
 *           <eslERANGE>:  <arg> is not within allowed <range>.
 *           <eslESYNTAX>: something wrong with <range> string.
 */
static int
verify_real_range(char *arg, char *range)
{
  double x;
  double upper, lower;		/* upper, lower bounds */
  char  *up, *lp;		
  int    geq, leq;	        /* use >=, <= instead of >, < */
  
  if (range == NULL) return eslOK;
  x = atof(arg);
  
  if (parse_rangestring(range, 'x', &lp, &geq, &up, &leq) != eslOK) 
    return eslESYNTAX;

  if (lp != NULL)
    {
      lower = atof(lp);
      if ((geq && ! (x >= lower)) || (!geq && !(x > lower)))
	return eslERANGE;
    }

  if (up != NULL) 
    {
      upper = atof(up);
      if ((leq && ! (x <= upper)) || (!leq && !(x < upper)))
	return eslERANGE;
    }
  return eslOK;
}


/* verify_char_range():
 * 
 * Verify that a string <arg>, when interpreted as a single
 * char argument, is a character that lies within the defined
 * <range>. If <range> is NULL, there is no range constraint,
 * and any <arg> is valid.
 *
 * Currently, <range> expression is limited to ASCII chars that can be
 * expressed as single chars. Could improve by allowing integer ASCII
 * codes, or backslash escapes.
 *
 * Returns:  <eslOK>:      <arg> is within allowed <range>.
 *           <eslERANGE>:  <arg> is not within allowed <range>.
 *           <eslESYNTAX>: something wrong with <range> string.
 */
static int
verify_char_range(char *arg, char *range)
{
  char   c;
  char  *upper, *lower;		
  int    geq, leq;	        /* use >=, <= instead of >, < */
  
  if (range == NULL) return eslOK;
  c = *arg;

  if (parse_rangestring(range, 'c', &lower, &geq, &upper, &leq) != eslOK) 
    return eslESYNTAX;

  if (lower != NULL)
    {
      if ((geq && ! (c >= *lower)) || (!geq && !(c > *lower)))
	return eslERANGE;
    }

  if (upper != NULL) 
    {
      if ((leq && ! (c <= *upper)) || (!leq && !(c < *upper)))
	return eslERANGE;
    }
  return eslOK;
}

/* parse_rangestring():
 * 
 * Given a range definition string in one of the following forms:
 *        c>=?(\d+)           lower bound 
 *        c<=?(\d+)           upper bound
 *        (\d+)<=?c<=?(\d+)   lower and upper bound
 * where <c> is a one-character marker expected for the 
 * argument type ('n' for integers, 'f' for floating-point values,
 * 'c' for characters).
 * 
 * Sets pointers to upper and lower bound strings, for parsing by
 * atoi() or atof() as appropriate.
 * Sets geq, leq flags to TRUE if bounds are supposed to be inclusive.
 * 
 * Returns <eslOK> on success, <eslESYNTAX> if the range string
 * is invalid. No errors are thrown here, so caller can format a
 * useful error message if range string is bogus.
 */
static int
parse_rangestring(char *range, char c, char **ret_lowerp, int *ret_geq, char **ret_upperp, int *ret_leq)
{
  char *ptr;

  *ret_geq    = *ret_leq    = FALSE;	/* 'til proven otherwise */
  *ret_lowerp = *ret_upperp = NULL;     /* 'til proven otherwise */

  if ((ptr = strchr(range, c)) == NULL) return eslESYNTAX;
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
      else return eslESYNTAX;
    }
  else	/* we're in a<=n<=b form; upper bound after n, lower before */
    {
      if (*(ptr+1) != '<') return eslESYNTAX;
      if (*(ptr+2) == '=') { *ret_leq = TRUE; *ret_upperp = ptr+3; } else *ret_upperp = ptr+2;

      ptr--;
      if (*ptr == '=') { *ret_geq = TRUE; ptr--; }
      if (*ptr != '<') return eslESYNTAX;
      *ret_lowerp = range;	/* start of string */
    }
  return eslOK;
}

/*-------------- end of private type/range-checking functions ----------------*/




/*****************************************************************
 * Private functions for checking optlists (toggles, required options, 
 * and incompatible options
 *****************************************************************/

/* process_optlist():
 *
 * Given a pointer <s> to the next option name in 
 * a comma-delimited list, figure out what option
 * this is; set <ret_opti> to its index. If another
 * option remains in the optlist, reset <s> to
 * the start of it, for the next call to process_optlist().
 * If no options remain after this one, reset <s> to NULL.
 * 
 * Returns: <eslOK> if an option has been successfully parsed
 *          out of the list and <ret_opti> is valid;
 *          <eslEOD> if no more option remains (<s> is NULL,
 *          or points to a \0).
 *          
 * Throws:  <eslEINVAL> if an option in the list isn't
 *          recognized.         
 */
static int 
process_optlist(ESL_GETOPTS *g, char **ret_s, int *ret_opti)
{
  char *s;
  int   i;
  int   n;
  
  if ((s = *ret_s) == NULL) return eslEOD;
  if (*s == '\0')           return eslEOD;

  n = strcspn(s, ",");

  /* a little weak here; we're only matching a n-long prefix, so we're
   * not going to catch the case where the optlist contains a
   * truncated, ambiguous optionname.  but optlists are not user
   * input, so the answer to this problem is "don't do that".
   */
  for (i = 0; i < g->nopts; i++)
    if (strncmp(g->opt[i].name, s, n) == 0) break;
  if (i == g->nopts) 
    ESL_ERROR(eslEINVAL, "no such option");

  *ret_opti = i;

  if (s[n] == ',') *ret_s = s+n+1; 
  else             *ret_s = NULL;

  return eslOK;
}

/*------- end of private functions for processing optlists -----------*/


/*****************************************************************
 * Code examples, and a test driver
 *****************************************************************/

/* The starting example of "standard" getopts behavior, without
 * any of the bells and whistles.
 *   gcc -g -Wall -o example -I. -DeslGETOPTS_EXAMPLE esl_getopts.c easel.c
 */
#ifdef eslGETOPTS_EXAMPLE
/*::cexcerpt::getopts_example::begin::*/
#include <stdio.h>
#include <easel.h>
#include <esl_getopts.h>

static ESL_OPTIONS options[] = {
  /* name        type       def   env  range toggles reqs incomp help                       docgroup*/
  { "-h",     eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "show help and usage",            0},
  { "-a",     eslARG_NONE,  FALSE, NULL, NULL, NULL, NULL, NULL, "a boolean switch",               0},
  { "-b",     eslARG_NONE,"default",NULL,NULL, NULL, NULL, NULL, "another boolean switch",         0},
  { "-n",     eslARG_INT,     "0", NULL, NULL, NULL, NULL, NULL, "an integer argument",            0},
  { "-x",     eslARG_REAL,  "1.0", NULL, NULL, NULL, NULL, NULL, "a real-valued argument",         0},
  { "--file", eslARG_STRING, NULL, NULL, NULL, NULL, NULL, NULL, "long option, with filename arg", 0},
  { "--char", eslARG_CHAR,     "", NULL, NULL, NULL, NULL, NULL, "long option, with character arg",0},
  {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
static char usage[] = "Usage: ./example [-options] <arg>";

int
main(int argc, char **argv)
{
  ESL_GETOPTS *go;
  int          show_help;
  int          opt_a, opt_b, opt_n;	
  float        opt_x;		
  char        *opt_file;
  char         opt_char;
  char        *arg;

  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);
  esl_opt_GetBooleanOption(go, "-h", &show_help);
  if (show_help) {
    puts(usage); 
    puts("\n  where options are:");
    esl_opt_DisplayHelp(stdout, go, 0, 2, 80); /* 0=all docgroups; 2=indentation; 80=width */
    return 0;
  }
  esl_opt_GetBooleanOption(go, "-a",     &opt_a);
  esl_opt_GetBooleanOption(go, "-b",     &opt_b);
  esl_opt_GetIntegerOption(go, "-n",     &opt_n);
  esl_opt_GetFloatOption(go,   "-x",     &opt_x);
  esl_opt_GetStringOption(go,  "--file", &opt_file);
  esl_opt_GetCharOption(go,    "--char", &opt_char);

  if (esl_opt_ArgNumber(go) != 1) {
    puts("Incorrect number of command line arguments.");
    puts(usage);
    return 1;
  }
  arg = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL);

  printf("Option -a:      %s\n", opt_a ? "on" : "off");
  printf("Option -b:      %s\n", opt_b ? "on" : "off");
  printf("Option -n:      %d\n", opt_n);
  printf("Option -x:      %f\n", opt_x);
  printf("Option --file:  %s\n", opt_file == NULL? "(null)" : opt_file);
  printf("Option --char:  %c\n", opt_char);
  printf("Cmdline arg:    %s\n", arg);

  esl_getopts_Destroy(go);
  return 0;
}
/*::cexcerpt::getopts_example::end::*/
#endif /*eslGETOPTS_EXAMPLE*/



#ifdef eslGETOPTS_TESTDRIVE 
/* gcc -g -Wall -o test -I. -DeslGETOPTS_TESTDRIVE esl_getopts.c easel.c
 */
#include <stdlib.h>
#include <stdio.h>

#include <easel.h>
#include <esl_getopts.h>

/*::cexcerpt::getopts_bigarray::begin::*/
static ESL_OPTIONS options[] = {
  /* name    type        default env_var  range toggles req  incompat help                  docgroup */
 { "-a",     eslARG_NONE, FALSE,"FOOTEST",NULL,  NULL,  NULL,  NULL,  "toggle b on",               1 },
 { "-b",     eslARG_NONE, FALSE,  NULL,   NULL,"--no-b",NULL,  NULL,  "toggle a on",               1 },
 { "--no-b", eslARG_NONE,"TRUE",  NULL,   NULL,   "-b", NULL,  NULL,  "toggle b off",              1 },
 { "-c",     eslARG_CHAR,   "x",  NULL,"a<=c<=z",NULL,  NULL,  NULL,  "character arg",             2 },
 { "-n",     eslARG_INT,    "0",  NULL,"0<=n<10",NULL,  NULL,  NULL,  "integer arg",               2 },
 { "-x",     eslARG_REAL, "0.8",  NULL, "0<x<1", NULL,  NULL,  NULL,  "real-value arg",            2 },
 { "--lowx", eslARG_REAL, "1.0",  NULL,   "x>0", NULL,  NULL,  NULL,  "real arg with lower bound", 2 },
 { "--hix",  eslARG_REAL, "0.9",  NULL,   "x<1", NULL,  NULL,  NULL,  "real arg with upper bound", 2 },
 { "--lown", eslARG_INT,   "42",  NULL,   "n>0", NULL,"-a,-b", NULL,  "int arg with lower bound",  2 },
 { "--hin",  eslARG_INT,   "-1",  NULL,   "n<0", NULL,  NULL,"--no-b","int arg with upper bound",  2 },
 { "--host", eslARG_STRING, "","HOSTTEST",NULL,  NULL,  NULL,  NULL,  "string arg with env var",   3 },
 {  0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
};
/*::cexcerpt::getopts_bigarray::end::*/

static char usage[] = "\
Usage: test [-options] <arg>\n\
";

    
int
main(void)
{
  ESL_GETOPTS *go;
  int   state;
  char  c;
  int   n;
  float x;
  char *s;
  char *file1 = "test.f1";
  char *file2 = "test.f2";
  FILE *f1, *f2;

  /* Declare a "command line" internally.
   */
  int   argc = 9;		/* progname; 5 options; 2 args */
  char *argv[] = { "progname", "-bc", "y", "-n9", "--hix=0.0", "--lown", "43", "arg1", "2005" };

  /* Create a config file #1.
   */
  if ((f1 = fopen(file1, "w")) == NULL) exit(1);
  fprintf(f1, "# Test config file #1\n");
  fprintf(f1, "#\n");
  fprintf(f1, "-b\n");
  fprintf(f1, "-n 3\n");
  fprintf(f1, "-x 0.5\n");
  fclose(f1);

  /* Create config file #2.
   */
  if ((f2 = fopen(file2, "w")) == NULL) exit(1);
  fprintf(f2, "# Test config file #2\n");
  fprintf(f2, "#\n");
  fprintf(f2, "--no-b\n");
  fprintf(f2, "--hin -33\n");
  fprintf(f2, "--host www.nytimes.com\n");
  fclose(f2);

  /* Put some test vars in the environment.
   */
  putenv("FOOTEST=");
  putenv("HOSTTEST=wasp.cryptogenomicon.org");

  /* Open the test config files for reading.
   */
  if ((f1 = fopen(file1, "r")) == NULL) abort();
  if ((f2 = fopen(file2, "r")) == NULL) abort();

  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessConfigfile(go, file1, f1);
  esl_opt_ProcessConfigfile(go, file2, f2);
  esl_opt_ProcessEnvironment(go);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);

  fclose(f1);
  fclose(f2);

  /* Option -a is set ON by an environment variable.
   */
  esl_opt_GetBooleanOption(go, "-a", &state);
  if (state != TRUE) abort();

  /* Option -b is overridden twice, and ends up being
   * ON because of the command line.
   */
  esl_opt_GetBooleanOption(go, "-b", &state);
  if (state != TRUE) abort();

  /* Option --no-b had better be off, therefore.
   */
  esl_opt_GetBooleanOption(go, "--no-b", &state);
  if (state != FALSE) abort();

  /* Option -c gets set to y by the command line, in
   * an optstring.
   */
  esl_opt_GetCharOption(go, "-c", &c);
  if (c != 'y') abort();

  /* Option -n gets set in a cfgfile, then overridden 
   * as a linked arg on the command line and set to 9.
   */
  esl_opt_GetIntegerOption(go, "-n", &n);
  if (n != 9) abort();

  /* Option -x is set from cfgfile #1 to 0.5.
   */
  esl_opt_GetFloatOption(go, "-x", &x);
  if (x != 0.5) abort();

  /* Option --lowx stays in default, 1.0.
   */
  esl_opt_GetFloatOption(go, "--lowx", &x);
  if (x != 1.0) abort();

  /* Option --hix is set to 0 on the command line in --arg=x format
   */
  esl_opt_GetFloatOption(go, "--hix", &x);
  if (x != 0.0) abort();

  /* Option --lown is set to 43 on the command line in --arg x format.
   * It requires -a and -b to be ON, which they should be.
   */
  esl_opt_GetIntegerOption(go, "--lown", &n);
  if (n != 43) abort();

  /* Option --hin is set to -33 in cfg file #2.
   * It requires --no-b to be OFF, which it should be.
   */
  esl_opt_GetIntegerOption(go, "--hin", &n);
  if (n != -33) abort();

  /* Option --host is set in cfg file #2,
   * then overridden in the environment.
   */
  esl_opt_GetStringOption(go, "--host", &s);
  if (strcmp(s, "wasp.cryptogenomicon.org") != 0) abort();

  /* Now the two remaining argv[] elements are the command line args
   */
  if (esl_opt_ArgNumber(go) != 2) abort();

  s = esl_opt_GetCmdlineArg(go, eslARG_STRING, NULL);
  if (strcmp(s, "arg1") != 0) abort();

  s = esl_opt_GetCmdlineArg(go, eslARG_INT, "2005<=n<=2005");
  if (strcmp(s, "2005") != 0) abort();

  esl_getopts_Destroy(go);
  exit(0);
}

#endif /*eslGETOPTS_TESTDRIVE*/
/*-------------- end of examples, test driver ********************/

/*****************************************************************  
 * @LICENSE@
 *****************************************************************/



