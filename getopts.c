/* getopts.c
 * Implements a somewhat more powerful command line getopt interface
 * than the standard UNIX/POSIX call.
 * 
 * SVN $Id$
 * SRE, Sat Jan  1 08:50:21 2005 [Panticosa, Spain]
 * xref STL8/p152; STL9/p5.
 */

#include <stdlib.h> 
#include <stdio.h> 
#include <string.h>
#include <ctype.h>

#include <easel/easel.h>
#include <easel/getopts.h>

/* Forward declarations of private functions.
 */
static int get_optidx_exactly(ESL_GETOPTS *g, char *optname, int *ret_opti);
static int get_optidx_abbrev(ESL_GETOPTS *g, char *optname, int n, int *ret_opti);
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
 *            an invalid <opts> structure.
 */
ESL_GETOPTS *
esl_getopts_Create(ESL_OPTIONS *opt, char *usage)
{
  ESL_GETOPTS *g;
  int i;

  if ((g = malloc(sizeof(ESL_GETOPTS))) == NULL) goto FAILURE;
  g->opt       = opt;
  g->argc      = 0;
  g->argv      = NULL;
  g->usage     = usage;
  g->optind    = 1;
  g->argi      = 1;		/* number cmdline arguments 1..n */
  g->nfiles    = 0;
  g->val       = NULL;
  g->setby     = NULL;
  g->optstring = NULL;

  /* Figure out the number of options.
   */
  g->nopts = 0;
  while (g->opt[g->nopts].name != NULL)
    g->nopts++;
  
  /* Set default values for all options.
   */
  if ((g->val   = malloc(sizeof(char *) * g->nopts)) == NULL) goto FAILURE;
  if ((g->setby = malloc(sizeof(int)    * g->nopts)) == NULL) goto FAILURE;
  for (i = 0; i < g->nopts; i++) 
    {
      g->val[i]   = g->opt[i].defval;
      g->setby[i] = eslARG_SETBY_DEFAULT;
    }

  /* Verify type/range of the defaults, even though it's
   * an application error (not user error) if they're invalid. 
   */
  for (i = 0; i < g->nopts; i++) 
    if (verify_type_and_range(g, i, g->val[i], eslARG_SETBY_DEFAULT) != ESL_OK)
      { esl_getopts_Destroy(g); return NULL; }

  /* Normal return.
   */
  return g;

  /* Abnormal return, on any alloc failure.
   */
 FAILURE:
  esl_getopts_Destroy(g);
  ESL_ERROR_NULL(ESL_EMEM, "allocation failed");
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
  if (g != NULL)
    {
      if (g->val   != NULL) free(g->val);
      if (g->setby != NULL) free(g->setby);
      free(g);
    }
}

/* Function:  esl_opt_ProcessConfigfile()
 * Incept:    SRE, Thu Jan 13 10:25:43 2005 [St. Louis]
 *
 * Purpose:   Given an open configuration file <fp> (and
 *            its name <filename>, for error reporting),
 *            parse it and set options in <g> accordingly.
 *            Anything following a <#> in the file is a
 *            comment. Blank (or all-comment) lines are
 *            ignored. Data lines contain one option and
 *            its optional argument: for example <--foo arg>
 *            or <-a>. All option arguments are type and
 *            range checked, as specified in <g>.
 *            
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_EFORMAT> on a parse or format error in the file.
 *            <ESL_EINVAL> if an option argument fails a type or range
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
  int   togi;
  int   status;

  line = 0;
  while (esl_fgets(&buf, &n, fp) == ESL_OK)
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
	esl_error(ESL_EFORMAT, __FILE__, __LINE__,  
		  "Parse failed at line %d of cfg file %s (saw %s, not an option)\n",
		  line, filename, optname);
	return ESL_EFORMAT;
      }
      
      /* Second token, if present, is the arg
       */
      esl_strtok(&s, " \t\n", &optarg, NULL);
      
      /* Anything else on the line had better be a comment
       */
      esl_strtok(&s, " \t\n", &comment, NULL);
      if (comment != NULL && *comment != '#') {
	esl_error(ESL_EFORMAT, __FILE__, __LINE__,  
		  "Parse failed at line %d of cfg file %s (saw %s, not a comment)\n",
		  line, filename, comment);
	return ESL_EFORMAT;
      }
	
      /* Now we've got an optname and an optional optarg;
       * process 'em.
       */
      if (get_optidx_exactly(g, optname, &opti) != ESL_OK) {
	esl_error(ESL_EFORMAT, __FILE__, __LINE__,  
		  "%s is not a recognized option (config file %s, line %d)\n",
		  optname, filename, line);
	return ESL_EFORMAT;
      }

      /* Have we already set this option in this file, even indirectly? 
       * Note the idiom for treating each configfile separately.
       */
      if (g->setby[opti] == eslARG_SETBY_CFGFILE + g->nfiles)
	  {
	    esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		      "Option %s was set more than once in cfg file %s.\n",
		      optname, filename);
	    return ESL_EINVAL;
	  }

      /* Type and range check the option argument.
       */
      if (verify_type_and_range(g, opti, optarg, eslARG_SETBY_CFGFILE+g->nfiles) != ESL_OK)
	return ESL_EINVAL;

      /* Set the option. 
       */
      g->setby[opti] = eslARG_SETBY_CFGFILE + g->nfiles;
      if (g->opt[opti].type == eslARG_NONE) /* booleans: anything non-NULL is true, so 0x1 is fine */
	g->val[opti] = (char *) TRUE;
      else
	g->val[opti] = optarg;

      /* Unset all options toggle-tied to this one.
       */
      s = g->opt[opti].toggle_opts;
      while ((status = process_optlist(g, &s, &togi)) == ESL_OK)
	{
	  if (g->setby[togi] == eslARG_SETBY_CFGFILE+g->nfiles)
	    {
	      esl_error(ESL_EINVAL, __FILE__, __LINE__,
			"Options %s and %s conflict in file %s, toggling each other.", 
			  g->opt[togi].name, g->opt[opti].name, filename);
	      return ESL_EINVAL;
	    }
	  g->setby[togi] = eslARG_SETBY_CFGFILE + g->nfiles; /* indirectly, but still */
	  g->val[togi] = NULL;	/* ok for false booleans too */
	}
      if (status != ESL_EOD) return status; /* not a normal end of optlist */
    }

  if (buf != NULL) free(buf);
  g->nfiles++;
  return ESL_OK;
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
 * Returns:   <ESL_OK> on success, and <g> is loaded with all
 *            options specified in the environment.
 *
 * Throws:    <ESL_EINVAL> on any failure, including type/range
 *            check failures.
 */
int
esl_opt_ProcessEnvironment(ESL_GETOPTS *g)
{
  int   i;
  char *optarg;
  char *s;
  int   togi;
  int   status;

  for (i = 0; i < g->nopts; i++)
    if (g->opt[i].envvar != NULL &&
	(optarg = getenv(g->opt[i].envvar)) != NULL)
      {
	/* Have we already set this option in the env, even indirectly? */
	if (g->setby[i] == eslARG_SETBY_ENV)
	  {
	    esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		      "Option %s was set more than once in the environment.\n\n%s", 
		      g->opt[i].name, g->usage);
	    return ESL_EINVAL;
	  }

	/* Type and range check the option argument.
	 */
	if (verify_type_and_range(g, i, optarg, eslARG_SETBY_ENV) != ESL_OK)
	  return ESL_EINVAL;

	/* Set the option. 
	 */
	g->setby[i] = eslARG_SETBY_ENV;
	if (g->opt[i].type == eslARG_NONE) /* booleans: anything non-NULL is true, so 0x1 is fine */
	  g->val[i] = (char *) TRUE;
	else
	  g->val[i] = optarg;

	/* Unset all options toggle-tied to this one.
	 */
	s = g->opt[i].toggle_opts;
	while ((status = process_optlist(g, &s, &togi)) == ESL_OK)
	  {
	    if (g->setby[togi] == eslARG_SETBY_ENV)
	      {
		esl_error(ESL_EINVAL, __FILE__, __LINE__,
			  "Options %s and %s conflict in environment, toggling each other.\n\n%s", 
			  g->opt[togi].name, g->opt[i].name, g->usage);
		return ESL_EINVAL;
	      }
	    g->setby[togi] = eslARG_SETBY_ENV; /* indirectly, but still */
	    g->val[togi] = NULL;	/* ok for false booleans too */
	  }
	if (status != ESL_EOD) return status; /* not a normal end of optlist */
      }
  return ESL_OK;
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
 *            pointing to the next argv[] element that is not an
 *            option; <esl_opt_GetArgument()> uses this to retrieves
 *            command line arguments in order of appearance.
 *            
 *            
 *            The parser starts with argv[1] and reads argv[] elements
 *            in order until it reaches an element that is not an option; 
 *            at this point, all subsequent argv[] elements are 
 *            interpreted as arguments to the application.
 *            
 *            Any argv[] element encountered in the command line that
 *            starts with "-" is an option, except "-" or "--" by
 *            themselves. "-" by itself is interpreted as a command
 *            line argument (usually meaning "read from stdin instead
 *            of a filename"). "--" by itself is interpreted as
 *            "end of options"; all subsequent argv[] elements are
 *            interpreted as command-line arguments even if they
 *            begin with "-". 
 *
 * Returns:   <ESL_OK> on success. <g> is loaded with
 *            all option settings specified on the cmdline.
 *
 * Throws:    <ESL_EINVAL> on any cmdline parsing problem,
 *            including option argument type/range check failures.
 */
int
esl_opt_ProcessCmdline(ESL_GETOPTS *g, int argc, char **argv)
{
  int   opti;
  char *optarg;
  char *s;		/* for walking thru toggle-tied optlist */
  int   togi;		/* index of a toggle-tied option        */
  int   status;

  g->argc      = argc;
  g->argv      = argv;
  g->optind    = 1;		/* start at argv[1]             */
  g->optstring = NULL;		/* not in a -abc optstring yet  */

  /* Walk through each option in the command line using esl_getopts(),
   * which advances g->optind as the index of the next argv element we need
   * to look at.
   */
  while (esl_getopts(g, &opti, &optarg) == ESL_OK)
    {
      /* Have we already set this option? */
      if (g->setby[opti] == eslARG_SETBY_CMDLINE)
	{
	  esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		    "Option %s has already been set on command line.\n\n%s", 
		    g->opt[opti].name, g->usage);
	  return ESL_EINVAL;
	}

      /* Type and range check the option argument.
       */
      if (verify_type_and_range(g, opti, optarg, eslARG_SETBY_CMDLINE) != ESL_OK)
	return ESL_EINVAL;

      /* Set the option. 
       */
      g->setby[opti] = eslARG_SETBY_CMDLINE;
      if (g->opt[opti].type == eslARG_NONE) /* booleans: anything non-NULL is true, so 0x1 is fine */
	g->val[opti] = (char *) TRUE;
      else
	g->val[opti] = optarg;

      /* Unset all options toggle-tied to this one.
       */
      s = g->opt[opti].toggle_opts;
      while ((status = process_optlist(g, &s, &togi)) == ESL_OK)
	{
	  if (g->setby[togi] == eslARG_SETBY_CMDLINE)
	    {
	      esl_error(ESL_EINVAL, __FILE__, __LINE__,
			"Options %s and %s conflict, toggling each other.\n\n%s", 
			g->opt[togi].name, g->opt[opti].name, g->usage);
	      return ESL_EINVAL;
	    }
	  
	  g->setby[togi] = eslARG_SETBY_CMDLINE; /* indirectly, but still */
	  g->val[togi] = NULL;	/* ok for false booleans too */
	}
      if (status != ESL_EOD) return status; /* not a normal end of optlist */
    }
  return ESL_OK;
}



/* Function:  esl_opt_VerifyConfig()
 * Incept:    SRE, Wed Jan 12 10:21:46 2005 [St. Louis]
 *
 * Purpose:   Given a <g> that we think is fully configured now --
 *            from config file(s), environment, and command line --
 *            verify that the configuration is self-consistent:
 *            for every option that is set, make sure that any
 *            required options are also set, and that no
 *            incompatible options are set. "Set" means
 *            the configured value is non-NULL (including booleans),
 *            and "not set" means the value is NULL. (That is,
 *            we don't go by <setby>, which refers to who
 *            determined the state of an option, even if
 *            it is turned off.)
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_EINVAL> if a required option is not set, or
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
	  while ((status = process_optlist(g, &s, &reqi)) == ESL_OK)
	    if (g->val[reqi] == NULL)
	      {
		esl_error(ESL_EINVAL, __FILE__, __LINE__,
			  "Option %s requires (or has no effect without) option(s) %s\n\n%s", 
			  g->opt[i].name, g->opt[i].required_opts, g->usage);
		return ESL_EINVAL;
	      }
	  if (status != ESL_EOD) return status;	/* non-normal end of optlist; throw status up */
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
	  while ((status = process_optlist(g, &s, &incompati)) == ESL_OK)
	    if (g->val[incompati] != NULL)
	      {
		esl_error(ESL_EINVAL, __FILE__, __LINE__,
			  "Option %s is incompatible with option(s) %s\n\n%s", 
			  g->opt[i].name, g->opt[i].incompat_opts, g->usage);
		return ESL_EINVAL;
	      }
	  /* non-normal end of optlist; throw status up */
	  if (status != ESL_EOD) return status;	
	}
    }

  return ESL_OK;
}

/* Function:  esl_opt_GetBooleanOption()
 * Incept:    SRE, Wed Jan 12 13:46:09 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured TRUE/FALSE value for option <optname>
 *            from <g>, leaving it in <ret_state>.
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_ENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetBooleanOption(ESL_GETOPTS *g, char *optname, int *ret_state)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == ESL_ENOTFOUND)
    ESL_ERROR(ESL_ENOTFOUND, "no such option");

  if (g->val[opti] == NULL) *ret_state = FALSE;
  else                      *ret_state = TRUE;
  return ESL_OK;
}

/* Function:  esl_opt_GetIntegerOption()
 * Incept:    SRE, Wed Jan 12 11:37:28 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured value for option <optname>
 *            from <g>, leaving it in <ret_n>.
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_ENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetIntegerOption(ESL_GETOPTS *g, char *optname, int *ret_n)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == ESL_ENOTFOUND)
    ESL_ERROR(ESL_ENOTFOUND, "no such option");
  *ret_n = atoi(g->val[opti]);
  return ESL_OK;
}
		
/* Function:  esl_opt_GetFloatOption()
 * Incept:    SRE, Wed Jan 12 13:46:27 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured value for option <optname>
 *            from <g>, leaving it in <ret_x>.
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_ENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetFloatOption(ESL_GETOPTS *g, char *optname, float *ret_x)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == ESL_ENOTFOUND)
    ESL_ERROR(ESL_ENOTFOUND, "no such option");
  *ret_x = atof(g->val[opti]);
  return ESL_OK;
}

/* Function:  esl_opt_GetDoubleOption()
 * Incept:    SRE, Wed Jan 12 13:46:27 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured value for option <optname>
 *            from <g>, leaving it in <ret_x>.
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_ENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetDoubleOption(ESL_GETOPTS *g, char *optname, double *ret_x)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == ESL_ENOTFOUND)
    ESL_ERROR(ESL_ENOTFOUND, "no such option");
  *ret_x = atof(g->val[opti]);
  return ESL_OK;
}

/* Function:  esl_opt_GetCharOption()
 * Incept:    SRE, Wed Jan 12 13:47:36 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured value for option <optname>
 *            from <g>, leaving it in <ret_c>.
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_ENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetCharOption(ESL_GETOPTS *g, char *optname, char *ret_c)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == ESL_ENOTFOUND)
    ESL_ERROR(ESL_ENOTFOUND, "no such option");
  *ret_c = *g->val[opti];
  return ESL_OK;
}

/* Function:  esl_opt_GetStringOption()
 * Incept:    SRE, Wed Jan 12 13:47:36 2005 [St. Louis]
 *
 * Purpose:   Retrieves the configured value for option <optname>
 *            from <g>, leaving it in <ret_s>.
 *
 * Returns:   <ESL_OK> on success.
 *
 * Throws:    <ESL_ENOTFOUND> if <optname> isn't a registered option.
 */
int
esl_opt_GetStringOption(ESL_GETOPTS *g, char *optname, char **ret_s)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == ESL_ENOTFOUND)
    ESL_ERROR(ESL_ENOTFOUND, "no such option");
  *ret_s = g->val[opti];
  return ESL_OK;
}


/* Function:  esl_opt_GetCmdlineArg()
 * Incept:    SRE, Thu Jan 13 09:21:34 2005 [St. Louis]
 *
 * Purpose:   Returns ptr to the next argv[] element in <g> that 
 *            is a command-line argument (as opposed to an
 *            option or an option's argument). Type check it
 *            with <type> (pass eslARG_NONE or eslARG_STRING to
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
      esl_error(ESL_EOD, __FILE__, __LINE__,
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
      status = ESL_OK;
      break;

    case eslARG_INT: 
      if (! is_integer(arg))
	{
	  esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d should be an integer; got %s\n\n%s",
		    g->argi, arg, g->usage);
	  return NULL;
	}
      status = verify_integer_range(arg, range);
      if (status == ESL_EINVAL)
	{
	  esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d should be integer in range %s; got %s\n\n%s", 
		    g->argi, range, arg, g->usage);
	  return NULL;
	}
      break;

    case eslARG_REAL:
      if (! is_real(arg))
	{
	  esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d should be a real-valued number; got %s\n\n%s",
		    g->argi, arg, g->usage);
	  return NULL;
	}
      status = verify_real_range(arg, range);
      if (status == ESL_EINVAL)
	{
	  esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d takes real number in range %s; got %s\n\n%s", 
		    g->argi, range, arg, g->usage);
	  return NULL;
	}
      break;

    case eslARG_CHAR:
      if (strlen(arg) > 1)
	{
	  esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d should be a single char; got %s\n\n%s",
		    g->argi, arg, g->usage);
	  return NULL;
	}
      status = verify_char_range(arg, range);
      if (status == ESL_EINVAL)
	{
	  esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		    "cmdline arg %d takes char in range %s; got %s\n\n%s", 
		    g->argi, range, arg, g->usage);
	  return NULL;
	}
      break;

    default: ESL_ERROR_NULL(ESL_EINCONCEIVABLE, "no such type");
    }

  /* We have some more possible generic errors to catch...
   */
  if (status == ESL_ESYNTAX)
    {
      esl_error(ESL_ESYNTAX, __FILE__, __LINE__, "range string %s for arg %d is corrupt",
		range, g->argi); 
      return NULL;
    }
  else if (status != ESL_OK)
    ESL_ERROR_NULL(ESL_EINCONCEIVABLE, "unexpected error code");

  /* Normal return. Bump the argi and optind counters.
   */
  g->optind++;
  g->argi++;
  return arg;
}


/*------------------ end of the public API -----------------------*/





/*****************************************************************
 * Private functions for retrieving option indices
 *****************************************************************/ 

/* get_optidx_exactly():
 * 
 * Find option named <optname> in <g>; set <ret_opti> to be
 * the index of the option, and return ESL_OK. <optname>
 * must exactly match one of the options in <g>.
 * 
 * If the option is not found, return ESL_ENOTFOUND.
 */
static int
get_optidx_exactly(ESL_GETOPTS *g, char *optname, int *ret_opti)
{
  int i;

  for (i = 0; i < g->nopts; i++)
    if (strcmp(optname, g->opt[i].name) == 0) { *ret_opti = i; return ESL_OK; }
  return ESL_ENOTFOUND;
}

/* get_optidx_abbrev():
 * 
 * Find option named <optname> in <g>; set <ret_opti> to be the index
 * of the option, and return ESL_OK. Allow <optname> to be an
 * abbreviation of one of the option names in <g>, so long as it is
 * unambiguous. If <n> is >0, the <optname> has an attached argument
 * (--foo=arg) and <n> is the # of characters before the = character
 * that we should match to find the option (5, in this example).
 * 
 * If the option is not found, return <ESL_ENOTFOUND>.
 * If <optname> ambiguously matches two or more options in <g>,
 * return <ESL_EAMBIGUOUS>.
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
  if (nmatch > 1)  return ESL_EAMBIGUOUS;
  if (nmatch == 0) return ESL_ENOTFOUND;
  return ESL_OK;
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
 * Return <ESL_OK> on success, <ESL_EOD> if we're out of
 * options. 
 * 
 * Throws <ESL_EINVAL> if something's wrong with the options.
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
    return ESL_EOD; 		/* normal end-of-data (end of options) return  */

  /* Check to see if we're being told that this is the end
   * of the options with the special "--" flag.
   */
  if (g->optstring == NULL &&
      strcmp(g->argv[g->optind], "--") == 0)
    { 
      g->optind++;
      return ESL_EOD; 		/* also a normal end-of-data return */
    }

  /* We have an option: an argv element that starts with -, but is
   * not "-" or "--".
   * 
   * We know the strncmp() test is ok for 2 chars, because if the option was
   * 1 char, we would've already caught it above (either it's a bare "-"
   * or it's a single non-option char, and in either case it's not an option
   * and we returned ESL_EOD.
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
  return ESL_OK;
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
 * Returns ESL_OK on success, returning the option number through
 * <ret_opti>, and a ptr to its argument through <ret_optarg> (or NULL
 * if this option takes no argument.) Internally, g->optind is
 * advanced to next argv element (+1, +2, +1, respectively, for --foo,
 * --foo arg, --foo=arg).
 *
 * Throws ESL_EINVAL and issues a useful error mesg if:
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
  if (status == ESL_EAMBIGUOUS)
    {
      esl_error(ESL_EINVAL, __FILE__, __LINE__,
		"Abbreviated option \"%s\" is ambiguous.\n\n%s", g->argv[g->optind], g->usage);
      return ESL_EINVAL;
    }
  if (status == ESL_ENOTFOUND)
    { 
      esl_error(ESL_EINVAL, __FILE__, __LINE__,
		"No such option \"%s\".\n\n%s", g->argv[g->optind], g->usage);
      return ESL_EINVAL;
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
	  esl_error(ESL_EINVAL, __FILE__, __LINE__,
		    "Option %s requires an argument\n\n%s", g->opt[opti].name, g->usage);
	  return ESL_EINVAL;
	}
      else			/* "--foo 666" style, with a space */
	*ret_optarg = g->argv[g->optind++];	/* assign optind as the arg, advance counter */
    }
  else  /* if there's not supposed to be an arg, but there is, then die */
    {
      if (argptr != NULL) 
	{
	  esl_error(ESL_EINVAL, __FILE__, __LINE__,
		    "Option %s does not take an argument\n\n%s", g->opt[opti].name, g->usage);
	  return ESL_EINVAL;
	}
      *ret_optarg = NULL;
    }

  return ESL_OK;
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
 * Process next optchar; return ESL_OK on success, returning option
 * number through <ret_opti> and a pointer to any argument through
 * <ret_optarg>. Internally, optind is advanced to the next argv element;
 * either 0, +1, or +2, depending on whether we're still processing an
 * optstring from a prev optind, starting a new optstring, or reading
 * a "-W arg" form, respectively.
 * 
 * Throws <ESL_EINVAL> and issues helpful error mesg if:
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
      esl_error(ESL_EINVAL, __FILE__, __LINE__,
		"No such option \"-%c\".\n\n%s", *(g->optstring), g->usage);
      return ESL_EINVAL;
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
	  esl_error(ESL_EINVAL, __FILE__, __LINE__,
		    "Option %s requires an argument\n\n%s", g->opt[opti].name, g->usage);
	  return ESL_EINVAL;
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
  return ESL_OK;
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
 * range.  If successful, return <ESL_OK>. 
 * 
 * The <setby> flag is used to help format useful error messages,
 * by saying who was responsible for a bad <val>.
 *
 * Returns: <ESL_OK> on success.
 *
 * Throws:  <ESL_EINVAL>:         <val> is not the right type.
 *          <ESL_ERANGE>:         <val> is out of allowed range.
 *          <ESL_ESYNTAX>:        a range string format was bogus.
 *          <ESL_EINCONCEIVABLE>: "can't happen" internal errors.
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
	esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		  "option %s takes integer arg; got %s %s\n\n%s", 
		  g->opt[i].name, val, where, g->usage);
	return ESL_EINVAL;
      }

    status = verify_integer_range(val, g->opt[i].range);
    if (status == ESL_ERANGE)
      {
	esl_error(ESL_ERANGE, __FILE__, __LINE__, 
		  "option %s takes integer arg in range %s; got %s %s\n\n%s", 
		  g->opt[i].name, g->opt[i].range, val, where, g->usage);
	return ESL_EINVAL;
      }
    else if (status == ESL_ESYNTAX) /* ESL_ESYNTAX, or anything else */
      {
	esl_error(ESL_ESYNTAX, __FILE__, __LINE__, "range string %s for option %s is corrupt",
		  g->opt[i].range, g->opt[i].name); 
	return ESL_ESYNTAX;
      }
    else if (status != ESL_OK) ESL_ERROR(ESL_EINCONCEIVABLE, "unexpected error code");
    break;

  case eslARG_REAL:
    if (! is_real(val))
      {
	esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		  "option %s takes real-valued arg; got %s %s\n\n%s", 
		  g->opt[i].name, val, where, g->usage);
	return ESL_EINVAL;
      }

    status = verify_real_range(val, g->opt[i].range);
    if (status == ESL_ERANGE)
      {
	esl_error(ESL_ERANGE, __FILE__, __LINE__, 
		  "option %s takes real-valued arg in range %s; got %s %s\n\n%s", 
		  g->opt[i].name, g->opt[i].range, val, where, g->usage);
	return ESL_ERANGE;
      }
    else if (status == ESL_ESYNTAX)
      {
	esl_error(ESL_ESYNTAX, __FILE__, __LINE__, "range string %s for option %s is corrupt",
		  g->opt[i].range, g->opt[i].name); 
	return ESL_ESYNTAX;
      }
    else if (status != ESL_OK) ESL_ERROR(ESL_EINCONCEIVABLE, "unexpected error code");
    break;

  case eslARG_CHAR:
    if (strlen(g->val[i]) > 1)
      {
	esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		  "option %s takes char arg; got %s %s\n\n%s", 
		  g->opt[i].name, val, where, g->usage);
	return ESL_EINVAL;
      }
    status = verify_char_range(val, g->opt[i].range);
    if (status == ESL_ERANGE)
      {
	esl_error(ESL_ERANGE, __FILE__, __LINE__, 
		  "option %s takes char arg in range %s; got %s %s\n\n%s", 
		  g->opt[i].name, g->opt[i].range, val, where, g->usage);
	return ESL_ERANGE;
      }
    else if (status == ESL_ESYNTAX)
      {
	esl_error(ESL_ESYNTAX, __FILE__, __LINE__, "range string %s for option %s is corrupt",
		  g->opt[i].range, g->opt[i].name); 
	return ESL_ESYNTAX;
      }
    else if (status != ESL_OK) ESL_ERROR(ESL_EINCONCEIVABLE, "unexpected error code");
    break;

  case eslARG_STRING: /* unchecked type. */
    if (g->opt[i].range != NULL)
      {
	esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		  "option %s takes a string arg that can't be range checked",  g->opt[i].name);
	return ESL_EINVAL;
      }
    break;			
    
  default: ESL_ERROR(ESL_EINVAL, "no such argument type");
  }

  return ESL_OK;
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


/* Function: is_real()
 * 
 * Purpose:  Returns TRUE if <s> is a string representation
 *           of a valid floating point number, convertable
 *           by atof().
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
 * Returns <ESL_OK> if the string <arg>, when converted 
 * to an integer by atoi(), gives a value that lies within
 * the given <range>, if <range> is non-NULL. (If
 * <range> is NULL, there is no constraint on the range
 * of this <arg>, so return TRUE.) Else, <arg> does
 * not lie in the <range>; return <ESL_ERANGE>. If
 * <range> is misformatted, return <ESL_ESYNTAX>, so caller
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
 * Returns:  <ESL_OK>:      <arg> is within allowed <range>.
 *           <ESL_ERANGE>:  <arg> is not within allowed <range>.
 *           <ESL_ESYNTAX>: something wrong with <range> string.
 */
static int
verify_integer_range(char *arg, char *range)
{
  int   n;
  int   upper, lower;		/* upper, lower bounds */
  char *up, *lp;		
  int   geq, leq;	        /* use >=, <= instead of >, < */
  
  if (range == NULL) return ESL_OK;
  n = atoi(arg);

  if (parse_rangestring(range, 'n', &lp, &geq, &up, &leq) != ESL_OK) 
    return ESL_ESYNTAX;

  if (lp != NULL)
    {
      lower = atoi(lp);
      if ((geq && ! (n >= lower)) || (!geq && !(n > lower)))
	return ESL_ERANGE;
    }

  if (up != NULL) 
    {
      upper = atoi(up);
      if ((leq && ! (n <= upper)) || (!leq && !(n < upper)))
	return ESL_ERANGE;
    }
  return ESL_OK;
}



/* verify_real_range():
 * 
 * Verify that a string <arg>, when converted to a
 * double-precision real by atof(), gives a value that lies
 * within the range defined by <range>. If <range> is NULL,
 * there is no range constraint, and any <arg> is valid.
 *
 * Returns:  <ESL_OK>:      <arg> is within allowed <range>.
 *           <ESL_ERANGE>:  <arg> is not within allowed <range>.
 *           <ESL_ESYNTAX>: something wrong with <range> string.
 */
static int
verify_real_range(char *arg, char *range)
{
  double x;
  double upper, lower;		/* upper, lower bounds */
  char  *up, *lp;		
  int    geq, leq;	        /* use >=, <= instead of >, < */
  
  if (range == NULL) return ESL_OK;
  x = atof(arg);
  
  if (parse_rangestring(range, 'x', &lp, &geq, &up, &leq) != ESL_OK) 
    return ESL_ESYNTAX;

  if (lp != NULL)
    {
      lower = atof(lp);
      if ((geq && ! (x >= lower)) || (!geq && !(x > lower)))
	return ESL_ERANGE;
    }

  if (up != NULL) 
    {
      upper = atof(up);
      if ((leq && ! (x <= upper)) || (!leq && !(x < upper)))
	return ESL_ERANGE;
    }
  return ESL_OK;
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
 * Returns:  <ESL_OK>:      <arg> is within allowed <range>.
 *           <ESL_ERANGE>:  <arg> is not within allowed <range>.
 *           <ESL_ESYNTAX>: something wrong with <range> string.
 */
static int
verify_char_range(char *arg, char *range)
{
  char   c;
  char  *upper, *lower;		
  int    geq, leq;	        /* use >=, <= instead of >, < */
  
  if (range == NULL) return ESL_OK;
  c = *arg;

  if (parse_rangestring(range, 'c', &lower, &geq, &upper, &leq) != ESL_OK) 
    return ESL_ESYNTAX;

  if (lower != NULL)
    {
      if ((geq && ! (c >= *lower)) || (!geq && !(c > *lower)))
	return ESL_ERANGE;
    }

  if (upper != NULL) 
    {
      if ((leq && ! (c <= *upper)) || (!leq && !(c < *upper)))
	return ESL_ERANGE;
    }
  return ESL_OK;
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
 * Returns <ESL_OK> on success, <ESL_ESYNTAX> if the range string
 * is invalid. No errors are thrown here, so caller can format a
 * useful error message if range string is bogus.
 */
static int
parse_rangestring(char *range, char c, char **ret_lowerp, int *ret_geq, char **ret_upperp, int *ret_leq)
{
  char *ptr;

  *ret_geq    = *ret_leq    = FALSE;	/* 'til proven otherwise */
  *ret_lowerp = *ret_upperp = NULL;     /* 'til proven otherwise */

  if ((ptr = strchr(range, c)) == NULL) return ESL_ESYNTAX;
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
      else return ESL_ESYNTAX;
    }
  else	/* we're in a<=n<=b form; upper bound after n, lower before */
    {
      if (*(ptr+1) != '<') return ESL_ESYNTAX;
      if (*(ptr+2) == '=') { *ret_leq = TRUE; *ret_upperp = ptr+3; } else *ret_upperp = ptr+2;

      ptr--;
      if (*ptr == '=') { *ret_geq = TRUE; ptr--; }
      if (*ptr != '<') return ESL_ESYNTAX;
      *ret_lowerp = range;	/* start of string */
    }
  return ESL_OK;
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
 * Returns: <ESL_OK> if an option has been successfully parsed
 *          out of the list and <ret_opti> is valid;
 *          <ESL_EOD> if no more option remains (<s> is NULL,
 *          or points to a \0).
 *          
 * Throws:  <ESL_EINVAL> if an option in the list isn't
 *          recognized.         
 */
static int 
process_optlist(ESL_GETOPTS *g, char **ret_s, int *ret_opti)
{
  char *s;
  int   i;
  int   n;
  
  if ((s = *ret_s) == NULL) return ESL_EOD;
  if (*s == '\0')           return ESL_EOD;

  n = strcspn(s, ",");

  /* a little weak here; we're only matching a n-long prefix, so we're
   * not going to catch the case where the optlist contains a
   * truncated, ambiguous optionname.  but optlists are not user
   * input, so the answer to this problem is "don't do that".
   */
  for (i = 0; i < g->nopts; i++)
    if (strncmp(g->opt[i].name, s, n) == 0) break;
  if (i == g->nopts) 
    ESL_ERROR(ESL_EINVAL, "no such option");

  *ret_opti = i;

  if (s[n] == ',') *ret_s = s+n+1; 
  else             *ret_s = NULL;

  return ESL_OK;
}

/*------- end of private functions for processing optlists -----------*/






#ifdef ESL_GETOPTS_TESTDRIVE 
/* gcc -g -Wall -o test -I. -DESL_GETOPTS_TESTDRIVE getopts.c easel.c
 */

#include <stdio.h>

#include <easel/easel.h>
#include <easel/getopts.h>

static ESL_OPTIONS options[] = {
  /* name          type         range   default   env_var  toggles  requires incompat_with */
  { "-a",     eslARG_NONE,       NULL,   FALSE, "FOOTEST",   NULL,    NULL,   NULL },
  { "-b",     eslARG_NONE,       NULL,   FALSE,     NULL, "--no-b",   NULL,   NULL },
  { "--no-b", eslARG_NONE,       NULL,   FALSE,     NULL,     "-b",   NULL,   NULL },
  { "-c",     eslARG_CHAR,  "a<=c<=z",    "x",      NULL,    NULL,    NULL,   NULL },
  { "-n",     eslARG_INT,   "0<=n<10",    "0",      NULL,    NULL,    NULL,   NULL },
  { "-x",     eslARG_REAL,    "0<x<1",  "0.5",      NULL,    NULL,    NULL,   NULL },
  { "--lowx", eslARG_REAL,      "x>0",  "1.0",      NULL,    NULL,    NULL,   NULL },
  { "--hix",  eslARG_REAL,      "x<1",  "0.9",      NULL,    NULL,    NULL,   NULL },
  { "--lown", eslARG_INT,       "n>0",   "42",      NULL,    NULL,  "-a,-b",  NULL },
  { "--hin",  eslARG_INT,       "n<0",   "-1",      NULL,    NULL,    NULL, "-a,-b" },
  { "--host", eslARG_STRING,     NULL,      "", "HOSTNAME",  NULL,    NULL,   NULL },
  {  0, 0, 0, 0, 0, 0, 0, 0 },
};

static char usage[] = "\
Usage: test [-options] <arg>\n\
";

    
int
main(int argc, char **argv)
{
  ESL_GETOPTS *go;
  int   state;
  char  c;
  int   n;
  float x;
  char *s;

  go = esl_getopts_Create(options, usage);
  esl_opt_ProcessEnvironment(go);
  esl_opt_ProcessCmdline(go, argc, argv);
  esl_opt_VerifyConfig(go);

  esl_opt_GetBooleanOption(go, "-a", &state);
  printf("Option -a:     %s\n", (state == FALSE)? "off" : "on");

  esl_opt_GetBooleanOption(go, "-b", &state);
  printf("Option -b:     %s\n", (state == FALSE)? "off" : "on");

  esl_opt_GetCharOption(go, "-c", &c);
  printf("Option -c:     %c\n", c);

  esl_opt_GetIntegerOption(go, "-n", &n);
  printf("Option -n:     %d\n", n);

  esl_opt_GetFloatOption(go, "-x", &x);
  printf("Option -x:     %.1f\n", x);

  esl_opt_GetStringOption(go, "--host", &s);
  printf("Option --host: %s\n", s);

  esl_getopts_Destroy(go);
  exit(0);
}

#endif /*ESL_GETOPTS_TESTDRIVE*/
/*-------------- end of examples, test driver ********************/

/*****************************************************************  
 * @LICENSE@
 *****************************************************************/
