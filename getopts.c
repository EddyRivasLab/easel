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
    if (verify_type_and_range(g, i, g->val[i], eslARG_SETBY_DEFAULT) != ESL_OK)
      {
	esl_getopts_Destroy(g); 
	return NULL;
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


/* Function:  esl_opt_ProcessCmdline()
 * Incept:    SRE, Wed Jan 12 10:12:43 2005 [St. Louis]
 *
 * Purpose:   
 *
 * Args:      
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
		    "Option %s appears more than once on command line.\n\n%s", 
		    g->opt[i].name, g->usage);
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
			g->opt[togi].name, g->opt[i].name, g->usage);
	      return ESL_EINVAL;
	    }
	  
	  g->setby[togi] = eslARG_SETBY_CMDLINE; /* indirectly, but still */
	  g->val[opti] = NULL;	/* ok for false booleans too */
	}
      if (status != ESL_EOD) return status; /* not a normal end of optlist */
    }
  return ESL_OK;
}



/* Function:  esl_opt_VerifyConfig()
 * Incept:    SRE, Wed Jan 12 10:21:46 2005 [St. Louis]
 *
 * Purpose:   For every option that is set, make sure any
 *            required options are also set, and no
 *            incompatible options are set. "Set" means
 *            <val> is non-NULL (including booleans),
 *            and "not set" means <val> is NULL. (That is,
 *            we don't go by <setby>, which refers to who
 *            determined the state of an option, even if
 *            it is off.)
 *
 * Args:      
 *
 * Returns:   
 *
 * Throws:    (no abnormal error conditions)
 *
 * Xref:      
 */
int
esl_opt_VerifyConfig(ESL_GETOPTS *g)
{
  int   i,reqi,incompati;
  char *s;
  int   status;

  /* For all options that are set (turned on with non-NULL vals), 
   * verify that all their required_opts are set.
   */
  for (i = 0; i < g->nopts; i++)
    {
      if (g->val[i] != NULL)
	{
	  s = g->opt[i].required_opts;
	  while ((status = process_optlist(g, &s, &reqi)) == ESL_OK)
	    if (g->val[reqi] == NULL)
	      {
		esl_error(ESL_EINVAL, __FILE__, __LINE__,
			  "Option %s requires (or has no effect without) option %s\n\n%s", 
			  g->opt[i].name, g->opt[reqi].name, g->usage);
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
      if (g->val[i] != NULL)
	{
	  s = g->opt[i].incompat_opts;
	  while ((status = process_optlist(g, &s, &incompati)) == ESL_OK)
	    if (g->val[incompati] != NULL)
	      {
		esl_error(ESL_EINVAL, __FILE__, __LINE__,
			  "Option %s is incompatible with option %s\n\n%s", 
			  g->opt[i].name, g->opt[incompati].name, g->usage);
		return ESL_EINVAL;
	      }
	  if (status != ESL_EOD) return status;	/* non-normal end of optlist; throw status up */
	}
    }

    return ESL_OK;
}



/* Function:  esl_opt_GetBooleanArg()
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
esl_opt_GetBooleanArg(ESL_GETOPTS *g, char *optname, int *ret_state)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == ESL_ENOTFOUND)
    ESL_ERROR(ESL_ENOTFOUND, "no such option");

  if (g->val[opti] == NULL) *ret_state = FALSE;
  else                      *ret_state = TRUE;
  return ESL_OK;
}

/* Function:  esl_opt_GetIntegerArg()
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
esl_opt_GetIntegerArg(ESL_GETOPTS *g, char *optname, int *ret_n)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == ESL_ENOTFOUND)
    ESL_ERROR(ESL_ENOTFOUND, "no such option");
  *ret_n = atoi(g->val[opti]);
  return ESL_OK;
}
		
/* Function:  esl_opt_GetRealArg()
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
esl_opt_GetRealArg(ESL_GETOPTS *g, char *optname, int *ret_x)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == ESL_ENOTFOUND)
    ESL_ERROR(ESL_ENOTFOUND, "no such option");
  *ret_x = atof(g->val[opti]);
  return ESL_OK;
}

/* Function:  esl_opt_GetCharArg()
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
esl_opt_GetCharArg(ESL_GETOPTS *g, char *optname, char *ret_c)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == ESL_ENOTFOUND)
    ESL_ERROR(ESL_ENOTFOUND, "no such option");
  *ret_c = *g->val[opti];
  return ESL_OK;
}

/* Function:  esl_opt_GetStringArg()
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
esl_opt_GetStringArg(ESL_GETOPTS *g, char *optname, char **ret_s)
{
  int opti;

  if (get_optidx_exactly(g, optname, &opti) == ESL_ENOTFOUND)
    ESL_ERROR(ESL_ENOTFOUND, "no such option");
  *ret_s = g->val[opti];
  return ESL_OK;
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
    if (strcmp(optname, g->opt[i].name) == 0) { *ret_opti = i; return ESL_OK}
  return ESL_ENOTFOUND;
}

/* get_optidx_abbrev():
 * 
 * Find option named <optname> in <g>; set <ret_opti> to be
 * the index of the option, and return ESL_OK. Allow <optname>
 * to be an abbreviation of one of the option names in <g>,
 * so long as it is unambiguous.
 * 
 * If the option is not found, return <ESL_ENOTFOUND>.
 * If <optname> ambiguously matches two or more options in <g>,
 * return <ESL_EAMBIGUOUS>.
 */
static int
get_optidx_abbrev(ESL_GETOPTS *g, char *optname, int *ret_opti)
{
  int nmatch = 0;
  int n;

  n = strlen(optname);		/* all of optname abbrev must match against the real name */
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
  char *argptr;
  int   opti;
  int   nmatch;

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
   * We know the strncmp() test is ok for 2 chars, because if the option was
   * 1 char, we would've already caught it above (either it's a bare "-"
   * or it's a single non-option char, and in either case it's not an option
   * and we returned ESL_EOD.
   * 
   * Watch out for the case where we're in the middle of a concatenated optstring
   * of single-letter options, a la -abc
   */
  if (g->optstring == NULL && strncmp(g->argv[optind], "--", 2) == 0)
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
  int   i;		/* counter for options                               */
  int   opti;		/* option number found                               */
  int   nmatch;		/* number of options in opt[] that match             */
  char *argptr;		/* ptr to arg in --foo=arg syntax                    */
  int   arglen;		/* length of argv elem's option part (up to = or \0) */
  int   status;

  /* Deal with options of syntax "--foo=arg" w/o modifying argv.
   */
  if ((argptr = strchr(g->argv[optind], '=')) != NULL)
    { arglen = argptr - g->argv[optind]; argptr++; } /* bump argptr off the = to the arg itself */
  else
    { arglen = strlen(g->argv[optind]); } /* and argptr == NULL from above. */

  /* Figure out which option this is.
   * The trick here is to allow abbreviations, and identify
   * ambiguities while we're doing it. (GNU getopt allows abbrevs.)
   */
  status = get_optidx_abbrev(g, g->argv[optind], &opti);
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
		    "Option %s requires an argument\n\n%s", g->opt[opti].name, usage);
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
		    "Option %s does not take an argument\n\n%s", g->opt[opti].name, usage);
	  return ESL_EINVAL;
	}
      *ret_optarg = NULL;
    }

  return;
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
      else if (g->optind+1 < g->argc)  /* unattached argument; assign optind, advance counter  */
	*ret_optarg = g->argv[g->optind++];
      else 
	{
	  esl_error(ESL_EINVAL, __FILE__, __LINE__,
		    "Option %s requires an argument\n\n%s", g->opt[opti].stdname, usage);
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
 *          <ESL_INCONCEIVABLE>:  "can't happen" internal errors.
 */
static int
verify_type_and_range(ESL_GETOPTS *g, int i, char *val, int setby)
{
  char *where;

  if       (setby == eslARG_SETBY_DEFAULT) where = "as default";
  else if  (setby == eslARG_SETBY_CMDLINE) where = "on cmdline";
  else if  (setby == eslARG_SETBY_ENV)     where = "in env";
  else if  (setby >= eslARG_SETBY_CFGFILE) where = "in cfgfile";

  switch (g->opt[i].type) {

  case eslARG_NONE:
    esl_error(ESL_EINVAL, __FILE__, __LINE__, 
	      "option %s takes no argument; got %s %s\n\n%s", 
	      g->opt[i].name, val, where, g->opt[i].usage);
    return ESL_EINVAL;

  case eslARG_INT:
    if (! is_integer(val))
      {
	esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		  "option %s takes integer arg; got %s %s\n\n%s", 
		  g->opt[i].name, val, where, g->opt[i].usage);
	return ESL_EINVAL;
      }

    status = verify_integer_range(val, range);
    if (status == ESL_EINVAL)
      {
	esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		  "option %s takes integer arg in range %s; got %s %s\n\n%s", 
		  g->opt[i].name, g->opt[i].range, val, where, g->opt[i].usage);
	return ESL_EINVAL;
      }
    else (status == ESL_ESYNTAX) /* ESL_ESYNTAX, or anything else */
      {
	esl_error(ESL_ESYNTAX, __FILE__, __LINE__, "range string %s for option %s is corrupt",
		  g->opt[i].range, g->opt[i].name); 
	return ESL_ESYNTAX;
      }
    else if (status != ESL_OK) ESL_ERROR(ESL_INCONCEIVABLE, "unexpected error code");
    break;

  case eslARG_REAL:
    if (! is_real(val))
      {
	esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		  "option %s takes real-valued arg; got %s %s\n\n%s", 
		  g->opt[i].name, val, where, g->opt[i].usage);
	return ESL_EINVAL;
      }

    status == verify_real_range(val, range);
    if (status == ESL_EINVAL)
      {
	esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		  "option %s takes real-valued arg in range %s; got %s %s\n\n%s", 
		  g->opt[i].name, g->opt[i].range, val, where, g->opt[i].usage);
	return ESL_EINVAL;
      }
    else if (status == ESL_ESYNTAX)
      {
	esl_error(ESL_ESYNTAX, __FILE__, __LINE__, "range string %s for option %s is corrupt",
		  g->opt[i].range, g->opt[i].name); 
	return ESL_ESYNTAX;
      }
    else if (status != ESL_OK) ESL_ERROR(ESL_INCONCEIVABLE, "unexpected error code");
    break;

  case eslARG_CHAR:
    if (strlen(g->val[i]) > 1)
      {
	esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		  "option %s takes char arg; got %s %s\n\n%s", 
		  g->opt[i].name, val, where, g->opt[i].usage);
	return ESL_EINVAL;
      }
    status = verify_char_range(val, range);
    if (status == ESL_EINVAL)
      {
	esl_error(ESL_EINVAL, __FILE__, __LINE__, 
		  "option %s takes char arg in range %s; got %s %s\n\n%s", 
		  g->opt[i].name, g->opt[i].range, val, where, g->opt[i].usage);
	return ESL_EINVAL;
      }
    else if (status == ESL_ESYNTAX)
      {
	esl_error(ESL_ESYNTAX, __FILE__, __LINE__, "range string %s for option %s is corrupt",
		  g->opt[i].range, g->opt[i].name); 
	return ESL_ESYNTAX;
      }
    else if (status != ESL_OK) ESL_ERROR(ESL_INCONCEIVABLE, "unexpected error code");
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
      upper = atoi(up);
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

  *ret_opti = opti;

  if (s[n] == ',') *ret_s = s+n+1; 
  else             *ret_s = NULL;

  return ESL_OK;
}

/*------- end of private functions for processing optlists -----------*/






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
