

struct esl_options_s {
  char *name;		/* either short "-a" or long "--foo" style */
  char *default;	/* default setting, or NULL */
  char *envvar;	        /* associated environ var ("BLASTDB"), or NULL */
  int   type;		/* arg type, for type checking: (eslARG_INT, etc.) */
  char *range;		/* for range checking arg: ("0<=x<=1", etc.) */
  char *toggle_opts;	/* comma-sep'd optlist: turn these off if this opt is on */
  char *required_opts;	/* comma-sep'd optlist: these must also be set */
  char *incompat_opts;	/* comma-sep'd optlist: these must not be set */
};

/* Argument types.
 */
#define eslARG_NONE      0	/* option takes no argument (so, is boolean) */
#define eslARG_INT       1	/* arg convertable by atoi()                 */
#define eslARG_REAL      2	/* arg convertable by atof()                 */
#define eslARG_CHAR      3	/* arg is a single character                 */
#define eslARG_STRING    4	/* unchecked arg type; includes filenames    */

typedef struct {
  struct esl_options_s *opt;    /* array of app-defined options        */
  int    nopts;                 /* number of options                   */

  int    argc;			/* argc from command line              */
  char **argv;			/* argv from command line              */
  char  *usage;			/* command-line usage                  */
  int    optind;		/* where we are in argc                */

  int   *islong;                /* TRUE if opt i is a long option (thus not concat'able)  */
  char **val;			/* configured value for each option (as a string) */
  int   *setby;			/* array [0..nopts-1] for who set each option     */

  char  *optstring;		/* internal state ptr into a string of single-char opts in some argv[] */
} ESL_GETOPTS;


/* Possible values of the <setby> variable in ESL_GETOPTS.
 * Additionally, values of >3 also indicate a config file, in order 
 * of _ProcessConfigFile() calls (that is, setby=3 is the first 
 * config file, setby=4 is the second, etc.)
 */
#define eslARG_SETBY_DEFAULT  0
#define eslARG_SETBY_CMDLINE  1
#define eslARG_SETBY_ENV      2
#define eslARG_SETBY_CFGFILE  3
