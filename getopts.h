

struct esl_options_s {
  char *stdname;	/* one char option ("-a"), or NULL  */
  char *longname;	/* long option ("--foo"), or NULL   */
  char *envname;	/* associated environ var ("BLASTDB"), or NULL */
  int   type;		/* for type checking: (eslARG_INT, etc.) */
  char *range;		/* for range checking: ("0<=x<=1", etc.) */
  char *required_opts;	/* comma-sep'd list of opts that must also be set */
  char *incompat_opts;	/* comma-sep'd list of opts that must not also be set */
};

/* Argument types.
 */
#define eslARG_NONE      0
#define eslARG_INT       1
#define eslARG_REAL      2
#define eslARG_CHAR      3
#define eslARG_STRING    4	/* unchecked type; includes filenames  */

struct esl_getopts_s {
  int    argc;			/* argc from command line              */
  char **argv;			/* argv from command line              */
  int    optind;		/* where we are in argc                */
  struct esl_options_s *opt;    /* array of app-defined options        */
  int    nopts;                 /* number of options                   */
  char  *optstring;		/* ptr into a string of single-char opts in some argv[] */

};
typedef struct esl_getopts_s ESL_GETOPTS;


