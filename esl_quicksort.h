/* Quicksort, reentrant.
 */
#ifndef eslQUICKSORT_INCLUDED
#define eslQUICKSORT_INCLUDED


extern int esl_quicksort(void *data, int n, int (*comparison)(void *data, int o1, int o2), int *sorted_at);


#endif /*eslQUICKSORT_INCLUDED*/

/*****************************************************************
 * @LICENSE@
 *****************************************************************/
