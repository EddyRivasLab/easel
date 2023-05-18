#ifndef ESL_QUEUELOCK_INCLUDED
#define ESL_QUEUELOCK_INCLUDED

#include "esl_threads.h" //This gets us pthreads.h and functions to determine how many cores we have
typedef struct{
  pthread_mutex_t lock;
  int locked; // 1 if locked, 0 otherwise
  int max_waiters;  // maximum number of simultaneous waiters the lock can handle
  // (size of wait_buffer) 
  int first_waiter; // index into waiting array of the first waiting thread
  int num_waiters; // number of threads currently waiting on the lock
  volatile int *wait_buffer; //buffer of locations threads will spin-wait on
} ESL_QUEUELOCK;


extern ESL_QUEUELOCK *esl_queuelock_Create(int max_lockers);
extern void esl_queuelock_Destroy(ESL_QUEUELOCK *the_lock);
extern void esl_queuelock_Lock(ESL_QUEUELOCK *the_lock);
extern void esl_queuelock_Unlock(ESL_QUEUELOCK *the_lock);
#endif //ESL_QUEUELOCK_INCLUDED