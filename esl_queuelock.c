/* User-space queueing lock for high-contention regions.  Should perform better
 * than base pthreads locks when a lock is often held at the time another thread
 * wants to acquire it, and osly trivially worse when the lock is available.
 * 
 * Contents:
 *    1. The ESL_QUEUELOCK API
 *    2. Unit tests
 */

#include "easel.h"
#include "esl_queuelock.h"
#include "esl_threads.h"

/*****************************************************************
 * 1. The ESL_QUEUELOCK API
 *****************************************************************/


/* Function:  esl_queuelock_Create()
 * Synopsis:  Creates an esl_queuelock object .
 * Incept:    NPC, Wed 5/17/2023 12:08 PM
 *
 * Purpose:   Allocates and returns a new ESL_QUEUELOCK object.  If max_lockers is 
 *            > 0, configures the queuelock to support at most that many threads
 *            waiting on the lock.  if max_lockers is <= 0, configures the queuelock 
 *            to support as many waiting threads as the system has CPU cores.
 *
 * Returns:   A pointer to the new ESL_QUEUELOCK object
 * 
 */
extern ESL_QUEUELOCK *esl_queuelock_Create(int max_lockers){
  ESL_QUEUELOCK *obj = NULL;
  int status;
  int real_lockers; 

  if(max_lockers >0){
    real_lockers = max_lockers;
  }
  else{
    real_lockers = esl_threads_GetCPUCount();
  }

  ESL_ALLOC(obj, sizeof(ESL_QUEUELOCK));
  if(pthread_mutex_init(&(obj->lock), NULL)){
    esl_fatal("Unable to initialize pthreads mutex in esl_queuelock_Create().\n");
  }
  obj->wait_buffer = NULL; // set this so we can tell if we need to free it on error
  obj->locked = 0;
  obj->max_waiters = real_lockers;
  obj->first_waiter = 0; 
  obj->num_waiters = 0;  
  ESL_ALLOC(obj->wait_buffer, real_lockers * sizeof(int));
  for(int i = 0; i < real_lockers; i++){
    obj->wait_buffer[i] = 0;
  }
  return obj; 
ERROR:
  if (obj != NULL){
    if(obj->wait_buffer != NULL){
      free((void *) obj->wait_buffer);
    }
    free(obj);
  }
  esl_fatal("Unable to allocate memory in esl_queuelock_Create()\n");
}

/* Function:  esl_queuelock_Destroy()
 * Synopsis:  Destroys an esl_queuelock object .
 * Incept:    NPC, Wed 5/17/2023 12:13 PM
 *
 * Purpose:   Destroys the provided ESL_QUEUELOCK object, freeing all of its internal
 *            storage
 *
 * Returns:   Nothing
 */
extern void esl_queuelock_Destroy(ESL_QUEUELOCK *the_lock){
  if(the_lock == NULL){
    return;
  }
  if(the_lock->wait_buffer != NULL){
    free((void *) the_lock->wait_buffer);
  }
  pthread_mutex_destroy(&(the_lock->lock)); 
  free(the_lock);
}

/* Function:  esl_queuelock_Lock()
 * Synopsis:  Locks an esl_queuelock object .
 * Incept:    NPC, Wed 5/17/2023 2:24 PM
 *
 * Purpose:   Acquires the provided ESL_QUEUELOCK.  If the lock is available
 *            marks the lock as locked and returns.  If not, adds itself to the 
 *            list of waiting threads and spin-waits on the appropriate location 
 *            in the(_lock->wait_buffer until notified that it has the lock
 *
 * Returns:   Nothing.  Returning indicates that the lock has been acquired
 * 
 * Throws:    Calls esl_fatal() if adding itself to the list wf waiters would mean
 *            that there are more threads waiting on the lock than it can support  
 */
extern void esl_queuelock_Lock(ESL_QUEUELOCK *the_lock){
  pthread_mutex_lock(&(the_lock->lock));
  if(!the_lock->locked){ // lock is available, so take it
    the_lock->locked = 1;
    the_lock->num_waiters = 0;
    the_lock->first_waiter = 0;
  }
  else{
    if(the_lock->num_waiters >= the_lock->max_waiters){
      // Can't add ourselves to the list of waiters because there isn't enough
      // space in the wait buffer
      esl_fatal("Attempted to add more waiters to an ESL_QUEUELOCK than it could support\n");
    }
    volatile int *my_wait_location;

    //Find the location in the wait buffer that we should wait on
    if(the_lock->first_waiter + the_lock->num_waiters < the_lock->max_waiters){
      my_wait_location = the_lock->wait_buffer+ (the_lock->first_waiter + the_lock->num_waiters);
    }
    else{
      my_wait_location = the_lock->wait_buffer + (the_lock->first_waiter + the_lock->num_waiters - the_lock->max_waiters);
    }
    *my_wait_location = 0;
    the_lock->num_waiters++;
    while(*my_wait_location ==0){} // spin here until signaled that I have the lock
  }

  // If we get here, we've acquired the queuelock, so release the lock on the data
  // structure
  pthread_mutex_unlock(&(the_lock->lock));
  return;
}

/* Function:  esl_queuelock_Unlock()
 * Synopsis:  Unlocks an esl_queuelock object .
 * Incept:    NPC, Wed 5/17/2023 2:51 PM
 *
 * Purpose:   Releases the provided ESL_QUEUELOCK.  If one or more threads are 
 *            waiting to acquire te lock, passes the lock to the first waiter in the 
 *            queue.  If the lock is not locked when this function is called, it 
 *            returns immediately.
 *
 * Returns:   Nothing.  Returning indicates that the lock has been released
 * 
 */
extern void esl_queuelock_Unlock(ESL_QUEUELOCK *the_lock){
  pthread_mutex_lock(&(the_lock->lock));
  
  if(!the_lock->locked){  // lock was already unlocked, so just return
    pthread_mutex_unlock(&(the_lock->lock));
    return;
  }
  if(the_lock->num_waiters ==0){ // No waiters, simple case
    the_lock->locked = 0;
  }
  else{  // at least one waiter, so need to signal
    volatile int *next_waiter = the_lock->wait_buffer + the_lock->first_waiter;
    the_lock->first_waiter++;
    if(the_lock->first_waiter == the_lock->max_waiters){ // wrap around to beginning
      the_lock->first_waiter = 0;
    }
    the_lock->num_waiters--;
    *next_waiter = 1; // signal the waiter
  }

  //getting here means we've released the lock
  pthread_mutex_unlock(&(the_lock->lock));
  return;
}

/*****************************************************************
 * 2. Unit tests
 *****************************************************************/

#ifdef eslQUEUELOCK_TESTDRIVE
#include <unistd.h>
typedef struct{
  pthread_mutex_t lock;
  ESL_QUEUELOCK *queuelock;
  volatile int counter;
} THREAD_ARGS;

static void *testthread_one(void *arg){
  THREAD_ARGS *args = (THREAD_ARGS *)arg;
  pthread_mutex_t lock; 
  esl_queuelock_Lock(args->queuelock); // Wait for master to release this lockj
  args->counter = 1;
  pthread_exit(NULL);  
}
int
main(int argc, char **argv)
{
  ESL_QUEUELOCK *the_lock;

  // First, test lock creation with specified number of max lockers
  the_lock = esl_queuelock_Create(4);
  esl_queuelock_Lock(the_lock);  // should go through because lock not locked at create
  esl_queuelock_Unlock(the_lock);
  esl_queuelock_Unlock(the_lock);  // test that unlocking an unlocked lock works
  esl_queuelock_Destroy(the_lock); 

  // Second, lock creation with the default number of max lockers
  int num_cores = esl_threads_GetCPUCount();
  the_lock = esl_queuelock_Create(0);
  if (the_lock->max_waiters != num_cores){
    esl_fatal("Queuelock created with default number of waiters had wrong number of waiters\n");
  }
  esl_queuelock_Lock(the_lock);  // should go through because lock not locked at create
  esl_queuelock_Unlock(the_lock);
  esl_queuelock_Unlock(the_lock);  // test that unlocking an unlocked lock works
  esl_queuelock_Destroy(the_lock); 


  // Third: Test that locking a locked lock stalls until the lock is released
  THREAD_ARGS args;

  args.queuelock = esl_queuelock_Create(2);
  args.counter = 0;

  //Plan:  lock the queuelock.  Start the worker thread,
  // which will try to acquire the queuelock and then change the value of args->counter.
  // Then, sleep for 1s and check the value of args->counter, which should be unchanged.
  // Next, release the queuelock, freeing the worker thread to change the counter and  terminate.
  // Finally, we check the counter to make sure that the worker thread set it.
  esl_queuelock_Lock(args.queuelock);
  
  pthread_attr_t attr;
  if(pthread_attr_init(&attr)){
    esl_fatal("Unable to create pthreads attribute\n");
  }

  pthread_t thread_ids[3];  // only need one now, but will want the others later
  pthread_create(&(thread_ids[0]), &attr, testthread_one, (void*) &args);
  sleep(1);
  if (args.counter !=0){
    esl_fatal("Worker thread proceeded even though lock was locked\n");
  }
  esl_queuelock_Unlock(args.queuelock);
  pthread_join(thread_ids[0], NULL);  // Wait for the worker to do its thing
  if(args.counter !=1){
    esl_fatal("Worker thread never acquired lock\n");
  }
  if(args.queuelock->locked == 0){
    esl_fatal("Worker thread should be holding lock, but isn't\n");
  }

  esl_queuelock_Destroy(args.queuelock);
  return eslOK;
}


#endif