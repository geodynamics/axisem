#include <stdio.h>
#include <pthread.h>
#include <unistd.h>

// compare https://computing.llnl.gov/tutorials/pthreads/

// careful with global variable, this only works if the calling programm is not threaded!
pthread_t thread;

// stub - so c knows what the function looks like
extern void __nc_routines_MOD_nc_dump_strain_to_disk();

// tread that calls the IO function from fortran
void *cfunc_thread(void* valp)
{
   int val;
   val = *((int*)valp);
   nc_dump_strain_to_disk();
   pthread_exit(NULL);
}

// create IO thread, to be called from fortran
void c_spawn_dumpthread(int* val){
   pthread_attr_t attr; 
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
   pthread_create(&thread, &attr, cfunc_thread, (void *)val);
}

// wait for the IO thread to finish, to be called from fortran
// global thead variable allows to come back to the thread
void c_wait_for_io() {
   void *status;
   pthread_join(thread, &status);
}

