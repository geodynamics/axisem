#include <stdio.h>
#include <pthread.h>
#include <unistd.h>

// compare https://computing.llnl.gov/tutorials/pthreads/

// careful with global variable, this only works if the calling programm is not threaded!
pthread_t thread;

// stub - so c knows what the function looks like
void __nc_routines_MOD_nc_dump_all_strain();
//void __nc_routines_MOD_nc_dump_all_strain(int* stepstodump);

// tread that calls the IO function from fortran
void *cfunc_thread(void* valp)
{
   int val;
   int one = 1;
   val = *((int*)valp);
   //printf("2nd value from fortran; %d \n", val);
#if defined(unc)
#if defined(__GNUC__)
   __nc_routines_MOD_nc_dump_strain_to_disk();
#elif defined(__INTEL_COMPILER)    
   nc_routines_mp_nc_dump_strain_to_disk_();
#else
   error "Compiler not supported";
#endif
#endif
   pthread_exit(NULL);
}

// create IO thread, to be called from fortran
void c_spawn_dumpthread_(int* val){
   //printf("In main: creating thread\n");
   pthread_attr_t attr; 
   pthread_attr_init(&attr);
   pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
   //printf("1st value from fortran; %d \n", val);
   pthread_create(&thread, &attr, cfunc_thread, (void *)val);
}

// wait for the IO thread to finish, to be called from fortran
// global thead variable allows to come back to the thread
void c_wait_for_io_() {
   void *status;
//   printf("Waiting in C\n");
   pthread_join(thread, &status);
}
