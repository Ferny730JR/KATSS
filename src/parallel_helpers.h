#ifndef PARALLELIZATION_HELPERS
#define PARALLELIZATION_HELPERS

#include <pthread.h>
#include "thpool.h"

pthread_mutex_t output_mutex;
unsigned int    max_threads;
threadpool      worker_pool;

/**
 *  @brief Specify a function to be atomic
*/
#define ATOMIC_BLOCK(a) { \
    if (max_threads > 1) { \
        pthread_mutex_lock(&output_mutex); \
        (a); \
        pthread_mutex_unlock(&output_mutex); \
    } else { \
        (a); \
    } \
}


#define THREADSAFE_STREAM_OUTPUT(a) { \
    if (max_threads > 1) { \
        pthread_mutex_lock(&output_mutex); \
        a; \
        pthread_mutex_unlock(&output_mutex); \
    } else { \
        a; \
    } \
}


#define INIT_PARALLELIZATION(a) { \
    max_threads = ((a) > 1) ? (unsigned int)(a) : 1; \
    /* initialize semaphores and thread pool */ \
    if (max_threads > 1) { \
        pthread_mutex_init(&output_mutex, NULL); \
        worker_pool = thpool_init(max_threads); \
    } \
}


#define UNINIT_PARALLELIZATION  { \
    if (max_threads > 1) { \
        thpool_wait(worker_pool); \
    } \
    pthread_mutex_destroy(&output_mutex); \
    if (max_threads > 1) { \
        thpool_destroy(worker_pool); \
    } \
}


#define RUN_IN_PARALLEL(fun, data)  { \
    if (max_threads > 1) { thpool_add_work(worker_pool, (void *)&fun, (void *)data); } \
    else { fun(data); } \
}


#define WAIT_FOR_FREE_SLOT(a) { \
    if (max_threads > 1) { \
        while (thpool_num_t_in_queue(worker_pool) >= (a)) { \
            usleep(1000); \
        } \
    } \
}


#define WAIT_FOR_THPOOL { \
    if (max_threads > 1) { \
        thpool_wait(worker_pool); \
    } \
}


int num_proc_cores(int  *num_cores,
               int  *num_cores_conf);


int max_user_threads(void);

#endif
