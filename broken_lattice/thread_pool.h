#ifndef THREAD_POOL_H_INC
#define THREAD_POOL_H_INC
#include <pthread.h>

struct pool_job {
	void * (*function) (void *);
	void *arg;
};

struct pool_work {
	struct pool_job *jobs;
	size_t n_jobs;
	/* current_job is ZERO INDEXED */
	size_t current_job;
	size_t n_collect;
	pthread_mutex_t job_mutex;
};



/* note: no need to set up the mutex before calling this */
void pool_run(struct pool_work *work, int n_threads);

struct test_arg{
	size_t n;
};
void *test_job(void *arg);


#endif
