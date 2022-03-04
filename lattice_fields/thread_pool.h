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


struct test_arg{
	size_t n;
};
void *test_job(void *arg);



/* note: no need to set up the mutex before calling this */
void pool_run(struct pool_work *work, int n_threads);

void thread_map(void (*function) (void *, void *, size_t, void *), void *general_args, void *in_array,
		size_t in_type_size, void *out_array, size_t out_type_size, size_t n_calls,
		size_t n_threads);

void map_test_function(void *in, void *out, size_t index, void *general_args);



#endif
