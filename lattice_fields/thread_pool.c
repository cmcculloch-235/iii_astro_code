#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#include "thread_pool.h"
#include "util.h"

struct pool_worker_arg {
	struct pool_work *work;
	int thread_n;
};

static void *pool_worker(void *arg);

void pool_run(struct pool_work *work, int n_threads)
{
	pthread_t *threads = calloc(n_threads, sizeof(pthread_t));
	struct pool_worker_arg *w_args = calloc(n_threads,
			sizeof(struct pool_worker_arg));
	/* set up mutex */
	pthread_mutex_init(&(work->job_mutex), NULL);

	/* start the jobs */
	for (int i = 0; i < n_threads; ++i) {
		w_args[i].thread_n = i;
		w_args[i].work = work;
		pthread_create(&(threads[i]), NULL, &pool_worker, (void *) &(w_args[i]));
	}


	/* wait for threads to finish */
	for (int i = 0; i < n_threads; ++i) {
		pthread_join(threads[i], NULL);
	}

	/* clean up */
	free(w_args);
	pthread_mutex_destroy(&(work->job_mutex));
}


static void *pool_worker(void *arg)
{
	/* unpack arguments */
	struct pool_worker_arg *arg_u = (struct pool_worker_arg *) arg;
	int thread_n = arg_u->thread_n;
	struct pool_work *work = arg_u->work;
	
	struct pool_job **job_buffer = calloc(work->n_collect, sizeof(struct pool_job *));

	eprintf("Thread %d started...", thread_n);

	/* acquire the mutex for work, take some tasks, and give it back */
	pthread_mutex_lock(&(work->job_mutex));
	while (work->current_job < work->n_jobs) {

		/* get some jobs, and be sure to update the work list */
		size_t remaining_jobs = work->n_jobs - work->current_job;
		size_t acquired_jobs = 0;
		for (size_t i = 0; i < work->n_collect && i < remaining_jobs; ++i) {
			job_buffer[i] = &(work->jobs[work->current_job]);
			++acquired_jobs;
			++(work->current_job);
		}
		eprintf("thread %d : %ld jobs left...", thread_n, remaining_jobs - acquired_jobs);
		pthread_mutex_unlock(&(work->job_mutex));

		/* do the jobs */
		for (size_t i = 0; i < acquired_jobs; ++i) {
			job_buffer[i]->function(job_buffer[i]->arg);
		}

		/* want to hold the mutex before we check how many jobs are left */
		pthread_mutex_lock(&(work->job_mutex));
	}
	/* we hold the lock still when we exit */
	pthread_mutex_unlock(&(work->job_mutex));

	free(job_buffer);

	return 0;
}


void *test_job(void *arg)
{
	struct test_arg *arg_u = (struct test_arg *) arg;
	printf("%ld\n", arg_u->n);
	
	return 0;
}


struct map_worker_arg {
	void (*function) (void *, void *, size_t, void *);
	void *general_args;
	size_t index;
	void *in_start;
	size_t in_type_size;
	void *out_start;
	size_t out_type_size;
	size_t n_calls;
};

static void *map_worker(void *arg_ptr_v);

void thread_map(void (*function) (void *, void *, size_t, void *), void *general_args,
	void *in_array, size_t in_type_size, void *out_array, size_t out_type_size,
	size_t n_calls, size_t n_threads)
{

	pthread_t *threads = calloc(n_threads, sizeof(pthread_t));
	struct map_worker_arg *w_args = calloc(n_threads,
			sizeof(struct map_worker_arg));

	size_t thread_calls = n_calls / n_threads;
	/* start the jobs */
	for (size_t i = 0; i < n_threads; ++i) {
		size_t calls = (i != n_threads - 1) ? thread_calls :
					n_calls - thread_calls * i;

		w_args[i].function = function;
		w_args[i].general_args = general_args;
		w_args[i].index = i * thread_calls;
		w_args[i].in_start = (void *) ((uint8_t *) in_array + i * in_type_size * thread_calls);
		w_args[i].in_type_size = in_type_size;
		w_args[i].out_start = (void *) ((uint8_t *) out_array + i * out_type_size * thread_calls);
		w_args[i].out_type_size = out_type_size;
		w_args[i].n_calls = calls; 
		pthread_create(&(threads[i]), NULL, &map_worker, (void *) &(w_args[i]));
	}


	/* wait for threads to finish */
	for (size_t i = 0; i < n_threads; ++i) {
		pthread_join(threads[i], NULL);
	}

	/* clean up */
	free(w_args);
	free(threads);

}

static void *map_worker(void *arg_ptr_v)
{
	struct map_worker_arg arg = *(struct map_worker_arg *) arg_ptr_v;
	for (size_t i = 0; i < arg.n_calls; ++i) {
		arg.function((void *)((uint8_t *)arg.in_start + i * arg.in_type_size),
				(void *)((uint8_t *)arg.out_start + i * arg.out_type_size),
				i + arg.index, arg.general_args);
	}

	return 0;
}


void map_test_function(void *in, void *out, size_t index, void *general_args)
{
	struct test_arg arg = *(struct test_arg *) general_args;
	*(int *) out = (int) (index + arg.n);
	eprintf("(%ld %ld) ", index, arg.n);
}
