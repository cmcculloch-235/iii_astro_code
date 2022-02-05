#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>

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
