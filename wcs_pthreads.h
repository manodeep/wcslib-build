/* Adapted by Manodeep Sinha from : https://nachtimwald.com/2019/04/05/cross-platform-thread-wrapper/
Copyright (c) 2019 John Schember <john@nachtimwald.com>

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the “Software”), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE. */
#pragma once


#if defined(_WIN32) || defined(_WIN64)
#define ON_WINDOWS
#endif

#ifdef __cplusplus
extern "C" {
#endif

#ifdef ON_WINDOWS
# include <stdbool.h>
# include <windows.h>

typedef CRITICAL_SECTION pthread_mutex_t;
typedef void pthread_mutexattr_t;
typedef void pthread_condattr_t;
// typedef void pthread_rwlockattr_t;
typedef HANDLE pthread_t;
typedef CONDITION_VARIABLE pthread_cond_t;

// int pthread_create(pthread_t *thread, pthread_attr_t *attr, void *(*start_routine)(void *), void *arg);
// int pthread_join(pthread_t thread, void **value_ptr);
// int pthread_detach(pthread_t);

int pthread_mutex_init(pthread_mutex_t *mutex, pthread_mutexattr_t *attr);
int pthread_mutex_destroy(pthread_mutex_t *mutex);
int pthread_mutex_lock(pthread_mutex_t *mutex);
int pthread_mutex_trylock(pthread_mutex_t *mutex);
int pthread_mutex_unlock(pthread_mutex_t *mutex);

int pthread_cond_init(pthread_cond_t *cond, pthread_condattr_t *attr);
int pthread_cond_destroy(pthread_cond_t *cond);
int pthread_cond_wait(pthread_cond_t *cond, pthread_mutex_t *mutex);
int pthread_cond_timedwait(pthread_cond_t *cond, pthread_mutex_t *mutex, const struct timespec *abstime);
int pthread_cond_signal(pthread_cond_t *cond);
int pthread_cond_broadcast(pthread_cond_t *cond);

#else
// on posix -> just use pthread.h
# include <pthread.h>
#endif

#define BEGIN_SINGLE_THREAD_REGION(wcs_ptr, mutex_ptr, control_var)     \
    if(pthread_mutex_trylock(mutex_ptr) == 0)                           \
    {                                                                   \
        if(abs(wcs->flag) != WCSSET) control_var = 0;                   \
        if(control_var == 0)                                            \
        {



#define END_SINGLE_THREAD_REGION(wcs_ptr, mutex_ptr, cond_ptr, control_var) \
        }    /* closes the control_var == 0 if cond */                      \
        if(control_var == 0 && wcs_ptr->naxis > 0) {                        \
            control_var = 1;                                                \
        }                                                                   \
        pthread_cond_broadcast(cond_ptr);                                   \
        pthread_mutex_unlock(mutex_ptr);                                    \
    } else {  /* the else for the mutex_trylock */                          \
        pthread_mutex_lock(mutex_ptr);                                      \
        while (control_var != 1) {                                          \
            pthread_cond_wait(cond_ptr, mutex_ptr);                         \
        }                                                                   \
        pthread_mutex_unlock(mutex_ptr);                                    \
    }


#ifdef __cplusplus
}
#endif
