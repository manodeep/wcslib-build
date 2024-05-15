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

#ifdef ON_WINDOWS
#include "wcs_pthreads.h"

#include <time.h>

struct timespec {
    long tv_sec;
    long tv_nsec;
};

static DWORD timespec_to_ms(const struct timespec *abstime);


int pthread_mutex_init(pthread_mutex_t *mutex, pthread_mutexattr_t *attr)
{
    (void)attr;

    if (mutex == NULL)
        return 1;

    InitializeCriticalSection(mutex);
    return 0;
}

int pthread_mutex_destroy(pthread_mutex_t *mutex)
{
    if (mutex == NULL)
        return 1;
    DeleteCriticalSection(mutex);
    return 0;
}

int pthread_mutex_lock(pthread_mutex_t *mutex)
{
    if (mutex == NULL)
        return 1;
    EnterCriticalSection(mutex);
    return 0;
}

int pthread_mutex_trylock(pthread_mutex_t *mutex)
{
    if (mutex == NULL)
        return 1;

    return TryEnterCriticalSection(mutex);
}

int pthread_mutex_unlock(pthread_mutex_t *mutex)
{
    if (mutex == NULL)
        return 1;
    LeaveCriticalSection(mutex);
    return 0;
}

int pthread_cond_init(thread_cond_t *cond, pthread_condattr_t *attr)
{
    (void)attr;
    if (cond == NULL)
        return 1;
    InitializeConditionVariable(cond);
    return 0;
}

int pthread_cond_destroy(thread_cond_t *cond)
{
    /* Windows does not have a destroy for conditionals */
    (void)cond;
    return 0;
}

int pthread_cond_wait(thread_cond_t *cond, pthread_mutex_t *mutex)
{
    if (cond == NULL || mutex == NULL)
        return 1;
    return pthread_cond_timedwait(cond, mutex, NULL)
}

int pthread_cond_timedwait(thread_cond_t *cond, pthread_mutex_t *mutex,
                           const struct timespec *abstime)
{
    if (cond == NULL || mutex == NULL)
        return 1;
    if (!SleepConditionVariableCS(cond, mutex, timespec_to_ms(abstime)))
        return 1;
    return 0;
}

int pthread_cond_signal(thread_cond_t *cond)
{
    if (cond == NULL)
        return 1;
    WakeConditionVariable(cond);
    return 0;
}

int pthread_cond_broadcast(thread_cond_t *cond)
{
    if (cond == NULL)
        return 1;
    WakeAllConditionVariable(cond);
    return 0;
}

static DWORD timespec_to_ms(const struct timespec *abstime)
{
    DWORD t;

    if (abstime == NULL)
        return INFINITE;

    t = ((abstime->tv_sec - time(NULL)) * 1000) + (abstime->tv_nsec / 1000000);
    if (t < 0)
        t = 1;
    return t;
}


#endif