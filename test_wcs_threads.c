#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <assert.h>
#include <string.h>
#include <time.h>

#include "wcslib.h"

typedef int (*roundtrip_coords_impls)(const int N, double *pixcrd, double *imgcrd, double *phi, double*theta, double *world, int *stat, double *dest, struct wcsprm *wcs, const int nthreads);


int set_wcs(struct wcsprm *wcs)
{
    //This replicates the wcs setup in astropy issue 16245
    assert(wcs != NULL);
    assert(wcs->naxis == 2);
    wcs->crpix[0] = -234.75;
    wcs->crpix[1] = 8.3393;
    wcs->cdelt[0] = -0.066667;
    wcs->cdelt[1] = 0.066667;
    wcs->crval[0] = 0.0;
    wcs->crval[1] = -90.0;
    sprintf(wcs->ctype[0], "RA---AIR");
    sprintf(wcs->ctype[1], "DEC--AIR");
#ifdef USE_FLAG_TO_BYPASS
#warning "Using the flag to bypass the wcsset() call in wcslib-8.3 onwards"
    wcs->flag = 1;
#else
    wcs->flag = 0;
#endif
    return wcsset(wcs);
}

int setup_random_pixcrd(const int N, const double min, const double max, double *pixcrd)
{
    unsigned int seed = 292364782;
    srand(seed);
    for(int i=0;i<N;i++) {
        double r = (double)rand()/(double)RAND_MAX;
        pixcrd[i] = r*(max - min) + min;
        r = (double)rand()/(double)RAND_MAX;
        pixcrd[i+N] = r*(max - min) + min;
    }
    return EXIT_SUCCESS;
}

int distribute_compute_over_ntasks(const int totn, const int NTasks, const int ThisTask, int *n_thistask, int *start)
{
    if(ThisTask > NTasks || ThisTask < 0 || NTasks < 1) {
        fprintf(stderr,"Error: ThisTask = %d and NTasks = %d must satisfy i) ThisTask < NTasks, ii) ThisTask > 0 and iii) NTasks >= 1\n",
                ThisTask, NTasks);
        return EXIT_FAILURE;
    }

    if(totn < 0) {
        fprintf(stderr,"Error: On ThisTask = %d: total number = %d must be >= 0\n", ThisTask, totn);
        return EXIT_FAILURE;
    }

    if(totn == 0) {
        fprintf(stderr,"Warning: Got 0 forests to distribute! Returning\n");
        *n_thistask = 0;
        *start = 0;
        return EXIT_SUCCESS;
    }


    // Assign each task an equal number of compute units. If we can't equally assign each task
    // the EXACT same number of compute, give each task an extra unit (if required).
    const int n_per_cpu =  totn/NTasks;
    const int rem_n = totn % NTasks;
    int n_this_task = n_per_cpu;
    if(ThisTask < rem_n) {
        n_this_task++;
    }

    int64_t start_num = n_per_cpu * ThisTask;
    if(ThisTask < rem_n) {
        start_num += ThisTask;
    } else {
        start_num += rem_n;  // All tasks that weren't given an extra forest will be offset by a constant amount.
    }

    /* Now fill up the destination */
    *n_thistask = n_this_task;
    *start = start_num;

    return EXIT_SUCCESS;
}


int roundtrip_coords_buggy(const int N, double *pixcrd, double *imgcrd, double *phi, double *theta, double *world, int *stat, double *dest, struct wcsprm *wcs, const int nthreads)
{
    const int naxis = wcs->naxis;
    fprintf(stderr,"Using %s with nthreads = %d\n", __FUNCTION__, nthreads);
    #pragma omp parallel num_threads(nthreads) shared(wcs, pixcrd, imgcrd, phi, theta, world, stat, dest)
    {
        if(omp_get_num_threads() != nthreads){
            fprintf(stderr,"Error: Expected %d threads, got %d threads\n", nthreads, omp_get_num_threads());
            exit(EXIT_FAILURE);
        }
        const int tid = omp_get_thread_num();
        int start, n_thistask;
        distribute_compute_over_ntasks(N, nthreads, tid, &n_thistask, &start);
        int end = start + n_thistask;

        const int offset = 2*start;
        double *local_pixcrd = pixcrd + offset;
        double *local_imgcrd = imgcrd + offset;
        double *local_phi = phi + offset;
        double *local_theta = theta + offset;
        double *local_world = world + offset;
        int *local_stat = stat + offset;
        double *local_dest = dest + offset;
        fprintf(stderr,"On thread %d: start = %d, end = %d n_thistask = %d\n", tid, start, end, n_thistask);
        int status = wcsp2s(wcs, n_thistask, naxis, local_pixcrd, local_imgcrd, local_phi, local_theta, local_world, local_stat);
        if(status != EXIT_SUCCESS) {
            fprintf(stderr, "wcsp2s failed with status %d\n", status);
        }

        //calling wcsset() causes the race condition (needs nowait to trigger it)
        wcsset(wcs);

        //convert back to pixel coordinates
        status = wcss2p(wcs, n_thistask, naxis, local_world, local_phi, local_theta, local_imgcrd, local_pixcrd, local_stat);
        if(status != EXIT_SUCCESS) {
            fprintf(stderr, "wcss2p failed with status %d\n", status);
        }

        memcpy(local_dest, local_pixcrd, 2*n_thistask*sizeof(*local_dest));
    }
    fprintf(stderr,"Using %s with nthreads = %d...done\n", __FUNCTION__, nthreads);
    return EXIT_SUCCESS;
}




int roundtrip_coords_omp_barrier(const int N, double *pixcrd, double *imgcrd, double *phi, double *theta, double *world, int *stat, double *dest, struct wcsprm *wcs, const int nthreads)
{
    const int naxis = wcs->naxis;
    fprintf(stderr,"Using %s with nthreads = %d\n", __FUNCTION__, nthreads);
    #pragma omp parallel num_threads(nthreads) shared(wcs, pixcrd, imgcrd, phi, theta, world, stat, dest)
    {
        assert(omp_get_num_threads() == nthreads);
        const int tid = omp_get_thread_num();
        int start, n_thistask;
        distribute_compute_over_ntasks(N, nthreads, tid, &n_thistask, &start);
        int end = start + n_thistask;

        const int offset = 2*start;
        double *local_pixcrd = pixcrd + offset;
        double *local_imgcrd = imgcrd + offset;
        double *local_phi = phi + offset;
        double *local_theta = theta + offset;
        double *local_world = world + offset;
        int *local_stat = stat + offset;
        double *local_dest = dest + offset;
        fprintf(stderr,"On thread %d: start = %d, end = %d n_thistask = %d\n", tid, start, end, n_thistask);

        int status = wcsp2s(wcs, n_thistask, naxis, local_pixcrd, local_imgcrd, local_phi, local_theta, local_world, local_stat);
        if(status != EXIT_SUCCESS) {
            fprintf(stderr, "wcsp2s failed with status %d\n", status);
        }

        //calling wcsset() causes the race condition (needs nowait to trigger it)
        #pragma omp barrier
        #pragma omp single
        wcsset(wcs);
        #pragma omp barrier

        //convert back to pixel coordinates
        status = wcss2p(wcs, n_thistask, naxis, local_world, local_phi, local_theta, local_imgcrd, local_pixcrd, local_stat);
        if(status != EXIT_SUCCESS) {
            fprintf(stderr, "wcss2p failed with status %d\n", status);
        }
        memcpy(local_dest, local_pixcrd, 2*n_thistask*sizeof(*local_dest));
    }
    fprintf(stderr,"Using %s with nthreads = %d...done\n", __FUNCTION__, nthreads);
    return EXIT_SUCCESS;
}


int count_n_mismatches(const int N, const double *input, const double *output)
{
    int n_mismatches = 0;
    #pragma omp simd reduction(+:n_mismatches)
    for(int i=0;i<2*N;i++) {
        if(fabs(input[i] - output[i]) > 1e-8) {
            n_mismatches++;
            // fprintf(stderr, "Mismatch at i=%d: input = %g, output = %g fabs(diff) = %0.10e\n",
            //         i, input[i], output[i], fabs(input[i] - output[i]));
        }
    }
    return n_mismatches/2;
}

int main(int argc, char **argv)
{
    int N = 100;
    int max_numthreads = 1;
    if(argc >= 2) {
        N = atoi(argv[1]);
    }

    if(argc >= 3) {
        max_numthreads = atoi(argv[2]);
    }

    struct wcsprm wcs;
    memset(&wcs, 0, sizeof(wcs));
    const int naxis = 2;
    wcsini(1, naxis, &wcs);
    set_wcs(&wcs);

    double *pixcrd = malloc(2*N*sizeof(*pixcrd));
    double *imgcrd = malloc(2*N*sizeof(*imgcrd));
    double *phi = malloc(2*N*sizeof(*phi));
    double *theta = malloc(2*N*sizeof(*theta));
    double *world = malloc(2*N*sizeof(*world));
    int *stat = malloc(2*N*sizeof(*stat));
    double *dest = malloc(2*N*sizeof(*dest));
    double *pixcrd_orig = malloc(2*N*sizeof(*pixcrd_orig));

    if(pixcrd == NULL || dest == NULL ||
       imgcrd == NULL || phi == NULL  ||
       theta == NULL || world == NULL ||
       stat == NULL || pixcrd_orig == NULL)
    {
        fprintf(stderr, "Memory allocation failed for N = %d\n", N);
        return EXIT_FAILURE;
    }
    // const roundtrip_coords_impls all_functions[] = { &roundtrip_coords_omp_barrier, &roundtrip_coords_buggy};
    // const char *function_names[72] = { "roundtrip_coords_omp_barrier", "roundtrip_coords_buggy"};
    // const int n_functions = sizeof(all_functions)/sizeof(all_functions[0]);

    const roundtrip_coords_impls all_functions[] = { &roundtrip_coords_buggy};
    const char *function_names[72] = { "roundtrip_coords_buggy"};
    const int n_functions = sizeof(all_functions)/sizeof(all_functions[0]);

    const double min = -1000.0;
    const double max = 1000.0;
    setup_random_pixcrd(N, min, max, pixcrd);
    memcpy(pixcrd_orig, pixcrd, 2*N*sizeof(*pixcrd));

    int64_t num_wrong = 0;
    for(int nthreads=1;nthreads<=max_numthreads;nthreads++) {
        for(int ifunc=0;ifunc<n_functions;ifunc++) {
            struct timespec t0, t1;
            clock_gettime(CLOCK_MONOTONIC, &t0);
            int status = all_functions[ifunc](N, pixcrd, imgcrd, phi, theta, world, stat, dest, &wcs, nthreads);
            if(status != EXIT_SUCCESS) {
                fprintf(stderr, "%s failed with status %d\n", function_names[ifunc], status);
                return status;
            }
            clock_gettime(CLOCK_MONOTONIC, &t1);
            const double func_time = (t1.tv_sec - t0.tv_sec) + 1e-9*(t1.tv_nsec - t0.tv_nsec);
            const int n_mismatches = count_n_mismatches(N, pixcrd_orig, dest);
            num_wrong += n_mismatches;
            fprintf(stderr,"Called function = '%s' with nthreads=%d, N=%d, n_mismatches = %d. Time taken = %g seconds\n",
                         function_names[ifunc], nthreads, N,  n_mismatches, func_time);
            memcpy(pixcrd, pixcrd_orig, 2*N*sizeof(*pixcrd));
        }
    }
    wcsfree(&wcs);

    free(pixcrd);
    free(world);
    free(theta);
    free(imgcrd);
    free(phi);
    free(dest);
    free(stat);
    free(pixcrd_orig);

    return num_wrong;
}
//I created the file in the cextern/wcslib/C directory and compiled with:
//clang -g -I.. -Wall -Wextra -fsanitize=leak -fsanitize=undefined -fsanitize=bounds -fsanitize=address -fsanitize-undefined-trap-on-error -fstack-protector-all test_wcs_threads.c -std=c99 -O3 -Wl,-no_compact_unwind -fopenmp=libomp libwcs-PIC.a -o test_wcs_threads

//Running on a single thread with 1 million points produces same results after roundtrip (tol = 1e-8):
/*
On thread 0: start = 0, end = 10000000 n_thistask = 10000000
Called function = 'roundtrip_coords_omp_barrier' with nthreads=1, N=10000000, n_mismatches = 0. Time taken = 2.9328 seconds
On thread 0: start = 0, end = 10000000 n_thistask = 10000000
Called function = 'roundtrip_coords_buggy' with nthreads=1, N=10000000, n_mismatches = 0. Time taken = 2.76603 seconds
*/

//Running on multiple threads sometimes produces a crash with double-free, and sometimes the n_mismatches > 0:
/*
(.venv) [~/sorsery-consulting/2024/reproject-astropy-wcs/astropy/cextern/wcslib/C @HAW10013056] ./test_wcs_threads 100000 2
On thread 0: start = 0, end = 100000 n_thistask = 100000
Called function = 'roundtrip_coords_omp_barrier' with nthreads=1, N=100000, n_mismatches = 0. Time taken = 0.029945 seconds
On thread 0: start = 0, end = 100000 n_thistask = 100000
Called function = 'roundtrip_coords_buggy' with nthreads=1, N=100000, n_mismatches = 0. Time taken = 0.028093 seconds
On thread 0: start = 0, end = 50000 n_thistask = 50000
On thread 1: start = 50000, end = 100000 n_thistask = 50000
Called function = 'roundtrip_coords_omp_barrier' with nthreads=2, N=100000, n_mismatches = 95872. Time taken = 0.016187 seconds
On thread 0: start = 0, end = 50000 n_thistask = 50000
On thread 1: start = 50000, end = 100000 n_thistask = 50000
=================================================================
==67903==ERROR: AddressSanitizer: attempting double-free on 0x0001034001b0 in thread T0:
    #0 0x100ce19f0 in wrap_free+0x8c (libclang_rt.asan_osx_dynamic.dylib:arm64+0x3d9f0) (BuildId: e7b5fb9640823e93be90c001748b8fc332000000200000000100000000000b00)
    #1 0x10055822c in linset lin.c:747
    #2 0x100574074 in wcsset.part.0 wcs.c:2807
    #3 0x10054bc70 in .omp_outlined. test_wcs_threads.c:89
    #4 0x10081daf8 in __kmp_invoke_microtask+0x98 (libomp.dylib:arm64+0x81af8) (BuildId: 583abbdbdb43371e8f4bf9100abf532b32000000200000000100000000000c00)

0x0001034001b0 is located 0 bytes inside of 16-byte region [0x0001034001b0,0x0001034001c0)
freed by thread T1 here:
    #0 0x100ce19f0 in wrap_free+0x8c (libclang_rt.asan_osx_dynamic.dylib:arm64+0x3d9f0) (BuildId: e7b5fb9640823e93be90c001748b8fc332000000200000000100000000000b00)
    #1 0x10055822c in linset lin.c:747
    #2 0x100574074 in wcsset.part.0 wcs.c:2807
    #3 0x10054bc70 in .omp_outlined. test_wcs_threads.c:89
    #4 0x10081daf8 in __kmp_invoke_microtask+0x98 (libomp.dylib:arm64+0x81af8) (BuildId: 583abbdbdb43371e8f4bf9100abf532b32000000200000000100000000000c00)

previously allocated by thread T0 here:
    #0 0x100ce1c68 in wrap_calloc+0x90 (libclang_rt.asan_osx_dynamic.dylib:arm64+0x3dc68) (BuildId: e7b5fb9640823e93be90c001748b8fc332000000200000000100000000000b00)
    #1 0x100558238 in linset lin.c:748
    #2 0x100574074 in wcsset.part.0 wcs.c:2807
    #3 0x10054c830 in .omp_outlined..13 test_wcs_threads.c:131
    #4 0x10081daf8 in __kmp_invoke_microtask+0x98 (libomp.dylib:arm64+0x81af8) (BuildId: 583abbdbdb43371e8f4bf9100abf532b32000000200000000100000000000c00)

Thread T1 created by T0 here:
    #0 0x100cdbc10 in wrap_pthread_create+0x50 (libclang_rt.asan_osx_dynamic.dylib:arm64+0x37c10) (BuildId: e7b5fb9640823e93be90c001748b8fc332000000200000000100000000000b00)
    #1 0x1007fff10 in __kmp_create_worker+0xcc (libomp.dylib:arm64+0x63f10) (BuildId: 583abbdbdb43371e8f4bf9100abf532b32000000200000000100000000000c00)
    #2 0x1007be35c in __kmp_allocate_thread+0x410 (libomp.dylib:arm64+0x2235c) (BuildId: 583abbdbdb43371e8f4bf9100abf532b32000000200000000100000000000c00)
    #3 0x1007b8f34 in __kmp_allocate_team+0x96c (libomp.dylib:arm64+0x1cf34) (BuildId: 583abbdbdb43371e8f4bf9100abf532b32000000200000000100000000000c00)
    #4 0x1007ba980 in __kmp_fork_call+0x13b4 (libomp.dylib:arm64+0x1e980) (BuildId: 583abbdbdb43371e8f4bf9100abf532b32000000200000000100000000000c00)
    #5 0x1007ade9c in __kmpc_fork_call+0xa8 (libomp.dylib:arm64+0x11e9c) (BuildId: 583abbdbdb43371e8f4bf9100abf532b32000000200000000100000000000c00)
    #6 0x10054c220 in roundtrip_coords_omp_barrier test_wcs_threads.c:131
    #7 0x10054d354 in main test_wcs_threads.c:233
    #8 0x1008b9088 in start+0x204 (dyld:arm64+0x5088) (BuildId: 75627683a78032adae34cf86dd23a26b32000000200000000100000000050c00)
    #9 0xf12efffffffffffc  (<unknown module>)

SUMMARY: AddressSanitizer: double-free (libclang_rt.asan_osx_dynamic.dylib:arm64+0x3d9f0) (BuildId: e7b5fb9640823e93be90c001748b8fc332000000200000000100000000000b00) in wrap_free+0x8c
==67903==ABORTING
Abort trap: 6
*/
