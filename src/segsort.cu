#include "kernels/segmentedsort.cuh"
#include "kernels/scan.cuh"
#include "device/ctasegsort.cuh"
#include "blast.h"
#include <inttypes.h>
#include "segsort.h"

using namespace mgpu;
template<typename Op, typename T1>
void Scan(T1 *source, T1 *dest, int count, typename Op::result_type identity,
        Op op, CudaContext& context) {

    typedef typename Op::result_type T2;
    Scan<MgpuScanTypeInc>(source, count, 
            identity, op, (T2*)0, (T2*)0, dest, context);
}


__global__ void filter_seq(
        uint64_t *hits, 
        int *binOffset, 
        int *numHits, 
        uint64_t *output, 
        int *numHits_out, 
        int num_bins, 
        int parameters_A, 
        int parameters_overlap)
{
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    for(int ii = tid; ii < num_bins; ii += gridDim.x * blockDim.x)
    {
        int num_extension = 0;
        int start = binOffset[ii] - numHits[ii]; // start_of_bin 
        int end = binOffset[ii]; // end_of_bin
        uint64_t prev = hits[start];
        uint64_t *output_w = output + start; 
        for(int jj = start + 1; jj < end; jj++)
        {
            uint64_t curr = hits[jj];
            uint64_t distance = curr - prev;
            if(distance >= parameters_A)
            {
                prev = curr;
            }
            else if(distance >= parameters_overlap)
            {
                output_w[num_extension] = prev;
                output_w[num_extension + 1] = curr;
                num_extension += 2;
                prev = curr;
            }
        }
        numHits_out[ii] = num_extension;
    }
}


int hit_sort_filter(
        int *blast_numHitsD, 
        int *blast_numExtD,
        int *binOffsetD,
        uint64_t *HitInfoD,
        int block_size, 
        int num_block, 
        int bin_x,
        int total_numHits,
        int parametersA,
        int parameters_overlap
        )
{

    ContextPtr context = CreateCudaDevice(0);
    struct timeval s0, s1;

    int num_bins = block_size * num_block * bin_x;
    cudaError_t cudaRes;

    gettimeofday(&s0, NULL);
    SegSortKeysFromIndices(HitInfoD, 
            total_numHits, 
            binOffsetD, 
            num_bins - 1, 
            *context);
    gettimeofday(&s1, NULL);

#ifdef VERBOSE
    printf("Hit SegSort Time: %f\n", (float)(1000000 * (s1.tv_sec - s0.tv_sec) + s1.tv_usec - s0.tv_usec)/1000000);
#endif
    cudaRes = cudaGetLastError(); 
    if( cudaRes != cudaSuccess)
    {
        printf("CUDA error: %s in %d\n", cudaGetErrorString(cudaRes), __LINE__);
        exit(-1);
    }

    gettimeofday(&s0, NULL);
    filter_seq<<<num_block, block_size>>>(
            HitInfoD, 
            binOffsetD, 
            blast_numHitsD, 
            HitInfoD, 
            blast_numExtD, 
            num_bins, 
            parametersA, 
            parameters_overlap);
    cudaRes = cudaGetLastError(); 
    if( cudaRes != cudaSuccess)
    {
        printf("CUDA error: %s in %d\n", cudaGetErrorString(cudaRes), __LINE__);
        exit(-1);
    }
    gettimeofday(&s1, NULL);

#ifdef VERBOSE
    printf("Hit Filter Time: %f\n", (float)(1000000 * (s1.tv_sec - s0.tv_sec) + s1.tv_usec - s0.tv_usec)/1000000);
#endif
    return 0;
}

int get_total_numHits(int *numHits, int num_bins)
{
    ContextPtr context = CreateCudaDevice(0);
    return Reduce(numHits, num_bins, *context);
}

void get_bin_offset(int *numHits, int *binOffset, int num_bins)
{
    ContextPtr context = CreateCudaDevice(0);
    typedef mgpu::plus<int> Op;
    Scan<Op, int>(numHits, binOffset, num_bins, (int)0, Op(),
            *context);
}


