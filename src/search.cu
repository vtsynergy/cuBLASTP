/*
 * cuBLASTP - Fine-Grained Parallelization of Protein Sequence Search on CPU+GPU 
 * Version 0.1 (beta)
 *
 * (c) 2015 Virginia Polytechnic Institute & State University (Virginia Tech) 
 * This version of cuBLASTP is licensed for non-commercial use only,
 * as specified in LICENSE files in licensing directory. For all other use
 * contact vtiplicensing@vtip.org
 * 
 * Developer: Jing Zhang 
 *
 */


#include "blast.h"
#include "wordLookupDFA.h"
#include <math.h>
#include <limits.h>
#include <sys/time.h>
#include <inttypes.h>
#include "segsort.h"

#ifdef SM_20
#define __ldg(x) (*(x))
#define NUM_BLOCK 112
#else
#define NUM_BLOCK 260
#endif

#define BLOCK_SIZE 128
#define BIN_X 4      // CHANGE BIN_SIZE
#define BIN_POWER 7  // BIN_X
#define BIN_MARK 127 // BIN_X

extern unsigned char *wordLookupDFA;
extern struct groupFP *wordLookupDFA_groupsFP;

TIMERECORD timeRecord;

struct parameters {
    char parameters_wordSize;
    unsigned char encoding_numCodes;
    char parameters_overlap;
    int4 wordLookupDFA_numCodes;
    uint4 additionalQueryPositionOffset;
    int4 statistics_ungappedNominalDropoff;
    int4 blast_ungappedNominalTrigger;
    int4 parameters_A;
    uint4 ungappedExtensionsPerThread;
    uint4 ungappedExtAdditionalStartLoc;
};

#define TARGET_THREAD 0
#define UNGAPEXT_PER_THREAD 150
#define TOTAL_UNGAPPED_EXT 500000

__device__ __constant__ int2 scoreMatrixC[1640];
__device__ __constant__ unsigned char querySequenceC[40000];

int compare_ungappedextension(const void *a, const void *b) {
    return ((struct ungappedExtension *)a)->sequenceCount -
        ((struct ungappedExtension *)b)->sequenceCount;
}
int findStartLoc(struct ungappedExtension *ungappedExtensionsPtr, int threadNo,
        int itemNum);
__device__ struct ungappedExtension *ungappedExtension_oneHitExtendD(
        unsigned char *, int4, unsigned char *, struct PSSMatrixFP, unsigned char *,
        uint4 *, unsigned char, int4, int4, struct ungappedExtension *, uint4 *,
        uint4, int4);

__device__ uint4 global_numAdditionalTriggerExtensions;

__device__ struct coordinate ungappedExtension_findProteinSeed(
        struct ungappedExtension *ungappedExtension, struct PSSMatrixFP PSSMatrixFP,
        unsigned char *subject, unsigned char encoding_numCodes);

__device__ struct coordinate ungappedExtension_findProteinSeed_sm(
        struct ungappedExtension *ungappedExtension,
        // struct PSSMatrixFP PSSMatrixFP,
        unsigned char *subject, unsigned char encoding_numCodes,
        unsigned char *querySequence, int2 *scoreMatrix);

__global__ void ungappedExtension_twoHitExtendG_findProteinSeed(
        struct parameters *parametersFP_g, struct PSSMatrixFP *PSSMatrixFP_g,
        unsigned char *sequence, struct sequenceDataFP *sequenceDataFP,
        struct ungappedExtension *ungappedExtension_extensionsp,
        uint4 *blast_numTriggerExtensions);

__global__ void ungappedExtension_twoHitExtendG_bin_sorted_sm_s(
        struct parameters *parametersFP_g,
        // char *matrixBody_g,
        uint64_t *HitInfo_g,
        // uint2 *num_hits,
        uint4 *blast_numUngappedExtensions, unsigned char *sequence,
        // unsigned char **hitMatrix_furthestp,
        // uint4 *hitMatrix_offsetp,
        struct PSSMatrixFP *PSSMatrixFP_g, struct sequenceDataFP *sequenceDataFP,
        struct ungappedExtension *ungappedExtension_extensionsp,
        uint4 *blast_numTriggerExtensions, int *numOneHitsD, int *numExtD,
        int *binOffset_g);

__device__ void ungappedExtension_twoHitExtendD_sm(
        unsigned char *sequenceStart, int4 queryOffset, unsigned char *subjectHit,
        unsigned char *lastHitFP,
        // struct PSSMatrixFP *PSSMatrixFP,
        unsigned char *subject, unsigned char **sequenceHitEnd,
        unsigned char encoding_numCodes, int4 statistics_ungappedNominalDropoff,
        int4 blast_ungappedNominalTrigger, int4 ungappedExtensionsPerThread,
        struct ungappedExtension *ungappedExtension_extensions,
        struct ungappedExtension *ungappedExtension_additonal,
        uint4 *numOfTriggerExtensions, uint4 sequenceCount, int2 *scoreMatrix,
        unsigned char *querySequence) {
    // int queryPosition;
    unsigned char *subjectPosition, *subjectStart, *subjectEnd;
    int4 changeSinceBest = 0;
    int4 dropoff, originalDropoff;
    int4 ungappedExtension_bestScore;

    originalDropoff = dropoff = -statistics_ungappedNominalDropoff;
    ungappedExtension_bestScore = 0;

    // Start at queryEnd,subjectEnd (right/last hit position)
    // queryPosition = matrix + queryOffset * encoding_numCodes;
    subjectPosition = subjectStart = subjectHit;
    int2 queryPosition = queryOffset + 1;

    while (changeSinceBest > dropoff) {
        // changeSinceBest += queryPosition[__ldg(subjectPosition)];
        changeSinceBest += scoreMatrix
            [querySequence[queryPosition] * encoding_numCodes + (*subjectPosition)];

        // If we have got a positive score
        if (changeSinceBest > 0) {
            // Keep updating best score and resetting change-since-best
            // whilst we are reading positive scores
            do {
                ungappedExtension_bestScore += changeSinceBest;
                // queryPosition = queryPosition - encoding_numCodes;
                queryPosition = queryPosition - 1;
                subjectPosition--;

                // changeSinceBest = queryPosition[__ldg(subjectPosition)];
                changeSinceBest =
                    scoreMatrix[querySequence[queryPosition] * encoding_numCodes +
                    (*subjectPosition)];
            } while (changeSinceBest > 0);

            subjectStart = subjectPosition;
        }

        // queryPosition = queryPosition - encoding_numCodes;
        queryPosition = queryPosition - 1;
        subjectPosition--;
    }

    // Correct for extra decrement
    subjectStart++;

    if (subjectStart > lastHitFP) {
        *sequenceHitEnd = subjectHit;
        return;
    }

    // Starting at right/last hit position again
    // queryPosition = matrix + (queryOffset + 1) * encoding_numCodes;
    queryPosition = (queryOffset + 2);
    subjectPosition = subjectHit + 1;
    subjectEnd = subjectHit;
    changeSinceBest = 0;

    // May need to alter dropoff so we also dropoff if below zero
    if (-ungappedExtension_bestScore > originalDropoff) {
        dropoff = -ungappedExtension_bestScore;
    }

    // Extend end of alignment until dropoff
    while (changeSinceBest > dropoff) {
        // Shucai
        // changeSinceBest += queryPosition[__ldg(subjectPosition)];
        changeSinceBest += scoreMatrix
            [querySequence[queryPosition] * encoding_numCodes + (*subjectPosition)];

        // If we have got a positive score
        if (changeSinceBest > 0) {
            // Keep updating best score and resetting change-since-best
            // whilst we are reading positive scores
            do {
                ungappedExtension_bestScore += changeSinceBest;
                // queryPosition = queryPosition + encoding_numCodes;
                queryPosition = queryPosition + 1;
                subjectPosition++;
                // changeSinceBest = queryPosition[__ldg(subjectPosition)];
                changeSinceBest =
                    scoreMatrix[querySequence[queryPosition] * encoding_numCodes +
                    (*subjectPosition)];
            } while (changeSinceBest > 0);

            subjectEnd = subjectPosition;

            // Check need for change in dropoff
            if ((dropoff = -ungappedExtension_bestScore) < originalDropoff) {
                dropoff = originalDropoff;
            }
        }
        // queryPosition = queryPosition + encoding_numCodes;
        queryPosition = queryPosition + 1;
        subjectPosition++;
    }

    subjectEnd--;

    *sequenceHitEnd = subjectEnd;

    if (ungappedExtension_bestScore >= blast_ungappedNominalTrigger) {
        int2 diagonal;
        struct ungappedExtension *newUngappedExtension = NULL;
        newUngappedExtension =
            *numOfTriggerExtensions >= ungappedExtensionsPerThread
            ? &ungappedExtension_additonal
            [atomicAdd(&global_numAdditionalTriggerExtensions, 1)]
            : &ungappedExtension_extensions
            [atomicAdd(numOfTriggerExtensions, 1)];
        // newUngappedExtension =
        // &ungappedExtension_extensions[*numOfTriggerExtensions];

        // Calculate diagonal
        diagonal = (subjectHit - subject) - queryOffset;
        // Determine offsets from pointers
        newUngappedExtension->start.subjectOffset = subjectStart - subject;
        newUngappedExtension->end.subjectOffset = subjectEnd - subject;
        newUngappedExtension->start.queryOffset =
            newUngappedExtension->start.subjectOffset - diagonal;
        newUngappedExtension->end.queryOffset =
            newUngappedExtension->end.subjectOffset - diagonal;
        newUngappedExtension->seed = ungappedExtension_findProteinSeed_sm(
                newUngappedExtension, subject, encoding_numCodes, querySequence,
                scoreMatrix);
        newUngappedExtension->next = NULL;
        newUngappedExtension->nominalScore = ungappedExtension_bestScore;
        newUngappedExtension->status = ungappedExtension_UNGAPPED;
        newUngappedExtension->sequenceCount = sequenceCount;
        // newUngappedExtension->tid = tid;

        // Shucai
        // Record the number of hits satisfying the next step
        //(*numOfTriggerExtensions)++;
    }
}

__global__ void search_protein_hit_detection_prescan(
        struct sequenceDataFP *sequenceDataFP, unsigned char *sequence,
        struct parameters *parametersFP_g, struct groupFP *wordLookupDFA_groupsFP,
        unsigned char *wordLookupDFAFP, int *blast_numHits_g,
        uint4 nTotalSequenceNum, unsigned int groupNum) {
    unsigned int tid = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int tt = threadIdx.x;
    unsigned int warpId = tid >> 5;
    unsigned char laneId = tt & 31;

    extern __shared__ unsigned char DFA_groupFP_array[];

    __shared__ unsigned int bin_numHits_share[BLOCK_SIZE * BIN_X];

    unsigned int *bin_numHits = bin_numHits_share + (tt >> 5 << BIN_POWER);

    for (int ii = 0; ii < BIN_X; ii++)
        bin_numHits[laneId + ii * 32] = 0;

    int *blast_numHits = blast_numHits_g + (warpId << BIN_POWER);
    unsigned char *subject, *sequenceEnd, *address, *start;
    int4 subjectOffset, count;
    unsigned char currentWord;
    const unsigned char *currentBlock;
    uint2 *queryOffsets, queryOffset;

    struct parameters parametersFP = *parametersFP_g;
    extern __shared__ struct groupFP DFA_groupFP_s[];

    for (unsigned int ii = tt; ii < groupNum; ii += blockDim.x) {
        DFA_groupFP_s[ii] = wordLookupDFA_groupsFP[ii];
    }
    __syncthreads();

    uint2 *wordLookupDFA_AddPositions =
        (uint2 *)((char *)wordLookupDFAFP +
                parametersFP.additionalQueryPositionOffset);

    uint4 sequenceCount = warpId;

    while (sequenceCount < nTotalSequenceNum) {
        if (sequenceDataFP[sequenceCount].sequenceLength >=
                parametersFP.parameters_wordSize) {
            start = subject = sequence + sequenceDataFP[sequenceCount].offset;
            sequenceEnd = subject + sequenceDataFP[sequenceCount].sequenceLength;
            address = start + laneId;

            while (address + parametersFP.parameters_wordSize - 1 < sequenceEnd) {
                struct groupFP *currentGroupFP = DFA_groupFP_s;
                unsigned char letter = *address;
                for (count = 1; count < parametersFP.parameters_wordSize; count++) {
                    currentGroupFP =
                        letter < parametersFP.wordLookupDFA_numCodes
                        ? &DFA_groupFP_s[currentGroupFP->nextGroups + letter]
                        : &DFA_groupFP_s[currentGroupFP->nextGroups];
                    address++;
                    letter = *address;
                }

                currentBlock = &wordLookupDFAFP[currentGroupFP->nextWords];

                // If current code is a regular letter
                currentWord = letter < parametersFP.wordLookupDFA_numCodes
                    ? __ldg(currentBlock + letter)
                    : __ldg(currentBlock);
                currentGroupFP =
                    letter < parametersFP.wordLookupDFA_numCodes
                    ? &DFA_groupFP_s[currentGroupFP->nextGroups + letter]
                    : &DFA_groupFP_s[currentGroupFP->nextGroups];

                if (currentWord) {
                    subjectOffset = address - subject;
                    // At least one query position, stored at an extenal address
                    queryOffsets = ((uint2 *)currentBlock) - currentWord;

                    if (!__ldg(queryOffsets)) {
                        // Go to an outside address for additional positions
                        queryOffsets = wordLookupDFA_AddPositions +
                            (__ldg(queryOffsets + 1) * constants_max_int2) +
                            __ldg(queryOffsets + 2);
                    }

                    do {
#ifndef NO_STAGE2
                        queryOffset = __ldg(queryOffsets);
                        int2 diagonal = subjectOffset - queryOffset;
                        unsigned char bin_id = (uint2)diagonal & BIN_MARK;
                        unsigned int nHits = atomicAdd(&(bin_numHits[bin_id]), 1);
#endif
                        queryOffsets++;
                    } while (__ldg(queryOffsets));
                }
                start += 32;
                address = start + laneId;
            }
        }
        sequenceCount += gridDim.x * 4;
    }

    for (int ii = 0; ii < BIN_X; ii++)
        blast_numHits[laneId + ii * 32] = bin_numHits[laneId + ii * 32];
}

__global__ void search_protein_hit_detection(
        struct sequenceDataFP *sequenceDataFP, unsigned char *sequence,
        struct parameters *parametersFP_g, struct groupFP *wordLookupDFA_groupsFP,
        const unsigned char *__restrict wordLookupDFAFP, int *blast_numHits_g,
        uint4 nTotalSequenceNum, uint64_t *HitInfo_g, unsigned int groupNum,
        int *binOffset_g) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int tt = threadIdx.x;
    unsigned int warpId = tid >> 5;
    unsigned char laneId = tt & 31;

    extern __shared__ unsigned char DFA_groupFP_array[];

    __shared__ unsigned int hitsDia_share[BLOCK_SIZE * BIN_X];
    __shared__ int binOffset_share[BLOCK_SIZE * BIN_X];

    unsigned int *hitsDia = hitsDia_share + (tt >> 5 << BIN_POWER);
    int *binOffset_t = binOffset_share + (tt >> 5 << BIN_POWER);

    // uint2 hitslast = 0;
    for (int ii = 0; ii < BIN_X; ii++) {
        int bin_id = (warpId << BIN_POWER) + laneId + ii * 32;
        binOffset_t[laneId + ii * 32] =
            binOffset_g[bin_id] - blast_numHits_g[bin_id];
        hitsDia[laneId + ii * 32] = 0;
    }
    unsigned char *subject, *sequenceEnd, *address, *start;
    int4 subjectOffset, count;
    unsigned char currentWord;
    const unsigned char *currentBlock;
    uint2 *queryOffsets, queryOffset;

    struct parameters parametersFP = *parametersFP_g;
    extern __shared__ struct groupFP DFA_groupFP_s[];

    for (unsigned int ii = tt; ii < groupNum; ii += blockDim.x) {
        DFA_groupFP_s[ii] = wordLookupDFA_groupsFP[ii];
    }
    __syncthreads();

    uint2 *wordLookupDFA_AddiPositions =
        (uint2 *)((char *)wordLookupDFAFP +
                parametersFP.additionalQueryPositionOffset);

    uint64_t sequenceCount = warpId;

    while (sequenceCount < nTotalSequenceNum) {
        if (sequenceDataFP[sequenceCount].sequenceLength >=
                parametersFP.parameters_wordSize) {
            start = subject = sequence + sequenceDataFP[sequenceCount].offset;
            sequenceEnd = subject + sequenceDataFP[sequenceCount].sequenceLength;
            address = start + laneId;

            while (address + parametersFP.parameters_wordSize - 1 < sequenceEnd) {
                struct groupFP *currentGroupFP = DFA_groupFP_s;
                unsigned char letter = *address;
                for (count = 1; count < parametersFP.parameters_wordSize; count++) {
                    currentGroupFP =
                        letter < parametersFP.wordLookupDFA_numCodes
                        ? &DFA_groupFP_s[currentGroupFP->nextGroups + letter]
                        : &DFA_groupFP_s[currentGroupFP->nextGroups];
                    address++;
                    letter = *address;
                }

                currentBlock = &wordLookupDFAFP[currentGroupFP->nextWords];

                // If current code is a regular letter
                currentWord = letter < parametersFP.wordLookupDFA_numCodes
                    ? __ldg(currentBlock + letter)
                    : __ldg(currentBlock);
                currentGroupFP =
                    letter < parametersFP.wordLookupDFA_numCodes
                    ? &DFA_groupFP_s[currentGroupFP->nextGroups + letter]
                    : &DFA_groupFP_s[currentGroupFP->nextGroups];

                if (currentWord) {
                    subjectOffset = address - subject;
                    // At least one query position, stored at an extenal address
                    queryOffsets = ((uint2 *)currentBlock) - currentWord;

                    if (!__ldg(queryOffsets)) {
                        // Go to an outside address for additional positions
                        queryOffsets = wordLookupDFA_AddiPositions +
                            (__ldg(queryOffsets + 1) * constants_max_int2) +
                            __ldg(queryOffsets + 2);
                    }

                    do {
#ifndef NO_STAGE2
                        queryOffset = __ldg(queryOffsets);
                        int2 diagonal = subjectOffset - queryOffset;
                        unsigned char bin_id = (uint2)diagonal & BIN_MARK;
                        unsigned int nHits = atomicAdd(&(hitsDia[bin_id]), 1);
                        unsigned int bin_p = binOffset_t[bin_id] + nHits;
                        HitInfo_g[bin_p] = (sequenceCount << 32) +
                            ((diagonal + 0x3fff) << 16) + subjectOffset;
                        queryOffsets++;
#endif
                    } while (__ldg(queryOffsets));
                }
                start += 32;
                address = start + laneId;
            }
        }
        sequenceCount += gridDim.x * 4;
    }
}

__global__ void search_protein_hit_detection_warp(
        struct sequenceDataFP *sequenceDataFP, unsigned char *sequence,
        struct parameters *parametersFP_g, struct groupFP *wordLookupDFA_groupsFP,
        const unsigned char *__restrict wordLookupDFAFP, int *blast_numHits_g,
        uint4 nTotalSequenceNum, uint64_t *HitInfo_g, unsigned int groupNum,
        unsigned int num_hits) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int tt = threadIdx.x;
    unsigned int warpId = tid >> 5;
    unsigned char laneId = tt & 31;

    extern __shared__ unsigned char DFA_groupFP_array[];

    __shared__ unsigned int hitsDia_share[BLOCK_SIZE * BIN_X];
    __shared__ unsigned int hitsLast_share[BLOCK_SIZE * BIN_X];

    unsigned int *hitsDia = hitsDia_share + (tt >> 5 << BIN_POWER);
    unsigned int *hitsLast = hitsLast_share + (tt >> 5 << BIN_POWER);

    uint64_t *HitInfo_t = HitInfo_g + num_hits * (warpId << BIN_POWER);
    // uint64_t *HitInfo_b = HitInfo_g + num_hits * (warpId << 7);

    // int *hitsoffset_warp = hitsoffset_g + NUM_SEQS * (warpId << 5);
    // int *hitsoffset_b = hitsoffset_g + NUM_SEQS * (warpId << 7);
    // int *hitsoffset = hitsoffset_warp + laneId * NUM_SEQS;
    // hitsoffset[0] = 0;
    // hitsoffset++;

    // uint2 hitslast = 0;
    for (int ii = 0; ii < BIN_X; ii++)
        hitsDia[laneId + ii * 32] = 0;

    int *blast_numHits = blast_numHits_g + (warpId << BIN_POWER);
    // int *blast_numSeqs = blast_numSeqs_g + (warpId << 5);
    unsigned char *subject, *sequenceEnd, *address, *start;
    int4 subjectOffset, count;
    unsigned char currentWord;
    const unsigned char *currentBlock;
    uint2 *queryOffsets, queryOffset;

    struct parameters parametersFP = *parametersFP_g;
    extern __shared__ struct groupFP DFA_groupFP_s[];

    for (unsigned int ii = tt; ii < groupNum; ii += blockDim.x) {
        DFA_groupFP_s[ii] = wordLookupDFA_groupsFP[ii];
    }
    __syncthreads();

    uint2 *wordLookupDFA_AddiPositions =
        (uint2 *)((char *)wordLookupDFAFP +
                parametersFP.additionalQueryPositionOffset);

    uint64_t sequenceCount = warpId;

    while (sequenceCount < nTotalSequenceNum) {
        if (sequenceDataFP[sequenceCount].sequenceLength >=
                parametersFP.parameters_wordSize) {
            start = subject = sequence + sequenceDataFP[sequenceCount].offset;
            sequenceEnd = subject + sequenceDataFP[sequenceCount].sequenceLength;
            address = start + laneId;

            while (address + parametersFP.parameters_wordSize - 1 < sequenceEnd) {
                struct groupFP *currentGroupFP = DFA_groupFP_s;
                unsigned char letter = *address;
                for (count = 1; count < parametersFP.parameters_wordSize; count++) {
                    currentGroupFP =
                        letter < parametersFP.wordLookupDFA_numCodes
                        ? &DFA_groupFP_s[currentGroupFP->nextGroups + letter]
                        : &DFA_groupFP_s[currentGroupFP->nextGroups];
                    address++;
                    letter = *address;
                }

                currentBlock = &wordLookupDFAFP[currentGroupFP->nextWords];

                // If current code is a regular letter
                currentWord = letter < parametersFP.wordLookupDFA_numCodes
                    ? __ldg(currentBlock + letter)
                    : __ldg(currentBlock);
                currentGroupFP =
                    letter < parametersFP.wordLookupDFA_numCodes
                    ? &DFA_groupFP_s[currentGroupFP->nextGroups + letter]
                    : &DFA_groupFP_s[currentGroupFP->nextGroups];

                if (currentWord) {
                    subjectOffset = address - subject;
                    // At least one query position, stored at an extenal address
                    queryOffsets = ((uint2 *)currentBlock) - currentWord;

                    if (!__ldg(queryOffsets)) {
                        // Go to an outside address for additional positions
                        queryOffsets = wordLookupDFA_AddiPositions +
                            (__ldg(queryOffsets + 1) * constants_max_int2) +
                            __ldg(queryOffsets + 2);
                    }

                    do {
#ifndef NO_STAGE2
                        queryOffset = __ldg(queryOffsets);
                        int2 diagonal = subjectOffset - queryOffset;
                        unsigned char bin_id = (uint2)diagonal & BIN_MARK;
                        unsigned int nHits = atomicAdd(&(hitsDia[bin_id]), 1);
                        unsigned int bin_p = bin_id * num_hits + nHits;
                        HitInfo_t[bin_p] = (sequenceCount << 32) +
                            ((diagonal + 0x3fff) << 16) + subjectOffset;
                        // uint32_t *hit2 = (uint32_t *)(&HitInfo_t[bin_p]);
                        // uint16_t *hit4 = (uint16_t *)(&HitInfo_t[bin_p]);
                        // hit2[1] = sequenceCount;
                        // hit4[1] = diagonal + 0x3fff;
                        // hit4[0] = subjectOffset;
                        queryOffsets++;
#endif
                    } while (__ldg(queryOffsets));
                }
                start += 32;
                address = start + laneId;
            }
        }

        sequenceCount += gridDim.x * 4;

        for (int ii = 0; ii < BIN_X; ii++) {
            if ((hitsDia[laneId + ii * 32] - hitsLast[laneId + ii * 32]) < 2) {
                hitsDia[laneId + ii * 32] = hitsLast[laneId + ii * 32];
            }
            hitsLast[laneId + ii * 32] = hitsDia[laneId + ii * 32];
        }
    }

    for (int ii = 0; ii < BIN_X; ii++)
        blast_numHits[laneId + ii * 32] = hitsDia[laneId + ii * 32];
}

void search_protein2hitParallel(struct scoreMatrix *scoreMatrixp,
        struct PSSMatrix PSSMatrix,
        struct PSSMatrixFP PSSMatrixFP,
        struct sequenceData *sequenceData_host,
        uint4 numSequences, uint4 tickFrequency) {
    // Shucai
    uint4 i, j, sequenceCount = 0;
    uint4 nRoundOffset;

    // PSSMatrix pointers
    struct PSSMatrixFP *PSSMatrixFPD;
    int2 *matrixBodyD;

    // Input database sequence
    struct sequenceDataFP *sequenceDataFP;
    struct sequenceDataFP *sequenceDataFPD;
    unsigned char *sequencesD;
    // unsigned char *sequencesH;
    unsigned char *roundStartAddress;

    // ungapped extension
    struct ungappedExtension *ungappedExtensionsD;
    struct ungappedExtension *ungappedExtension;
    struct ungappedExtension *ungappedExtensionCur, *newUngappedExtension,
                             *additionalUngappedExtension;

    struct sequenceData *sequenceData;
    cudaMallocHost((void **)&sequenceData,
            sizeof(struct sequenceData) * numSequences);
    memcpy(sequenceData, sequenceData_host,
            sizeof(struct sequenceData) * numSequences);
    // ungapped extension numbers
    uint4 *blast_numUngappedExtensionsD, *blast_numUngappedExtensionsH;
    uint4 *blast_numTriggerExtensionsD, *blast_numTriggerExtensionsH;
    uint4 numAdditionalTriggerExtensions, numExtensions;
    int *blast_numHitsD, *blast_numHitsH;
    int *blast_numExtD;
    int *binOffsetD;
    int4 preSequenceCount;

    // For time record
    struct timeval t0, t1, t2, t3, t4, t5, t6, t7, t8, t9;

    int4 wordNum, groupNum;

    // parameters
    struct parameters strParameters;
    struct parameters *parametersD;

    // word lookup table
    struct groupFP *wordLookupDFA_groupD;
    unsigned char *wordLookupDFAD;
    uint4 wordLookupDFA_size;

    // grid and block dimensions
    int nBlockNum = NUM_BLOCK;
    int nBlockSize = BLOCK_SIZE;
    int nTotalThreadNum = nBlockNum * nBlockSize;

    // get t0
    gettimeofday(&t0, NULL);

    wordNum = wordLookupDFA_numWords;
    groupNum = wordLookupDFA_numGroups;

    // printf("\n");
    // Allocate GPU buffer for PSSMatrix
    cudaMalloc((void **)&PSSMatrixFPD, sizeof(struct PSSMatrixFP));
    cudaMalloc((void **)&matrixBodyD,
            sizeof(int2) * (PSSMatrixFP.length + 2) * encoding_numCodes);
    // printf("matrixBody: %d\n", (PSSMatrixFP.length + 2) * encoding_numCodes);
    char *matrixBodyH = (char *)malloc(sizeof(char) * (PSSMatrixFP.length + 2) *
            encoding_numCodes);

    for (int ii = 0; ii < PSSMatrixFP.length + 2; ii++) {
        for (int jj = 0; jj < encoding_numCodes; jj++) {
            // printf(".%d.", (PSSMatrixFP.matrix - encoding_numCodes)[ii *
            // encoding_numCodes + jj]);
            matrixBodyH[ii * encoding_numCodes + jj] =
                (PSSMatrixFP.matrix - encoding_numCodes)[ii * encoding_numCodes + jj];
            if (ii == 0 || ii == PSSMatrixFP.length + 1 ||
                    jj == encoding_numCodes - 1)
                matrixBodyH[ii * encoding_numCodes + jj] = -127;
            // printf(".%d.", matrixBodyH[ii * encoding_numCodes + jj]);
        }

        // printf("\n");
    }

    // Copy PSSMatrix to device memory
    cudaMemcpy(PSSMatrixFPD, &PSSMatrixFP, sizeof(struct PSSMatrixFP),
            cudaMemcpyHostToDevice);
    cudaMemcpy(matrixBodyD, matrixBodyH,
            sizeof(int2) * (PSSMatrixFP.length + 2) * encoding_numCodes,
            cudaMemcpyHostToDevice);
    free(matrixBodyH);
    // cudaMemcpyToSymbol(matrixBody_c, (PSSMatrixFP.matrix - encoding_numCodes),
    // sizeof(int2) * (PSSMatrixFP.length + 2) * encoding_numCodes);

    // Each thread is for align of one database sequence
    // sequenceDataFP = (struct sequenceDataFP *)global_malloc(numSequences *
    // sizeof(struct sequenceDataFP));
    cudaMallocHost((void **)&sequenceDataFP,
            sizeof(struct sequenceDataFP) * numSequences);
    cudaMalloc((void **)&sequenceDataFPD,
            numSequences * sizeof(struct sequenceDataFP));

    // Allocate buffer for hit matrix offset
    // hitMatrix_offsetH = (uint4 *)global_malloc((nBlockNum + 1) *
    // sizeof(uint4));
    // cudaMalloc((void **)&hitMatrix_offsetD, (nBlockNum + 1) * sizeof(uint4));

    // Allocate ungapped extension buffer on device
    // int4 nUngappedExtensionNum = UNGAPEXT_PER_THREAD * nTotalThreadNum;
    int4 nUngappedExtensionNum = TOTAL_UNGAPPED_EXT;
    strParameters.ungappedExtensionsPerThread =
        nUngappedExtensionNum / nTotalThreadNum - 1;
    strParameters.ungappedExtAdditionalStartLoc =
        strParameters.ungappedExtensionsPerThread * nTotalThreadNum;

    cudaMalloc((void **)&ungappedExtensionsD,
            nUngappedExtensionNum * sizeof(struct ungappedExtension));
    ungappedExtension = (struct ungappedExtension *)global_malloc(
            nUngappedExtensionNum * sizeof(struct ungappedExtension));

    // Allocate numbers for ungapped extensions
    blast_numUngappedExtensionsH =
        (uint4 *)global_malloc(sizeof(uint4) * nTotalThreadNum);
    blast_numTriggerExtensionsH =
        (uint4 *)global_malloc(sizeof(uint4) * nTotalThreadNum);
    blast_numHitsH = (int *)global_malloc(sizeof(int) * nTotalThreadNum * BIN_X);

    cudaMalloc((void **)&blast_numUngappedExtensionsD,
            sizeof(uint4) * nTotalThreadNum);
    cudaMalloc((void **)&blast_numTriggerExtensionsD,
            sizeof(uint4) * nTotalThreadNum);
    cudaMalloc((void **)&blast_numHitsD, sizeof(int) * nTotalThreadNum * BIN_X);
    cudaMalloc((void **)&blast_numExtD, sizeof(int) * nTotalThreadNum * BIN_X);
    cudaMalloc((void **)&binOffsetD, sizeof(int) * nTotalThreadNum * BIN_X);

    // Allocate device memory, about 132Mbytes (according to texture limit)
    cudaMalloc((void **)&sequencesD,
            sizeof(unsigned char) * (parameters_batchSize + 50000));

    // Allocate parameters buffer on device
    cudaMalloc((void **)&parametersD, sizeof(struct parameters));
    strParameters.parameters_wordSize = parameters_wordSize;
    strParameters.encoding_numCodes = encoding_numCodes;
    strParameters.wordLookupDFA_numCodes = wordLookupDFA_numCodes;
    strParameters.additionalQueryPositionOffset =
        wordNum * sizeof(char) + sizeof(int2) * wordLookupDFA_numExtPositions;
    strParameters.blast_ungappedNominalTrigger = blast_ungappedNominalTrigger;
    strParameters.statistics_ungappedNominalDropoff =
        statistics_ungappedNominalDropoff;
    strParameters.parameters_A = parameters_A;
    strParameters.parameters_overlap = parameters_overlap;
    // printf("parameters_A: %d parameters_overlap: %d\n", parameters_A,
    // parameters_overlap);
    cudaMemcpy(parametersD, &strParameters, sizeof(struct parameters),
            cudaMemcpyHostToDevice);
    // printf("parameters_Size: %d\n", sizeof(struct parameters));

    // Allocate word lookup table
    wordLookupDFA_size =
        sizeof(char) * wordNum + 2 * sizeof(int2) * wordLookupDFA_numExtPositions;
    cudaMalloc((void **)&wordLookupDFA_groupD, sizeof(struct groupFP) * groupNum);
    cudaMalloc((void **)&wordLookupDFAD, wordLookupDFA_size);
    cudaMemcpy(wordLookupDFAD, wordLookupDFA, wordLookupDFA_size,
            cudaMemcpyHostToDevice);
    cudaMemcpy(wordLookupDFA_groupD, wordLookupDFA_groupsFP,
            sizeof(struct groupFP) * groupNum, cudaMemcpyHostToDevice);
    // printf("numDFAGroup: %d DFA_group_size: %u wordLookupDFA_size: %u\n",
    // groupNum, sizeof(struct groupFP) * groupNum, wordLookupDFA_size);

    cudaMemset(blast_numUngappedExtensionsD, 0, sizeof(uint4) * nTotalThreadNum);
    cudaMemset(blast_numHitsD, 0, sizeof(uint4) * nBlockNum);

    //	//Use constant memory for the word lookup table group
    //	cudaMemcpyToSymbol(wordLookupDFA_groupsC, wordLookupDFA_groupsFP,
    //sizeof(struct groupFP) * groupNum);
    //
    // Use constant memory to store score matrix
    int scoreMatrixSize = encoding_numCodes * encoding_numCodes;
    cudaMemcpyToSymbol(scoreMatrixC, ((char *)scoreMatrixp->matrix) +
            sizeof(int2 *) * encoding_numCodes,
            sizeof(int2) * scoreMatrixSize);

    // Use constant memory to store query sequence
    unsigned char *tempQueryCode;
    tempQueryCode = (unsigned char *)global_malloc(sizeof(unsigned char) *
            (PSSMatrixFP.length + 2));
    memcpy(&tempQueryCode[1], PSSMatrixFP.queryCodes,
            sizeof(unsigned char) * PSSMatrixFP.length);
    tempQueryCode[0] = encoding_sentinalCode;
    tempQueryCode[PSSMatrixFP.length + 1] = encoding_sentinalCode;
    cudaMemcpyToSymbol(querySequenceC, tempQueryCode,
            sizeof(unsigned char) * (PSSMatrixFP.length + 2));
    free(tempQueryCode);

    // uint4 iniVal = nTotalThreadNum;
    // printf("PSSMatrixSize: %d scoreMatrixSize: %d querySize: %d DFA_group: %d
    // DFA_qp: %d\n", sizeof(char) * (PSSMatrixFP.length + 2) * encoding_numCodes
    // >> 10, sizeof(int2) * scoreMatrixSize >> 10, (PSSMatrixFP.length + 2) *
    // sizeof(unsigned char) >> 10, (sizeof(struct groupFP) * groupNum) >> 10,
    // wordLookupDFA_size >> 10);
    // get t1
    gettimeofday(&t1, NULL);
    int4 numSequencesRound, numSequenceProcessed;
    numSequenceProcessed = 0;
    // int totalSeqLength = 0;
    size_t dmem_tot = 0, dmem_free = 0;
    cudaMemGetInfo(&dmem_free, &dmem_tot);
    printf("Dmem total: %d (MB) Batch size: %d (MB) Dmem free: %d (MB)\n", dmem_tot >> 20, parameters_batchSize >> 20, dmem_free >> 20);

    // int *HitNumD;
    while (numSequenceProcessed < numSequences) {

        // get t2
        gettimeofday(&t2, NULL);

        cudaMemset(blast_numTriggerExtensionsD, 0, sizeof(uint4) * nTotalThreadNum);

        roundStartAddress = sequenceData[sequenceCount].sequence - 1;
        for (i = 0; sequenceCount < numSequences; i++, sequenceCount++) {
            sequenceDataFP[i].descriptionLength =
                sequenceData[sequenceCount].descriptionLength;
            sequenceDataFP[i].descriptionStart =
                sequenceData[sequenceCount].descriptionStart;
            sequenceDataFP[i].sequenceLength =
                sequenceData[sequenceCount].sequenceLength;
            sequenceDataFP[i].encodedLength =
                sequenceData[sequenceCount].encodedLength;
            sequenceDataFP[i].offset =
                sequenceData[sequenceCount].sequence - roundStartAddress;

            // about 130MB
            if (sequenceDataFP[i].offset + sequenceData[sequenceCount].encodedLength >
                    parameters_batchSize) {
                i++;
                sequenceCount++;
                break;
            }
        }

        nRoundOffset =
            sequenceDataFP[i - 1].offset + sequenceDataFP[i - 1].encodedLength;
        numSequencesRound = i;
        cudaMemcpy(sequencesD, roundStartAddress,
                sizeof(unsigned char) * (nRoundOffset + 2),
                cudaMemcpyHostToDevice);
        cudaMemcpy(sequenceDataFPD, sequenceDataFP,
                sizeof(struct sequenceDataFP) * numSequencesRound,
                cudaMemcpyHostToDevice);

        gettimeofday(&t3, NULL);
        // cudaMemcpyToSymbol(global_sequenceCount, &iniVal, sizeof(uint4));
        numAdditionalTriggerExtensions = 0;
        cudaMemcpyToSymbol(global_numAdditionalTriggerExtensions,
                &numAdditionalTriggerExtensions, sizeof(uint4));

        // get t4
        gettimeofday(&t4, NULL);

        struct timeval s0, s1;
        cudaError_t cudaRes;

        dim3 dimGrid(nBlockNum, 1);
        dim3 dimBlock(nBlockSize, 1);

        int kernel_time = 0;
        gettimeofday(&s0, NULL);
        search_protein_hit_detection_prescan
            <<<dimGrid, dimBlock, sizeof(struct groupFP) * groupNum>>>
            (sequenceDataFPD, sequencesD, parametersD, wordLookupDFA_groupD,
             wordLookupDFAD, blast_numHitsD, numSequencesRound, groupNum);
        cudaThreadSynchronize();
        cudaRes = cudaGetLastError();
        if (cudaRes != cudaSuccess) {
            printf("CUDA error: %s in %d\n", cudaGetErrorString(cudaRes), __LINE__);
            exit(-1);
        }

        int total_numHits =
            get_total_numHits(blast_numHitsD, nTotalThreadNum * BIN_X);
        get_bin_offset(blast_numHitsD, binOffsetD, nTotalThreadNum * BIN_X);


        cudaMemGetInfo(&dmem_free, &dmem_tot);

        size_t est_mem_usage = total_numHits * sizeof(uint64_t) * 2 +
            (parameters_batchSize + 50000) * sizeof(unsigned char);
        printf("Est mem usage: %d (MB) numSeqProc: %d percentProc: %d\n", est_mem_usage >> 20, numSequenceProcessed, 100 * numSequenceProcessed/numSequences);

        // printf("total_numHits: %d hit buffer memory size: %d (MB)\n",
        // total_numHits, sizeof(uint64_t) * total_numHits >> 20);
        // printf("Estimated memory usage: %d (MB)\n", (sizeof(uint64_t) *
        // total_numHits * 2 + BATCH_SIZE + 2000000) >> 20);

        uint64_t *HitInfoD;
        cudaMalloc((void **)&HitInfoD, sizeof(uint64_t) * total_numHits);

        gettimeofday(&s1, NULL);
        kernel_time +=
            (1000000 * (s1.tv_sec - s0.tv_sec) + s1.tv_usec - s0.tv_usec);
#ifdef VERBOSE
        printf(
                "Hit Detection Prescan Time: %f\n",
                (float)(1000000 * (s1.tv_sec - s0.tv_sec) + s1.tv_usec - s0.tv_usec) /
                1000000);
#endif

        gettimeofday(&s0, NULL);
        search_protein_hit_detection
            <<<dimGrid, dimBlock, sizeof(struct groupFP) * groupNum>>>
            (sequenceDataFPD, sequencesD, parametersD, wordLookupDFA_groupD,
             wordLookupDFAD, blast_numHitsD, numSequencesRound, HitInfoD, groupNum,
             binOffsetD);
        cudaThreadSynchronize();
        cudaRes = cudaGetLastError();
        if (cudaRes != cudaSuccess) {
            printf("CUDA error: %s in %d\n", cudaGetErrorString(cudaRes), __LINE__);
            exit(-1);
        }

        gettimeofday(&s1, NULL);
        kernel_time +=
            (1000000 * (s1.tv_sec - s0.tv_sec) + s1.tv_usec - s0.tv_usec);

#ifdef VERBOSE
        printf(
                "Hit Detection Time: %f\n",
                (float)(1000000 * (s1.tv_sec - s0.tv_sec) + s1.tv_usec - s0.tv_usec) /
                1000000);
#endif

        hit_sort_filter(blast_numHitsD, blast_numExtD, binOffsetD, HitInfoD,
                BLOCK_SIZE, NUM_BLOCK, BIN_X, total_numHits,
                strParameters.parameters_A,
                strParameters.parameters_overlap);

        gettimeofday(&s1, NULL);

        kernel_time +=
            (1000000 * (s1.tv_sec - s0.tv_sec) + s1.tv_usec - s0.tv_usec);
        // printf("Hit Sort Time: %d\n", (1000000 * (s1.tv_sec - s0.tv_sec) +
        // s1.tv_usec - s0.tv_usec));

        gettimeofday(&s0, NULL);
        ungappedExtension_twoHitExtendG_bin_sorted_sm_s
            <<<dimGrid, dimBlock,
            sizeof(unsigned char) * (PSSMatrixFP.length + 2)>>>
                (parametersD, HitInfoD, blast_numUngappedExtensionsD, sequencesD,
                 PSSMatrixFPD, sequenceDataFPD, ungappedExtensionsD,
                 blast_numTriggerExtensionsD, blast_numHitsD, blast_numExtD,
                 binOffsetD);
        cudaThreadSynchronize();
        gettimeofday(&s1, NULL);
        cudaRes = cudaGetLastError();
        if (cudaRes != cudaSuccess) {
            printf("CUDA error: %s in %d\n", cudaGetErrorString(cudaRes), __LINE__);
            exit(-1);
        }

        kernel_time +=
            (1000000 * (s1.tv_sec - s0.tv_sec) + s1.tv_usec - s0.tv_usec);

#ifdef VERBOSE
        printf(
                "Hit Extension Time: %f\n",
                (float)(1000000 * (s1.tv_sec - s0.tv_sec) + s1.tv_usec - s0.tv_usec) /
                1000000);
#endif

        // get t5
        gettimeofday(&t5, NULL);

        // Post processing
        // copy hit results back
        cudaMemcpy(blast_numTriggerExtensionsH, blast_numTriggerExtensionsD,
                sizeof(uint4) * nTotalThreadNum, cudaMemcpyDeviceToHost);

        cudaMemcpy(ungappedExtension, ungappedExtensionsD,
                nUngappedExtensionNum * sizeof(struct ungappedExtension),
                cudaMemcpyDeviceToHost);

        cudaMemcpyFromSymbol(&numAdditionalTriggerExtensions,
                global_numAdditionalTriggerExtensions, sizeof(uint4));

        // get t6
        gettimeofday(&t6, NULL);

        // Add hits to the alignment list
        // Additional buffer is used, sort ungapped extensions
        // according to sequence index

        long qsort_time = 0;
        if (numAdditionalTriggerExtensions > 0) {
            additionalUngappedExtension =
                ungappedExtension + strParameters.ungappedExtAdditionalStartLoc;
            gettimeofday(&s0, NULL);
            qsort(additionalUngappedExtension, numAdditionalTriggerExtensions,
                    sizeof(struct ungappedExtension), compare_ungappedextension);
            gettimeofday(&s1, NULL);
            qsort_time +=
                (1000000 * (s1.tv_sec - s0.tv_sec) + s1.tv_usec - s0.tv_usec);
        }

        for (i = 0; i < nTotalThreadNum; i++) {
            // printf("%d %d\n", i, blast_numTriggerExtensionsH[i]);
            if (blast_numTriggerExtensionsH[i] > 0) {
                // ungappedExtensionCur =  ungappedExtension + i * UNGAPEXT_PER_THREAD;
                ungappedExtensionCur =
                    ungappedExtension + i * strParameters.ungappedExtensionsPerThread;
                preSequenceCount = -1;
                numExtensions = (blast_numTriggerExtensionsH[i] >
                        strParameters.ungappedExtensionsPerThread)
                    ? strParameters.ungappedExtensionsPerThread
                    : blast_numTriggerExtensionsH[i];

                gettimeofday(&s0, NULL);
                qsort(ungappedExtensionCur, numExtensions,
                        sizeof(struct ungappedExtension), compare_ungappedextension);
                gettimeofday(&s1, NULL);
                qsort_time +=
                    (1000000 * (s1.tv_sec - s0.tv_sec) + s1.tv_usec - s0.tv_usec);

                for (j = 0; j < numExtensions; j++) {

                    // printf("seq: %d sub: %d - %d qury: %d - %d seed: %d - %d\n",
                    ////i,
                    // ungappedExtensionCur[j].sequenceCount + numSequenceProcessed,
                    // ungappedExtensionCur[j].start.subjectOffset,
                    // ungappedExtensionCur[j].end.subjectOffset,
                    // ungappedExtensionCur[j].start.queryOffset,
                    // ungappedExtensionCur[j].end.queryOffset,
                    // ungappedExtensionCur[j].seed.queryOffset,
                    // ungappedExtensionCur[j].seed.subjectOffset
                    //);

                    if (ungappedExtensionCur[j].sequenceCount != preSequenceCount) {

                        alignments_createNew(
                                sequenceData[ungappedExtensionCur[j].sequenceCount +
                                numSequenceProcessed].descriptionStart,
                                sequenceData[ungappedExtensionCur[j].sequenceCount +
                                numSequenceProcessed].descriptionLength,
                                sequenceData[ungappedExtensionCur[j].sequenceCount +
                                numSequenceProcessed].sequence,
                                sequenceData[ungappedExtensionCur[j].sequenceCount +
                                numSequenceProcessed].sequenceLength,
                                sequenceData[ungappedExtensionCur[j].sequenceCount +
                                numSequenceProcessed].encodedLength);
                        preSequenceCount = ungappedExtensionCur[j].sequenceCount;
                    }

                    newUngappedExtension = (struct ungappedExtension *)memBlocks_newEntry(
                            ungappedExtension_extensions);
                    memcpy(newUngappedExtension, &ungappedExtensionCur[j],
                            sizeof(struct ungappedExtension));
                    alignments_addUngappedExtension(newUngappedExtension);
                }

                // Add additional extensions
                if (blast_numTriggerExtensionsH[i] >
                        strParameters.ungappedExtensionsPerThread) {
                    int tempStartLoc = findStartLoc(additionalUngappedExtension, i,
                            numAdditionalTriggerExtensions);
                    numExtensions = blast_numTriggerExtensionsH[i] -
                        strParameters.ungappedExtensionsPerThread;

                    for (j = tempStartLoc; j < numExtensions + tempStartLoc; j++) {
                        if (additionalUngappedExtension[j].sequenceCount !=
                                preSequenceCount) {
                            alignments_createNew(
                                    sequenceData[additionalUngappedExtension[j].sequenceCount +
                                    numSequenceProcessed].descriptionStart,
                                    sequenceData[additionalUngappedExtension[j].sequenceCount +
                                    numSequenceProcessed].descriptionLength,
                                    sequenceData[additionalUngappedExtension[j].sequenceCount +
                                    numSequenceProcessed].sequence,
                                    sequenceData[additionalUngappedExtension[j].sequenceCount +
                                    numSequenceProcessed].sequenceLength,
                                    sequenceData[additionalUngappedExtension[j].sequenceCount +
                                    numSequenceProcessed].encodedLength);
                            preSequenceCount = additionalUngappedExtension[j].sequenceCount;
                        }

                        newUngappedExtension =
                            (struct ungappedExtension *)memBlocks_newEntry(
                                    ungappedExtension_extensions);
                        memcpy(newUngappedExtension, &additionalUngappedExtension[j],
                                sizeof(struct ungappedExtension));
                        alignments_addUngappedExtension(newUngappedExtension);
                    }
                }

                blast_numTriggerExtensions += blast_numTriggerExtensionsH[i];
            }
        }

#ifdef VERBOSE
        printf("Extension Sort Time: %d\n", (float)qsort_time / 1000000);
#endif
        numSequenceProcessed += numSequencesRound;

        // cudaFree(hitMatrix_furthestD);

        // get t7
        gettimeofday(&t7, NULL);

        // gapped extension for the current chunk of sequences on the GPU
        // alignments_fingGoodAlignmentsGPU(&PSSMatrixFPD,	//GPU buffer
        // PSSMatrixFP,
        // scoreMatrixp,
        //&matrixBodyD,		//GPU buffer
        //&sequenceDataFPD[w], //GPU buffer
        //&sequencesD[w],		//GPU buffer
        // nRoundOffset);
        // use cpu for gapped extension
        // alignments_findGoodAlignments(PSSMatrix, PSSMatrixFP);
        // get t9
        gettimeofday(&t9, NULL);
        timeRecord.gappedAlignmentTime +=
            (1000000 * (t9.tv_sec - t7.tv_sec) + t9.tv_usec - t7.tv_usec);
        // aggregate execution time
        timeRecord.preProcessTime +=
            (1000000 * (t3.tv_sec - t2.tv_sec) + t3.tv_usec - t2.tv_usec);
        timeRecord.dataCopyTimeH2D +=
            (1000000 * (t4.tv_sec - t3.tv_sec) + t4.tv_usec - t3.tv_usec);
        timeRecord.searchTime +=
            (1000000 * (t5.tv_sec - t4.tv_sec) + t5.tv_usec - t4.tv_usec);
        timeRecord.dataCopyTimeD2H +=
            (1000000 * (t6.tv_sec - t5.tv_sec) + t6.tv_usec - t5.tv_usec);
        timeRecord.addUngappedExtensionTime +=
            (1000000 * (t7.tv_sec - t6.tv_sec) + t7.tv_usec - t6.tv_usec);
        blast_numHits += total_numHits;
        cudaFree(HitInfoD);
    }

    cudaMemcpy(blast_numUngappedExtensionsH, blast_numUngappedExtensionsD,
            sizeof(uint4) * nTotalThreadNum, cudaMemcpyDeviceToHost);
    for (j = 0; j < nTotalThreadNum; j++)
        blast_numUngappedExtensions += blast_numUngappedExtensionsH[j];

    cudaFree(PSSMatrixFPD);
    cudaFree(matrixBodyD);
    cudaFree(ungappedExtensionsD);
    cudaFree(blast_numUngappedExtensionsD);
    cudaFree(blast_numTriggerExtensionsD);
    cudaFree(blast_numHitsD);
    cudaFree(blast_numExtD);
    cudaFree(binOffsetD);
    cudaFree(parametersD);
    cudaFree(wordLookupDFA_groupD);
    cudaFree(wordLookupDFAD);
    cudaFree(sequenceDataFPD);
    cudaFree(sequencesD);
    cudaFree(sequenceDataFP);

    free(ungappedExtension);
    free(blast_numUngappedExtensionsH);
    free(blast_numTriggerExtensionsH);
    free(blast_numHitsH);
    cudaFreeHost(sequenceData);

    // get t8
    gettimeofday(&t8, NULL);

    // Record time
    timeRecord.iniTime =
        1000000 * (t1.tv_sec - t0.tv_sec) + t1.tv_usec - t0.tv_usec;
    timeRecord.postProcessTime =
        1000000 * (t8.tv_sec - t7.tv_sec) + t8.tv_usec - t7.tv_usec;
    timeRecord.hitUngappedExtTime =
        1000000 * (t8.tv_sec - t1.tv_sec) + t8.tv_usec - t1.tv_usec;
}
// Embarrassingly parallel approach is used. One thread is used for
// the hit detection of one sequence
__global__ void search_protein1hitKernel(
        struct PSSMatrixFP *PSSMatrixFP, int2 *matrixBody,
        struct sequenceDataFP *sequenceDataFP, unsigned char *sequence,
        struct parameters *parametersFP, struct groupFP *wordLookupDFA_groupsFP,
        unsigned char *wordLookupDFAFP, uint4 *blast_numUngappedExtensions,
        uint4 *blast_numTriggerExtensions, uint4 *blast_numHits,
        uint4 *hitMatrix_furthestp, uint4 *hitMatrix_offsetp,
        struct ungappedExtension *ungappedExtension_extensionsp,
        uint4 nTotalSequenceNum) {
    int bid = blockIdx.x * gridDim.y + blockIdx.y;
    int tid =
        bid * blockDim.x * blockDim.y + threadIdx.x * blockDim.y + threadIdx.y;

    unsigned char *subject, *sequenceEnd, *address;
    int4 subjectOffset, count;
    unsigned char currentWord, *currentBlock;
    struct groupFP *currentGroupFP;
    uint2 *wordLookupDFA_AddiPositions;
    uint4 numOfTriggerExtensions = 0;
    uint2 *queryOffsets, queryOffset;
    struct ungappedExtension *ungappedExtension_current;
    int4 diagonal;
    uint4 *lastHitFP;
    uint4 ungappedExtension_subjectEndReachedFP;
    uint4 *hitMatrix_Local;
    uint4 sequenceCount;

    hitMatrix_Local =
        hitMatrix_furthestp + hitMatrix_offsetp[tid] + PSSMatrixFP->length;
    ungappedExtension_extensionsp->start.subjectOffset = 0;
    ungappedExtension_current =
        ungappedExtension_extensionsp + tid * UNGAPEXT_PER_THREAD;
    wordLookupDFA_AddiPositions =
        (uint2 *)((char *)wordLookupDFAFP +
                parametersFP->additionalQueryPositionOffset);

    // Set the PSSMatrix body
    PSSMatrixFP->matrix = matrixBody + parametersFP->encoding_numCodes;

    sequenceCount = tid;
    while (sequenceCount < nTotalSequenceNum) {
        subject = address = sequence + sequenceDataFP[sequenceCount].offset;

        if (sequenceDataFP[sequenceCount].sequenceLength >=
                parametersFP->parameters_wordSize) {
            currentGroupFP = wordLookupDFA_groupsFP;
            // currentGroupFP = wordLookupDFA_groupsC;

            count = 1;
            while (count < parametersFP->parameters_wordSize) {
                if (*address < parametersFP->wordLookupDFA_numCodes) {
                    currentGroupFP =
                        &wordLookupDFA_groupsFP[currentGroupFP->nextGroups + *address];
                    // currentGroupFP = &wordLookupDFA_groupsC[currentGroupFP->nextGroups
                    // + *address];
                } else {
                    currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups];
                    // currentGroupFP =
                    // &wordLookupDFA_groupsC[currentGroupFP->nextGroups];
                }
                address++;
                count++;
            }

            sequenceEnd = subject + sequenceDataFP[sequenceCount].sequenceLength;

            while (address < sequenceEnd) {
                currentBlock = &wordLookupDFAFP[currentGroupFP->nextWords];

                // If current code is a regular letter
                if (*address < parametersFP->wordLookupDFA_numCodes) {

                    currentWord = currentBlock[*address];
                    currentGroupFP =
                        &wordLookupDFA_groupsFP[currentGroupFP->nextGroups + *address];
                    // currentGroupFP = &wordLookupDFA_groupsC[currentGroupFP->nextGroups
                    // + *address];
                } else {
                    if (address >= sequenceEnd)
                        break;

                    currentWord = *currentBlock;
                    currentGroupFP = &wordLookupDFA_groupsFP[currentGroupFP->nextGroups];
                    // currentGroupFP =
                    // &wordLookupDFA_groupsC[currentGroupFP->nextGroups];
                }

                if (currentWord) {
                    subjectOffset = address - subject;
                    // At least one query position, stored at an extenal address
                    queryOffsets = ((uint2 *)currentBlock) - currentWord;

                    if (!*queryOffsets) {
                        // Go to an outside address for additional positions
                        queryOffsets = wordLookupDFA_AddiPositions +
                            (*(queryOffsets + 1) * constants_max_int2) +
                            *(queryOffsets + 2);
                    }

                    do {
                        queryOffset = *queryOffsets;

#ifndef NO_STAGE2
                        // Calculate the diagonal this hit is on
                        diagonal = subjectOffset - queryOffset;

                        // If we have not extended past this point on this diagonal
                        lastHitFP = hitMatrix_Local + diagonal;

                        if (*lastHitFP < address - sequence) {
                            // Number of extensions for each subject sequence
                            blast_numUngappedExtensions[tid]++;

                            // If only one hit triggered this extension
                            ungappedExtension_oneHitExtendD(
                                    sequence, queryOffset, address, *PSSMatrixFP, subject,
                                    &ungappedExtension_subjectEndReachedFP,
                                    parametersFP->encoding_numCodes,
                                    parametersFP->statistics_ungappedNominalDropoff,
                                    parametersFP->blast_ungappedNominalTrigger,
                                    ungappedExtension_current, &numOfTriggerExtensions,
                                    sequenceCount, tid);

                            // Update furthest reached value for the diagonal
                            *lastHitFP = ungappedExtension_subjectEndReachedFP;
                        }
#endif

                        queryOffsets++;
                        blast_numHits[tid]++;
                    } while (*queryOffsets);
                }
                address++;
            }
        }
        // option=======================================================
        // sequenceCount = atomicAdd(&global_sequenceCount, 1);
        sequenceCount += gridDim.x * blockDim.x;
        //============================================================
    }

    blast_numTriggerExtensions[tid] = (uint4)numOfTriggerExtensions;
    return;
}

__device__ struct ungappedExtension *ungappedExtension_oneHitExtendD(
        unsigned char *sequenceStart, int4 queryOffset, unsigned char *subjectHit,
        struct PSSMatrixFP PSSMatrixFP, unsigned char *subject,
        uint4 *sequenceHitEnd, unsigned char encoding_numCodes,
        int4 statistics_ungappedNominalDropoff, int4 blast_ungappedNominalTrigger,
        struct ungappedExtension *ungappedExtension_extensions,
        uint4 *numOfTriggerExtensions, uint4 sequenceCount, int4 tid) {
    int2 *queryPosition;
    // int4 queryPosition;
    unsigned char *subjectPosition, *subjectStart, *subjectEnd;
    int4 changeSinceBest = 0;
    int4 dropoff, originalDropoff;
    int4 ungappedExtension_bestScore;

    originalDropoff = dropoff = -statistics_ungappedNominalDropoff;
    ungappedExtension_bestScore = 0;

    // Start at queryEnd,subjectEnd (right/last hit position)
    queryPosition = PSSMatrixFP.matrix + queryOffset * encoding_numCodes;
    // queryPosition = queryOffset + 1;
    subjectPosition = subjectStart = subjectHit;

    while (changeSinceBest > dropoff) {
        changeSinceBest += queryPosition[*subjectPosition];
        // changeSinceBest += scoreMatrixC[querySequenceC[queryPosition] *
        // encoding_numCodes + (*subjectPosition)];

        // If we have got a positive score
        if (changeSinceBest > 0) {
            // Keep updating best score and resetting change-since-best
            // whilst we are reading positive scores
            do {
                ungappedExtension_bestScore += changeSinceBest;
                queryPosition = queryPosition - encoding_numCodes;
                // queryPosition = queryPosition - 1;
                subjectPosition--;

                changeSinceBest = queryPosition[*subjectPosition];
                // changeSinceBest = scoreMatrixC[querySequenceC[queryPosition] *
                // encoding_numCodes + (*subjectPosition)];
            } while (changeSinceBest > 0);

            subjectStart = subjectPosition;
        }

        queryPosition = queryPosition - encoding_numCodes;
        // queryPosition = queryPosition - 1;
        subjectPosition--;
    }

    // Correct for extra decrement
    subjectStart++;

    // Starting at right/last hit position again
    queryPosition = PSSMatrixFP.matrix + (queryOffset + 1) * encoding_numCodes;
    // queryPosition = (queryOffset + 2);
    subjectPosition = subjectEnd = subjectHit + 1;
    changeSinceBest = 0;

    // May need to alter dropoff so we also dropoff if below zero
    if (-ungappedExtension_bestScore > originalDropoff) {
        dropoff = -ungappedExtension_bestScore;
    }

    // Extend end of alignment until dropoff
    while (changeSinceBest > dropoff) {
        // Shucai
        changeSinceBest += queryPosition[*subjectPosition];
        // changeSinceBest += scoreMatrixC[querySequenceC[queryPosition] *
        // encoding_numCodes + (*subjectPosition)];

        // If we have got a positive score
        if (changeSinceBest > 0) {
            // Keep updating best score and resetting change-since-best
            // whilst we are reading positive scores
            do {
                ungappedExtension_bestScore += changeSinceBest;
                queryPosition = queryPosition + encoding_numCodes;
                // queryPosition = queryPosition + 1;

                subjectPosition++;
                changeSinceBest = queryPosition[*subjectPosition];
                // changeSinceBest = scoreMatrixC[querySequenceC[queryPosition] *
                // encoding_numCodes + (*subjectPosition)];
            } while (changeSinceBest > 0);

            subjectEnd = subjectPosition;

            // Check need for change in dropoff
            if ((dropoff = -ungappedExtension_bestScore) < originalDropoff) {
                dropoff = originalDropoff;
            }
        }
        queryPosition = queryPosition + encoding_numCodes;
        // queryPosition = queryPosition + 1;

        subjectPosition++;
    }

    subjectEnd--;

    //*sequenceHitEnd = subjectPosition - subject;
    *sequenceHitEnd = subjectPosition - sequenceStart;

    if (ungappedExtension_bestScore >= blast_ungappedNominalTrigger) {
        int4 diagonal;
        struct ungappedExtension *newUngappedExtension = NULL;

        newUngappedExtension =
            &ungappedExtension_extensions[*numOfTriggerExtensions];
        // Calculate diagonal
        diagonal = (subjectHit - subject) - queryOffset;

        // Determine offsets from pointers
        newUngappedExtension->start.subjectOffset = subjectStart - subject;
        newUngappedExtension->end.subjectOffset = subjectEnd - subject;
        newUngappedExtension->start.queryOffset =
            newUngappedExtension->start.subjectOffset - diagonal;
        newUngappedExtension->end.queryOffset =
            newUngappedExtension->end.subjectOffset - diagonal;
        // newUngappedExtension->seed =
        // ungappedExtension_findProteinSeed(newUngappedExtension,
        //        PSSMatrixFP, subject, encoding_numCodes);

        newUngappedExtension->next = NULL;
        newUngappedExtension->nominalScore = ungappedExtension_bestScore;
        newUngappedExtension->status = ungappedExtension_UNGAPPED;
        newUngappedExtension->sequenceCount = sequenceCount;

        // Shucai
        // Record the number of hits satisfying the next step
        (*numOfTriggerExtensions)++;

        return newUngappedExtension;
    } else {
        return NULL;
    }
}

// Shucai
// Search a protein database using 1-hit extension mode
void search_protein1hitParallel(struct scoreMatrix *scoreMatrixp,
        struct PSSMatrixFP PSSMatrixFP,
        struct sequenceData *sequenceData,
        uint4 numSequences, uint4 tickFrequency) {
    // Shucai
    uint4 i, j, sequenceCount = 0;
    uint4 nRoundOffset;

    // PSSMatrix pointers
    struct PSSMatrixFP *PSSMatrixFPD;
    int2 *matrixBodyD;

    // Input database sequence
    struct sequenceDataFP *sequenceDataFP;
    struct sequenceDataFP *sequenceDataFPD;
    unsigned char *sequencesD;
    unsigned char *roundStartAddress;

    // ungapped extension
    struct ungappedExtension *ungappedExtensionsD;
    struct ungappedExtension *ungappedExtension;
    struct ungappedExtension *ungappedExtensionCur, *newUngappedExtension;

    // ungapped extension numbers
    uint4 *blast_numUngappedExtensionsD, *blast_numUngappedExtensionsH;
    uint4 *blast_numTriggerExtensionsD, *blast_numTriggerExtensionsH;
    uint4 *blast_numHitsD, *blast_numHitsH;
    uint4 *hitMatrix_furthestD;
    uint4 *hitMatrix_offsetH;
    uint4 *hitMatrix_offsetD;
    uint4 preSequenceCount;

    // For time record
    struct timeval t0, t1, t2, t3, t4, t5, t6, t7, t8;

    int4 wordNum, groupNum;

    // parameters
    struct parameters strParameters;
    struct parameters *parametersD;

    // word lookup table
    struct groupFP *wordLookupDFA_groupD;
    unsigned char *wordLookupDFAD;
    uint4 wordLookupDFA_size;

    // grid and block dimensions
    int nBlockNum = NUM_BLOCK;
    int nBlockSize = BLOCK_SIZE;
    int nTotalThreadNum = nBlockNum * nBlockSize;
    dim3 dimGrid(nBlockNum, 1);
    dim3 dimBlock(nBlockSize, 1);

    // get t0
    gettimeofday(&t0, NULL);

    wordNum = wordLookupDFA_numWords;
    groupNum = wordLookupDFA_numGroups;

    // Allocate GPU buffer for PSSMatrix
    cudaMalloc((void **)&PSSMatrixFPD, sizeof(struct PSSMatrixFP));
    cudaMalloc((void **)&matrixBodyD,
            sizeof(int2) * (PSSMatrixFP.length + 2) * encoding_numCodes);

    // Copy PSSMatrix to device memory
    cudaMemcpy(PSSMatrixFPD, &PSSMatrixFP, sizeof(struct PSSMatrixFP),
            cudaMemcpyHostToDevice);
    cudaMemcpy(matrixBodyD, (PSSMatrixFP.matrix - encoding_numCodes),
            sizeof(int2) * (PSSMatrixFP.length + 2) * encoding_numCodes,
            cudaMemcpyHostToDevice);

    // Each thread is for align of one database sequence
    sequenceDataFP = (struct sequenceDataFP *)global_malloc(
            numSequences * sizeof(struct sequenceDataFP));
    cudaMalloc((void **)&sequenceDataFPD,
            numSequences * sizeof(struct sequenceDataFP));

    // Allocate buffer for hit matrix offset
    hitMatrix_offsetH =
        (uint4 *)global_malloc((nTotalThreadNum + 1) * sizeof(uint4));
    cudaMalloc((void **)&hitMatrix_offsetD,
            (nTotalThreadNum + 1) * sizeof(uint4));

    // Allocate ungapped extension buffer on device
    int4 nUngappedExtensionNum = UNGAPEXT_PER_THREAD * nTotalThreadNum;

    cudaMalloc((void **)&ungappedExtensionsD,
            nUngappedExtensionNum * sizeof(struct ungappedExtension));
    ungappedExtension = (struct ungappedExtension *)global_malloc(
            nUngappedExtensionNum * sizeof(struct ungappedExtension));

    // Allocate numbers for ungapped extensions
    blast_numUngappedExtensionsH =
        (uint4 *)global_malloc(sizeof(uint4) * nTotalThreadNum);
    blast_numTriggerExtensionsH =
        (uint4 *)global_malloc(sizeof(uint4) * nTotalThreadNum);
    blast_numHitsH = (uint4 *)global_malloc(sizeof(uint4) * nTotalThreadNum);

    cudaMalloc((void **)&blast_numUngappedExtensionsD,
            sizeof(uint4) * nTotalThreadNum);
    cudaMalloc((void **)&blast_numTriggerExtensionsD,
            sizeof(uint4) * nTotalThreadNum);
    cudaMalloc((void **)&blast_numHitsD, sizeof(uint4) * nTotalThreadNum);

    // Allocate device memory, about 132Mbytes (according to texture limit)
    cudaMalloc((void **)&sequencesD, sizeof(unsigned char) * 132000000);

    // Allocate parameters buffer on device
    cudaMalloc((void **)&parametersD, sizeof(struct parameters));
    strParameters.parameters_wordSize = parameters_wordSize;
    strParameters.encoding_numCodes = encoding_numCodes;
    strParameters.wordLookupDFA_numCodes = wordLookupDFA_numCodes;
    strParameters.additionalQueryPositionOffset =
        wordNum * sizeof(char) + sizeof(int2) * wordLookupDFA_numExtPositions;
    strParameters.blast_ungappedNominalTrigger = blast_ungappedNominalTrigger;
    strParameters.statistics_ungappedNominalDropoff =
        statistics_ungappedNominalDropoff;
    cudaMemcpy(parametersD, &strParameters, sizeof(struct parameters),
            cudaMemcpyHostToDevice);

    // Allocate word lookup table
    wordLookupDFA_size =
        sizeof(char) * wordNum + 2 * sizeof(int2) * wordLookupDFA_numExtPositions;
    cudaMalloc((void **)&wordLookupDFA_groupD, sizeof(struct groupFP) * groupNum);
    cudaMalloc((void **)&wordLookupDFAD, wordLookupDFA_size);

    cudaMemset(blast_numUngappedExtensionsD, 0, sizeof(uint4) * nTotalThreadNum);
    cudaMemset(blast_numHitsD, 0, sizeof(uint4) * nTotalThreadNum);

    cudaMemcpy(wordLookupDFA_groupD, wordLookupDFA_groupsFP,
            sizeof(struct groupFP) * groupNum, cudaMemcpyHostToDevice);
    //	//Use constant memory for the word lookup table group
    //	cudaMemcpyToSymbol(wordLookupDFA_groupsC, wordLookupDFA_groupsFP,
    //sizeof(struct groupFP) * groupNum);
    //
    //	//Use constant memory to store score matrix
    //	int scoreMatrixSize = encoding_numCodes * encoding_numCodes;
    //	cudaMemcpyToSymbol(scoreMatrixC,
    //					  ((char *)scoreMatrixp->matrix) + sizeof(int2 *) *
    //encoding_numCodes,
    //					  sizeof(int2) * scoreMatrixSize);

    //	//Use constant memory to store query sequence
    //	unsigned char *tempQueryCode;
    //	tempQueryCode = (unsigned char *)global_malloc(sizeof(unsigned char) *
    //(PSSMatrixFP.length + 2));
    //	memcpy(&tempQueryCode[1], PSSMatrixFP.queryCodes, sizeof(unsigned char)
    //* PSSMatrixFP.length);
    //	tempQueryCode[0] = encoding_sentinalCode;
    //	tempQueryCode[PSSMatrixFP.length + 1] = encoding_sentinalCode;
    //	cudaMemcpyToSymbol(querySequenceC, tempQueryCode, sizeof(unsigned char)
    //* (PSSMatrixFP.length + 2));
    //	free(tempQueryCode);

    cudaMemcpy(wordLookupDFAD, wordLookupDFA, wordLookupDFA_size,
            cudaMemcpyHostToDevice);
    // uint4 iniVal = nTotalThreadNum;

    // get t1
    gettimeofday(&t1, NULL);
    int4 numSequencesRound, numSequenceProcessed;
    numSequenceProcessed = 0;
    while (sequenceCount < numSequences) {
        // get t2
        gettimeofday(&t2, NULL);

        memset(hitMatrix_offsetH, 0, sizeof(int4) * (nTotalThreadNum + 1));
        roundStartAddress = sequenceData[sequenceCount].sequence - 1;
        for (i = 0; sequenceCount < numSequences; i++, sequenceCount++) {
            sequenceDataFP[i].descriptionLength =
                sequenceData[sequenceCount].descriptionLength;
            sequenceDataFP[i].descriptionStart =
                sequenceData[sequenceCount].descriptionStart;
            sequenceDataFP[i].sequenceLength =
                sequenceData[sequenceCount].sequenceLength;
            sequenceDataFP[i].encodedLength =
                sequenceData[sequenceCount].encodedLength;
            sequenceDataFP[i].offset =
                sequenceData[sequenceCount].sequence - roundStartAddress;

            // Calculate the longest sequence size aligned by the current thread
            if (sequenceDataFP[i].sequenceLength >
                    hitMatrix_offsetH[(i % nTotalThreadNum) + 1]) {
                hitMatrix_offsetH[(i % nTotalThreadNum) + 1] =
                    sequenceDataFP[i].sequenceLength;
            }

            // about 130MB
            if (sequenceDataFP[i].offset + sequenceData[sequenceCount].encodedLength >
                    130000000) {
                i++;
                sequenceCount++;
                break;
            }
        }
        nRoundOffset =
            sequenceDataFP[i - 1].offset + sequenceDataFP[i - 1].encodedLength;
        numSequencesRound = i;

        // Calculate the offset of each thread
        for (i = 1; i < nTotalThreadNum + 1; i++) {
            hitMatrix_offsetH[i] += hitMatrix_offsetH[i - 1] +
                (PSSMatrixFP.length - parameters_wordSize + 1);
        }

        // copy offset info to device
        cudaMemcpy(hitMatrix_offsetD, hitMatrix_offsetH,
                (nTotalThreadNum + 1) * sizeof(int4), cudaMemcpyHostToDevice);

        // get t3
        gettimeofday(&t3, NULL);

        // Allocate device memory
        //		cudaMalloc((void **)&sequencesD, sizeof(unsigned char) * (nRoundOffset
        //+ 2));

        // Allocate diagonal buffers
        int nElemNum = hitMatrix_offsetH[nTotalThreadNum];

        cudaMalloc((void **)&hitMatrix_furthestD, sizeof(uint4) * nElemNum);
        cudaMemset(hitMatrix_furthestD, 0, sizeof(uint4) * nElemNum);
        cudaMemset(blast_numTriggerExtensionsD, 0, sizeof(uint4) * nTotalThreadNum);

        // Copy data to device
        cudaMemcpy(sequenceDataFPD, sequenceDataFP,
                sizeof(struct sequenceDataFP) * numSequencesRound,
                cudaMemcpyHostToDevice);
        cudaMemcpy(sequencesD, roundStartAddress,
                sizeof(unsigned char) * (nRoundOffset + 2),
                cudaMemcpyHostToDevice);

        // cudaMemcpyToSymbol(global_sequenceCount, &iniVal, sizeof(uint4));

        // get t4
        gettimeofday(&t4, NULL);

        // all the required data are copied to device, launch the kernel
        search_protein1hitKernel << <dimGrid, dimBlock>>>
            (PSSMatrixFPD, matrixBodyD, sequenceDataFPD, sequencesD, parametersD,
             wordLookupDFA_groupD, wordLookupDFAD, blast_numUngappedExtensionsD,
             blast_numTriggerExtensionsD, blast_numHitsD, hitMatrix_furthestD,
             hitMatrix_offsetD, ungappedExtensionsD, numSequencesRound);
        cudaThreadSynchronize();

        // get t5
        gettimeofday(&t5, NULL);

        // Post processing
        // copy hit results back
        cudaMemcpy(blast_numTriggerExtensionsH, blast_numTriggerExtensionsD,
                sizeof(uint4) * nTotalThreadNum, cudaMemcpyDeviceToHost);

        cudaMemcpy(ungappedExtension, ungappedExtensionsD,
                nUngappedExtensionNum * sizeof(struct ungappedExtension),
                cudaMemcpyDeviceToHost);

        cudaMemcpy(blast_numUngappedExtensionsH, blast_numUngappedExtensionsD,
                sizeof(uint4) * nTotalThreadNum, cudaMemcpyDeviceToHost);

        // get t6
        gettimeofday(&t6, NULL);
        // Add hits to the alignment list
        for (i = 0; i < nTotalThreadNum; i++) {
            if (blast_numTriggerExtensionsH[i] > 0) {
                ungappedExtensionCur = ungappedExtension + i * UNGAPEXT_PER_THREAD;
                preSequenceCount = INT_MAX;
                for (j = 0; j < blast_numTriggerExtensionsH[i]; j++) {
                    if (ungappedExtensionCur[j].sequenceCount != preSequenceCount) {
                        alignments_createNew(
                                sequenceData[ungappedExtensionCur[j].sequenceCount +
                                numSequenceProcessed].descriptionStart,
                                sequenceData[ungappedExtensionCur[j].sequenceCount +
                                numSequenceProcessed].descriptionLength,
                                sequenceData[ungappedExtensionCur[j].sequenceCount +
                                numSequenceProcessed].sequence,
                                sequenceData[ungappedExtensionCur[j].sequenceCount +
                                numSequenceProcessed].sequenceLength,
                                sequenceData[ungappedExtensionCur[j].sequenceCount +
                                numSequenceProcessed].encodedLength);
                        preSequenceCount = ungappedExtensionCur[j].sequenceCount;
                    }

                    newUngappedExtension = (struct ungappedExtension *)memBlocks_newEntry(
                            ungappedExtension_extensions);
                    memcpy(newUngappedExtension, &ungappedExtensionCur[j],
                            sizeof(struct ungappedExtension));
                    alignments_addUngappedExtension(newUngappedExtension);
                }

                blast_numTriggerExtensions += blast_numTriggerExtensionsH[i];
            }
        }

        numSequenceProcessed += numSequencesRound;

        cudaFree(hitMatrix_furthestD);
        //		cudaFree(sequencesD);

        // get t7
        gettimeofday(&t7, NULL);

        // aggregate execution time
        timeRecord.preProcessTime +=
            (1000000 * (t3.tv_sec - t2.tv_sec) + t3.tv_usec - t2.tv_usec);
        timeRecord.dataCopyTimeH2D +=
            (1000000 * (t4.tv_sec - t3.tv_sec) + t4.tv_usec - t3.tv_usec);
        timeRecord.searchTime +=
            (1000000 * (t5.tv_sec - t4.tv_sec) + t5.tv_usec - t4.tv_usec);
        timeRecord.dataCopyTimeD2H +=
            (1000000 * (t6.tv_sec - t5.tv_sec) + t6.tv_usec - t5.tv_usec);
        timeRecord.addUngappedExtensionTime +=
            (1000000 * (t7.tv_sec - t6.tv_sec) + t7.tv_usec - t6.tv_usec);
    }

    // After all sequences are processed
    cudaMemcpy(blast_numUngappedExtensionsH, blast_numUngappedExtensionsD,
            sizeof(uint4) * nTotalThreadNum, cudaMemcpyDeviceToHost);

    cudaMemcpy(blast_numHitsH, blast_numHitsD, sizeof(uint4) * nTotalThreadNum,
            cudaMemcpyDeviceToHost);

    for (j = 0; j < nTotalThreadNum; j++) {
        blast_numUngappedExtensions += blast_numUngappedExtensionsH[j];
        blast_numHits += blast_numHitsH[j];
    }

    cudaFree(PSSMatrixFPD);
    cudaFree(matrixBodyD);
    cudaFree(sequenceDataFPD);
    cudaFree(ungappedExtensionsD);
    cudaFree(blast_numUngappedExtensionsD);
    cudaFree(blast_numTriggerExtensionsD);
    cudaFree(blast_numHitsD);
    cudaFree(parametersD);
    cudaFree(wordLookupDFA_groupD);
    cudaFree(wordLookupDFAD);
    cudaFree(hitMatrix_offsetD);
    cudaFree(sequencesD);

    free(sequenceDataFP);
    free(ungappedExtension);
    free(blast_numUngappedExtensionsH);
    free(blast_numTriggerExtensionsH);
    free(blast_numHitsH);
    free(hitMatrix_offsetH);

    // get t8
    gettimeofday(&t8, NULL);

    // Record time
    timeRecord.iniTime =
        1000000 * (t1.tv_sec - t0.tv_sec) + t1.tv_usec - t0.tv_usec;
    timeRecord.postProcessTime =
        1000000 * (t8.tv_sec - t7.tv_sec) + t8.tv_usec - t7.tv_usec;
    timeRecord.hitUngappedExtTime =
        1000000 * (t8.tv_sec - t1.tv_sec) + t8.tv_usec - t1.tv_usec;
}

int findStartLoc(struct ungappedExtension *ungappedExtensionsPtr, int threadNo,
        int itemNum) {
    int i;

    for (i = 0; i < itemNum; i++) {
        if (ungappedExtensionsPtr[i].tid == threadNo) {
            return i;
        }
    }

    return -1;
}

__device__ struct coordinate ungappedExtension_findProteinSeed_sm(
        struct ungappedExtension *ungappedExtension,
        // struct PSSMatrixFP PSSMatrixFP,
        unsigned char *subject, unsigned char encoding_numCodes,
        // char *matrix
        unsigned char *querySequence, int2 *scoreMatrix) {
    // char *queryWindowStart, *queryWindowEnd;
    int2 queryWindowStart, queryWindowEnd;
    unsigned char *subjectWindowStart, *subjectWindowEnd;

    int2 bestQueryPosition;
    unsigned char *bestSubjectPosition;
    int4 bestSegmentScore;
    int4 nominalScore, count;
    struct coordinate seed;

    if (ungappedExtension->end.queryOffset -
            ungappedExtension->start.queryOffset <
            11) {
        // The seed point is the middle of the extension
        seed.queryOffset = (ungappedExtension->end.queryOffset +
                ungappedExtension->start.queryOffset) /
            2;
        seed.subjectOffset = (ungappedExtension->end.subjectOffset +
                ungappedExtension->start.subjectOffset) /
            2;
    } else {
        // Else find the highest scoring length-11 segment of the ungapped extension
        // queryWindowStart = queryWindowEnd = matrix +
        // ungappedExtension->start.queryOffset * encoding_numCodes;
        queryWindowStart = queryWindowEnd =
            ungappedExtension->start.queryOffset + 1;
        // subjectWindowStart = subjectWindowEnd = subject +
        // ungappedExtension->start.subjectOffset;
        subjectWindowStart = subjectWindowEnd =
            subject + ungappedExtension->start.subjectOffset;

        // Find initial score for first 11 positions
        nominalScore = 0;
        count = 0;
        while (count < 11) {
            // nominalScore += queryWindowEnd[*subjectWindowEnd];
            // queryWindowEnd += encoding_numCodes;
            nominalScore +=
                scoreMatrix[querySequence[queryWindowEnd] * encoding_numCodes +
                (*subjectWindowEnd)];
            queryWindowEnd++;
            subjectWindowEnd++;
            count++;
        }

        // queryWindowEnd -= encoding_numCodes;
        queryWindowEnd--;
        subjectWindowEnd--;

        // By default first-11 positions gives best position and score
        bestQueryPosition = queryWindowStart;
        bestSubjectPosition = subjectWindowStart;
        bestSegmentScore = nominalScore;

        // Now slide the window across and record the better scores/positions
        // while (queryWindowEnd < matrix + ungappedExtension->end.queryOffset *
        // encoding_numCodes)
        while (queryWindowEnd < ungappedExtension->end.queryOffset + 1) {
            // Advance window end, add new position value
            // queryWindowEnd += encoding_numCodes;
            queryWindowEnd++;
            subjectWindowEnd++;

            // nominalScore += queryWindowEnd[*subjectWindowEnd];
            nominalScore +=
                scoreMatrix[querySequence[queryWindowEnd] * encoding_numCodes +
                (*subjectWindowEnd)];

            // Remove position that we will leave behind
            // nominalScore -= queryWindowStart[*subjectWindowStart];
            nominalScore -=
                scoreMatrix[querySequence[queryWindowStart] * encoding_numCodes +
                (*subjectWindowStart)];

            // queryWindowStart += encoding_numCodes;
            queryWindowStart++;
            subjectWindowStart++;

            // Check if best window position yet
            if (nominalScore > bestSegmentScore) {
                bestSegmentScore = nominalScore;
                bestQueryPosition = queryWindowStart;
                bestSubjectPosition = subjectWindowStart;
            }
        }

        // Middle of the best window is the seed position
        seed.queryOffset = (bestQueryPosition - 1) + 5;
        seed.subjectOffset = bestSubjectPosition + 5 - subject;
    }
    return seed;
}

__global__ void ungappedExtension_twoHitExtendG_findProteinSeed(
        struct parameters *parametersFP_g, struct PSSMatrixFP *PSSMatrixFP_g,
        unsigned char *sequence, struct sequenceDataFP *sequenceDataFP,
        struct ungappedExtension *ungappedExtension_extensionsp,
        uint4 *blast_numTriggerExtensions) {

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int tt = threadIdx.x;
    struct parameters *parametersFP = parametersFP_g;
    __shared__ int2 scoreMatrix_s[1024];
    extern __shared__ unsigned char querySequence_s[];
    unsigned encoding_numCodes = parametersFP->encoding_numCodes;
    for (unsigned int ii = tt; ii < encoding_numCodes * encoding_numCodes;
            ii += blockDim.x) {
        scoreMatrix_s[ii] = scoreMatrixC[ii];
    }

    for (unsigned int ii = tt; ii < PSSMatrixFP_g->length + 2; ii += blockDim.x) {
        querySequence_s[ii] = querySequenceC[ii];
    }

    __syncthreads();

    struct ungappedExtension *ungappedExtension_current =
        ungappedExtension_extensionsp +
        tid * parametersFP->ungappedExtensionsPerThread;
    // unsigned int tBins = BLOCK_SIZE * NUM_BLOCK;
    // for( unsigned int bb = tid; bb < tBins; bb += gridDim.x * blockDim.x)
    {
        unsigned int numExtensions = blast_numTriggerExtensions[tid];
        for (unsigned int w_id = 0; w_id < numExtensions; w_id++) {
            int4 sequenceCount = ungappedExtension_current[w_id].sequenceCount;
            unsigned char *subject = sequence + sequenceDataFP[sequenceCount].offset;
            ungappedExtension_current[w_id].seed =
                ungappedExtension_findProteinSeed_sm(
                        &(ungappedExtension_current[w_id]), subject, encoding_numCodes,
                        querySequence_s, scoreMatrix_s);
        }
    }
}

__global__ void ungappedExtension_twoHitExtendG_bin_sorted_sm_s(
        struct parameters *parametersFP, uint64_t *HitInfo_g,
        uint4 *blast_numUngappedExtensions, unsigned char *sequence,
        struct PSSMatrixFP *PSSMatrixFP, struct sequenceDataFP *sequenceDataFP,
        struct ungappedExtension *ungappedExtension_extensionsp,
        uint4 *blast_numTriggerExtensions, int *numOneHitsD, int *numExtD,
        int *binOffset_g) {

    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int tt = threadIdx.x;
    int laneId = tt & 31;
    // int warpId = tid >> 5;
    struct ungappedExtension *ungappedExtension_additional =
        ungappedExtension_extensionsp +
        parametersFP->ungappedExtAdditionalStartLoc;
    unsigned char *ungappedExtension_subjectEndReachedFP = 0;
    __shared__ int2 scoreMatrix_s[1024];
    extern __shared__ unsigned char querySequence_s[];
    unsigned encoding_numCodes = parametersFP->encoding_numCodes;

    for (unsigned int ii = tt; ii < encoding_numCodes * encoding_numCodes;
            ii += blockDim.x) {
        scoreMatrix_s[ii] = scoreMatrixC[ii];
    }

    for (unsigned int ii = tt; ii < PSSMatrixFP->length + 2; ii += blockDim.x) {
        querySequence_s[ii] = querySequenceC[ii];
    }

    __syncthreads();

    // uint4 numOfTriggerExtensions_s = 0;
    __shared__ uint4 numOfTriggerExtensions_s[BLOCK_SIZE];
    uint4 *numOfTriggerExtensions_w = numOfTriggerExtensions_s + (tt >> 5 << 5);
    // uint4 *numOfTriggerExtensions_w = numOfTriggerExtensions_s;
    numOfTriggerExtensions_w[laneId] = 0;
    // numOfTriggerExtensions_w[tt] = 0;

    uint4 blast_numUngappedExtensions_s = 0;
    // ungappedExtension_extensionsp->start.subjectOffset = 0;

    struct ungappedExtension *ungappedExtension_w =
        ungappedExtension_extensionsp +
        (tid >> 5 << 5) * parametersFP->ungappedExtensionsPerThread;
    // struct ungappedExtension *ungappedExtension_w =
    // ungappedExtension_extensionsp +
    // blockIdx.x * BLOCK_SIZE * parametersFP->ungappedExtensionsPerThread;

    unsigned int num_bins = BLOCK_SIZE * NUM_BLOCK * BIN_X;
    for (unsigned int bb = tid; bb < num_bins; bb += gridDim.x * blockDim.x)
        // unsigned int b_start = warpId << BIN_POWER;
        // unsigned int b_end   = (warpId + 1) << BIN_POWER;
        // for(unsigned int bb  =  b_start + laneId; bb < b_end; bb += 32)
    {
        uint64_t *HitInfo_t = HitInfo_g + binOffset_g[bb] - numOneHitsD[bb];
        unsigned int numHits_t = numExtD[bb];
        uint64_t prev_ext = 0;
        for (unsigned int w_id = 0; w_id < numHits_t; w_id += 2) {
            uint64_t prev = HitInfo_t[w_id] > prev_ext ? HitInfo_t[w_id] : prev_ext;
            uint64_t curr = HitInfo_t[w_id + 1];

            if (prev < curr) {

                blast_numUngappedExtensions_s++;
                uint4 sequenceCount = (uint4)(curr >> 32);
                int2 diagonal = (int2)((curr >> 16) & 0xffff) - 0x3fff;
                uint2 subjectOffset = (uint2)(curr & 0xffff);
                uint2 queryOffset = subjectOffset - diagonal;

                unsigned char *subject =
                    sequence + sequenceDataFP[sequenceCount].offset;
                unsigned char *address = subject + subjectOffset;
                unsigned char *lastHit_addr = subject + (uint2)(prev & 0xffff);

                // int bin_id = sequenceCount & 127;
                int bin_id = sequenceCount & 31;

                struct ungappedExtension *ungappedExtension_current =
                    ungappedExtension_w +
                    bin_id * parametersFP->ungappedExtensionsPerThread;

                // If only one hit triggered this extension
                ungappedExtension_twoHitExtendD_sm(
                        sequence, queryOffset, address, lastHit_addr, subject,
                        &ungappedExtension_subjectEndReachedFP,
                        parametersFP->encoding_numCodes,
                        parametersFP->statistics_ungappedNominalDropoff,
                        parametersFP->blast_ungappedNominalTrigger,
                        parametersFP->ungappedExtensionsPerThread,
                        ungappedExtension_current, ungappedExtension_additional,
                        &(numOfTriggerExtensions_w[bin_id]), sequenceCount, scoreMatrix_s,
                        querySequence_s);

                prev_ext = (curr & 0xffffffffffff0000) +
                    (ungappedExtension_subjectEndReachedFP - subject);
            }
        }
    }
    blast_numTriggerExtensions[tid] = numOfTriggerExtensions_w[laneId];
    blast_numUngappedExtensions[tid] += blast_numUngappedExtensions_s;
}

__device__ struct coordinate ungappedExtension_findProteinSeed(
        struct ungappedExtension *ungappedExtension, struct PSSMatrixFP PSSMatrixFP,
        unsigned char *subject, unsigned char encoding_numCodes) {
    int2 *queryWindowStart, *queryWindowEnd;
    unsigned char *subjectWindowStart, *subjectWindowEnd;

    int2 *bestQueryPosition;
    unsigned char *bestSubjectPosition;
    int4 bestSegmentScore;
    int4 nominalScore, count;
    struct coordinate seed;

    if (ungappedExtension->end.queryOffset -
            ungappedExtension->start.queryOffset <
            11) {
        // The seed point is the middle of the extension
        seed.queryOffset = (ungappedExtension->end.queryOffset +
                ungappedExtension->start.queryOffset) /
            2;
        seed.subjectOffset = (ungappedExtension->end.subjectOffset +
                ungappedExtension->start.subjectOffset) /
            2;
    } else {
        // Else find the highest scoring length-11 segment of the ungapped extension
        queryWindowStart = queryWindowEnd =
            PSSMatrixFP.matrix +
            ungappedExtension->start.queryOffset * encoding_numCodes;
        subjectWindowStart = subjectWindowEnd =
            subject + ungappedExtension->start.subjectOffset;

        // Find initial score for first 11 positions
        nominalScore = 0;
        count = 0;
        while (count < 11) {
            nominalScore += queryWindowEnd[*subjectWindowEnd];
            queryWindowEnd += encoding_numCodes;
            subjectWindowEnd++;
            count++;
        }

        queryWindowEnd -= encoding_numCodes;
        subjectWindowEnd--;

        // By default first-11 positions gives best position and score
        bestQueryPosition = queryWindowStart;
        bestSubjectPosition = subjectWindowStart;
        bestSegmentScore = nominalScore;

        // Now slide the window across and record the better scores/positions
        while (queryWindowEnd <
                PSSMatrixFP.matrix +
                ungappedExtension->end.queryOffset * encoding_numCodes) {
            // Advance window end, add new position value
            queryWindowEnd += encoding_numCodes;
            subjectWindowEnd++;

            nominalScore += queryWindowEnd[*subjectWindowEnd];
            // Remove position that we will leave behind
            nominalScore -= queryWindowStart[*subjectWindowStart];

            queryWindowStart += encoding_numCodes;
            subjectWindowStart++;

            // Check if best window position yet
            if (nominalScore > bestSegmentScore) {
                bestSegmentScore = nominalScore;
                bestQueryPosition = queryWindowStart;
                bestSubjectPosition = subjectWindowStart;
            }
        }

        // Middle of the best window is the seed position
        seed.queryOffset =
            (bestQueryPosition - PSSMatrixFP.matrix) / encoding_numCodes + 5;
        seed.subjectOffset = bestSubjectPosition + 5 - subject;
    }
    return seed;
}
