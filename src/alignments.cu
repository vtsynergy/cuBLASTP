/*
 * cuBLASTP - Fine-Grained Parallelization of Protein Sequence Search on CPU+GPU 
 * Version 0.1 (beta)
 *
 * (c) 2015 Virginia Polytechnic Institute & State University (Virginia Tech) 
 *
 * This version of cuBLASTP is licensed for non-commercial use only,
 * as specified in LICENSE files in licensing directory. For all other use
 * contact vtiplicensing@vtip.org
 * 
 * Developer: Jing Zhang, Shucai Xiao 
 *
 */




#include "blast.h"
#include <semiGappedScoringGPU.cu>
#include <sys/time.h>

struct memBlocks* alignments_alignments;

void getAlignmentsNum(int *alignmentsNum,
        int *ungappedExtensionNum)
{
    struct alignment *alignment;
    struct ungappedExtension *ungappedExtension;
    int alignmentsNo, ungappedExtensionNo, totalUngappedNum = 0;
    alignmentsNo = 0;
    memBlocks_resetCurrent(alignments_alignments);
    while ((alignment = (struct alignment *)memBlocks_getCurrent(alignments_alignments)) != NULL)
    {
        alignmentsNo++;
        ungappedExtensionNo = 0;
        ungappedExtension = alignment->ungappedExtensions;
        while (ungappedExtension != NULL)
        {
            ungappedExtensionNo++;
            ungappedExtension = ungappedExtension->next;
        }
        totalUngappedNum += ungappedExtensionNo;
    }

    *alignmentsNum = alignmentsNo;
    *ungappedExtensionNum = totalUngappedNum;

    return;
}


//void copyAlignmentsToBuffer(struct alignment *alignmentArray,
void copyAlignmentsToBuffer(struct ungappedExtension *ungappedExtensionArray,
        int *ungappedExtensionStartLoc,
        int *ungappedExtensionNum)
{
    struct alignment *alignment;
    struct ungappedExtension *ungappedExtension;
    struct ungappedExtension *curUngappedExtension;

    int alignmentNo, ungappedExtensionNo;
    int ungappedStartLoc = 0;
    alignmentNo = 0;
    memBlocks_resetCurrent(alignments_alignments);
    while ((alignment = (struct alignment *)memBlocks_getCurrent(alignments_alignments)) != NULL)
    {
        ungappedExtensionStartLoc[alignmentNo] = ungappedStartLoc;
        curUngappedExtension = &ungappedExtensionArray[ungappedStartLoc];
        ungappedExtension = alignment->ungappedExtensions;
        ungappedExtensionNo = 0;
        while (ungappedExtension != NULL)
        {
            memcpy(curUngappedExtension + ungappedExtensionNo, ungappedExtension, sizeof(struct ungappedExtension));
            ungappedExtensionNo++;
            ungappedExtension = ungappedExtension->next;
        }

        ungappedExtensionNum[alignmentNo] = ungappedExtensionNo;
        ungappedStartLoc += ungappedExtensionNo;
        alignmentNo++;
    }
}

void updataUngappedExtensions(struct ungappedExtension *ungappedExtensionArray)
{
    struct alignment *alignment;
    struct ungappedExtension *ungappedExtension;
    struct ungappedExtension *curUngappedExtension;

    int alignmentNo, ungappedExtensionNo = 0;
    memBlocks_resetCurrent(alignments_alignments);
    while ((alignment = (struct alignment *)memBlocks_getCurrent(alignments_alignments)) != NULL)
    {
        ungappedExtension = alignment->ungappedExtensions;
        while (ungappedExtension != NULL)
        {
            memcpy(ungappedExtension, ungappedExtensionArray + ungappedExtensionNo, sizeof(struct ungappedExtension));
            ungappedExtensionNo++;
            ungappedExtension = ungappedExtension->next;
        }
    }

    return;
}

void addGoodAlignment(int4 *bestScoreArray, 
        int4 *numGoodExtensions, 
        int4 *numSemiGapping)
{
    struct alignment *alignment;
    int alignmentNo = 0;
    memBlocks_resetCurrent(alignments_alignments);
    while ((alignment = (struct alignment *)memBlocks_getCurrent(alignments_alignments)) != NULL)
    {
        if (bestScoreArray[alignmentNo] >= blast_nominalR1cutoff)
        {
            alignments_addGoodAlignment(bestScoreArray[alignmentNo], alignment);
            blast_numGoodExtensions += numGoodExtensions[alignmentNo];
            blast_numGoodAlignments++;
        }

        blast_numSemiGapped += numSemiGapping[alignmentNo];

        alignmentNo++;
    }

    return;
}

//sort the alignment according to the num of ungapped extensions contained.
void sortAlignments(int *ungappedExtensionNum, 
        int *orderArray,
        int nBegin, int nEnd)
{
    int piv, l, r, temp;

    if (nEnd > nBegin + 1)
    {
        piv = ungappedExtensionNum[nBegin];
        l = nBegin + 1;
        r = nEnd;

        while (l < r)
        {
            if (ungappedExtensionNum[l] <= piv)
            {
                l++;
            }
            else
            {
                --r;
                temp = ungappedExtensionNum[r];
                ungappedExtensionNum[r] = ungappedExtensionNum[l];
                ungappedExtensionNum[l] = temp;

                temp = orderArray[r];
                orderArray[r] = orderArray[l];
                orderArray[l] = temp;
            }
        }
        --l;
        temp = ungappedExtensionNum[l];
        ungappedExtensionNum[l] = ungappedExtensionNum[nBegin];
        ungappedExtensionNum[nBegin] = temp;

        temp = orderArray[l];
        orderArray[l] = orderArray[nBegin];
        orderArray[nBegin] = temp;

        sortAlignments(ungappedExtensionNum, orderArray, nBegin, l);
        sortAlignments(ungappedExtensionNum, orderArray, r, nEnd);
    }

    return;
}

void alignments_fingGoodAlignmentsGPU(struct PSSMatrixFP **PSSMatrixFP,		//GPU buffer
        struct PSSMatrixFP PSSMatrixFPH,
        struct scoreMatrix *scoreMatrixp,
        char **matrixBody,					//GPU buffer
        struct sequenceDataFP **sequenceDataFP, //GPU buffer
        unsigned char **subjectSeqs,		//GPU buffer
        int4 totalSubSize)
{
    int alignmentsNum, ungappedExtensionsNum;
    int *startLocArray, *ungappedExtensionsNumArray, *tempUngappedNumArray;
    int *startLocArrayD, *ungappedExtensionsNumArrayD;
    int4 	*orderArray, *orderArrayD;
    int4	*bestScoreArray, *bestScoreArrayD;
    struct ungappedExtension *ungappedExtensions;
    struct ungappedExtension *ungappedExtensionsD;
    int4	*numGoodExtensions, *numGoodExtensionsD;
    int4	*numSemiGapping, *numSemiGappingD;
    struct gappedExtensionParameters parameters, *parametersD;
    int4 i;

    //struct timeval t1, t2;
    timerStart();

    //get the number of alignments and ungappedExtensions
    getAlignmentsNum(&alignmentsNum, &ungappedExtensionsNum);
    startLocArray = (int *)global_malloc(sizeof(int) * alignmentsNum);
    orderArray = (int *)global_malloc(sizeof(int) * alignmentsNum);
    ungappedExtensionsNumArray = (int *)global_malloc(sizeof(int) * alignmentsNum);
    tempUngappedNumArray = (int *)global_malloc(sizeof(int) * alignmentsNum);
    ungappedExtensions = (struct ungappedExtension *)global_malloc(sizeof(struct ungappedExtension) * ungappedExtensionsNum);
    if (startLocArray == NULL ||
            orderArray == NULL ||
            ungappedExtensionsNumArray == NULL ||
            ungappedExtensions == NULL)
    {
        printf("In findGoodAlignmentGPU, allocate buffers error!\n");
        exit (1);
    }

    bestScoreArray = (int4 *)global_malloc(sizeof(int4) * alignmentsNum);
    numGoodExtensions = (int4 *)global_malloc(sizeof(int4) * alignmentsNum);
    numSemiGapping = (int4 *)global_malloc(sizeof(int4) * alignmentsNum);
    if (bestScoreArray == NULL ||
            numGoodExtensions == NULL ||
            numSemiGapping == NULL)
    {
        printf("Allocate buffer 1 error!\n");
        exit(1);
    }

    //allocate gpu memory for alignments
    cudaError_t cudaRes0, cudaRes1, cudaRes2, cudaRes3, cudaRes4;
    cudaRes1 = cudaMalloc((void **)&bestScoreArrayD, sizeof(int4) * alignmentsNum);
    cudaRes2 = cudaMalloc((void **)&orderArrayD, sizeof(int4) * alignmentsNum);
    cudaRes3 = cudaMalloc((void **)&ungappedExtensionsD, sizeof(struct ungappedExtension) * ungappedExtensionsNum);
    if (cudaRes1 != cudaSuccess ||
            cudaRes2 != cudaSuccess ||
            cudaRes3 != cudaSuccess)
    {
        printf("Allocate miscellaneous GPU memory error, %d, %d, %d\n", 
                cudaRes1,
                cudaRes2,
                cudaRes3);
        exit(1);
    }

    cudaRes0 = cudaMalloc((void **)&startLocArrayD, sizeof(int4) * alignmentsNum);
    cudaRes1 = cudaMalloc((void **)&ungappedExtensionsNumArrayD, sizeof(int4) * alignmentsNum);
    cudaRes2 = cudaMalloc((void **)&parametersD, sizeof(struct gappedExtensionParameters));
    cudaRes3 = cudaMalloc((void **)&numGoodExtensionsD, sizeof(int4) * alignmentsNum);
    cudaRes4 = cudaMalloc((void **)&numSemiGappingD, sizeof(int4) * alignmentsNum);
    if (cudaRes0 != cudaSuccess ||
            cudaRes1 != cudaSuccess ||
            cudaRes2 != cudaSuccess ||
            cudaRes3 != cudaSuccess ||
            cudaRes4 != cudaSuccess)
    {
        printf("Allocate miscellaneous1 GPU memory error, %d, %d, %d, %d, %d\n", 
                cudaRes0,
                cudaRes1,
                cudaRes2,
                cudaRes3,
                cudaRes4);
        exit(1);
    }

    //buffers for bestRow and insertRow
    int4 *bestRowD, *insertQRowD;
    cudaRes0 = cudaMalloc((void **)&bestRowD, sizeof(int4) * totalSubSize);
    cudaRes1 = cudaMalloc((void **)&insertQRowD, sizeof(int4) * totalSubSize);
    if (cudaRes0 != cudaSuccess ||
            cudaRes1 != cudaSuccess)
    {
        printf("Allocate row memory buffers error! %d, %d\n",
                cudaRes0,
                cudaRes1);
        exit(1);
    }

    //time
    timerEnd();
    strTime.iniTime += elapsedTime();

    //time
    timerStart();

    //copy alignments and ungaped extensions to the allocated buffer
    copyAlignmentsToBuffer(ungappedExtensions, startLocArray, ungappedExtensionsNumArray);

    //sort the alignments according to the number of ungapped extensions
    for (i = 0; i < alignmentsNum; i++)
    {
        orderArray[i] = i;
    }

    //time
    timerEnd();
    strTime.preprocessingTime += elapsedTime();

    //gettimeofday(&t1, NULL);
    //memcpy(tempUngappedNumArray, ungappedExtensionsNumArray, sizeof(int) * alignmentsNum);
    //sortAlignments(tempUngappedNumArray, orderArray, 0, alignmentsNum);
    //gettimeofday(&t2, NULL);

    //printf("sort time = %.3f\n", 1000.0 * (t2.tv_sec - t1.tv_sec) + (t2.tv_usec - t1.tv_usec)/1000.0);

    //time
    timerStart();
    //copy alignments info to GPU memory
    cudaMemcpy(ungappedExtensionsD, ungappedExtensions, sizeof(struct ungappedExtension) * ungappedExtensionsNum, cudaMemcpyHostToDevice);
    cudaMemset(numGoodExtensionsD, 0, sizeof(int4) * alignmentsNum);
    cudaMemset(numSemiGappingD, 0, sizeof(int4) * alignmentsNum);

    parameters.semiGappedOpenGap = parameters_semiGappedOpenGap;
    parameters.semiGappedExtendGap = parameters_semiGappedExtendGap;
    parameters.semiGappedExtensionN = parameters_semiGappedExtensionN;
    parameters.encoding_numCodes = encoding_numCodes;
    parameters.semiGappedDropoffIncrease = parameters_semiGappedDropoffIncrease;

    //copy parameters to GPU
    cudaMemcpy(parametersD, &parameters, sizeof(struct gappedExtensionParameters), 
            cudaMemcpyHostToDevice);

    cudaMemcpy(startLocArrayD, startLocArray, sizeof(int4) * alignmentsNum, cudaMemcpyHostToDevice);
    cudaMemcpy(orderArrayD, orderArray, sizeof(int4) * alignmentsNum, cudaMemcpyHostToDevice);
    cudaMemcpy(ungappedExtensionsNumArrayD, ungappedExtensionsNumArray,
            sizeof(int4) * alignmentsNum, cudaMemcpyHostToDevice);
    //debug=================================
    //Use constant memory to store score matrix
    int scoreMatrixSize = encoding_numCodes * encoding_numCodes;
    cudaMemcpyToSymbol(scoreMatrixC, 
            ((char *)scoreMatrixp->matrix) + sizeof(int2 *) * encoding_numCodes, 
            sizeof(int2) * scoreMatrixSize);

    //Use constant memory to store query sequence
    unsigned char *tempQueryCode;
    tempQueryCode = (unsigned char *)global_malloc(sizeof(unsigned char) * (PSSMatrixFPH.length + 2));
    memcpy(&tempQueryCode[1], PSSMatrixFPH.queryCodes, sizeof(unsigned char) * PSSMatrixFPH.length);
    tempQueryCode[0] = encoding_sentinalCode;
    tempQueryCode[PSSMatrixFPH.length + 1] = encoding_sentinalCode;
    cudaMemcpyToSymbol(querySequenceC, tempQueryCode, sizeof(unsigned char) * (PSSMatrixFPH.length + 2));
    free(tempQueryCode);
    //==============================================
    //timer
    timerEnd();
    strTime.copyTimeHostToDevice += elapsedTime();

    int4 blockSize = parameters_threadNum;
    int blockNum = (alignmentsNum - 1)/blockSize + 1;
    //	if (blockNum < 14) //14 is the number of SMs on the Tesla C2050
    //	{
    //		blockSize = 64;
    //		blockNum = (alignmentsNum - 1)/blockSize + 1;
    //	}
    dim3 dimGrid(blockNum, 1);
    dim3 dimBlock(blockSize, 1);

    //configure L1 cache as 48KB
    cudaFuncSetCacheConfig(semiGappedScoring_kernel, cudaFuncCachePreferL1);

    //time
    timerStart();

    semiGappedScoring_kernel<<<dimGrid, dimBlock>>>(*sequenceDataFP,
            *subjectSeqs,
            *PSSMatrixFP,
            *matrixBody,
            parametersD,
            startLocArrayD,
            ungappedExtensionsNumArrayD,
            alignmentsNum,
            ungappedExtensionsD,
            bestScoreArrayD,
            numGoodExtensionsD,
            numSemiGappingD,
            orderArrayD,
            bestRowD,
            insertQRowD,
            statistics_gappedNominalDropoff,
            blast_nominalR1cutoff);
    cudaThreadSynchronize();
    //time
    timerEnd();
    strTime.kernelTime += elapsedTime();

    //debug=================================
    printf("%s\n", cudaGetErrorString(cudaGetLastError()));
    //========================================

    //time
    timerStart();

    cudaMemcpy(numGoodExtensions, numGoodExtensionsD, sizeof(int4) * alignmentsNum, cudaMemcpyDeviceToHost);
    cudaMemcpy(numSemiGapping, numSemiGappingD, sizeof(int4) * alignmentsNum, cudaMemcpyDeviceToHost);
    cudaMemcpy(bestScoreArray, bestScoreArrayD, sizeof(int4) * alignmentsNum, cudaMemcpyDeviceToHost);
    cudaMemcpy(ungappedExtensions, ungappedExtensionsD, sizeof(struct ungappedExtension) * ungappedExtensionsNum, cudaMemcpyDeviceToHost);
    //time
    timerEnd();
    strTime.copyTimeDeviceToHost += elapsedTime();

    //time
    timerStart();

    updataUngappedExtensions(ungappedExtensions);
    addGoodAlignment(bestScoreArray, numGoodExtensions, numSemiGapping);

    alignments_volumeAlignments[alignments_numVolumes] = alignments_alignments;
    alignments_numVolumes++;
    alignments_alignments = memBlocks_initialize(sizeof(struct alignment),
            constants_initialAllocAlignments);

    //free memory
    free(startLocArray);
    free(ungappedExtensionsNumArray);
    free(tempUngappedNumArray);
    free(bestScoreArray);
    free(ungappedExtensions);
    free(numGoodExtensions);
    free(numSemiGapping);
    free(orderArray);


    cudaFree(startLocArrayD);
    cudaFree(ungappedExtensionsNumArrayD);
    cudaFree(bestScoreArrayD);
    cudaFree(ungappedExtensionsD);
    cudaFree(numGoodExtensionsD);
    cudaFree(numSemiGappingD);
    cudaFree(orderArrayD);
    cudaFree(parametersD);
    cudaFree(bestRowD);
    cudaFree(insertQRowD);

    //time
    timerEnd();
    strTime.postprocessingTime += elapsedTime();
}
