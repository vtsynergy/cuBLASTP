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
 * Developer: Jing Zhang 
 *
 */


#define _GNU_SOURCE
#include "blast.h"
#include <pthread.h>

int numAboveCutoff_g = 0;

int get_affinity(int *cpu_map) {
#pragma omp parallel
  { cpu_map[omp_get_thread_num()] = sched_getcpu(); }
  return 0;
}

typedef struct {
  int tid;
  struct alignment **alignment_list;
  int blast_numTriggerSequences;
  struct PSSMatrix PSSMatrix;
  struct PSSMatrixFP PSSMatrixFP;
  struct alignment **good_alignment_list;
  int *bestScore_list;
} thread_aux_findGoodAlignment;

typedef struct {
  int tid;
  struct finalAlignment **finalAlignment_list;
  struct PSSMatrix PSSMatrix;
  struct PSSMatrixFP PSSMatrixFP;
  int numOfAlignment;
  int4 *numProcessed;
  int4 *numAboveCutoff;
} thread_aux_getTracebacks;

void findGoodAlignments_core(int tid, struct alignment **alignment_list,
                             int blast_numTriggerSequences,
                             struct PSSMatrix PSSMatrix,
                             struct PSSMatrixFP PSSMatrixFP,
                             struct alignment **good_alignment_list,
                             int *bestScore_list) {

  int ii = 0;
  struct alignment *alignment;
  struct ungappedExtension *ungappedExtension;
  int4 bestScore, numExtensions, hasChildren;
  for (ii = 0; ii < blast_numTriggerSequences; ii++) {
    alignment = alignment_list[ii];
    bestScore = 0;
    //blast_dloc = alignment->descriptionLocation;

    // debug=======================
    // printf("count = %d\n", count++);
    //===============================

    // Record if subject has children
    if (encoding_alphabetType == encoding_protein &&
        alignment->encodedLength > alignment->subjectLength + 2)
      hasChildren = 1;
    else
      hasChildren = 0;

    // For each ungapped extension (in descending order of score)
    numExtensions = 0;
    ungappedExtension = alignment->ungappedExtensions;
    while (ungappedExtension != NULL) {
      if (ungappedExtension->status != ungappedExtension_DELETED) {
        // Find the seed
        ungappedExtension_findSeed(ungappedExtension, PSSMatrix, PSSMatrixFP,
                                   alignment->subject);

        // Semi-gapped scoring
        if (parameters_semiGappedScoring) {
          ungappedExtension->nominalScore = semiGappedScoring_score_multi(
              ungappedExtension, PSSMatrix, alignment->subjectLength,
              alignment->subject, statistics_gappedNominalDropoff, tid);
          // ungappedExtension->nominalScore
          //= semiGappedScoring_score(ungappedExtension, PSSMatrix,
          //alignment->subjectLength,
          // alignment->subject, statistics_gappedNominalDropoff);

          // Mark as semigapped
          ungappedExtension->status = ungappedExtension_SEMIGAPPED;
        } else {
          printf("Error: multithreading goodAlignment only supports "
                 "semiGappedAlignment!\n");
          exit(0);
        }

        // If alignment scores above R1 cutoff
        if (ungappedExtension->nominalScore >= blast_nominalR1cutoff) {
          if (hasChildren) {
            // Subject has children so perform stage1 and 2 on children
            alignments_expandCluster(alignment, PSSMatrix, PSSMatrixFP);
            bestScore = 0;
            break;
          } else if (ungappedExtension->nominalScore > bestScore) {
            // Update best score for the alignment
            bestScore = ungappedExtension->nominalScore;
          }
        } else {
          // Else mark it as deleted
          ungappedExtension->status = ungappedExtension_DELETED;
        }

        // Remove any ungapped extensions in this alignment that are in the area
        // covered
        // by the gapped scoring just performed
        alignments_pruneRegion(alignment, ungappedExtension);

        numExtensions++;
      }

      ungappedExtension = ungappedExtension->next;
    }

    // If this alignment contains gapped extensions that could score above
    // cutoff
    if (bestScore >= blast_nominalR1cutoff) {
      // If a single sequence add to list of "good" alignments
      // alignments_addGoodAlignment(bestScore, alignment);
      // int pos_good_alignment = __sync_fetch_and_add(&num_good_alignment, 1);
      int pos_good_alignment =
          __sync_fetch_and_add(&blast_numGoodAlignments, 1);
      bestScore_list[pos_good_alignment] = bestScore;
      good_alignment_list[pos_good_alignment] = alignment;
      __sync_fetch_and_add(&blast_numGoodExtensions, numExtensions);
    }
  }
}

static void *findGoodAlignments_worker(void *data) {
  thread_aux_findGoodAlignment *d = (thread_aux_findGoodAlignment *)data;
  findGoodAlignments_core(d->tid, d->alignment_list,
                          d->blast_numTriggerSequences, d->PSSMatrix,
                          d->PSSMatrixFP,
                          // d->num_good_alignment,
                          d->good_alignment_list, d->bestScore_list);
}

void alignments_findGoodAlignments_multithread(struct PSSMatrix PSSMatrix,
                                               struct PSSMatrixFP PSSMatrixFP,
                                               int n_threads) {
  struct alignment *alignment;
  struct alignment **alignment_list = (struct alignment **)malloc(
      sizeof(struct alignment *) * blast_numTriggerSequences);
  struct alignment **good_alignment_list = (struct alignment **)malloc(
      sizeof(struct alignment *) * blast_numTriggerSequences);
  int4 *bestScore_list =
      (int4 *)malloc(sizeof(int4) * blast_numTriggerSequences);

  int ii = 0;
  memBlocks_resetCurrent(alignments_alignments);
  int max_sub_len = 0;
  while ((alignment = (struct alignment *)memBlocks_getCurrent(
              alignments_alignments)) != NULL) {
    max_sub_len = max(alignment->subjectLength, max_sub_len);
    alignment_list[ii] = alignment;
    ii++;
  }
  // printf("Max Subject Length: %d\n", max_sub_len);
  int num_alignment = ii;
  // debug=========================
  int count = 0;
  //==============================

  // For each alignment
  volatile int num_good_alignment = 0;

  pthread_t *tid;
  pthread_attr_t attr;
  thread_aux_findGoodAlignment *data;
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
  data = (thread_aux_findGoodAlignment *)calloc(
      n_threads, sizeof(thread_aux_findGoodAlignment));
  tid = (pthread_t *)calloc(n_threads, sizeof(pthread_t));
  int block_size = num_alignment / n_threads;
  int cpu_map[MAX_NUM_THREADS];
  get_affinity(cpu_map);
  for (ii = 0; ii < n_threads; ii++) {

    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(ii, &cpuset);
    data[ii].tid = ii;
    data[ii].alignment_list = alignment_list + block_size * ii;
    if (ii == n_threads - 1)
      data[ii].blast_numTriggerSequences =
          block_size + num_alignment % n_threads;
    else
      data[ii].blast_numTriggerSequences = block_size;
    data[ii].PSSMatrix = PSSMatrix;
    data[ii].PSSMatrixFP = PSSMatrixFP;
    data[ii].good_alignment_list = good_alignment_list;
    data[ii].bestScore_list = bestScore_list;
    pthread_create(&tid[ii], &attr, findGoodAlignments_worker, data + ii);
    pthread_setaffinity_np(tid[ii], sizeof(cpu_set_t), &cpuset);
  }

  for (ii = 0; ii < n_threads; ++ii)
    pthread_join(tid[ii], 0);
  free(data);
  free(tid);
  // while ((alignment = (struct
  // alignment*)memBlocks_getCurrent(alignments_alignments)) != NULL)
  //#pragma omp parallel for shared(alignment_list, PSSMatrix, PSSMatrixFP)
  for (ii = 0; ii < blast_numGoodAlignments; ii++)
    alignments_addGoodAlignment(bestScore_list[ii], good_alignment_list[ii]);

  // Record point to list of alignments for this volume
  alignments_volumeAlignments[alignments_numVolumes] = alignments_alignments;
  alignments_numVolumes++;

  // Construct new list for next volume (if there is one)
  alignments_alignments = memBlocks_initialize(
      sizeof(struct alignment), constants_initialAllocAlignments);
  free(alignment_list);
  free(good_alignment_list);
  free(bestScore_list);
}

void getTracebacks_core(int tid, struct finalAlignment **finalAlignment_list,
                        struct PSSMatrix PSSMatrix,
                        struct PSSMatrixFP PSSMatrixFP, int numOfAlignment,
                        int4 *numProcessed, int4 *numAboveCutoff) {
  struct finalAlignment *finalAlignment;
  struct gappedExtension *gappedExtension;
  struct ungappedExtension *ungappedExtension, *highestScoringExtension,
      oldUngappedExtension;
  struct alignment *alignment;
  struct unpackRegion *unpackRegion;
  int4 alignmentCount;
  unsigned char *subject;
  int ii = 0;

  // printf("numOfAlignment: %d\n", numOfAlignment);
  // for(ii = 0; ii < numOfAlignment; ii++)
  while ((ii < numOfAlignment) &&
         (numAboveCutoff_g < parameters_numDisplayTracebacks ||
          parameters_numDisplayTracebacks == 0)) {
    // printf("numAboveCutoff: %d\n", numAboveCutoff_g);
    finalAlignment = finalAlignment_list[ii];
    alignment = finalAlignment_list[ii]->alignment;
    //blast_dloc = alignment->descriptionLocation;

    // If traceback haven't been computed for this alignment
    if (alignment->gappedExtensions == NULL) {
      // Unpack part or all of the subject (for nucleotide)
      unpack_unpackSubject(PSSMatrix, alignment);

      // For each ungapped extension that hasn't been deleted
      ungappedExtension = alignment->ungappedExtensions;
      while (ungappedExtension != NULL) {
        // If extension scores above cutoff
        if (ungappedExtension->status != ungappedExtension_DELETED) {
          // Make copy of ungapped extension
          oldUngappedExtension = *ungappedExtension;

          // If subject and query are short enough and sequence does not have
          // multiple
          // unpack regions, use faster but less memory efficient gapped
          // alignment with traceback
          if (((uint8)PSSMatrix.length * (uint8)alignment->subjectLength <
               (uint8)constants_maximumTracebackSize) &&
              unpack_entireSubjectUnpacked(alignment)) {
            gappedExtension = fasterGappedExtension_build_multi(
                ungappedExtension, PSSMatrix, alignment->subjectLength,
                alignment->unpackRegions[0].unpackedSubject,
                statistics_gappedFinalNominalDropoff, tid);
          }
          // Otherwise use slower but more memory-efficient gapped alignment
          else {
            unpackRegion = unpack_selectRegion(
                alignment->unpackRegions, alignment->numUnpackRegions,
                ungappedExtension->seed.subjectOffset);

            gappedExtension = gappedExtension_build_multi(
                ungappedExtension, PSSMatrix, alignment->subjectLength,
                alignment->subject, unpackRegion,
                statistics_gappedFinalNominalDropoff, tid);
          }

          // Calculate normalized score and e-value
          gappedExtension_score(gappedExtension);

          // Add it to the current alignment
          if (gappedExtension->nominalScore >= blast_gappedNominalCutoff)
            alignments_addGappedExtension(alignment, gappedExtension);
          else {
            free(gappedExtension->trace.traceCodes);
            free(gappedExtension);
          }

          // Check for ungapped extensions that were mistakenly pruned
          if (parameters_tableScoring)
            alignments_unpruneRegion(alignment, &oldUngappedExtension,
                                     ungappedExtension, PSSMatrix, PSSMatrixFP);

          // Remove any ungapped extensions in this alignment that are in the
          // area covered
          // by the gapped scoring just performed
          alignments_pruneRegion(alignment, ungappedExtension);
        }
        ungappedExtension = ungappedExtension->next;
      }

      // Finally prune extensions that share a start or end point
      alignments_pruneOverlappingExtensions(alignment);

      //                printf("Was %d Now %d\n",
      // finalAlignment->highestNominalScore,
      // alignment->gappedExtensions->nominalScore);

      // Update final alignment's high-score to that of the first (and
      // highest-scoring)
      // gapped extension in the list
      if (alignment->gappedExtensions != NULL)
        finalAlignment->highestNominalScore =
            alignment->gappedExtensions->nominalScore;
      else
        finalAlignment->highestNominalScore = 0;

      // numProcessed++;
      //__sync_fetch_and_add(numProcessed, 1);
      //	            printf("Computed alignment score=%d\n",
      //finalAlignment->highestNominalScore);
    }

    // Tally number of final alignments above cutoff
    if (finalAlignment->highestNominalScore > 0)
      __sync_fetch_and_add(&numAboveCutoff_g, 1);
    ii++;
  }
}

static void *getTracebacks_worker(void *data) {
  thread_aux_getTracebacks *d = (thread_aux_getTracebacks *)data;
  getTracebacks_core(d->tid, d->finalAlignment_list, d->PSSMatrix,
                     d->PSSMatrixFP, d->numOfAlignment, d->numProcessed,
                     d->numAboveCutoff);
}
// Get the tracebacks for all of the final alignments
void alignments_getTracebacks_multithread(struct PSSMatrix PSSMatrix,
                                          struct PSSMatrixFP PSSMatrixFP,
                                          int n_threads) {
  struct finalAlignment *finalAlignment;
  struct gappedExtension *gappedExtension;
  struct ungappedExtension *ungappedExtension, *highestScoringExtension,
      oldUngappedExtension;
  struct alignment *alignment;
  struct unpackRegion *unpackRegion;
  int4 numProcessed = 0, numAboveCutoff = 0, repeatComputeTracebacks = 1;
  int4 alignmentCount;
  unsigned char *subject;

  // Sort final alignments by score
  alignments_sortFinalAlignments();

  //    alignments_printFinalAlignments();

  // Only keep alignments above cutoff
  while (alignments_finalAlignments->numEntries > 0) {
    // Get last alignment in list
    finalAlignment = (struct finalAlignment *)memSingleBlock_getLastEntry(
        alignments_finalAlignments);

    // Stop if above cutoff
    if (finalAlignment->highestNominalScore != 0)
      break;

    // Otherwise remove it from the list
    alignments_finalAlignments->numEntries--;
  }

  // For each alignment that is in the top numDisplayAlignments but not the top
  // numDisplayTracebacks
  alignmentCount = parameters_numDisplayTracebacks;
  while (alignmentCount < parameters_numDisplayAlignments &&
         alignmentCount < alignments_finalAlignments->numEntries) {
    finalAlignment = (struct finalAlignment *)memSingleBlock_getEntry(
        alignments_finalAlignments, alignmentCount);
    alignment = finalAlignment->alignment;
    //blast_dloc = alignment->descriptionLocation;

    // Get the highest scoring ungapped extension
    ungappedExtension = alignment->ungappedExtensions;
    highestScoringExtension = ungappedExtension;
    while (ungappedExtension != NULL) {
      if (ungappedExtension->nominalScore >
              highestScoringExtension->nominalScore &&
          ungappedExtension->status != ungappedExtension_DELETED) {
        highestScoringExtension = ungappedExtension;
      }
      ungappedExtension = ungappedExtension->next;
    }

    if (highestScoringExtension != NULL) {
      subject = alignments_selectRegion(alignment, highestScoringExtension);

      // Perform gapped scoring with higher dropoff
      highestScoringExtension->nominalScore = gappedScoring_score(
          highestScoringExtension, PSSMatrix, alignment->subjectLength, subject,
          statistics_gappedFinalNominalDropoff);

      finalAlignment->highestNominalScore =
          highestScoringExtension->nominalScore;
    }

    //        printf("Rescore with larger dropoff num %d: Score=%d\n",
    // alignmentCount, highestScoringExtension->nominalScore);
    alignmentCount++;
  }

  while (repeatComputeTracebacks) {
    numAboveCutoff = 0;
    numProcessed = 0;
    memSingleBlock_resetCurrent(alignments_finalAlignments);

    struct finalAlignment **finalAlignment_list =
        (struct finalAlignment **)malloc(
            sizeof(struct finalAlignment *) *
            alignments_finalAlignments->numEntries);

    int ii = 0;
    while ((finalAlignment = (struct finalAlignment *)memSingleBlock_getCurrent(
                alignments_finalAlignments)) != NULL) {
      finalAlignment_list[ii] = finalAlignment;
      ii++;
    }

    int cpu_map[MAX_NUM_THREADS];
    int num_finalAlignment = ii;
    // printf("num_finalAlignment: %d\n", num_finalAlignment);
    pthread_t *tid;
    pthread_attr_t attr;
    thread_aux_getTracebacks *data;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
    data = (thread_aux_getTracebacks *)calloc(n_threads,
                                              sizeof(thread_aux_getTracebacks));
    tid = (pthread_t *)calloc(n_threads, sizeof(pthread_t));
    int block_size = num_finalAlignment / n_threads;
    numAboveCutoff_g = 0;
    get_affinity(cpu_map);
    for (ii = 0; ii < n_threads; ii++) {
      cpu_set_t cpuset;
      CPU_ZERO(&cpuset);
      CPU_SET(cpu_map[ii], &cpuset);
      data[ii].tid = ii;
      data[ii].finalAlignment_list = finalAlignment_list + block_size * ii;
      if (ii == n_threads - 1)
        data[ii].numOfAlignment = block_size + num_finalAlignment % n_threads;
      else
        data[ii].numOfAlignment = block_size;
      data[ii].PSSMatrix = PSSMatrix;
      data[ii].PSSMatrixFP = PSSMatrixFP;
      data[ii].numProcessed = &numProcessed;
      data[ii].numAboveCutoff = &numAboveCutoff;
      pthread_create(&tid[ii], &attr, getTracebacks_worker, data + ii);
      pthread_setaffinity_np(tid[ii], sizeof(cpu_set_t), &cpuset);
    }

    for (ii = 0; ii < n_threads; ++ii)
      pthread_join(tid[ii], 0);
    free(data);
    free(tid);

    free(finalAlignment_list);

    //    	printf("Traceback alignments performed=%d\n", numProcessed);

    // Sort final alignments by score
    alignments_sortFinalAlignments();

    //        printf("repeatComputeTracebacks:");

    // If the first numDisplayTracebacks alignments have traceback computed,
    // stop
    numProcessed = 0;
    repeatComputeTracebacks = 0;
    memSingleBlock_resetCurrent(alignments_finalAlignments);
    while ((finalAlignment = (struct finalAlignment *)memSingleBlock_getCurrent(
                alignments_finalAlignments)) != NULL &&
           numProcessed < parameters_numDisplayTracebacks) {
      alignment = finalAlignment->alignment;
      if (alignment->gappedExtensions == NULL &&
          finalAlignment->highestNominalScore != 0) {
        //            	printf("1");
        repeatComputeTracebacks = 1;
        break;
      }
      //            printf("0");

      numProcessed++;
    }
  }

  // Only keep top N alignments
  if (parameters_numDisplayAlignments != 0 &&
      alignments_finalAlignments->numEntries >
          parameters_numDisplayAlignments) {
    alignments_finalAlignments->numEntries = parameters_numDisplayAlignments;
  }

  // Only keep alignments above cutoff
  while (alignments_finalAlignments->numEntries > 0) {
    // Get last alignment in list
    finalAlignment = (struct finalAlignment *)memSingleBlock_getLastEntry(
        alignments_finalAlignments);

    // Stop if above cutoff
    if (finalAlignment->highestNominalScore != 0)
      break;

    // Otherwise remove it from the list
    alignments_finalAlignments->numEntries--;
  }
}

