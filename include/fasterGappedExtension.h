#ifndef _fasterGappedExtension_
#define _fasterGappedExtension_

#ifdef __cplusplus
extern "C" {
#endif

// Build a gapped extension with a trace and nominal score from the seed point4 of an ungapped
// extension using dynamic programming
struct gappedExtension* fasterGappedExtension_build(struct ungappedExtension* ungappedExtension,
                                  struct PSSMatrix PSSMatrix, int4 subjectSize,
                                  unsigned char* subject, int4 dropoff);
struct gappedExtension* fasterGappedExtension_build_multi(struct ungappedExtension* ungappedExtension,
                                  struct PSSMatrix PSSMatrix, int4 subjectSize,
                                  unsigned char* subject, int4 dropoff, int tid);


// Debugging routine
void fasterGappedExtension_printBeforeRow(int4* row, unsigned char* subject, unsigned char* rowDropoff,
                                    unsigned char* columnDropoff);
// Debugging routine
void fasterGappedExtension_printAfterRow(int4* row, unsigned char* subject, unsigned char* rowDropoff,
                                   unsigned char* columnDropoff);

// Given a gapped extension with a nominal score, calculate the normalized score
// and E-Value
void fasterGappedExtension_score(struct gappedExtension* gappedExtension);

// Given a gappedExtension and list of ungappedExtensions, prune the latter to
// remove those which overlap/int4ersect with the gappedExtension
void fasterGappedExtension_prune(struct gappedExtension* gappedExtension,
                           struct ungappedExtension* ungappedExtension);

void fasterGappedExtension_free();
void fasterGappedExtension_free_multi(int n_thread);


struct trace fasterGappedExtension_traceBeforeSeed(struct dpResults beforeDpResults, struct coordinate seed);
struct trace fasterGappedExtension_traceAfterSeed(struct dpResults beforeDpResults, int4 queryLength);

struct trace fasterGappedExtension_joinTraces(struct trace beforeTrace, struct trace afterTrace);

#ifdef __cplusplus
}
#endif


#endif

