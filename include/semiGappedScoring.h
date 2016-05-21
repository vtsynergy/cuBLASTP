#ifndef _semiGappedScoring_
#define _semiGappedScoring_

#ifdef __cplusplus
extern "C" {
#endif

// Perform semi-gapped alignment with restricted insertion
int4 semiGappedScoring_score(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                        int4 subjectSize, unsigned char* subject, int4 dropoff);
int4 semiGappedScoring_score_multi(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                        int4 subjectSize, unsigned char* subject, int4 dropoff, int TID);

void semiGappedScoring_free();
void semiGappedScoring_free_multi(int n_threads);

#ifdef __cplusplus
}
#endif

#endif
