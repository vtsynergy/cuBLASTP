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
        );
 
int get_total_numHits(int *numHits, int num_bins);

void get_bin_offset(int *numHits, int *binOffset, int num_bins);
