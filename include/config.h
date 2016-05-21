#ifdef USE_SM_20
#define __ldg(x) (*(x)) 
#define NUM_BLOCK 112
#else
#define NUM_BLOCK 260
#endif

#define BLOCK_SIZE 128
#define BIN_X 4 // CHANGE BIN_SIZE 
#define BIN_POWER 7 // BIN_X
#define BIN_MARK 127 // BIN_X
