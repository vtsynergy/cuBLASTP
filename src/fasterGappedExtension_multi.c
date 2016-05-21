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


#include "blast.h"


int4* fasterGappedExtension_matchRow_t[MAX_NUM_THREADS] = {NULL};
int4* fasterGappedExtension_insertQrow_t[MAX_NUM_THREADS] = {NULL};
int4* fasterGappedExtension_insertSrow_t[MAX_NUM_THREADS] = {NULL};
int4 fasterGappedExtension_rowSizes_t[MAX_NUM_THREADS] = {0};
unsigned char **fasterGappedExtension_traceback_t[MAX_NUM_THREADS] = {NULL};
int4 fasterGappedExtension_numRows_t[MAX_NUM_THREADS] = {0};
int4 fasterGappedExtension_numColumns_t[MAX_NUM_THREADS] = {0};

struct dpResults fasterGappedExtension_dpAfterSeed_multi(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                             int4 dropoff, int4 subjectLength, int tid)
{
	int2 **queryPosition, **bestQueryPosition, **queryEnd;
	int2 *matrixColumn;
	unsigned char *rowDropoff, *columnDropoff;
	unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd;
	unsigned char **tracebackRow, *tracebackColumn;
	int4 bestScore = 0, dropoffThreshold;
	int4 *matchRow, *insertQrow, *insertSrow, rowOffset;
	int4 queryDistance, subjectDistance;
	int4 oldMatch, match, previousOldMatch, previousOldInsertS, previousOldInsertQ;
	int4 previousMatch, previousInsertS;
	struct dpResults dpResults;
	unsigned char leftOfDropoff;
	int4 queryLength;

	queryLength = PSSMatrix.length;
	subjectEnd = subject + subjectLength;
	queryEnd = PSSMatrix.matrix + queryLength;

	// Declare processing rows for storing match, insert-subject and insert-query values
	// If current malloced rows aren't big enough
	if (subjectLength >= fasterGappedExtension_rowSizes_t[tid])
	{
		// Free existing rows
		free(fasterGappedExtension_matchRow_t[tid]);
		free(fasterGappedExtension_insertQrow_t[tid]);
		free(fasterGappedExtension_insertSrow_t[tid]);
		// Set size to double current needed length
		fasterGappedExtension_rowSizes_t[tid] = subjectLength * 2;
		// Malloc new rows
		fasterGappedExtension_matchRow_t[tid] = (int4*)malloc(sizeof(int4) * fasterGappedExtension_rowSizes_t[tid]);
		fasterGappedExtension_insertQrow_t[tid] = (int4*)malloc(sizeof(int4) * fasterGappedExtension_rowSizes_t[tid]);
		fasterGappedExtension_insertSrow_t[tid] = (int4*)malloc(sizeof(int4) * fasterGappedExtension_rowSizes_t[tid]);
	}

	// Determine lowest score before dropoff
	dropoffThreshold = -dropoff;

//    printf("%d,%d (%d,%d) --\n", queryLength, subjectLength,
//                                 fasterGappedExtension_numRows_t[tid], fasterGappedExtension_numColumns_t[tid]); fflush(stdout);

    // Increase number of columns in traceback array if neccessary
    if (subjectLength > fasterGappedExtension_numColumns_t[tid])
    {
    	// For each existing row
        queryDistance = 0;
        while (queryDistance < fasterGappedExtension_numRows_t[tid])
        {
            // Increase number of columns
            fasterGappedExtension_traceback_t[tid][queryDistance]
                = (unsigned char*)realloc(fasterGappedExtension_traceback_t[tid][queryDistance],
                                          sizeof(unsigned char) * (subjectLength));

            queryDistance++;
        }

        // Update number of columns
        fasterGappedExtension_numColumns_t[tid] = subjectLength;
	}

    // If more rows are required
    if (queryLength > fasterGappedExtension_numRows_t[tid])
    {
        // Increase number of row pointers
        fasterGappedExtension_traceback_t[tid] = (unsigned char**)realloc(fasterGappedExtension_traceback_t[tid],
                                    sizeof(unsigned char*) * queryLength);

        // Declare new rows
        while (fasterGappedExtension_numRows_t[tid] < queryLength)
        {
			fasterGappedExtension_traceback_t[tid][fasterGappedExtension_numRows_t[tid]]
            	= (unsigned char*)malloc(sizeof(unsigned char) * fasterGappedExtension_numColumns_t[tid]);

            fasterGappedExtension_numRows_t[tid]++;
        }
    }

    bestSubjectPosition = subjectPosition = subject + 1;
	bestQueryPosition = queryPosition = PSSMatrix.matrix + 1;

	// Initialize rows
	matchRow = fasterGappedExtension_matchRow_t[tid] + 1;
	insertQrow = fasterGappedExtension_insertQrow_t[tid] + 1;
	insertSrow = fasterGappedExtension_insertSrow_t[tid] + 1;

	// Set initial row dropoff and column dropoff
	rowDropoff = subject + subjectLength - 1;
	columnDropoff = subject;

	// Initialize traceback pointers
	tracebackRow = fasterGappedExtension_traceback_t[tid] + (queryPosition - PSSMatrix.matrix);
	tracebackColumn = *tracebackRow + (subjectPosition - subject);

	// -----FIRST ROW-----

	// Using first column of the query matrix
	matrixColumn = (*queryPosition);

	// -----FIRST CELL-----
	// Set M value for top-left cell
	match = matrixColumn[*subjectPosition];
	*matchRow = match;
	// Set DUMMY Ix and Iy values, which should never be used
	*insertSrow = constants_gappedExtensionDummyValue;
	*insertQrow = constants_gappedExtensionDummyValue;
	// M came from M
	*tracebackColumn = 0;

	// If this is the best-yet scoring cell
	if (match > bestScore)
	{
		// Update best start cell data
		bestScore = match;
		dropoffThreshold = bestScore - dropoff;
		bestQueryPosition = queryPosition;
		bestSubjectPosition = subjectPosition;
	}

	// Record match and insertS for this about-to-be-previous cell
	previousMatch = match;
	previousInsertS = *insertSrow;

	subjectDistance = 0;
	subjectPosition++; matchRow++; insertQrow++; insertSrow++; tracebackColumn++;

	// ----- REMAINING CELLS -----
	// For each remaining columns in the top row, scanning from left-to-right
	while (subjectPosition < subjectEnd)
	{
		// Set value for M
		match = matrixColumn[*subjectPosition]
		      - parameters_openGap - subjectDistance * parameters_extendGap;
		*matchRow = match;

		// Set value for Ix
		if (previousInsertS - parameters_extendGap >
			previousMatch - parameters_openGap)
		{
			*insertSrow = previousInsertS - parameters_extendGap;
			// M came from Ix and Ix came from Ix
			*tracebackColumn = 5;
		}
		else
		{
			*insertSrow = previousMatch - parameters_openGap;
			// M came from Ix and Ix came from M
			*tracebackColumn = 1;
		}

		// Set DUMMY Iy value, which should never be used
		*insertQrow = constants_gappedExtensionDummyValue;

		// If this is the best-yet scoring cell
		if (match > bestScore)
		{
			// Update best start cell data
			bestScore = match;
			dropoffThreshold = bestScore - dropoff;
			bestQueryPosition = queryPosition;
			bestSubjectPosition = subjectPosition;
		}

		// If score at current cell is below dropoff
		if (dropoffThreshold > match &&
		    dropoffThreshold > *insertSrow)
		{
			// Record dropoff position
			rowDropoff = subjectPosition;
			// And stop processing row
			break;
		}

		// Record match and insertS for this about-to-be-previous cell
		previousMatch = match;
		previousInsertS = *insertSrow;

		subjectPosition++; matchRow++; insertQrow++; insertSrow++; tracebackColumn++;
		subjectDistance++;
	}

//    if (dloc==88197331)
//    print2(fasterGappedExtension_matchRow_t[tid], subject, rowDropoff, columnDropoff);

	queryDistance = 0;
	queryPosition++; tracebackRow++;

	// -----REMAINING ROWS-----
	while (queryPosition < queryEnd && rowDropoff > columnDropoff)
	{
		subjectPosition = columnDropoff + 1;
		tracebackColumn = *tracebackRow + (subjectPosition - subject);

		// Reset rows
		rowOffset = (subjectPosition - subject);
		matchRow = fasterGappedExtension_matchRow_t[tid] + rowOffset;
		insertQrow = fasterGappedExtension_insertQrow_t[tid] + rowOffset;
		insertSrow = fasterGappedExtension_insertSrow_t[tid] + rowOffset;

		// Using next column of the query matrix
		matrixColumn = (*queryPosition);

		// -----FAR LEFT CELL-----
		// Record some old values
		previousOldMatch = *matchRow;
		previousOldInsertQ = *insertQrow;
		previousOldInsertS = *insertSrow;

		// Set Iy value
		if (*insertQrow - parameters_extendGap > *matchRow - parameters_openGap)
		{
			*insertQrow = *insertQrow - parameters_extendGap;
			// Iy is derived from Iy, M is derived from Iy
			*tracebackColumn = 34;
		}
		else
		{
			*insertQrow = *matchRow - parameters_openGap;
			// Iy is derived from M, M is derived from Iy
			*tracebackColumn = 2;
		}

		// Set DUMMY values for M and Iy, which should never be used
		match = *matchRow = constants_gappedExtensionDummyValue;
		*insertSrow = constants_gappedExtensionDummyValue;

		// If score at current cell is below dropoff
		if (dropoffThreshold > *insertQrow)
		{
			// Record dropoff position
			columnDropoff = subjectPosition;
			leftOfDropoff = 1;
		}
		else
		{
			// We are left of the column dropoff for this row
			leftOfDropoff = 0;
		}

		// Record match and insertS for this about-to-be-previous cell
		previousMatch = match;
		previousInsertS = *insertSrow;

		subjectPosition++; matchRow++; insertQrow++; insertSrow++; tracebackColumn++;

		// -----CELLS LEFT OF ROW DROPOFF-----
		while (subjectPosition <= rowDropoff)
		{
			// Remember old M value (for cell below this one)
			oldMatch = *matchRow;

			// Calculate new M value
			if (previousOldMatch >= previousOldInsertQ)
			{
				if (previousOldMatch >= previousOldInsertS)
				{
					match = matrixColumn[*subjectPosition] + previousOldMatch;
					// M is derived from M
					*tracebackColumn = 0;
				}
				else
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertS;
					// M is derived from Ix
					*tracebackColumn = 1;
				}
			}
			else
			{
				if (previousOldInsertQ >= previousOldInsertS)
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertQ;
					// M is derived from Iy
					*tracebackColumn = 2;
				}
				else
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertS;
					// M is derived from Ix
					*tracebackColumn = 1;
				}
			}
			
			*matchRow = match;

			// Record some old values
			previousOldMatch = oldMatch;
			previousOldInsertQ = *insertQrow;
			previousOldInsertS = *insertSrow;

			// Set new Iy value
			if (oldMatch - parameters_openGap >= *insertQrow - parameters_extendGap)
			{
				*insertQrow = oldMatch - parameters_openGap;
				// Iy is derived from M
				// No change to traceback
			}
			else
			{
				*insertQrow = *insertQrow - parameters_extendGap;
				// Iy is derived from Iy
				*tracebackColumn |= 32;
			}
			// Calculate new Ix
			if (previousMatch - parameters_openGap >= previousInsertS - parameters_extendGap)
			{
				*insertSrow = previousMatch - parameters_openGap;
				// Ix is derived from M
				// No change to traceback
			}
			else
			{
				*insertSrow = previousInsertS - parameters_extendGap;
				// Ix is derived from Ix
				*tracebackColumn |= 4;
			}

			// If this is the best-yet scoring cell
			if (match > bestScore)
			{
				// Update best start cell data
				bestScore = match;
				dropoffThreshold = bestScore - dropoff;
				bestQueryPosition = queryPosition;
				bestSubjectPosition = subjectPosition;
			}

			// If score at current cell (and cells to its left) are below dropoff
			if (leftOfDropoff)
			{
				if (dropoffThreshold > match &&
					dropoffThreshold > *insertSrow &&
					dropoffThreshold > *insertQrow)
				{
					// Record dropoff position
					columnDropoff = subjectPosition;
				}
				else
				{
					// We are left of the column dropoff for this row
					leftOfDropoff = 0;
				}
			}

			// Record match and insertS for this about-to-be-previous cell
			previousMatch = match;
			previousInsertS = *insertSrow;

			subjectPosition++; matchRow++; insertQrow++; insertSrow++; tracebackColumn++;
		}

		// -----CELLS RIGHT OF ROW DROPOFF -----
		if (!(dropoffThreshold > previousMatch &&
			dropoffThreshold > previousInsertS &&
			dropoffThreshold > *(insertQrow - 1)))
		{
			while (subjectPosition < subjectEnd)
			{
				// Set value for Ix
				*insertSrow = previousInsertS - parameters_extendGap;
				// Ix came from Ix
				*tracebackColumn = 4;

				// Set DUMMY values for M and Ix, which should never be used
				*matchRow = constants_gappedExtensionDummyValue;
				*insertQrow = constants_gappedExtensionDummyValue;

				// If score at current cell is below dropoff
				if (dropoffThreshold > *insertSrow)
				{
					// And stop processing row
					subjectPosition++;
					break;
				}

				// Record insertS for this about-to-be-previous cell
				previousInsertS = *insertSrow;

				subjectPosition++; matchRow++; insertQrow++; insertSrow++; tracebackColumn++;
				subjectDistance++;
			}
		}

		// Record dropoff position
		rowDropoff = subjectPosition - 1;

		queryDistance++;
		queryPosition++; tracebackRow++;

//        if (dloc==88197331)
//		print2(fasterGappedExtension_matchRow_t[tid], subject, rowDropoff, columnDropoff);
	}

	dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
	dpResults.best.subjectOffset = bestSubjectPosition - subject;
	dpResults.bestScore = bestScore;
	dpResults.traceback = fasterGappedExtension_traceback_t[tid];

	return dpResults;
}
struct dpResults fasterGappedExtension_dpBeforeSeed_multi(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                              struct coordinate seed, int4 dropoff, int tid)
{
	int2 **queryPosition, **bestQueryPosition;
	int2 *matrixColumn;
	unsigned char *rowDropoff, *columnDropoff;
	unsigned char *subjectPosition, *bestSubjectPosition;
	unsigned char **tracebackRow, *tracebackColumn;
	int4 bestScore = 0, dropoffThreshold;
	int4 *matchRow, *insertQrow, *insertSrow, rowOffset;
	int4 queryDistance, subjectDistance;
	int4 oldMatch, match, previousOldMatch, previousOldInsertS, previousOldInsertQ;
	int4 previousMatch, previousInsertS;
	struct dpResults dpResults;
	unsigned char rightOfDropoff;

	// Declare processing rows for storing match, insert-subject and insert-query values
	// If current malloced rows aren't big enough
	if (seed.subjectOffset >= fasterGappedExtension_rowSizes_t[tid])
	{
		// Free existing rows
		free(fasterGappedExtension_matchRow_t[tid]);
		free(fasterGappedExtension_insertQrow_t[tid]);
		free(fasterGappedExtension_insertSrow_t[tid]);
		// Set size to double current needed length
		fasterGappedExtension_rowSizes_t[tid] = (seed.subjectOffset) * 2;
		// Malloc new rows
		fasterGappedExtension_matchRow_t[tid] = (int4*)malloc(sizeof(int4) * fasterGappedExtension_rowSizes_t[tid]);
		fasterGappedExtension_insertQrow_t[tid] = (int4*)malloc(sizeof(int4) * fasterGappedExtension_rowSizes_t[tid]);
		fasterGappedExtension_insertSrow_t[tid] = (int4*)malloc(sizeof(int4) * fasterGappedExtension_rowSizes_t[tid]);
	}

	// Determine lowest score before dropoff
	dropoffThreshold = -dropoff;

//    printf("%d,%d (%d,%d)\n", seed.queryOffset, seed.subjectOffset,
//                              fasterGappedExtension_numRows_t[tid], fasterGappedExtension_numColumns_t[tid]); fflush(stdout);

    // Increase number of columns in traceback array if neccessary
    if (seed.subjectOffset > fasterGappedExtension_numColumns_t[tid])
    {
    	// For each existing row
        queryDistance = 0;
        while (queryDistance < fasterGappedExtension_numRows_t[tid])
        {
            // Increase number of columns
            fasterGappedExtension_traceback_t[tid][queryDistance]
                = (unsigned char*)realloc(fasterGappedExtension_traceback_t[tid][queryDistance],
                                          sizeof(unsigned char) * (seed.subjectOffset));

            queryDistance++;
        }

        // Update number of columns
        fasterGappedExtension_numColumns_t[tid] = seed.subjectOffset;
	}

    // If more rows are required
    if (seed.queryOffset > fasterGappedExtension_numRows_t[tid])
    {
        // Increase number of row pointers
        fasterGappedExtension_traceback_t[tid] = (unsigned char**)realloc(fasterGappedExtension_traceback_t[tid],
                                    sizeof(unsigned char*) * (seed.queryOffset));

        // Declare new rows
        while (fasterGappedExtension_numRows_t[tid] < seed.queryOffset)
        {
			fasterGappedExtension_traceback_t[tid][fasterGappedExtension_numRows_t[tid]]
            	= (unsigned char*)malloc(sizeof(unsigned char) * (fasterGappedExtension_numColumns_t[tid]));

            fasterGappedExtension_numRows_t[tid]++;
        }
    }

//    printf("%d,%d (%d,%d) AFTER\n", seed.queryOffset, seed.subjectOffset,
//                              fasterGappedExtension_numRows_t[tid], fasterGappedExtension_numColumns_t[tid]); fflush(stdout);

    bestSubjectPosition = subjectPosition = subject + seed.subjectOffset - 1;
	bestQueryPosition = queryPosition = PSSMatrix.matrix + seed.queryOffset - 1;

	// Initialize row pointers
	rowOffset = (subjectPosition - subject);
	matchRow = fasterGappedExtension_matchRow_t[tid] + rowOffset;
	insertQrow = fasterGappedExtension_insertQrow_t[tid] + rowOffset;
	insertSrow = fasterGappedExtension_insertSrow_t[tid] + rowOffset;

	// Set initial row dropoff and column dropoff
	rowDropoff = subject;
	columnDropoff = subject + seed.subjectOffset;

	// Initialize traceback pointers
	tracebackRow = fasterGappedExtension_traceback_t[tid] + (queryPosition - PSSMatrix.matrix);
	tracebackColumn = *tracebackRow + (subjectPosition - subject);

//    printf("[%d,%d]\n", seed.subjectOffset, (subjectPosition - subject)); fflush(stdout);
//    printf("[%d,%d]\n", seed.queryOffset, (queryPosition - PSSMatrix.matrix)); fflush(stdout);

	// -----FIRST ROW-----

	// Using first column of query matrix
	matrixColumn = *queryPosition;

//    printf("[%d]", matrixColumn); fflush(stdout);
//    printf("[%d]", subjectPosition); fflush(stdout);

	// -----FIRST CELL-----
	// Set M value for bottom-right cell
	match = matrixColumn[*subjectPosition];
	*matchRow = match;

	// Set DUMMY Ix and Iy values, which should never be used
	*insertSrow = constants_gappedExtensionDummyValue;
	*insertQrow = constants_gappedExtensionDummyValue;

	// M came from M
	*tracebackColumn = 0;

	// If this is the best-yet scoring cell
	if (match > bestScore)
	{
		// Update best start cell data
		bestScore = *matchRow;
		dropoffThreshold = bestScore - dropoff;
		bestQueryPosition = queryPosition;
		bestSubjectPosition = subjectPosition;
	}

	// Record match and insertS for this about-to-be-previous cell
	previousMatch = match;
	previousInsertS = *insertSrow;

	subjectDistance = 0;
	subjectPosition--; matchRow--; insertSrow--; insertQrow--; tracebackColumn--;

	// ----- REMAINING CELLS -----
	// For each remaining column in the bottom row, scanning from right-to-left
	while (subjectPosition >= subject)
	{
		// Set value for M
		match = matrixColumn[*subjectPosition]
		      - parameters_openGap - subjectDistance * parameters_extendGap;
		*matchRow = match;

		// Set value for Ix
		if (previousInsertS - parameters_extendGap > previousMatch - parameters_openGap)
		{
			*insertSrow = previousInsertS - parameters_extendGap;
			// M came from Ix and Ix came from Ix
			*tracebackColumn = 5;
		}
		else
		{
			*insertSrow = previousMatch - parameters_openGap;
			// M came from Ix and Ix came from M
			*tracebackColumn = 1;
		}

		// Set DUMMY Iy value, which should never be used
		*insertQrow = constants_gappedExtensionDummyValue;

		// If this is the best-yet scoring cell
		if (match > bestScore)
		{
			// Update best start cell data
			bestScore = match;
			dropoffThreshold = bestScore - dropoff;
			bestQueryPosition = queryPosition;
			bestSubjectPosition = subjectPosition;
		}

		// If score at current cell is below dropoff
		if (dropoffThreshold > match &&
		    dropoffThreshold > *insertSrow)
		{
			// Record dropoff position
			rowDropoff = subjectPosition;
			// And stop processing row
			break;
		}

		// Record match and insertS for this about-to-be-previous cell
		previousMatch = match;
		previousInsertS = *insertSrow;

		subjectPosition--; matchRow--; insertSrow--; insertQrow--; tracebackColumn--;
		subjectDistance++;
	}

//	print(fasterGappedExtension_matchRow_t[tid], subject, rowDropoff, columnDropoff);

	queryDistance = 0;

	// -----REMAINING ROWS-----
	while (queryPosition > PSSMatrix.matrix && rowDropoff < columnDropoff)
	{
		queryPosition--; tracebackRow--;
		subjectPosition = columnDropoff - 1;
		tracebackColumn = *tracebackRow + (subjectPosition - subject);

		// Reset row pointers to start of rows
		rowOffset = (subjectPosition - subject);
		matchRow = fasterGappedExtension_matchRow_t[tid] + rowOffset;
		insertQrow = fasterGappedExtension_insertQrow_t[tid] + rowOffset;
		insertSrow = fasterGappedExtension_insertSrow_t[tid] + rowOffset;

		// Using next column of query matrix
		matrixColumn = *queryPosition;

		// -----FAR RIGHT CELL-----
		// Record some old values
		previousOldMatch = *matchRow;
		previousOldInsertQ = *insertQrow;
		previousOldInsertS = *insertSrow;

		// Set Iy value
		if (*insertQrow - parameters_extendGap > *matchRow - parameters_openGap)
		{
			*insertQrow = *insertQrow - parameters_extendGap;
			// Iy is derived from Iy, M is derived from Iy
			*tracebackColumn = 34;
		}
		else
		{
			*insertQrow = *matchRow - parameters_openGap;
			// Iy is derived from M, M is derived from Iy
			*tracebackColumn = 2;
		}

		// Set DUMMY values for M and Iy, which should never be used
		match = *matchRow = constants_gappedExtensionDummyValue;
		*insertSrow = constants_gappedExtensionDummyValue;

		// If score at current cell is below dropoff
		if (dropoffThreshold > *insertQrow)
		{
			// Record dropoff position
			columnDropoff = subjectPosition;
			rightOfDropoff = 1;
		}
		else
		{
			// We are left of the column dropoff for this row
			rightOfDropoff = 0;
		}

		// Record match and insertS for this about-to-be-previous cell
		previousMatch = match;
		previousInsertS = *insertSrow;

		subjectPosition--; matchRow--; insertSrow--; insertQrow--; tracebackColumn--;

		// -----CELLS RIGHT OF ROW DROPOFF-----
		while (subjectPosition >= rowDropoff)
		{
			// Remember old M value (for cell below this one)
			oldMatch = *matchRow;

			// Calculate new M value
			if (previousOldMatch >= previousOldInsertQ)
			{
				if (previousOldMatch >= previousOldInsertS)
				{
					match = matrixColumn[*subjectPosition] + previousOldMatch;
					// M is derived from M
					*tracebackColumn = 0;
				}
				else
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertS;
					// M is derived from Ix
					*tracebackColumn = 1;
				}
			}
			else
			{
				if (previousOldInsertQ >= previousOldInsertS)
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertQ;
					// M is derived from Iy
					*tracebackColumn = 2;
				}
				else
				{
					match = matrixColumn[*subjectPosition] + previousOldInsertS;
					// M is derived from Ix
					*tracebackColumn = 1;
				}
			}

			*matchRow = match;

			// Record some old values
			previousOldMatch = oldMatch;
			previousOldInsertQ = *insertQrow;
			previousOldInsertS = *insertSrow;

			// Set new Iy value
			if (oldMatch - parameters_openGap >= *insertQrow - parameters_extendGap)
			{
				*insertQrow = oldMatch - parameters_openGap;
				// Iy is derived from M
				// No change to traceback
			}
			else
			{
				*insertQrow = *insertQrow - parameters_extendGap;
				// Iy is derived from Iy
				*tracebackColumn |= 32;
			}
			// Calculate new Ix
			if (previousMatch - parameters_openGap >= previousInsertS - parameters_extendGap)
			{
				*insertSrow = previousMatch - parameters_openGap;
				// Ix is derived from M
				// No change to traceback
			}
			else
			{
				*insertSrow = previousInsertS - parameters_extendGap;
				// Ix is derived from Ix
				*tracebackColumn |= 4;
			}

			// If this is the best-yet scoring cell
			if (match > bestScore)
			{
				// Update best start cell data
				bestScore = match;
				dropoffThreshold = bestScore - dropoff;
				bestQueryPosition = queryPosition;
				bestSubjectPosition = subjectPosition;
			}

			// If score at current cell (and cells to its right) are below dropoff
			if (rightOfDropoff)
			{
				if (dropoffThreshold > match &&
					dropoffThreshold > *insertSrow &&
					dropoffThreshold > *insertQrow)
				{
					// Record dropoff position
					columnDropoff = subjectPosition;
				}
				else
				{
					// We are left of the column dropoff for this row
					rightOfDropoff = 0;
				}
			}

			// Record match and insertS for this about-to-be-previous cell
			previousMatch = match;
			previousInsertS = *insertSrow;

			subjectPosition--; matchRow--; insertSrow--; insertQrow--; tracebackColumn--;
		}

		// -----CELLS LEFT OF ROW DROPOFF -----
		if (!(dropoffThreshold > previousMatch &&
			dropoffThreshold > previousInsertS &&
			dropoffThreshold > *(insertQrow + 1)))
		{
			while (subjectPosition >= subject)
			{
				// Set value for Ix
				*insertSrow = previousInsertS - parameters_extendGap;
				// Ix came from Ix
				*tracebackColumn = 4;

				// Set DUMMY values for M and Ix, which should never be used
				*matchRow = constants_gappedExtensionDummyValue;
				*insertQrow = constants_gappedExtensionDummyValue;

				// If score at current cell is below dropoff
				if (dropoffThreshold > *insertSrow)
				{
					// Stop processing row
					subjectPosition--;
					break;
				}

				// Record match and insertS for this about-to-be-previous cell
				previousInsertS = *insertSrow;

				subjectPosition--; matchRow--; insertSrow--; insertQrow--; tracebackColumn--;
				subjectDistance++;
			}
		}

		// Record dropoff position
		rowDropoff = subjectPosition + 1;

//		print(fasterGappedExtension_matchRow_t[tid], subject, rowDropoff, columnDropoff);

		queryDistance++;
	}

	dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
	dpResults.best.subjectOffset = bestSubjectPosition - subject;
	dpResults.bestScore = bestScore;
	dpResults.traceback = fasterGappedExtension_traceback_t[tid];

	return dpResults;
}

struct gappedExtension* fasterGappedExtension_build_multi(struct ungappedExtension* ungappedExtension,
                                  struct PSSMatrix PSSMatrix, int4 subjectSize,
                                  unsigned char* subject, int4 dropoff, int tid)
{
	struct coordinate seed;
	unsigned char *choppedSubject;
	struct dpResults beforeDpResults, afterDpResults;
	struct trace beforeTrace, afterTrace, trace;
	struct PSSMatrix choppedPSSMatrix;
	int4 choppedSubjectSize;
	struct gappedExtension* gappedExtension;
    int4 strandOffset = 0;

	// Perform dynamic programming for points before the seed
	seed = ungappedExtension->seed;
    if (seed.queryOffset > PSSMatrix.strandLength)
    {
	    // If query position is in the second strand, remove first strand from PSSM
        strandOffset = PSSMatrix.strandLength;
		seed.queryOffset -= PSSMatrix.strandLength;
		PSSMatrix = PSSMatrix_chop(PSSMatrix, PSSMatrix.strandLength);
    }
    else
    {
    	// Otherwise remove second strand
    	PSSMatrix.length = PSSMatrix.strandLength;
    }

//    printf("Seed=%d,%d Length=%d,%d\n", seed.queryOffset, seed.subjectOffset, PSSMatrix.length, subjectSize);

	beforeDpResults = fasterGappedExtension_dpBeforeSeed_multi(subject, PSSMatrix, seed, dropoff, tid);

	// Trace back and create the trace which specifies the first half of the alignment
	beforeTrace = fasterGappedExtension_traceBeforeSeed(beforeDpResults, seed);

	// Chop the start off the query and subject so they begin at the seed
	choppedPSSMatrix = PSSMatrix_chop(PSSMatrix, seed.queryOffset);
	choppedSubject = subject + seed.subjectOffset;
	choppedSubjectSize = subjectSize - (seed.subjectOffset);

	// Perform dynamic programming for points after the seed
	afterDpResults = fasterGappedExtension_dpAfterSeed_multi(choppedSubject, choppedPSSMatrix,
	                                              dropoff, choppedSubjectSize, tid);

	// Trace back to get the trace for the seed onwards
	afterTrace = fasterGappedExtension_traceAfterSeed(afterDpResults, choppedPSSMatrix.length);

	// Join afterTrace to the end of beforeTrace
	trace = fasterGappedExtension_joinTraces(beforeTrace, afterTrace);
	free(afterTrace.traceCodes);

    // Adjust coordinates if extension was performed in the second strand
    afterDpResults.best.queryOffset += strandOffset;
    beforeDpResults.best.queryOffset += strandOffset;
	trace.queryStart += strandOffset;

//	printf("Final trace length=%d\n", trace.length);

	// Create gapped extension
	gappedExtension = (struct gappedExtension*)malloc(sizeof(struct gappedExtension));
	gappedExtension->trace = trace;
	gappedExtension->next = NULL;

	// Start of afterTrace is end of the gapped extension, but we need to add seed position
	// to get correct offset
	gappedExtension->queryEnd = seed.queryOffset + afterTrace.queryStart + strandOffset;
	gappedExtension->subjectEnd = seed.subjectOffset + afterTrace.subjectStart;

//	if (dloc == 88197331)
//		printf("final[%d,%d,%d](%d)\n", beforeDpResults.bestScore, afterDpResults.bestScore,
//		choppedPSSMatrix.matrix[0][choppedSubject[0]], seed.queryOffset);

    // Determine score by combining score from the two traces, and the match score at
	// the seed position
	gappedExtension->nominalScore = beforeDpResults.bestScore + afterDpResults.bestScore
	                             + choppedPSSMatrix.matrix[0][choppedSubject[0]];

    // Update ungappedExtension start/end
    ungappedExtension->start.queryOffset = trace.queryStart;
    ungappedExtension->end.queryOffset = gappedExtension->queryEnd;
    ungappedExtension->start.subjectOffset = trace.subjectStart;
    ungappedExtension->end.subjectOffset = gappedExtension->subjectEnd;
	ungappedExtension->nominalScore = gappedExtension->nominalScore;

    return gappedExtension;
}




void fasterGappedExtension_free_multi(int n_thread)
{
	// Free memory used by traceback array
    int ii = 0;
    for(ii = 0; ii < n_thread; ii++)
    {
        while (fasterGappedExtension_numRows_t[ii] > 0)
        {
            fasterGappedExtension_numRows_t[ii]--;
            free(fasterGappedExtension_traceback_t[ii][fasterGappedExtension_numRows_t[ii]]);
        }
        free(fasterGappedExtension_traceback_t[ii]);

        free(fasterGappedExtension_matchRow_t[ii]);
        free(fasterGappedExtension_insertQrow_t[ii]);
        free(fasterGappedExtension_insertSrow_t[ii]);
    }
}
