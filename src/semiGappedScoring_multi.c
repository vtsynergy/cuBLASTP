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
#include <pthread.h>
int4* semiGappedScoring_bestRow_t[MAX_NUM_THREADS] = {NULL};
int4* semiGappedScoring_insertQrow_t[MAX_NUM_THREADS] = {NULL};
int4 semiGappedScoring_rowSizes_t[MAX_NUM_THREADS] = {0};

// Perform dynamic programming to explore possible END points and alignments that start at
// the given seed and find the best score
struct dpResults semiGappedScoring_dpAfterSeed_multi(unsigned char* subject, struct PSSMatrix PSSMatrix,
                              int4 dropoff, int4 subjectLength, int TID)
{
    int2 **queryPosition, **bestQueryPosition, **queryEnd;
    int2* matrixColumn;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd, *startSubjectPosition;
    int4 bestScore = 0;
    int4 *bestRow, *insertQrow, insertS, rowOffset;
    int4 subjectDistance;
    int4 oldBest, match, previousOldBest;
    unsigned char leftOfDropoff;
    int4 queryLength;
    int4 queryCount, subjectCount;
    struct dpResults dpResults;

    queryLength = PSSMatrix.length;
    subjectEnd = subject + subjectLength;
    queryEnd = PSSMatrix.matrix + queryLength;

    // Declare processing rows for storing match, insert-subject and insert-query values
    // If current malloced rows aren't big enough
    if (subjectLength >= semiGappedScoring_rowSizes_t[TID])
    {
        // Free existing rows
        if(semiGappedScoring_rowSizes_t[TID] > 0)
        {
            free(semiGappedScoring_bestRow_t[TID]);
            free(semiGappedScoring_insertQrow_t[TID]);
        }
        // Set size to double current needed length
        semiGappedScoring_rowSizes_t[TID] = subjectLength * 2;
        // Malloc new rows
        semiGappedScoring_bestRow_t[TID] = (int4*)malloc(sizeof(int4) * semiGappedScoring_rowSizes_t[TID]);
        semiGappedScoring_insertQrow_t[TID] = (int4*)malloc(sizeof(int4) * semiGappedScoring_rowSizes_t[TID]);
    }

    bestSubjectPosition = subjectPosition = startSubjectPosition = subject + 1;
    bestQueryPosition = queryPosition = PSSMatrix.matrix + 1;

    // Initialize rows
    bestRow = semiGappedScoring_bestRow_t[TID] + 1;
    insertQrow = semiGappedScoring_insertQrow_t[TID] + 1;

    // Set initial row dropoff and column dropoff
    rowDropoff = subject + subjectLength - 1;
    columnDropoff = subject;

    // -----FIRST ROW-----

    // Using first column of query matrix
    matrixColumn = *queryPosition;

    // -----FIRST CELL-----
    // Set M value for top-left cell
    match = matrixColumn[*subjectPosition];

    // M must be the best
    *bestRow = match;

    // Only gap opens possible
    *insertQrow = insertS = match - parameters_semiGappedOpenGap;

    // If this is the best-yet scoring cell
    if (match > bestScore)
    {
        // Update best start cell data
        bestScore = match;
        bestQueryPosition = queryPosition;
        bestSubjectPosition = subjectPosition;
    }

    subjectDistance = 0;
    subjectPosition++; bestRow++; insertQrow++;

    // ----- REMAINING CELLS -----
    // For each remaining columns in the top row, scanning from left-to-right
    while (subjectPosition < subjectEnd)
    {
        // Set value for M
        match = matrixColumn[*subjectPosition]
              - parameters_semiGappedOpenGap - subjectDistance * parameters_semiGappedExtendGap;

        // Determine the best of M and Iy
        if (match > insertS)
        {
            *bestRow = match;

            // Calculate new Iy
            insertS = maximum(match - parameters_semiGappedOpenGap,
                              insertS - parameters_semiGappedExtendGap);
        }
        else
        {
            *bestRow = insertS;

            // Since M <= Iy, new Iy must derive from Iy
            insertS -= parameters_semiGappedExtendGap;
        }

        // Set DUMMY Ix value, which should never be used
        *insertQrow = constants_gappedExtensionDummyValue;

        // If this is the best-yet scoring cell
        if (match > bestScore)
        {
            // Update best start cell data
            bestScore = match;
            bestQueryPosition = queryPosition;
            bestSubjectPosition = subjectPosition;
        }

        // If score at current cell is below dropoff
        if (bestScore > *bestRow + dropoff)
        {
            // Record dropoff position
            rowDropoff = subjectPosition;
            // And stop processing row
            break;
        }

        subjectPosition++; bestRow++; insertQrow++;
        subjectDistance++;
    }

//    if (dloc == 88197331)
//        print2(semiGappedScoring_bestRow_t[TID] + 1, subject, rowDropoff, columnDropoff);

    // Start queryCount at N. Only allow insertS for every Nth row when queryCount
    // reaches 0
    queryCount = parameters_semiGappedExtensionN;

    queryPosition++; queryCount--;

    // -----REMAINING ROWS-----
    while (queryPosition < queryEnd && rowDropoff > columnDropoff)
    {
        subjectPosition = columnDropoff + 1;

        // Determine subjectCount for initial subjectPosition. Is used to only allow
        // insertQ when subjectOffset % parameters_semiGappedExtensionN == 0
        subjectCount = ((int4)(subjectPosition - startSubjectPosition) % parameters_semiGappedExtensionN);
        if (subjectCount)
            subjectCount = parameters_semiGappedExtensionN - subjectCount;

        // Reset rows
        rowOffset = (subjectPosition - subject);
        bestRow = semiGappedScoring_bestRow_t[TID] + rowOffset;
        insertQrow = semiGappedScoring_insertQrow_t[TID] + rowOffset;

        // Using next column of query matrix
        matrixColumn = *queryPosition;


        // ************ All rows we are not allowing insertS
        if (queryCount)
        {
            // ** No insertQ allowed this column, this cell will only get a DUMMY score
            if (subjectCount)
            {
                previousOldBest = *bestRow;
                *bestRow = constants_gappedExtensionDummyValue;

                // Score at this cell is below dropoff
                columnDropoff = subjectPosition;
                leftOfDropoff = 1;
            }
            // ** We are allowing insertQ this column
            else
            {
                // -----FAR LEFT CELL-----
                // Record some old values
                previousOldBest = *bestRow;

                // Set Ix value
                *bestRow = *insertQrow;
                *insertQrow -= parameters_semiGappedExtendGap;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
                {
                    // Record dropoff position
                    columnDropoff = subjectPosition;
                    leftOfDropoff = 1;
                }
                else
                {
                    // We are right of the column dropoff for this row
                    leftOfDropoff = 0;
                }

                // Reset subjectCount
                subjectCount = parameters_semiGappedExtensionN;
            }

            subjectPosition++; bestRow++; insertQrow++; subjectCount--;

            // -----CELLS LEFT OF ROW DROPOFF-----
            while (subjectPosition <= rowDropoff)
            {
                // ** We are not allowing insertQ this column
                if (subjectCount)
                {
                    // Calculate new M value, which is also the best
                    oldBest = *bestRow;
                    match = *bestRow = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }
                }
                // We are allowing insertQ this column
                else
                {
                    // Calculate new M value
                    oldBest = *bestRow;
                    match = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // Determine the best of M and Ix
					if (match > *insertQrow)
                    {
						*bestRow = match;

                        // Calculate new Ix
                        *insertQrow = maximum(match - parameters_semiGappedOpenGap,
                                              *insertQrow - parameters_semiGappedExtendGap);
                    }
                    else
                    {
						*bestRow = *insertQrow;

                        // Since M <= Ix, new Ix must derive from Ix
						*insertQrow -= parameters_semiGappedExtendGap;
                    }

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }

                    // Reset subjectCount
                    subjectCount = parameters_semiGappedExtensionN;
                }

                subjectPosition++; bestRow++; insertQrow++; subjectCount--;
            }

            // -----SINGLE CELL RIGHT OF ROW DROPOFF -----
            if (!(bestScore > previousOldBest + dropoff) && (subjectPosition < subjectEnd))
            {
                // Set value for best
                *bestRow = match = previousOldBest + matrixColumn[*subjectPosition];

                // Set DUMMY values for Ix
                *insertQrow = constants_gappedExtensionDummyValue;

                if (match + dropoff >= bestScore)
                {
                    // Record dropoff position
                    rowDropoff = subjectPosition;
                }
            }
        }

        // ************ Every Nth row we allow insertS
        else
        {
            // -----FAR LEFT CELL-----

            // ** No insertQ allowed this column, this cell will only get a DUMMY score
            if (subjectCount)
            {
                previousOldBest = *bestRow;
                *bestRow = constants_gappedExtensionDummyValue;

                // Score at this cell is below dropoff
                columnDropoff = subjectPosition;
                leftOfDropoff = 1;
            }
            // ** We are allowing insertQ this column
            else
            {
                // Record some old values
                previousOldBest = *bestRow;

                 // Set Ix value
                *bestRow = *insertQrow;
                *insertQrow -= parameters_semiGappedExtendGap;

                // Set DUMMY value for Iy, which should never be used
                insertS = constants_gappedExtensionDummyValue;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
                {
                    // Record dropoff position
                    columnDropoff = subjectPosition;
                    leftOfDropoff = 1;
                }
                else
                {
                    // We are right of the column dropoff for this row
                    leftOfDropoff = 0;
                }

                // Reset subjectCount
                subjectCount = parameters_semiGappedExtensionN;
            }

            subjectPosition++; bestRow++; insertQrow++; subjectCount--;

            // -----CELLS LEFT OF ROW DROPOFF-----
            while (subjectPosition <= rowDropoff)
            {
                // ** We are not allowing insertQ this column
                if (subjectCount)
                {
                    // Remember old M value (for cell below this one)
                    oldBest = *bestRow;
                    match = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // Determine the best of M and Iy
					if (match > insertS)
                    {
						*bestRow = match;

                        // Calculate new Iy
                        insertS = maximum(match - parameters_semiGappedOpenGap,
                                          insertS - parameters_semiGappedExtendGap);
                    }
                    else
                    {
						*bestRow = insertS;

                        // Since M <= Iy, new Iy must derive from Iy
						insertS -= parameters_semiGappedExtendGap;
                    }

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }

                    // If score at current cell (and cells to its left) are below dropoff
                    if (leftOfDropoff)
                    {
                        if (bestScore > *bestRow + dropoff)
                        {
                            // Record dropoff position
                            columnDropoff = subjectPosition;
                        }
                        else
                        {
                            // We are right of the column dropoff for this row
                            leftOfDropoff = 0;
                        }
                    }
                }
                // ** We are allowing insertQ this column
                else
                {
                    // Remember old M value (for cell below this one)
                    oldBest = *bestRow;
                    match = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // Determine the best of M, Ix and Iy
                    if (match > insertS)
                    {
                        if (match > *insertQrow)
                        {
                            // Match is largest
                            *bestRow = match;

                            // Calculate new Ix
                            *insertQrow = maximum(match - parameters_semiGappedOpenGap,
                                                  *insertQrow - parameters_semiGappedExtendGap);

                            // Calculate new Iy
                            insertS = maximum(match - parameters_semiGappedOpenGap,
                                              insertS - parameters_semiGappedExtendGap);

                            // If this is the best-yet scoring cell
                            if (match > bestScore)
                            {
                                // Update best start cell data
                                bestScore = match;
                                bestQueryPosition = queryPosition;
                                bestSubjectPosition = subjectPosition;
                            }
                        }
                        else
                        {
                            // insertQ is largest
                            *bestRow = *insertQrow;

                            // Calculate new Ix
                            *insertQrow -= parameters_semiGappedExtendGap;

							// Dummy Iy
                            insertS = constants_gappedExtensionDummyValue;
                        }
                    }
                    else
                    {
                        if (insertS > *insertQrow)
                        {
                            // insertS is largest
                            *bestRow = insertS;

                            // Dummy Ix
                            *insertQrow = constants_gappedExtensionDummyValue;

                            // Calculate new Iy
                            insertS -= parameters_semiGappedExtendGap;

                        }
                        else
                        {
                            // insertQ is largest
                            *bestRow = *insertQrow;

                            // Calculate new Ix
                            *insertQrow -= parameters_semiGappedExtendGap;

							// Dummy Iy
                            insertS = constants_gappedExtensionDummyValue;
                        }
                    }

                    // If score at current cell (and cells to its left) are below dropoff
                    if (leftOfDropoff)
                    {
                        if (bestScore > *bestRow + dropoff)
                        {
                            // Record dropoff position
                            columnDropoff = subjectPosition;
                        }
                        else
                        {
                            // We are right of the column dropoff for this row
                            leftOfDropoff = 0;
                        }
                    }

                    // Reset subjectCount
                    subjectCount = parameters_semiGappedExtensionN;
                }

                subjectPosition++; bestRow++; insertQrow++; subjectCount--;
            }

            // -----SINGLE CELL RIGHT OF ROW DROPOFF -----
            if (!(bestScore > previousOldBest + dropoff) && (subjectPosition < subjectEnd))
            {
                // Calculate match value
                match = previousOldBest + matrixColumn[*subjectPosition];

                // Set value for best
                *bestRow = maximum(match, insertS);

                // Calculate new Iy
                insertS = maximum(match - parameters_semiGappedOpenGap,
                                  insertS - parameters_semiGappedExtendGap);

                // Set DUMMY values for Ix
                *insertQrow = constants_gappedExtensionDummyValue;

                subjectPosition++; bestRow++; insertQrow++;
            }

            // -----CELLS RIGHT OF ROW DROPOFF -----
            if (!(bestScore > *(bestRow - 1) + dropoff))
            {
                while (subjectPosition < subjectEnd)
                {
                    // Set value for Iy and best
                    *bestRow = insertS;
                    insertS = insertS - parameters_semiGappedExtendGap;

                    // Set DUMMY values for Ix
                    *insertQrow = constants_gappedExtensionDummyValue;

                    // If score at current cell is below dropoff
                    if (bestScore > *bestRow + dropoff)
                    {
                        // Stop processing row
                        subjectPosition++;
                        break;
                    }

                    subjectPosition++; bestRow++; insertQrow++;
                }
            }

            // Record dropoff position
            rowDropoff = subjectPosition - 1;

            // Clear insertS for next row
            insertS = constants_gappedExtensionDummyValue;

            // Reset queryCount
            queryCount = parameters_semiGappedExtensionN;
        }
//    if (dloc == 88197331)
//        print2(semiGappedScoring_bestRow_t[TID] + 1, subject, rowDropoff, columnDropoff);

        queryPosition++; queryCount--;
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}

struct dpResults semiGappedScoring_dpBeforeSeed_multi(unsigned char* subject, struct PSSMatrix PSSMatrix,
                                 struct coordinate seed, int4 dropoff, int TID)
{
    int2 **queryPosition, **bestQueryPosition;
    int2* matrixColumn;
    unsigned char *rowDropoff, *columnDropoff;
    unsigned char *subjectPosition, *bestSubjectPosition, *startSubjectPosition;
    int4 bestScore = 0;
    int4 *bestRow, *insertQrow, insertS, rowOffset;
    int4 subjectDistance;
    int4 oldBest, match, previousOldBest;
    unsigned char rightOfDropoff;
    int4 queryCount, subjectCount;
    struct dpResults dpResults;

    // Declare processing rows for storing match, insert-subject and insert-query values
    // If current malloced rows aren't big enough
    if (seed.subjectOffset >= semiGappedScoring_rowSizes_t[TID])
    {
        // Free existing rows
        if(semiGappedScoring_rowSizes_t[TID] > 0)
        {
            free(semiGappedScoring_bestRow_t[TID]);
            free(semiGappedScoring_insertQrow_t[TID]);
        }
        // Set size to double current needed length
        semiGappedScoring_rowSizes_t[TID] = (seed.subjectOffset) * 2;
        // Malloc new rows
        semiGappedScoring_bestRow_t[TID] = (int4*)malloc(sizeof(int4) * semiGappedScoring_rowSizes_t[TID]);
        semiGappedScoring_insertQrow_t[TID] = (int4*)malloc(sizeof(int4) * semiGappedScoring_rowSizes_t[TID]);
    }

    bestSubjectPosition = subjectPosition = startSubjectPosition = subject + seed.subjectOffset - 1;
    bestQueryPosition = queryPosition = PSSMatrix.matrix + seed.queryOffset - 1;

    // Initialize row pointers
    rowOffset = (subjectPosition - subject);
//    printf("rowOffset=%d Dloc=%d\n", rowOffset, dloc); fflush(stdout);
    bestRow = semiGappedScoring_bestRow_t[TID] + rowOffset;
    insertQrow = semiGappedScoring_insertQrow_t[TID] + rowOffset;

    // Set initial row dropoff and column dropoff
    rowDropoff = subject;
    columnDropoff = subject + seed.subjectOffset;

    // Using first column of query matrix
    matrixColumn = *queryPosition;

    // -----FIRST ROW-----

    // -----FIRST CELL-----
    // Set M value for bottom-right cell
    match = matrixColumn[*subjectPosition];

    // M must be the best
    *bestRow = match;

    // Only gap opens possible
    *insertQrow = insertS = match - parameters_semiGappedOpenGap;

    // If this is the best-yet scoring cell
    if (match > bestScore)
    {
        // Update best start cell data
        bestScore = match;
        bestQueryPosition = queryPosition;
        bestSubjectPosition = subjectPosition;
    }

    subjectDistance = 0;
    subjectPosition--; bestRow--; insertQrow--;

    // ----- REMAINING CELLS -----
    // For each remaining column in the bottom row, scanning from right-to-left
    while (subjectPosition >= subject)
    {
        // Set value for M
        match = matrixColumn[*subjectPosition]
              - parameters_semiGappedOpenGap - subjectDistance * parameters_semiGappedExtendGap;

        // Determine the best of M and Iy
        if (match > insertS)
        {
            *bestRow = match;

            // Calculate new Iy
            insertS = maximum(match - parameters_semiGappedOpenGap,
                              insertS - parameters_semiGappedExtendGap);
        }
        else
        {
            *bestRow = insertS;

            // Since M <= Iy, new Iy must derive from Iy
            insertS -= parameters_semiGappedExtendGap;
        }

        // Set DUMMY Ix value, which should never be used
        *insertQrow = constants_gappedExtensionDummyValue;

        // If this is the best-yet scoring cell
        if (match > bestScore)
        {
            // Update best start cell data
            bestScore = match;
            bestQueryPosition = queryPosition;
            bestSubjectPosition = subjectPosition;
        }

        // If score at current cell is below dropoff
        if (bestScore > *bestRow + dropoff)
        {
            // Record dropoff position
            rowDropoff = subjectPosition;
            // And stop processing row
            break;
        }

        subjectPosition--; bestRow--; insertQrow--;
        subjectDistance++;
    }

//    if (dloc == 746829265)
//    print(semiGappedScoring_bestRow_t[TID], subject, rowDropoff, columnDropoff);

    // Start queryCount at N. Only allow insertS for every Nth row when queryCount
    // reaches 0
    queryCount = parameters_semiGappedExtensionN;

    // -----REMAINING ROWS-----
    while (queryPosition > PSSMatrix.matrix && rowDropoff < columnDropoff)
    {
        queryPosition--;
        queryCount--;
        subjectPosition = columnDropoff - 1;

        // Determine subjectCount for initial subjectPosition. Is used to only allow
        // insertQ when subjectOffset % parameters_semiGappedExtensionN == 0
        subjectCount = (int4)(startSubjectPosition - subjectPosition) % parameters_semiGappedExtensionN;
        if (subjectCount)
            subjectCount = parameters_semiGappedExtensionN - subjectCount;

        // Reset row pointers to start of rows
        rowOffset = (subjectPosition - subject);
        bestRow = semiGappedScoring_bestRow_t[TID] + rowOffset;
        insertQrow = semiGappedScoring_insertQrow_t[TID] + rowOffset;

        // Using next column of query matrix
        matrixColumn = *queryPosition;

        // ************ All rows we are not allowing insertS
        if (queryCount)
        {
            // ** No insertQ allowed this column, this cell will only get a DUMMY score
            if (subjectCount)
            {
                previousOldBest = *bestRow;
                *bestRow = constants_gappedExtensionDummyValue;

                // Score at this cell is below dropoff
                columnDropoff = subjectPosition;
                rightOfDropoff = 1;
            }
            // ** We are allowing insertQ this column
            else
            {
                // -----FAR RIGHT CELL-----
                // Record some old values
                previousOldBest = *bestRow;

                // Set Ix value
                *bestRow = *insertQrow;
                *insertQrow -= parameters_semiGappedExtendGap;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
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

                // Reset subjectCount
                subjectCount = parameters_semiGappedExtensionN;
            }

            subjectPosition--; bestRow--; insertQrow--; subjectCount--;

            // -----CELLS RIGHT OF ROW DROPOFF-----
            while (subjectPosition >= rowDropoff)
            {
                // ** We are not allowing insertQ this column
                if (subjectCount)
                {
                    // Calculate new M value, which is also the best
                    oldBest = *bestRow;
                    match = *bestRow = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }
                }
                // We are allowing insertQ this column
                else
                {
                    // Calculate new M value
                    oldBest = *bestRow;
                    match = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // Determine the best of M and Ix
					if (match > *insertQrow)
                    {
						*bestRow = match;

                        // Calculate new Ix
                        *insertQrow = maximum(match - parameters_semiGappedOpenGap,
                                              *insertQrow - parameters_semiGappedExtendGap);
                    }
                    else
                    {
						*bestRow = *insertQrow;

                        // Since M <= Ix, new Ix must derive from Ix
						*insertQrow -= parameters_semiGappedExtendGap;
                    }

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }

                    // Reset subjectCount
                    subjectCount = parameters_semiGappedExtensionN;
                }

                subjectPosition--; bestRow--; insertQrow--; subjectCount--;
            }

            // -----SINGLE CELL LEFT OF ROW DROPOFF -----
            if (!(bestScore > previousOldBest + dropoff) && (subjectPosition >= subject))
            {
                // Set value for best
                *bestRow = match = previousOldBest + matrixColumn[*subjectPosition];

                // Set DUMMY values for Ix
                *insertQrow = constants_gappedExtensionDummyValue;

                if (match + dropoff >= bestScore)
                {
                    // Record dropoff position
                    rowDropoff = subjectPosition;
                }
            }
        }

        // ************ Every Nth row we allow insertS
        else
        {
            // -----FAR RIGHT CELL-----

            // ** No insertQ allowed this column, this cell will only get a DUMMY score
            if (subjectCount)
            {
                previousOldBest = *bestRow;
                *bestRow = constants_gappedExtensionDummyValue;

                // Score at this cell is below dropoff
                columnDropoff = subjectPosition;
                rightOfDropoff = 1;
            }
            // ** We are allowing insertQ this column
            else
            {
                // Record some old values
                previousOldBest = *bestRow;

                // Set Ix value
                *bestRow = *insertQrow;
                *insertQrow -= parameters_semiGappedExtendGap;

                // Set DUMMY value for Iy, which should never be used
                insertS = constants_gappedExtensionDummyValue;

                // If score at current cell is below dropoff
                if (bestScore > *bestRow + dropoff)
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

                // Reset subjectCount
                subjectCount = parameters_semiGappedExtensionN;
            }

            subjectPosition--; bestRow--; insertQrow--; subjectCount--;

            // -----CELLS RIGHT OF ROW DROPOFF-----
            while (subjectPosition >= rowDropoff)
            {
                // ** We are not allowing insertQ this column
                if (subjectCount)
                {
                    // Remember old M value (for cell below this one)
                    oldBest = *bestRow;
                    match = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // Determine the best of M and Iy
					if (match > insertS)
                    {
						*bestRow = match;

                        // Calculate new Iy
                        insertS = maximum(match - parameters_semiGappedOpenGap,
                                          insertS - parameters_semiGappedExtendGap);
                    }
                    else
                    {
						*bestRow = insertS;

                        // Since M <= Iy, new Iy must derive from Iy
						insertS -= parameters_semiGappedExtendGap;
                    }

                    // If this is the best-yet scoring cell
                    if (match > bestScore)
                    {
                        // Update best start cell data
                        bestScore = match;
                        bestQueryPosition = queryPosition;
                        bestSubjectPosition = subjectPosition;
                    }

                    // If score at current cell (and cells to its right) are below dropoff
                    if (rightOfDropoff)
                    {
                        if (bestScore > *bestRow + dropoff)
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
                }
                // ** We are allowing insertQ this column
                else
                {
                    // Remember old M value (for cell below this one)
                    oldBest = *bestRow;
                    match = matrixColumn[*subjectPosition] + previousOldBest;
                    previousOldBest = oldBest;

                    // Determine the best of M, Ix and Iy
                    if (match > insertS)
                    {
                        if (match > *insertQrow)
                        {
                            // Match is largest
                            *bestRow = match;

                            // Calculate new Ix
                            *insertQrow = maximum(match - parameters_semiGappedOpenGap,
                                                  *insertQrow - parameters_semiGappedExtendGap);

                            // Calculate new Iy
                            insertS = maximum(match - parameters_semiGappedOpenGap,
                                              insertS - parameters_semiGappedExtendGap);

                            // If this is the best-yet scoring cell
                            if (match > bestScore)
                            {
                                // Update best start cell data
                                bestScore = match;
                                bestQueryPosition = queryPosition;
                                bestSubjectPosition = subjectPosition;
                            }
                        }
                        else
                        {
                            // insertQ is largest
                            *bestRow = *insertQrow;

                            // Calculate new Ix
                            *insertQrow -= parameters_semiGappedExtendGap;

							// Dummy Iy
                            insertS = constants_gappedExtensionDummyValue;
                        }
                    }
                    else
                    {
                        if (insertS > *insertQrow)
                        {
                            // insertS is largest
                            *bestRow = insertS;

                            // Dummy Ix
                            *insertQrow = constants_gappedExtensionDummyValue;

                            // Calculate new Iy
                            insertS -= parameters_semiGappedExtendGap;

                        }
                        else
                        {
                            // insertQ is largest
                            *bestRow = *insertQrow;

                            // Calculate new Ix
                            *insertQrow -= parameters_semiGappedExtendGap;

							// Dummy Iy
                            insertS = constants_gappedExtensionDummyValue;
                        }
                    }

                    // If score at current cell (and cells to its right) are below dropoff
                    if (rightOfDropoff)
                    {
                        if (bestScore > *bestRow + dropoff)
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

                    // Reset subjectCount
                    subjectCount = parameters_semiGappedExtensionN;
                }

                subjectPosition--; bestRow--; insertQrow--; subjectCount--;
            }

            // -----SINGLE CELL LEFT OF ROW DROPOFF -----
            if (!(bestScore > previousOldBest + dropoff) && (subjectPosition >= subject))
            {
                // Calculate match value
                match = previousOldBest + matrixColumn[*subjectPosition];

                // Set value for best
                *bestRow = maximum(match, insertS);

                // Calculate new Iy
                insertS = maximum(match - parameters_semiGappedOpenGap,
                                  insertS - parameters_semiGappedExtendGap);
 
                // Set DUMMY values for Ix
                *insertQrow = constants_gappedExtensionDummyValue;

                subjectPosition--; bestRow--; insertQrow--;
            }

            // -----CELLS LEFT OF ROW DROPOFF -----
            if (!(bestScore > *(bestRow + 1) + dropoff))
            {
                while (subjectPosition >= subject)
                {
                    // Set value for Iy and best
                    *bestRow = insertS;
                    insertS = insertS - parameters_semiGappedExtendGap;

                    // Set DUMMY values for Ix
                    *insertQrow = constants_gappedExtensionDummyValue;

                    // If score at current cell is below dropoff
                    if (bestScore > *bestRow + dropoff)
                    {
                        // Stop processing row
                        subjectPosition--;
                        break;
                    }

                    subjectPosition--; bestRow--; insertQrow--;
                }
            }

            // Record dropoff position
            rowDropoff = subjectPosition + 1;

            // Clear insertS for next row
            insertS = constants_gappedExtensionDummyValue;

            // Reset queryCount
            queryCount = parameters_semiGappedExtensionN;
        }
//        if (dloc == 746829265)
//                print(semiGappedScoring_bestRow_t[TID], subject, rowDropoff, columnDropoff);
    }

    dpResults.best.queryOffset = bestQueryPosition - PSSMatrix.matrix;
    dpResults.best.subjectOffset = bestSubjectPosition - subject;
    dpResults.bestScore = bestScore;
    dpResults.traceback = NULL;
    return dpResults;
}

int4 semiGappedScoring_score_multi(struct ungappedExtension* ungappedExtension, struct PSSMatrix PSSMatrix,
                        int4 subjectSize, unsigned char* subject, int4 dropoff, int TID)
{
    struct coordinate seed;
    unsigned char *choppedSubject;
    struct PSSMatrix choppedPSSMatrix;
    int4 choppedSubjectSize;
    struct dpResults beforeDpResults, afterDpResults;
    int4 strandOffset = 0;

    // Increment counter number of HSP's semi-gapped
    //blast_numSemiGapped++;
    
    __sync_fetch_and_add(&blast_numSemiGapped, 1);

    // Use old code if restricted insertion has been disabled
    if (!parameters_restrictedInsertionScoring)
    {
        printf("multithreading semiGappedScoring_score_pthread does not support oldSemiGappedScoring_score!\n");
        exit(0);
        //return oldSemiGappedScoring_score(ungappedExtension, PSSMatrix, subjectSize, subject, dropoff);
    }

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

    beforeDpResults = semiGappedScoring_dpBeforeSeed_multi(subject, PSSMatrix,
                    seed, dropoff + parameters_semiGappedDropoffIncrease, TID);

    // Chop the start off the query and subject so they begin at the seed
    choppedPSSMatrix = PSSMatrix_chop(PSSMatrix, seed.queryOffset);
    choppedSubject = subject + seed.subjectOffset;
    choppedSubjectSize = subjectSize - seed.subjectOffset;

    // Perform dynamic programming for points after the seed
    afterDpResults = semiGappedScoring_dpAfterSeed_multi(choppedSubject, choppedPSSMatrix,
                    dropoff + parameters_semiGappedDropoffIncrease, choppedSubjectSize, TID);

    // Re-adjust result change due to chopping subject/query and strand adjustment
    afterDpResults.best.queryOffset += seed.queryOffset + strandOffset;
    afterDpResults.best.subjectOffset += seed.subjectOffset;
    beforeDpResults.best.queryOffset += strandOffset;

    // Associate best scoring start and end points with the ungapped extension
    ungappedExtension->start = beforeDpResults.best;
    ungappedExtension->end = afterDpResults.best;

//    if (dloc == 88197331)
//        printf("semiGapped[%d,%d,%d]\n", beforeDpResults.bestScore, afterDpResults.bestScore,
//        choppedPSSMatrix.matrix[0][choppedSubject[0]]);

    // Determine score by combining score from the two traces, and the match score at
    // the seed position
    return beforeDpResults.bestScore + afterDpResults.bestScore +
           choppedPSSMatrix.matrix[0][choppedSubject[0]];
}





void semiGappedScoring_free_multi(int n_threads)
{
    int ii = 0;
    for(ii = 0; ii < n_threads; ii++)
    {
        free(semiGappedScoring_bestRow_t[ii]);
        free(semiGappedScoring_insertQrow_t[ii]);
    }
}
