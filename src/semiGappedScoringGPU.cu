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
 * Developer: Shucai Xiao 
 *
 */



//#define maximum(a,b) (((a) > (b)) ? (a) : (b))
#define TARGET 0

__device__ __constant__ int2 scoreMatrixC[1640];
__device__ __constant__ unsigned char querySequenceC[40000];

struct gappedExtensionParameters {
	int2 semiGappedOpenGap;
	int2 semiGappedExtendGap;
	int4 semiGappedExtensionN;
	int4 semiGappedDropoffIncrease;
	unsigned char encoding_numCodes;
};

__device__ struct dpResults semiGappedScoring_dpBeforeSeedGPU(unsigned char* subject, struct PSSMatrixFP PSSMatrixFP,
                                 struct coordinate seed, int4 dropoff, int4 *bestRow, int4 *insertQrow, 
								 unsigned char encoding_numCodes, int4 openGapPenalty, int4 extensionGapPenalty,
								 int4 semiGappedExtensionN);

__device__ struct dpResults semiGappedScoring_dpAfterSeedGPU(unsigned char* subject, struct PSSMatrixFP PSSMatrixFP,
                              int4 dropoff, int4 subjectLength, int4 *bestRow, int4 *insertQrow,
							  unsigned char encoding_numCodes, int4 openGapPenalty, int4 extensionGapPenalty,
							  int4 semiGappedExtensionN, int4 queryOffset);

__device__  int4 semiGappedScoring_scoreGPU(struct ungappedExtension *ungappedExtension,
									  struct PSSMatrixFP PSSMatrix,
									  int4	 subjectSize,
									  unsigned char *subject,
									  int4	 dropoff,
									  int4	 semiGappedDropoffIncrease,
									  int4	 *bestRow,
									  int4	 *insertQRow,
									  unsigned char encoding_numCodes,
									  int4	 openGapPenalty,
									  int4	 extensionGapPenalty,
									  int4	 semiGappedExtensionN);

__device__ struct PSSMatrixFP PSSMatrixFP_chop(struct PSSMatrixFP PSSMatrixFP, int4 amount, unsigned char encoding_numCodes)
{
	struct PSSMatrixFP chopped;
	chopped.matrix = PSSMatrixFP.matrix + amount * encoding_numCodes;
	chopped.queryCodes = PSSMatrixFP.queryCodes + amount;
	chopped.bytePackedCodes = PSSMatrixFP.bytePackedCodes + amount;
	chopped.xorCodes = PSSMatrixFP.xorCodes + amount;
	chopped.length = PSSMatrixFP.length - amount;
	chopped.highestValue = PSSMatrixFP.highestValue;
	chopped.lowestValue = PSSMatrixFP.lowestValue;
	chopped.strandLength = PSSMatrixFP.strandLength - amount;

	if (chopped.strandLength < 0)
		chopped.strandLength = 0;
	
	return chopped;
}

__device__ struct coordinate ungappedExtension_findProteinSeed1(
					struct ungappedExtension* ungappedExtension,
					struct PSSMatrixFP PSSMatrixFP, 
					unsigned char* subject, 
					unsigned char encoding_numCodes)
{
	char *queryWindowStart, *queryWindowEnd;
	unsigned char *subjectWindowStart, *subjectWindowEnd;

	char* bestQueryPosition;
	unsigned char* bestSubjectPosition;
	int4 bestSegmentScore;
	int4 nominalScore, count;
	struct coordinate seed;

	if (ungappedExtension->end.queryOffset - ungappedExtension->start.queryOffset < 11)
	{
		// The seed point is the middle of the extension
		seed.queryOffset = (ungappedExtension->end.queryOffset +
							ungappedExtension->start.queryOffset) / 2;
		seed.subjectOffset = (ungappedExtension->end.subjectOffset +
							  ungappedExtension->start.subjectOffset) / 2;
	}
	else
	{
		// Else find the highest scoring length-11 segment of the ungapped extension
		queryWindowStart = queryWindowEnd = PSSMatrixFP.matrix + ungappedExtension->start.queryOffset * encoding_numCodes;
		subjectWindowStart = subjectWindowEnd = subject + ungappedExtension->start.subjectOffset;

		// Find initial score for first 11 positions
		nominalScore = 0;
		count = 0;
		while (count < 11)
		{
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
		while (queryWindowEnd < PSSMatrixFP.matrix + ungappedExtension->end.queryOffset * encoding_numCodes)
		{
			// Advance window end, add new position value
			queryWindowEnd += encoding_numCodes;
			subjectWindowEnd++;

			nominalScore += queryWindowEnd[*subjectWindowEnd];
			// Remove position that we will leave behind
			nominalScore -= queryWindowStart[*subjectWindowStart];

			queryWindowStart += encoding_numCodes;
			subjectWindowStart++;

			// Check if best window position yet
			if (nominalScore > bestSegmentScore)
			{
				bestSegmentScore = nominalScore;
				bestQueryPosition = queryWindowStart;
				bestSubjectPosition = subjectWindowStart;
			}
		}

		// Middle of the best window is the seed position
		seed.queryOffset = (bestQueryPosition - PSSMatrixFP.matrix) / encoding_numCodes + 5;
		seed.subjectOffset = bestSubjectPosition + 5 - subject;
	}

	return seed;
}

__device__ void alignments_pruneRegion(struct ungappedExtension *ungappedExtension, 
									   struct ungappedExtension *curExtension,
									   int ungappedExtensionNum)
{
	int i;
	for (i = 0; i < ungappedExtensionNum; i++)
	{
		if (ungappedExtension[i].status != ungappedExtension_DELETED)
		{
			if (ungappedExtension[i].start.queryOffset >= curExtension->start.queryOffset &&
				ungappedExtension[i].end.queryOffset <= curExtension->end.queryOffset &&
				ungappedExtension[i].start.subjectOffset >= curExtension->start.subjectOffset &&
				ungappedExtension[i].end.subjectOffset <= curExtension->end.subjectOffset &&
				ungappedExtension[i].nominalScore <= curExtension->nominalScore &&
				ungappedExtension + i != curExtension)
			{
				ungappedExtension[i].status = ungappedExtension_DELETED;
			}
		}
	}

	return;
}

// Perform semi-gapped alignment with restricted insertion
__global__ void semiGappedScoring_kernel(struct sequenceDataFP 	*sequenceDataFP, 
										 unsigned char			*sequences,
										 struct PSSMatrixFP		*PSSMatrix,
										 char					*PSSMatrixBody,
										 struct gappedExtensionParameters *parameters,
										 int					*startLocArray,
										 int					*ungappedExtensionNumArray,
										 int					alignmentNum,
										 struct ungappedExtension *ungappedExtensions,
										 int4					*bestScores,
										 int4					*numGoodExtensions,
										 int4					*numSemiGapping,
										 int4					*orderArray,
										 int4					*bestRowAll,
										 int4					*insertQRowAll,
                        				 int4 					dropoff,
										 int4					nominalR1cutoff)
{
	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= alignmentNum)
	{
		return;
	}

	int2 openGapPenalty = parameters->semiGappedOpenGap;
	int2 extensionGapPenalty = parameters->semiGappedExtendGap;
	int4 semiGappedExtensionN = parameters->semiGappedExtensionN;
	int4 semiGappedDropoffIncrease = parameters->semiGappedDropoffIncrease;
	unsigned char encoding_numCodes = parameters->encoding_numCodes;
	PSSMatrix->matrix = PSSMatrixBody + encoding_numCodes;
	
	//int4 ungappedStartLoc = startLocArray[orderArray[tid]];
	int4 ungappedStartLoc = startLocArray[tid];
	int4 sequenceCount;
	struct ungappedExtension *ungappedExtension;

	sequenceCount = ungappedExtensions[ungappedStartLoc].sequenceCount;
	int4		   subjectOffset = sequenceDataFP[sequenceCount].offset;
	unsigned char *subject = sequences + subjectOffset;
	int4		  *bestRow = bestRowAll + subjectOffset;
	int4		  *insertQRow = insertQRowAll + subjectOffset;
	uint4 subjectLength = sequenceDataFP[sequenceCount].sequenceLength;
	int4 bestScore = 0;
	int numSemiGappingDevice = 0;
	uint4 goodExtensionNo, ungappedExtensionNo;

	goodExtensionNo = 0;
	//for (ungappedExtensionNo = 0; ungappedExtensionNo < ungappedExtensionNumArray[orderArray[tid]]; ungappedExtensionNo++)
	for (ungappedExtensionNo = 0; ungappedExtensionNo < ungappedExtensionNumArray[tid]; ungappedExtensionNo++)
	{
		ungappedExtension = &ungappedExtensions[ungappedStartLoc + ungappedExtensionNo];
		if (ungappedExtension->status != ungappedExtension_DELETED)
		{
			if (ungappedExtension->seed.queryOffset == -1 &&
				ungappedExtension->seed.subjectOffset == -1)
			{
				ungappedExtension->seed = ungappedExtension_findProteinSeed1(ungappedExtension, *PSSMatrix, subject, encoding_numCodes);
			}

			numSemiGappingDevice++;
			ungappedExtension->nominalScore = semiGappedScoring_scoreGPU(ungappedExtension, *PSSMatrix, 
											  subjectLength, subject, dropoff, semiGappedDropoffIncrease, 
											  bestRow, insertQRow, encoding_numCodes, openGapPenalty, 
											  extensionGapPenalty, semiGappedExtensionN);
			ungappedExtension->status = ungappedExtension_SEMIGAPPED;

			if (ungappedExtension->nominalScore >= nominalR1cutoff)
			{
				if (ungappedExtension->nominalScore > bestScore)
				{
					bestScore = ungappedExtension->nominalScore;
				}
			}
			else
			{
				ungappedExtension->status = ungappedExtension_DELETED;
			}
			goodExtensionNo++;
			//alignments_pruneRegion(&ungappedExtensions[ungappedStartLoc], ungappedExtension, ungappedExtensionNumArray[orderArray[tid]]);
			alignments_pruneRegion(&ungappedExtensions[ungappedStartLoc], ungappedExtension, ungappedExtensionNumArray[tid]);
		}
	}

	numSemiGapping[tid] = numSemiGappingDevice;
	//bestScores[orderArray[tid]] = bestScore;
	bestScores[tid] = bestScore;
	if (bestScore >= nominalR1cutoff)
	{
		numGoodExtensions[tid] += goodExtensionNo;
	}

	return;
}

__device__ int4 semiGappedScoring_scoreGPU(struct ungappedExtension *ungappedExtension,
									  struct PSSMatrixFP PSSMatrix,
									  int4	 subjectSize,
									  unsigned char *subject,
									  int4	 dropoff,
									  int4	 semiGappedDropoffIncrease,
									  int4	 *bestRow,
									  int4	 *insertQRow,
									  unsigned char encoding_numCodes,
									  int4	 openGapPenalty,
									  int4	 extensionGapPenalty,
									  int4	 semiGappedExtensionN)
{
	struct coordinate seed;
	unsigned char *choppedSubject;
	struct PSSMatrixFP choppedPSSMatrix;
	int4 choppedSubjectSize;
	struct dpResults beforeDpResults, afterDpResults;
	int4 strandOffset = 0;

	// Perform dynamic programming for points before the seed
	seed = ungappedExtension->seed;
	if (seed.queryOffset > PSSMatrix.strandLength)
	{
		// If query position is in the second strand, remove first strand from PSSM
		strandOffset = PSSMatrix.strandLength;
		seed.queryOffset -= PSSMatrix.strandLength;
		PSSMatrix = PSSMatrixFP_chop(PSSMatrix, PSSMatrix.strandLength, encoding_numCodes);
	}
	else
	{
		// Otherwise remove second strand
		PSSMatrix.length = PSSMatrix.strandLength;
	}

	beforeDpResults = semiGappedScoring_dpBeforeSeedGPU(subject, PSSMatrix,
					seed, dropoff + semiGappedDropoffIncrease, bestRow, 
					insertQRow, encoding_numCodes, openGapPenalty, 
					extensionGapPenalty, semiGappedExtensionN);

	// Chop the start off the query and subject so they begin at the seed
	choppedPSSMatrix = PSSMatrixFP_chop(PSSMatrix, seed.queryOffset, encoding_numCodes);
	choppedSubject = subject + seed.subjectOffset;
	choppedSubjectSize = subjectSize - seed.subjectOffset;

	// Perform dynamic programming for points after the seed
	afterDpResults = semiGappedScoring_dpAfterSeedGPU(choppedSubject, choppedPSSMatrix,
					dropoff + semiGappedDropoffIncrease, choppedSubjectSize, bestRow, 
					insertQRow, encoding_numCodes, openGapPenalty, 
					extensionGapPenalty, semiGappedExtensionN, seed.queryOffset);

	// Re-adjust result change due to chopping subject/query and strand adjustment
	afterDpResults.best.queryOffset += seed.queryOffset + strandOffset;
	afterDpResults.best.subjectOffset += seed.subjectOffset;
	beforeDpResults.best.queryOffset += strandOffset;

	// Associate best scoring start and end points with the ungapped extension
	ungappedExtension->start = beforeDpResults.best;
	ungappedExtension->end = afterDpResults.best;

    // Determine score by combining score from the two traces, and the match score at
    // the seed position
    return beforeDpResults.bestScore + afterDpResults.bestScore +
           choppedPSSMatrix.matrix[choppedSubject[0]];
}

// Perform dynamic programming to explore possible start points and alignments that end at
// the given seed and find the best score
__device__ struct dpResults semiGappedScoring_dpBeforeSeedGPU(unsigned char* subject, struct PSSMatrixFP PSSMatrixFP,
                                 struct coordinate seed, int4 dropoff, int4 *bestRowCur, int4 *insertQrowCur, 
								 unsigned char encoding_numCodes, int4 openGapPenalty, int4 extensionGapPenalty,
								 int4 semiGappedExtensionN)
{
	int4	queryPosition, bestQueryPosition;
	int4	matrixColumn;
	unsigned char *rowDropoff, *columnDropoff;
	
	int4 *bestRow, *insertQrow;
	unsigned char *subjectPosition, *bestSubjectPosition, *startSubjectPosition;
	int4 bestScore = 0;
	int4 insertS, rowOffset;
	int4 subjectDistance;
	int4 oldBest, match, previousOldBest;
	unsigned char rightOfDropoff;
	int4 queryCount, subjectCount;
	struct dpResults dpResults;

	bestSubjectPosition = subjectPosition = startSubjectPosition = subject + seed.subjectOffset - 1;
	bestQueryPosition = queryPosition = seed.queryOffset - 1;

	// Initialize row pointers
	rowOffset = (subjectPosition - subject);

	bestRow = bestRowCur + rowOffset;
	insertQrow = insertQrowCur + rowOffset;

	// Set initial row dropoff and column dropoff
	rowDropoff = subject;
	columnDropoff = subject + seed.subjectOffset;

	// Using first column of query matrix
	matrixColumn = queryPosition + 1;

	// -----FIRST ROW-----

	// -----FIRST CELL-----
	// Set M value for bottom-right cell
	match = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)];

	// M must be the best
	*bestRow = match;

	// Only gap opens possible
	//*insertQrow = insertS = match - parameters_semiGappedOpenGap;
	*insertQrow = insertS = match - openGapPenalty;

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
		match = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)]
			  - openGapPenalty - subjectDistance * extensionGapPenalty;

		// Determine the best of M and Iy
		if (match > insertS)
		{
			*bestRow = match;

			// Calculate new Iy
			insertS = maximum(match - openGapPenalty,
							  insertS - extensionGapPenalty);
		}
		else
		{
			*bestRow = insertS;

			// Since M <= Iy, new Iy must derive from Iy
			insertS -= extensionGapPenalty;
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

	// Start queryCount at N. Only allow insertS for every Nth row when queryCount
	// reaches 0
	queryCount = semiGappedExtensionN;

	// -----REMAINING ROWS-----
	//while (queryPosition > PSSMatrixFP.matrix && rowDropoff < columnDropoff)
	while (queryPosition > 0 && rowDropoff < columnDropoff)
	{
		queryPosition--;
		queryCount--;
		subjectPosition = columnDropoff - 1;

		// Determine subjectCount for initial subjectPosition. Is used to only allow
		// insertQ when subjectOffset % parameters_semiGappedExtensionN == 0
		subjectCount = (int4)(startSubjectPosition - subjectPosition) % semiGappedExtensionN;
		if (subjectCount)
			subjectCount = semiGappedExtensionN - subjectCount;

		// Reset row pointers to start of rows
		rowOffset = (subjectPosition - subject);
		bestRow = bestRowCur + rowOffset;
		insertQrow = insertQrowCur + rowOffset;

		// Using next column of query matrix
		matrixColumn = queryPosition + 1;

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
				*insertQrow -= extensionGapPenalty;

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
				subjectCount = semiGappedExtensionN;
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
					match = *bestRow = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)] + previousOldBest;
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
					match = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)] + previousOldBest;
					previousOldBest = oldBest;

					// Determine the best of M and Ix
					if (match > *insertQrow)
					{
						*bestRow = match;

						// Calculate new Ix
						*insertQrow = maximum(match - openGapPenalty,
											  *insertQrow - extensionGapPenalty);
					}
					else
					{
						*bestRow = *insertQrow;

						// Since M <= Ix, new Ix must derive from Ix
						*insertQrow -= extensionGapPenalty;
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
					subjectCount = semiGappedExtensionN;
				}

				subjectPosition--; bestRow--; insertQrow--; subjectCount--;
			}

			// -----SINGLE CELL LEFT OF ROW DROPOFF -----
			if (!(bestScore > previousOldBest + dropoff) && (subjectPosition >= subject))
			{
				// Set value for best
				*bestRow = match = previousOldBest + scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)];

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
				*insertQrow -= extensionGapPenalty;

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
				subjectCount = semiGappedExtensionN;
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
					match = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)] + previousOldBest;
					previousOldBest = oldBest;

					// Determine the best of M and Iy
					if (match > insertS)
					{
						*bestRow = match;

						// Calculate new Iy
						insertS = maximum(match - openGapPenalty,
										  insertS - extensionGapPenalty);
					}
					else
					{
						*bestRow = insertS;

						// Since M <= Iy, new Iy must derive from Iy
						insertS -= extensionGapPenalty;
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
					match = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)] + previousOldBest;
					previousOldBest = oldBest;

					// Determine the best of M, Ix and Iy
					if (match > insertS)
					{
						if (match > *insertQrow)
						{
							// Match is largest
							*bestRow = match;

							// Calculate new Ix
							*insertQrow = maximum(match - openGapPenalty,
												  *insertQrow - extensionGapPenalty);

							// Calculate new Iy
							insertS = maximum(match - openGapPenalty,
											  insertS - extensionGapPenalty);

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
							*insertQrow -= extensionGapPenalty;

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
							insertS -= extensionGapPenalty;

						}
						else
						{
							// insertQ is largest
							*bestRow = *insertQrow;

							// Calculate new Ix
							*insertQrow -= extensionGapPenalty;

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
					subjectCount = semiGappedExtensionN;
				}

				subjectPosition--; bestRow--; insertQrow--; subjectCount--;
			}

			// -----SINGLE CELL LEFT OF ROW DROPOFF -----
			if (!(bestScore > previousOldBest + dropoff) && (subjectPosition >= subject))
			{
				// Calculate match value
				match = previousOldBest + scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)];

				// Set value for best
				*bestRow = maximum(match, insertS);

				// Calculate new Iy
				insertS = maximum(match - openGapPenalty,
								  insertS - extensionGapPenalty);

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
					insertS = insertS - extensionGapPenalty;

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
			queryCount = semiGappedExtensionN;
		}
	}

	dpResults.best.queryOffset = bestQueryPosition;
	dpResults.best.subjectOffset = bestSubjectPosition - subject;
	dpResults.bestScore = bestScore;
	dpResults.traceback = NULL;
	return dpResults;
}

// Perform dynamic programming to explore possible END points and alignments that start at
// the given seed and find the best score
__device__ struct dpResults semiGappedScoring_dpAfterSeedGPU(unsigned char* subject, struct PSSMatrixFP PSSMatrixFP,
                              int4 dropoff, int4 subjectLength, int4 *bestRowCur, int4 *insertQrowCur,
							  unsigned char encoding_numCodes, int4 openGapPenalty, int4 extensionGapPenalty,
							  int4 semiGappedExtensionN, int4 queryOffset)
{
	int4 queryPosition, bestQueryPosition, queryEnd;
	int4 matrixColumn;

	unsigned char *rowDropoff, *columnDropoff;
	unsigned char *subjectPosition, *bestSubjectPosition, *subjectEnd, *startSubjectPosition;
	int4 *bestRow, *insertQrow;
	int4 bestScore = 0;
	int4 insertS, rowOffset;
	int4 subjectDistance;
	int4 oldBest, match, previousOldBest;
	unsigned char leftOfDropoff;
	int4 queryLength;
	int4 queryCount, subjectCount;
	struct dpResults dpResults;

	queryLength = PSSMatrixFP.length;
	subjectEnd = subject + subjectLength;
	queryEnd = queryLength;

	bestSubjectPosition = subjectPosition = startSubjectPosition = subject + 1;
	bestQueryPosition = queryPosition = 1;

	// Initialize rows
	bestRow = bestRowCur + 1;
	insertQrow = insertQrowCur + 1;

	// Set initial row dropoff and column dropoff
	rowDropoff = subject + subjectLength - 1;
	columnDropoff = subject;

	// -----FIRST ROW-----

	// Using first column of query matrix
	matrixColumn = queryPosition + 1 + queryOffset;

	// -----FIRST CELL-----
	// Set M value for top-left cell
	match = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)];

	// M must be the best
	*bestRow = match;

	// Only gap opens possible
	*insertQrow = insertS = match - openGapPenalty;

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
		match = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)]
			  - openGapPenalty - subjectDistance * extensionGapPenalty;

		// Determine the best of M and Iy
		if (match > insertS)
		{
			*bestRow = match;

			// Calculate new Iy
			insertS = maximum(match - openGapPenalty,
							  insertS - extensionGapPenalty);
		}
		else
		{
			*bestRow = insertS;

			// Since M <= Iy, new Iy must derive from Iy
			insertS -= extensionGapPenalty;
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

	// Start queryCount at N. Only allow insertS for every Nth row when queryCount
	// reaches 0
	queryCount = semiGappedExtensionN;

	queryPosition += 1; 
	queryCount--;

	// -----REMAINING ROWS-----
	while (queryPosition < queryEnd && rowDropoff > columnDropoff)
	{
		subjectPosition = columnDropoff + 1;

		// Determine subjectCount for initial subjectPosition. Is used to only allow
		// insertQ when subjectOffset % parameters_semiGappedExtensionN == 0
		subjectCount = ((int4)(subjectPosition - startSubjectPosition) % semiGappedExtensionN);
		if (subjectCount)
			subjectCount = semiGappedExtensionN - subjectCount;

		// Reset rows
		rowOffset = (subjectPosition - subject);
		bestRow = bestRowCur + rowOffset;
		insertQrow = insertQrowCur + rowOffset;

		// Using next column of query matrix
		matrixColumn = queryPosition + 1 + queryOffset;


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
				*insertQrow -= extensionGapPenalty;

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
				subjectCount = semiGappedExtensionN;
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
					match = *bestRow = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)] + previousOldBest;
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
					match = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)] + previousOldBest;
					previousOldBest = oldBest;

					// Determine the best of M and Ix
					if (match > *insertQrow)
					{
						*bestRow = match;

						// Calculate new Ix
						*insertQrow = maximum(match - openGapPenalty,
											  *insertQrow - extensionGapPenalty);
					}
					else
					{
						*bestRow = *insertQrow;

						// Since M <= Ix, new Ix must derive from Ix
						*insertQrow -= extensionGapPenalty;
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
					subjectCount = semiGappedExtensionN;
				}

				subjectPosition++; bestRow++; insertQrow++; subjectCount--;
			}

			// -----SINGLE CELL RIGHT OF ROW DROPOFF -----
			if (!(bestScore > previousOldBest + dropoff) && (subjectPosition < subjectEnd))
			{
				// Set value for best
				*bestRow = match = previousOldBest + scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)];

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
				*insertQrow -= extensionGapPenalty;

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
				subjectCount = semiGappedExtensionN;
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
					match = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)] + previousOldBest;
					previousOldBest = oldBest;

					// Determine the best of M and Iy
					if (match > insertS)
					{
						*bestRow = match;

						// Calculate new Iy
						insertS = maximum(match - openGapPenalty,
										  insertS - extensionGapPenalty);
					}
					else
					{
						*bestRow = insertS;

						// Since M <= Iy, new Iy must derive from Iy
						insertS -= extensionGapPenalty;
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
					match = scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)] + previousOldBest;
					previousOldBest = oldBest;

					// Determine the best of M, Ix and Iy
					if (match > insertS)
					{
						if (match > *insertQrow)
						{
							// Match is largest
							*bestRow = match;

							// Calculate new Ix
							*insertQrow = maximum(match - openGapPenalty,
												  *insertQrow - extensionGapPenalty);

							// Calculate new Iy
							insertS = maximum(match - openGapPenalty,
											  insertS - extensionGapPenalty);

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
							*insertQrow -= extensionGapPenalty;

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
							insertS -= extensionGapPenalty;

						}
						else
						{
							// insertQ is largest
							*bestRow = *insertQrow;

							// Calculate new Ix
							*insertQrow -= extensionGapPenalty;

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
					subjectCount = semiGappedExtensionN;
				}

				subjectPosition++; bestRow++; insertQrow++; subjectCount--;
			}

			// -----SINGLE CELL RIGHT OF ROW DROPOFF -----
			if (!(bestScore > previousOldBest + dropoff) && (subjectPosition < subjectEnd))
			{
				// Calculate match value
				match = previousOldBest + scoreMatrixC[querySequenceC[matrixColumn] * encoding_numCodes + (*subjectPosition)];

				// Set value for best
				*bestRow = maximum(match, insertS);

				// Calculate new Iy
				insertS = maximum(match - openGapPenalty,
								  insertS - extensionGapPenalty);

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
					insertS = insertS - extensionGapPenalty;

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
			queryCount = semiGappedExtensionN;
		}

		//queryPosition += encoding_numCodes; queryCount--;
		queryPosition += 1; queryCount--;
	}

	dpResults.best.queryOffset = bestQueryPosition;
	dpResults.best.subjectOffset = bestSubjectPosition - subject;
	dpResults.bestScore = bestScore;
	dpResults.traceback = NULL;
	return dpResults;
}


