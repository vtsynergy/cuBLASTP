INC=include
SRC=src
OBJ=obj

MGPU_PATH=/home/zhang/moderngpu-1.1/
CUDA_PATH=/usr/local/cuda-6.0/

CC=gcc -g
CCPLUS=g++ -g 
NVCC=nvcc -g -gencode arch=compute_35,code=sm_35 
CFLAGS= -I$(CUDA_PATH)/include -I$(INC) -I$(MGPU_PATH)/include -O3 -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
NVCCFLAGS= -I$(CUDA_PATH)/include -I$(MGPU_PATH)/include -O3 -I$(INC) -Isrc -D_FILE_OFFSET_BITS=64 -D_LARGEFILE64_SOURCE
LDFLAGS= -lpthread -lgomp -lm

OBJFILES= $(OBJ)/alignments_multi.o  $(OBJ)/alignments.o $(OBJ)/mgpucontextcu.o $(OBJ)/mgpuutil.o $(OBJ)/bytepackGappedScoring.o $(OBJ)/descriptions.o $(OBJ)/encoding.o $(OBJ)/fasterGappedExtension_multi.o $(OBJ)/fasterGappedExtension.o $(OBJ)/gappedExtension.o $(OBJ)/gappedScoring.o $(OBJ)/nuGappedScoring.o $(OBJ)/global.o $(OBJ)/hitMatrix.o $(OBJ)/karlin.o $(OBJ)/memBlocks.o $(OBJ)/memSingleBlock.o $(OBJ)/nucleotideLookup.o $(OBJ)/oldGappedScoring.o $(OBJ)/oldSemiGappedScoring.o $(OBJ)/parameters.o $(OBJ)/print.o $(OBJ)/PSSMatrix.o $(OBJ)/qPosList.o $(OBJ)/readFasta.o $(OBJ)/readFile.o $(OBJ)/scoreMatrix.o $(OBJ)/semiGappedScoring_multi.o $(OBJ)/semiGappedScoring.o $(OBJ)/statistics.o $(OBJ)/ungappedExtension.o $(OBJ)/wordLookupDFA.o $(OBJ)/writeFile.o $(OBJ)/constants.o $(OBJ)/smithWatermanTraceback.o $(OBJ)/smithWatermanScoring.o $(OBJ)/tableGappedScoring.o $(OBJ)/vbyte.o $(OBJ)/unpack.o $(OBJ)/index.o $(OBJ)/postings.o $(OBJ)/hashcounter.o $(OBJ)/writedb.o $(OBJ)/readdb.o $(OBJ)/searchcu.o $(OBJ)/search.o $(OBJ)/wildcards.o $(OBJ)/dust.o $(OBJ)/seg.o $(OBJ)/timeRec.o $(OBJ)/segsortcu.o

HEADERFILES= $(INC)/segsort.h $(INC)/alignments.h $(INC)/blast.h $(INC)/bytepackGappedScoring.h $(INC)/descriptions.h $(INC)/encoding.h $(INC)/fasterGappedExtension.h $(INC)/gappedExtension.h $(INC)/nuGappedScoring.h $(INC)/gappedScoring.h $(INC)/global.h $(INC)/hitMatrix.h $(INC)/karlin.h $(INC)/memBlocks.h $(INC)/memSingleBlock.h $(INC)/nucleotideLookup.h $(INC)/oldGappedScoring.h $(INC)/oldSemiGappedScoring.h $(INC)/parameters.h $(INC)/print.h $(INC)/PSSMatrix.h $(INC)/qPosList.h $(INC)/readFasta.h $(INC)/readFile.h $(INC)/scoreMatrix.h $(INC)/semiGappedScoring.h $(INC)/statistics.h $(INC)/ungappedExtension.h $(INC)/wordLookupDFA.h $(INC)/writeFile.h $(INC)/constants.h $(INC)/smithWatermanTraceback.h $(INC)/smithWatermanScoring.h $(INC)/tableGappedScoring.h $(INC)/vbyte.h $(INC)/unpack.h $(INC)/index.h $(INC)/postings.h $(INC)/hashcounter.h $(INC)/writedb.h $(INC)/readdb.h $(INC)/search.h $(INC)/wildcards.h $(INC)/dust.h $(INC)/seg.h $(INC)/timeRec.h

SRCFILES=$(SRC)/mgpucontext.cu $(SRC)/mgpuutil.cpp $(SRC)/alignments_multi.c  $(SRC)/alignments.c $(SRC)/cuBLASTP.c $(SRC)/bytepackGappedScoring.c $(SRC)/descriptions.c $(SRC)/encoding.c $(SRC)/fasterGappedExtension_multi.c $(SRC)/fasterGappedExtension.c $(SRC)/gappedExtension.c $(SRC)/gappedScoring.c $(SRC)/nuGappedScoring.c $(SRC)/global.c $(SRC)/hitMatrix.c $(SRC)/karlin.c $(SRC)/memBlocks.c $(SRC)/memSingleBlock.c $(SRC)/nucleotideLookup.c $(SRC)/oldGappedScoring.c $(SRC)/oldSemiGappedScoring.c $(SRC)/parameters.c $(SRC)/print.c $(SRC)/PSSMatrix.c $(SRC)/qPosList.c $(SRC)/readFasta.c $(SRC)/readFile.c $(SRC)/scoreMatrix.c $(SRC)/semiGappedScoring_multi.c $(SRC)/semiGappedScoring.c $(SRC)/statistics.c $(SRC)/ungappedExtension.c $(SRC)/wordLookupDFA.c $(SRC)/writeFile.c $(SRC)/constants.c $(SRC)/smithWatermanTraceback.c $(SRC)/smithWatermanScoring.c $(SRC)/tableGappedScoring.c $(SRC)/vbyte.c $(SRC)/unpack.c $(SRC)/index.c $(SRC)/postings.c $(SRC)/hashcounter.c $(SRC)/writedb.c $(SRC)/readdb.c $(SRC)/search.cu $(SRC)/search.c $(SRC)/wildcards.c $(SRC)/dust.c $(SRC)/seg.c $(SRC)/timeRec.c $(SRC)/segsort.cu

all: cuBLASTP formatdb #readdb ssearch cluster

cuBLASTP: $(OBJ)/cuBLASTP.o $(OBJFILES)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) -o cuBLASTP $(OBJ)/cuBLASTP.o $(OBJFILES)
ssearch: $(OBJ)/ssearch.o $(OBJFILES)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) -o ssearch $(OBJ)/ssearch.o $(OBJFILES)
verboseBlast: $(SRC)/cuBLASTP.c $(SRCFILES)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) -DVERBOSE -o verboseBlast $(SRCFILES)
formatdb: $(OBJ)/formatdb.o $(OBJFILES)
	$(NVCC) $(CFLAGS) $(LDFLAGS) -o formatdb $(OBJ)/formatdb.o $(OBJFILES)
createindex: $(OBJ)/createindex.o $(OBJFILES)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) -o createindex $(OBJ)/createindex.o $(OBJFILES)
readdb: $(OBJ)/readdbApp.o $(OBJFILES)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) -o readdb $(OBJ)/readdbApp.o $(OBJFILES)
chooseWilds: $(OBJ)/chooseWilds.o $(OBJFILES)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) -o chooseWilds $(OBJ)/chooseWilds.o $(OBJFILES)
dust: $(OBJ)/dustApp.o $(OBJFILES)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) -o dust $(OBJ)/dustApp.o $(OBJFILES)
printDescription: $(OBJ)/printDescription.o $(OBJFILES)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) -o printDescription $(OBJ)/printDescription.o $(OBJFILES)
cluster: $(OBJ)/cluster.o $(OBJFILES)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) -o cluster $(OBJ)/cluster.o $(OBJFILES)
rsdb: $(OBJ)/rsdb.o $(OBJFILES)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) -o rsdb $(OBJ)/rsdb.o $(OBJFILES)
readNcbidb: $(OBJ)/readNcbidb.o $(OBJFILES)
	$(NVCC) $(NVCCFLAGS) $(LDFLAGS) -o readNcbidb $(OBJ)/readNcbidb.o $(OBJFILES)


$(OBJ)/chooseWilds.o: $(SRC)/chooseWilds.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/chooseWilds.o $(SRC)/chooseWilds.c
$(OBJ)/readdbApp.o: $(SRC)/readdbApp.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/readdbApp.o $(SRC)/readdbApp.c
$(OBJ)/dustApp.o: $(SRC)/dustApp.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/dustApp.o $(SRC)/dustApp.c
$(OBJ)/printDescription.o: $(SRC)/printDescription.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/printDescription.o $(SRC)/printDescription.c
$(OBJ)/formatdb.o: $(SRC)/formatdb.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/formatdb.o $(SRC)/formatdb.c
$(OBJ)/createindex.o: $(SRC)/createindex.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/createindex.o $(SRC)/createindex.c
$(OBJ)/cuBLASTP.o: $(HEADERFILES) $(SRC)/cuBLASTP.c
	$(CCPLUS) $(CFLAGS) -c -o $(OBJ)/cuBLASTP.o $(SRC)/cuBLASTP.c
$(OBJ)/cluster.o: $(SRC)/cluster.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/cluster.o $(SRC)/cluster.c
$(OBJ)/rsdb.o: $(SRC)/rsdb.c $(SRC)/identityAlign.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/rsdb.o $(SRC)/rsdb.c
$(OBJ)/readNcbidb.o: $(SRC)/readNcbidb.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/readNcbidb.o $(SRC)/readNcbidb.c

$(OBJ)/alignments.o: $(SRC)/alignments.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/alignments.o $(SRC)/alignments.c
$(OBJ)/alignmentscu.o: $(SRC)/alignments.cu $(HEADERFILES)
	$(NVCC) $(NVCCFLAGS) -c -o $(OBJ)/alignmentscu.o $(SRC)/alignments.cu

$(OBJ)/mgpuutil.o: $(SRC)/mgpuutil.cpp $(HEADERFILES)
	$(NVCC) $(CFLAGS) -c -o $(OBJ)/mgpuutil.o $(SRC)/mgpuutil.cpp
$(OBJ)/mgpucontextcu.o: $(SRC)/mgpucontext.cu $(HEADERFILES)
	$(NVCC) $(NVCCFLAGS) -c -o $(OBJ)/mgpucontextcu.o $(SRC)/mgpucontext.cu

$(OBJ)/bytepackGappedScoring.o: $(SRC)/bytepackGappedScoring.c $(SRC)/fasterBytepackGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/bytepackGappedScoring.o $(SRC)/bytepackGappedScoring.c
$(OBJ)/constants.o: $(SRC)/constants.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/constants.o $(SRC)/constants.c
$(OBJ)/descriptions.o: $(SRC)/descriptions.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/descriptions.o $(SRC)/descriptions.c
$(OBJ)/encoding.o: $(SRC)/encoding.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/encoding.o $(SRC)/encoding.c
$(OBJ)/gappedExtension.o: $(SRC)/gappedExtension.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/gappedExtension.o $(SRC)/gappedExtension.c
$(OBJ)/fasterGappedExtension.o: $(SRC)/fasterGappedExtension.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/fasterGappedExtension.o $(SRC)/fasterGappedExtension.c
$(OBJ)/gappedScoring.o: $(SRC)/gappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/gappedScoring.o $(SRC)/gappedScoring.c
$(OBJ)/nuGappedScoring.o: $(SRC)/nuGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/nuGappedScoring.o $(SRC)/nuGappedScoring.c
$(OBJ)/global.o: $(SRC)/global.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/global.o $(SRC)/global.c
$(OBJ)/hitMatrix.o: $(SRC)/hitMatrix.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/hitMatrix.o $(SRC)/hitMatrix.c
$(OBJ)/karlin.o: $(SRC)/karlin.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/karlin.o $(SRC)/karlin.c
$(OBJ)/memBlocks.o: $(SRC)/memBlocks.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/memBlocks.o $(SRC)/memBlocks.c
$(OBJ)/memSingleBlock.o: $(SRC)/memSingleBlock.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/memSingleBlock.o $(SRC)/memSingleBlock.c
$(OBJ)/nucleotideLookup.o: $(SRC)/nucleotideLookup.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/nucleotideLookup.o $(SRC)/nucleotideLookup.c
$(OBJ)/oldGappedScoring.o: $(SRC)/oldGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/oldGappedScoring.o $(SRC)/oldGappedScoring.c
$(OBJ)/oldSemiGappedScoring.o: $(SRC)/oldSemiGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/oldSemiGappedScoring.o $(SRC)/oldSemiGappedScoring.c
$(OBJ)/parameters.o: $(SRC)/parameters.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/parameters.o $(SRC)/parameters.c
$(OBJ)/print.o: $(SRC)/print.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/print.o $(SRC)/print.c
$(OBJ)/PSSMatrix.o: $(SRC)/PSSMatrix.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/PSSMatrix.o $(SRC)/PSSMatrix.c
$(OBJ)/qPosList.o: $(SRC)/qPosList.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/qPosList.o $(SRC)/qPosList.c
$(OBJ)/readFasta.o: $(SRC)/readFasta.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/readFasta.o $(SRC)/readFasta.c
$(OBJ)/readFile.o: $(SRC)/readFile.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/readFile.o $(SRC)/readFile.c
$(OBJ)/scoreMatrix.o: $(SRC)/scoreMatrix.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/scoreMatrix.o $(SRC)/scoreMatrix.c
$(OBJ)/semiGappedScoring.o: $(SRC)/semiGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/semiGappedScoring.o $(SRC)/semiGappedScoring.c
$(OBJ)/statistics.o: $(SRC)/statistics.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/statistics.o $(SRC)/statistics.c
$(OBJ)/ungappedExtension.o: $(SRC)/ungappedExtension.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/ungappedExtension.o $(SRC)/ungappedExtension.c
$(OBJ)/wordLookupDFA.o: $(SRC)/wordLookupDFA.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/wordLookupDFA.o $(SRC)/wordLookupDFA.c
$(OBJ)/writeFile.o: $(SRC)/writeFile.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/writeFile.o $(SRC)/writeFile.c
$(OBJ)/smithWatermanTraceback.o: $(SRC)/smithWatermanTraceback.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/smithWatermanTraceback.o $(SRC)/smithWatermanTraceback.c
$(OBJ)/smithWatermanScoring.o: $(SRC)/smithWatermanScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/smithWatermanScoring.o $(SRC)/smithWatermanScoring.c
$(OBJ)/tableGappedScoring.o: $(SRC)/tableGappedScoring.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/tableGappedScoring.o $(SRC)/tableGappedScoring.c
$(OBJ)/vbyte.o: $(SRC)/vbyte.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/vbyte.o $(SRC)/vbyte.c
$(OBJ)/unpack.o: $(SRC)/unpack.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/unpack.o $(SRC)/unpack.c
$(OBJ)/index.o: $(SRC)/index.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/index.o $(SRC)/index.c
$(OBJ)/hashcounter.o: $(SRC)/hashcounter.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/hashcounter.o $(SRC)/hashcounter.c
$(OBJ)/postings.o: $(SRC)/postings.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/postings.o $(SRC)/postings.c
$(OBJ)/writedb.o: $(SRC)/writedb.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/writedb.o $(SRC)/writedb.c
$(OBJ)/readdb.o: $(SRC)/readdb.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/readdb.o $(SRC)/readdb.c
$(OBJ)/searchcu.o: $(SRC)/search.cu $(HEADERFILES)
	$(NVCC) $(NVCCFLAGS) -c -o $(OBJ)/searchcu.o $(SRC)/search.cu
$(OBJ)/segsortcu.o: $(SRC)/segsort.cu $(HEADERFILES)
	$(NVCC) $(NVCCFLAGS) -c -o $(OBJ)/segsortcu.o $(SRC)/segsort.cu
$(OBJ)/search.o: $(SRC)/search.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/search.o $(SRC)/search.c
$(OBJ)/wildcards.o: $(SRC)/wildcards.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/wildcards.o $(SRC)/wildcards.c
$(OBJ)/ssearch.o: $(SRC)/cuBLASTP.c $(HEADERFILES)
	$(CC) $(CFLAGS) -DSSEARCH -c -o $(OBJ)/ssearch.o $(SRC)/cuBLASTP.c
$(OBJ)/dust.o: $(SRC)/dust.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/dust.o $(SRC)/dust.c
$(OBJ)/seg.o: $(SRC)/seg.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/seg.o $(SRC)/seg.c
$(OBJ)/timeRec.o: $(SRC)/timeRec.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/timeRec.o $(SRC)/timeRec.c

$(OBJ)/alignments_multi.o: $(SRC)/alignments_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/alignments_multi.o $(SRC)/alignments_multi.c
$(OBJ)/fasterGappedExtension_multi.o: $(SRC)/fasterGappedExtension_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/fasterGappedExtension_multi.o $(SRC)/fasterGappedExtension_multi.c
$(OBJ)/semiGappedScoring_multi.o: $(SRC)/semiGappedScoring_multi.c $(HEADERFILES)
	$(CC) $(CFLAGS) -c -o $(OBJ)/semiGappedScoring_multi.o $(SRC)/semiGappedScoring_multi.c

kclean:
	rm -rf cuBLASTP $(OBJ)/searchcu.o $(OBJ)/segsortcu.o

clean:
	rm -f cuBLASTP ssearch formatdb dust printDescription readdb verboseBlast $(OBJ)/* cluster rsdb chooseWilds createindex readNcbidb simpleCluster

