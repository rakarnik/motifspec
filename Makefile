#
# Macros
#
CC = /usr/bin/g++
CC_OPTIONS = -O3 -DNDEBUG -Wall -Wextra
CC_DEBUG_OPTIONS = -O0 -g -pg -Wall -Wextra
LNK_OPTIONS = -pg
LNK_DEBUG_OPTIONS = -pg

#
# INCLUDE directories for im
#
INCLUDE = -I.\
					-Isrc


#
# Build im
#
all: im im-debug

im: \
		bin/archivesites.o\
		bin/bgmodel.o\
		bin/im.o\
		bin/motif.o\
		bin/motifcompare.o\
		bin/motifsearch.o\
		bin/motifsearchexpr.o\
		bin/motifsearchscore.o\
		bin/motifsearchsubset.o\
		bin/seqset.o\
		bin/site.o\
		bin/standard.o
	$(CC) $(LNK_OPTIONS) \
		bin/archivesites.o\
		bin/bgmodel.o\
		bin/im.o\
		bin/motif.o\
		bin/motifcompare.o\
		bin/motifsearch.o\
		bin/motifsearchexpr.o\
		bin/motifsearchscore.o\
		bin/motifsearchsubset.o\
		bin/seqset.o\
		bin/site.o\
		bin/standard.o\
		-o bin/im

im-debug: \
		debug/archivesites.o\
		debug/bgmodel.o\
		debug/im.o\
		debug/motif.o\
		debug/motifcompare.o\
		debug/motifsearch.o\
		debug/motifsearchexpr.o\
		debug/motifsearchscore.o\
		debug/motifsearchsubset.o\
		debug/seqset.o\
		debug/site.o\
		debug/standard.o
	$(CC) $(LNK_DEBUG_OPTIONS) \
		debug/archivesites.o\
		debug/bgmodel.o\
		debug/im.o\
		debug/motif.o\
		debug/motifcompare.o\
		debug/motifsearch.o\
		debug/motifsearchexpr.o\
		debug/motifsearchscore.o\
		debug/motifsearchsubset.o\
		debug/seqset.o\
		debug/site.o\
		debug/standard.o\
		-o debug/im-debug

clean: 
	rm -f bin/*.o debug/*.o bin/im debug/im-debug

bin/archivesites.o: src/archivesites.h src/archivesites.cpp
	$(CC) $(CC_OPTIONS) src/archivesites.cpp -c $(INCLUDE) -o bin/archivesites.o

bin/bgmodel.o: src/bgmodel.h src/bgmodel.cpp
	$(CC) $(CC_OPTIONS) src/bgmodel.cpp -c $(INCLUDE) -o bin/bgmodel.o

bin/im.o: src/im.h src/im.cpp
	$(CC) $(CC_OPTIONS) src/im.cpp -c $(INCLUDE) -o bin/im.o

bin/motif.o: src/motif.h src/motif.cpp
	$(CC) $(CC_OPTIONS) src/motif.cpp -c $(INCLUDE) -o bin/motif.o

bin/motifcompare.o: src/motifcompare.h src/motifcompare.cpp
	$(CC) $(CC_OPTIONS) src/motifcompare.cpp -c $(INCLUDE) -o bin/motifcompare.o

bin/motifsearch.o: src/motifsearch.h src/motifsearch.cpp
	$(CC) $(CC_OPTIONS) src/motifsearch.cpp -c $(INCLUDE) -o bin/motifsearch.o

bin/motifsearchexpr.o: src/motifsearchexpr.h src/motifsearchexpr.cpp
	$(CC) $(CC_OPTIONS) src/motifsearchexpr.cpp -c $(INCLUDE) -o bin/motifsearchexpr.o

bin/motifsearchscore.o: src/motifsearchscore.h src/motifsearchscore.cpp
	$(CC) $(CC_OPTIONS) src/motifsearchscore.cpp -c $(INCLUDE) -o bin/motifsearchscore.o

bin/motifsearchsubset.o: src/motifsearchsubset.h src/motifsearchsubset.cpp
	$(CC) $(CC_OPTIONS) src/motifsearchsubset.cpp -c $(INCLUDE) -o bin/motifsearchsubset.o

bin/seqset.o: src/seqset.h src/seqset.cpp
	$(CC) $(CC_OPTIONS) src/seqset.cpp -c $(INCLUDE) -o bin/seqset.o

bin/site.o: src/site.h src/site.cpp
	$(CC) $(CC_OPTIONS) src/site.cpp -c $(INCLUDE) -o bin/site.o

bin/standard.o: src/standard.h src/standard.cpp
	$(CC) $(CC_OPTIONS) src/standard.cpp -c $(INCLUDE) -o bin/standard.o

debug/archivesites.o: src/archivesites.h src/archivesites.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/archivesites.cpp -c $(INCLUDE) -o debug/archivesites.o

debug/bgmodel.o: src/bgmodel.h src/bgmodel.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/bgmodel.cpp -c $(INCLUDE) -o debug/bgmodel.o

debug/compareace.o: src/compareace.h src/compareace.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/compareace.cpp -c $(INCLUDE) -o debug/compareace.o

debug/im.o: src/im.h src/im.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/im.cpp -c $(INCLUDE) -o debug/im.o

debug/motif.o: src/motif.h src/motif.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/motif.cpp -c $(INCLUDE) -o debug/motif.o

debug/motifcompare.o: src/motifcompare.h src/motifcompare.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/motifcompare.cpp -c $(INCLUDE) -o debug/motifcompare.o

debug/motifsearch.o: src/motifsearch.h src/motifsearch.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/motifsearch.cpp -c $(INCLUDE) -o debug/motifsearch.o

debug/motifsearchexpr.o: src/motifsearchexpr.h src/motifsearchexpr.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/motifsearchexpr.cpp -c $(INCLUDE) -o debug/motifsearchexpr.o

debug/motifsearchscore.o: src/motifsearchscore.h src/motifsearchscore.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/motifsearchscore.cpp -c $(INCLUDE) -o debug/motifsearchscore.o

debug/motifsearchsubset.o: src/motifsearchsubset.h src/motifsearchsubset.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/motifsearchsubset.cpp -c $(INCLUDE) -o debug/motifsearchsubset.o

debug/seqset.o: src/seqset.h src/seqset.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/seqset.cpp -c $(INCLUDE) -o debug/seqset.o

debug/site.o: src/site.h src/site.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/site.cpp -c $(INCLUDE) -o debug/site.o

debug/standard.o: src/standard.h src/standard.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/standard.cpp -c $(INCLUDE) -o debug/standard.o
