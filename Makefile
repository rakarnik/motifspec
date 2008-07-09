###################################################
#
# Makefile for im
# Creator [Xcode -> Makefile Ver: Feb 14 2007 09:18:41]
# Created: [Tue Aug 21 16:03:42 2007]
#
###################################################

#
# Macros
#

CC = /usr/bin/g++
CC_OPTIONS = -Os
CC_DEBUG_OPTIONS = -O0 -g -pg
LNK_OPTIONS =
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

im : \
		bin/archivesites.o\
		bin/compareace.o\
		bin/im.o\
		bin/motif.o\
		bin/myheap.o\
		bin/semodel.o\
		bin/seqset.o\
		bin/site.o\
		bin/standard.o
	$(CC) $(LNK_OPTIONS) \
		bin/archivesites.o\
		bin/compareace.o\
		bin/im.o\
		bin/motif.o\
		bin/myheap.o\
		bin/semodel.o\
		bin/seqset.o\
		bin/site.o\
		bin/standard.o\
		-o im

im-debug : \
		debug/archivesites.o\
		debug/compareace.o\
		debug/im.o\
		debug/motif.o\
		debug/myheap.o\
		debug/semodel.o\
		debug/seqset.o\
		debug/site.o\
		debug/standard.o
	$(CC) $(LNK_DEBUG_OPTIONS) \
		debug/archivesites.o\
		debug/compareace.o\
		debug/im.o\
		debug/motif.o\
		debug/myheap.o\
		debug/semodel.o\
		debug/seqset.o\
		debug/site.o\
		debug/standard.o\
		-o im-debug

clean : 
	rm -f\
	bin/*.o\
	debug/*.o\
	im\
	im-debug

bin/archivesites.o : src/archivesites.h src/archivesites.cpp
	$(CC) $(CC_OPTIONS) src/archivesites.cpp -c $(INCLUDE) -o bin/archivesites.o

bin/compareace.o : src/compareace.h src/compareace.cpp
	$(CC) $(CC_OPTIONS) src/compareace.cpp -c $(INCLUDE) -o bin/compareace.o

bin/im.o : src/im.h src/im.cpp
	$(CC) $(CC_OPTIONS) src/im.cpp -c $(INCLUDE) -o bin/im.o

bin/motif.o : src/motif.h src/motif.cpp
	$(CC) $(CC_OPTIONS) src/motif.cpp -c $(INCLUDE) -o bin/motif.o

bin/myheap.o : src/myheap.cpp
	$(CC) $(CC_OPTIONS) src/myheap.cpp -c $(INCLUDE) -o bin/myheap.o

bin/seqset.o : src/seqset.h src/seqset.cpp
	$(CC) $(CC_OPTIONS) src/seqset.cpp -c $(INCLUDE) -o bin/seqset.o

bin/semodel.o : src/semodel.h src/semodel.cpp
	$(CC) $(CC_OPTIONS) src/semodel.cpp -c $(INCLUDE) -o bin/semodel.o

bin/site.o : src/site.h src/site.cpp
	$(CC) $(CC_OPTIONS) src/site.cpp -c $(INCLUDE) -o bin/site.o

bin/standard.o : src/standard.h src/standard.cpp
	$(CC) $(CC_OPTIONS) src/standard.cpp -c $(INCLUDE) -o bin/standard.o

debug/archivesites.o : src/archivesites.h src/archivesites.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/archivesites.cpp -c $(INCLUDE) -o debug/archivesites.o

debug/compareace.o : src/compareace.h src/compareace.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/compareace.cpp -c $(INCLUDE) -o debug/compareace.o

debug/im.o : src/im.h src/im.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/im.cpp -c $(INCLUDE) -o debug/im.o

debug/motif.o : src/motif.h src/motif.cpp
	$(CC) $(CC_OPTIONS) src/motif.cpp -c $(INCLUDE) -o debug/motif.o

debug/myheap.o : src/myheap.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/myheap.cpp -c $(INCLUDE) -o debug/myheap.o

debug/seqset.o : src/seqset.h src/seqset.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/seqset.cpp -c $(INCLUDE) -o debug/seqset.o

debug/semodel.o : src/semodel.h src/semodel.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/semodel.cpp -c $(INCLUDE) -o debug/semodel.o

debug/site.o : src/site.h src/site.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/site.cpp -c $(INCLUDE) -o debug/site.o

debug/standard.o : src/standard.h src/standard.cpp
	$(CC) $(CC_DEBUG_OPTIONS) src/standard.cpp -c $(INCLUDE) -o debug/standard.o
