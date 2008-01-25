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
CC_OPTIONS = -O3
LNK_OPTIONS =


#
# INCLUDE directories for im
#

INCLUDE = -I.\
		-Isrc


#
# Build im
#

im : \
		bin/archivesites.o\
		bin/compareace.o\
		bin/im.o\
		bin/myheap.o\
		bin/semodel.o\
		bin/seqset.o\
		bin/sites.o\
		bin/standard.o
	$(CC) $(LNK_OPTIONS) \
		bin/archivesites.o\
		bin/compareace.o\
		bin/im.o\
		bin/myheap.o\
		bin/semodel.o\
		bin/seqset.o\
		bin/sites.o\
		bin/standard.o\
		-o im

clean : 
		rm -f\
		bin/*.o\
		im

install : im
		cp im im

#
# Build the parts of im
#


bin/archivesites.o : src/archivesites.cpp
	$(CC) $(CC_OPTIONS) src/archivesites.cpp -c $(INCLUDE) -o bin/archivesites.o

bin/compareace.o : src/compareace.cpp
	$(CC) $(CC_OPTIONS) src/compareace.cpp -c $(INCLUDE) -o bin/compareace.o

bin/im.o : src/im.cpp
	$(CC) $(CC_OPTIONS) src/im.cpp -c $(INCLUDE) -o bin/im.o

bin/myheap.o : src/myheap.cpp
	$(CC) $(CC_OPTIONS) src/myheap.cpp -c $(INCLUDE) -o bin/myheap.o

bin/seqset.o : src/seqset.cpp
	$(CC) $(CC_OPTIONS) src/seqset.cpp -c $(INCLUDE) -o bin/seqset.o

bin/semodel.o : src/semodel.cpp
	$(CC) $(CC_OPTIONS) src/semodel.cpp -c $(INCLUDE) -o bin/semodel.o

bin/sites.o : src/sites.cpp
	$(CC) $(CC_OPTIONS) src/sites.cpp -c $(INCLUDE) -o bin/sites.o

bin/standard.o : src/standard.cpp
	$(CC) $(CC_OPTIONS) src/standard.cpp -c $(INCLUDE) -o bin/standard.o

##### END RUN ####
