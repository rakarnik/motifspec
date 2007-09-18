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
CC_OPTIONS = -g
LNK_OPTIONS = -lpthread


#
# INCLUDE directories for im
#

INCLUDE = -I.\
		-Isrc


#
# Build im
#

im : \
		bin/alignace.o\
		bin/archivesites.o\
		bin/cluster.o\
		bin/compareace.o\
		bin/im.o\
		bin/myheap.o\
		bin/seqset.o\
		bin/sites.o\
		bin/standard.o
	$(CC) $(LNK_OPTIONS) \
		bin/alignace.o\
		bin/archivesites.o\
		bin/cluster.o\
		bin/compareace.o\
		bin/im.o\
		bin/myheap.o\
		bin/seqset.o\
		bin/sites.o\
		bin/standard.o\
		-o im

clean : 
		rm -f\
		bin/alignace.o\
		bin/archivesites.o\
		bin/cluster.o\
		bin/compareace.o\
		bin/im.o\
		bin/main-AlignACE.o\
		bin/main-CompareACE.o\
		bin/myheap.o\
		bin/seqset.o\
		bin/sites.o\
		bin/standard.o\
		im

install : im
		cp im im

#
# Build the parts of im
#


# Item # 1 -- alignace --
bin/alignace.o : src/alignace.cpp
	$(CC) $(CC_OPTIONS) src/alignace.cpp -c $(INCLUDE) -o bin/alignace.o


# Item # 2 -- archivesites --
bin/archivesites.o : src/archivesites.cpp
	$(CC) $(CC_OPTIONS) src/archivesites.cpp -c $(INCLUDE) -o bin/archivesites.o


# Item # 3 -- cluster --
bin/cluster.o : src/cluster.cpp
	$(CC) $(CC_OPTIONS) src/cluster.cpp -c $(INCLUDE) -o bin/cluster.o


# Item # 4 -- compareace --
bin/compareace.o : src/compareace.cpp
	$(CC) $(CC_OPTIONS) src/compareace.cpp -c $(INCLUDE) -o bin/compareace.o


# Item # 5 -- im --
bin/im.o : src/im.cpp
	$(CC) $(CC_OPTIONS) src/im.cpp -c $(INCLUDE) -o bin/im.o


# Item # 6 -- main-AlignACE --
bin/main-AlignACE.o : src/main-AlignACE.cpp
	$(CC) $(CC_OPTIONS) src/main-AlignACE.cpp -c $(INCLUDE) -o bin/main-AlignACE.o


# Item # 7 -- main-CompareACE --
bin/main-CompareACE.o : src/main-CompareACE.cpp
	$(CC) $(CC_OPTIONS) src/main-CompareACE.cpp -c $(INCLUDE) -o bin/main-CompareACE.o


# Item # 8 -- myheap --
bin/myheap.o : src/myheap.cpp
	$(CC) $(CC_OPTIONS) src/myheap.cpp -c $(INCLUDE) -o bin/myheap.o


# Item # 9 -- seqset --
bin/seqset.o : src/seqset.cpp
	$(CC) $(CC_OPTIONS) src/seqset.cpp -c $(INCLUDE) -o bin/seqset.o


# Item # 10 -- sites --
bin/sites.o : src/sites.cpp
	$(CC) $(CC_OPTIONS) src/sites.cpp -c $(INCLUDE) -o bin/sites.o


# Item # 11 -- standard --
bin/standard.o : src/standard.cpp
	$(CC) $(CC_OPTIONS) src/standard.cpp -c $(INCLUDE) -o bin/standard.o


##### END RUN ####
