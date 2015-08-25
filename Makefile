#
# Macros
#
CC = /usr/bin/g++
CC_OPTIONS = -O3 -g -DNDEBUG -Wall -Wextra
CC_DEBUG_OPTIONS = -O0 -g -pg -Wall -Wextra
LNK_OPTIONS =
LNK_DEBUG_OPTIONS = -pg
BIN_DIR = bin
DEBUG_DIR = debug

#
# INCLUDE directories for motifspec
#
INCLUDE = -I.\
					-Isrc

#
# Build motifspec
#
all: motifspec motifspec-debug

motifspec: \
		bin/archivesites.o\
		bin/bgmodel.o\
		bin/motifspec.o\
		bin/fastmath.o\
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
		bin/motifspec.o\
		bin/fastmath.o\
		bin/motif.o\
		bin/motifcompare.o\
		bin/motifsearch.o\
		bin/motifsearchexpr.o\
		bin/motifsearchscore.o\
		bin/motifsearchsubset.o\
		bin/seqset.o\
		bin/site.o\
		bin/standard.o\
		-o bin/motifspec

motifspec-debug: \
		debug/archivesites.o\
		debug/bgmodel.o\
		debug/motifspec.o\
		debug/fastmath.o\
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
		debug/motifspec.o\
		debug/fastmath.o\
		debug/motif.o\
		debug/motifcompare.o\
		debug/motifsearch.o\
		debug/motifsearchexpr.o\
		debug/motifsearchscore.o\
		debug/motifsearchsubset.o\
		debug/seqset.o\
		debug/site.o\
		debug/standard.o\
		-o debug/motifspec-debug

clean: 
	rm -f $(BIN_DIR)/*.o $(BIN_DIR)/motifspec $(DEBUG_DIR)/*.o $(DEBUG_DIR)/motifspec-debug

dir_guard=@mkdir -p $(@D)

$(BIN_DIR)/%.o: src/%.cpp
	$(dir_guard)
	$(CC) $(CC_OPTIONS) -c $(INCLUDE) -o $@ $<

$(DEBUG_DIR)/%.o: src/%.cpp
	$(dir_guard)
	$(CC) $(CC_DEBUG_OPTIONS) -c $(INCLUDE) -o $@ $<
