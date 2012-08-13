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
# INCLUDE directories for inspector
#
INCLUDE = -I.\
					-Isrc

#
# Build inspector
#
all: inspector inspector-debug

inspector: \
		bin/archivesites.o\
		bin/bgmodel.o\
		bin/inspector.o\
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
		bin/inspector.o\
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
		-o bin/inspector

inspector-debug: \
		debug/archivesites.o\
		debug/bgmodel.o\
		debug/inspector.o\
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
		debug/inspector.o\
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
		-o debug/inspector-debug

clean: 
	rm -f $(BIN_DIR)/*.o $(BIN_DIR)/inspector $(DEBUG_DIR)/*.o $(DEBUG_DIR)/inspector-debug

dir_guard=@mkdir -p $(@D)

$(BIN_DIR)/%.o: src/%.cpp
	$(dir_guard)
	$(CC) $(CC_OPTIONS) -c $(INCLUDE) -o $@ $<

$(DEBUG_DIR)/%.o: src/%.cpp
	$(dir_guard)
	$(CC) $(CC_DEBUG_OPTIONS) -c $(INCLUDE) -o $@ $<
