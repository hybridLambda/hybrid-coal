EXTRA_DIST = trees/*tax* bootstrap
#DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST) example/eg*
distdir = $(PACKAGE)-$(VERSION)

bin_PROGRAMS = hybrid-coal
PROG = HYBRID-COAL
check_PROGRAMS = hybrid-coal_dbg hybrid-coal_prof unittest
TESTS = unittest

common_flags = -O3 -g -std=c++0x -Isrc/graph/
hybrid_coal_CXXFLAGS = -DNDEBUG $(common_flags)#options for not showing debug messages.
hybrid_coal_dbg_CXXFLAGS = $(common_flags)
hybrid_coal_prof_CXXFLAGS = $(common_flags) -DNDEBUG -pg #options for not showing debug messages. -fpermissive

shared_src = src/graph/node.cpp \
			 src/graph/nodeContainer.cpp \
			 src/graph/graph.cpp \
			 src/plot/figure.cpp \
			 src/tree_topo/all_gene_topo.cpp \
			 src/coal.cpp \
			 src/hybridcoal.cpp


hybrid_coal_SOURCES      = $(shared_src) src/main.cpp
hybrid_coal_dbg_SOURCES  = $(shared_src) src/main.cpp
hybrid_coal_prof_SOURCES = $(shared_src) src/main.cpp

test_src = tests/unittest/test_runner.cpp tests/unittest/test_coal.cpp tests/unittest/test_Rm.cpp tests/unittest/test_4taxa_dist.cpp
unittest_CXXFLAGS = -DUNITTEST -DNDEBUG -std=c++0x -Isrc/graph/ -Isrc/
unittest_LDADD = -lcppunit -ldl
unittest_SOURCES = ${shared_src} ${test_src}

include_HEADERS = src/graph/exceptions.hpp \
                  src/graph/nodeContainer.hpp \
                  src/graph/graph.hpp \
                  src/graph/node.hpp \
                  src/plot/figure.hpp \
                  src/tree_topo/all_gene_topo.hpp \
                  src/coal.hpp \
                  src/hybridcoal.hpp
