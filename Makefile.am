EXTRA_DIST = trees/*tax* bootstrap example/eg* 
#DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
distdir = $(PACKAGE)-$(VERSION)

bin_PROGRAMS = hybrid-coal
PROG = HYBRID-COAL
#check_PROGRAMS = hybrid-coal_dbg hybrid-coal_prof unit_tests
#TESTS = unit_tests 

common_flags = -O3 -g -std=c++0x -I src/ 
hybrid_coal_CXXFLAGS = -DNDEBUG $(common_flags)#options for not showing debug messages.
hybrid_coal_dbg_CXXFLAGS = $(common_flags)
hybrid_coal_prof_CXXFLAGS = $(common_flags) -DNDEBUG -pg #options for not showing debug messages. -fpermissive

shared_src = src/node.cpp \
			 src/net.cpp \
			 src/plot/figure.cpp  \
			 src/tree_topo/all_gene_topo.cpp \
			 src/hybridcoal.cpp 
			  
hybrid_coal_SOURCES = $(shared_src) src/main.cpp  
hybrid_coal_dbg_SOURCES = $(shared_src) src/main.cpp src/debug/node_debug.cpp src/debug/net_debug.cpp
hybrid_coal_prof_SOURCES = $(shared_src) src/main.cpp  

#include_HEADERS = src/net.hpp src/node.hpp src/utility.hpp src/export.hpp
