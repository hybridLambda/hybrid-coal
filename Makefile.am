EXTRA_DIST = trees/*tax* bootstrap example/eg* 
#DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
distdir = $(PACKAGE)-$(VERSION)

bin_PROGRAMS = hybrid-coal
PROG = HYBRID-COAL
#check_PROGRAMS = hybrid-coal_dbg hybrid-coal_prof unit_tests
#TESTS = unit_tests 

common_flags = -O3 -g -std=c++0x -I src/ 
hybrid_coal_CXXFLAGS = -DNDEBUG $(common_flags)#options for not showing debug messages.
#hybrid_coal_dbg_CXXFLAGS = $(common_flags)
#hybrid_coal_prof_CXXFLAGS = $(common_flags) -DNDEBUG -pg #options for not showing debug messages. -fpermissive

shared_SOURCES = src/hybridcoal.cpp \
				 src/node.cpp \
				 src/net.cpp \
				 src/sim_gt.cpp \
				 src/plot/figure.cpp  

hybrid_Lambda_SOURCES = $(shared_SOURCES) src/main.cpp  

test_hybrid_coal_SOURCES = all_gene_topo.cpp  all_gene_topo.hpp \
coal.cpp  coal.hpp  \
node.cpp node.hpp \
newmain.cpp  \
net.cpp net.hpp \
utility.cpp  utility.hpp \
export.cpp export.hpp



include_HEADERS = net.hpp node.hpp utility.hpp export.hpp
