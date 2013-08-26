bin_PROGRAMS = test_hybrid_coal
PROG = TEST_HYBRID_COAL
test_hybrid_coal_SOURCES = all_gene_topo.cpp  all_gene_topo.hpp \
coal.cpp  coal.hpp  \
node.cpp node.hpp \
newmain.cpp  \
net.cpp net.hpp \
utility.cpp  utility.hpp \
export.cpp export.hpp


test_hybrid_coal_CXXFLAGS= -g -O3 -Wall  #options for showing debug messages
#test_hybrid_coal_CXXFLAGS= -O3 -Wall -DNDEBUG   #options for not showing debug messages

include_HEADERS = net.hpp node.hpp utility.hpp export.hpp


##################### old make from here
#bin_PROGRAMS = hybrid_coal #test_hybrid_coal
#PROG = HYBRID_COAL
#hybrid_coal_SOURCES = all_gene_topo.cpp  all_gene_topo.hpp \
#node.cpp node.hpp \
#coal.cpp  coal.hpp  \
#main.cpp  \
#rm_hybrid.cpp  rm_hybrid.hpp \
#utility.cpp  utility.hpp \
#export.cpp export.hpp
##stells.cpp  stells.hpp 

#include_HEADERS = all_gene_topo.hpp coal.hpp rm_hybrid.hpp utility.hpp export.hpp
################################### old make til here



#DISTFILES = $(DIST_COMMON) $(DIST_SOURCES) $(TEXINFOS) $(EXTRA_DIST)
#PROG = TEST_HYBRID_COAL
#test_hybrid_coal_SOURCES = all_gene_topo.cpp  all_gene_topo.hpp \  
#coal.cpp  coal.hpp  \ 
#test.cpp  \ 
#rm_hybrid.cpp  rm_hybrid.hpp \
#utility.cpp  utility.hpp \
#export.cpp export.hpp

#include_HEADERS = all_gene_topo.hpp coal.hpp rm_hybrid.hpp utility.hpp export.hpp




#if DX_COND_man

# You'd probably want to post-process man pages and installed the patched
# versions.
#dist_man1_MANS = @DX_DOCDIR@/man/man1/doxample.c.1

#$(dist_man1_MANS): doxygen-doc

#endif

#MOSTLYCLEANFILES = $(DX_CLEANFILES)

#EXTRA_DIST = cleanup bootstrap $(TESTS) $(DX_CONFIG)
