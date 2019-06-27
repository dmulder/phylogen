AUTOMAKE_OPTIONS = subdir-objects

python_PYTHON = phylogen.py

iqtree_wrap.cxx: iqtree.i iqtree.h
	$(SWIG) $(AX_SWIG_PYTHON_OPT) -c++ -o $@ $<

git_checkout:
	cd $(top_srcdir) && git submodule init && git submodule update && cd $(srcdir)

astral.jar: git_checkout
	pushd $(srcdir)/ASTRAL; ./make.sh; popd
	cp $(srcdir)/ASTRAL/astral*.jar $(srcdir)/astral.jar
	cp -r $(srcdir)/ASTRAL/lib $(srcdir)/

$(srcdir)/IQ-TREE/iqtree_config.h: git_checkout
	pushd $(srcdir)/IQ-TREE && cmake CMakeLists.txt && popd

lib_LTLIBRARIES = libiqtree.la
libiqtree_la_DEPENDENCIES = $(srcdir)/IQ-TREE/iqtree_config.h
libiqtree_la_LIBADD = -L$(srcdir)/IQ-TREE/lib -lm -lz -lgomp
libiqtree_la_LDFLAGS = -Wl,-rpath,$(srcdir)/IQ-TREE/lib -Wl,--gc-sections -rdynamic ${PYTHON_LDFLAGS}
libiqtree_la_CPPFLAGS = -DIQ_TREE -D_USE_PTHREADS -I$(srcdir)/IQ-TREE -fPIC -std=c++11 -fopenmp -pthread -O3 -g -ffunction-sections -fdata-sections ${PYTHON_CPPFLAGS} -I/usr/include/eigen3 -DUSE_EIGEN3 -D__SSE3 -msse3 -Dmain=main_entry -D__AVX -mavx -fabi-version=0 -mfma
libiqtree_la_SOURCES =	iqtree_wrap.cxx \
			$(srcdir)/IQ-TREE/tree/phylokernelsse.cpp \
			$(srcdir)/IQ-TREE/gsl/binomial_tpe.cpp \
			$(srcdir)/IQ-TREE/gsl/multinomial.cpp \
			$(srcdir)/IQ-TREE/gsl/pow_int.cpp \
			$(srcdir)/IQ-TREE/gsl/gauss.cpp \
			$(srcdir)/IQ-TREE/gsl/gaussinv.cpp \
			$(srcdir)/IQ-TREE/gsl/gausspdf.cpp \
			$(srcdir)/IQ-TREE/ncl/nxsassumptionsblock.cpp \
			$(srcdir)/IQ-TREE/ncl/nxsblock.cpp \
			$(srcdir)/IQ-TREE/ncl/nxscharactersblock.cpp \
			$(srcdir)/IQ-TREE/ncl/nxsdatablock.cpp \
			$(srcdir)/IQ-TREE/ncl/nxsdiscretedatum.cpp \
			$(srcdir)/IQ-TREE/ncl/nxsdiscretematrix.cpp \
			$(srcdir)/IQ-TREE/ncl/nxsdistancedatum.cpp \
			$(srcdir)/IQ-TREE/ncl/nxsdistancesblock.cpp \
			$(srcdir)/IQ-TREE/ncl/nxsemptyblock.cpp \
			$(srcdir)/IQ-TREE/ncl/nxsexception.cpp \
			$(srcdir)/IQ-TREE/ncl/nxsreader.cpp \
			$(srcdir)/IQ-TREE/ncl/nxssetreader.cpp \
			$(srcdir)/IQ-TREE/ncl/nxsstring.cpp \
			$(srcdir)/IQ-TREE/ncl/nxstaxablock.cpp \
			$(srcdir)/IQ-TREE/ncl/nxstoken.cpp \
			$(srcdir)/IQ-TREE/ncl/nxstreesblock.cpp \
			$(srcdir)/IQ-TREE/alignment/pattern.cpp \
			$(srcdir)/IQ-TREE/alignment/alignment.cpp \
			$(srcdir)/IQ-TREE/alignment/alignmentpairwise.cpp \
			$(srcdir)/IQ-TREE/alignment/maalignment.cpp \
			$(srcdir)/IQ-TREE/alignment/superalignment.cpp \
			$(srcdir)/IQ-TREE/alignment/superalignmentpairwise.cpp \
			$(srcdir)/IQ-TREE/alignment/superalignmentpairwiseplen.cpp \
			$(srcdir)/IQ-TREE/pll/alignment.c \
			$(srcdir)/IQ-TREE/pll/bipartitionList.c \
			$(srcdir)/IQ-TREE/pll/evaluateGenericSpecial.c \
			$(srcdir)/IQ-TREE/pll/evaluatePartialGenericSpecial.c \
			$(srcdir)/IQ-TREE/pll/fastDNAparsimony.c	\
			$(srcdir)/IQ-TREE/pll/hardware.c \
			$(srcdir)/IQ-TREE/pll/hash.c \
			$(srcdir)/IQ-TREE/pll/lexer.c \
			$(srcdir)/IQ-TREE/pll/makenewzGenericSpecial.c \
			$(srcdir)/IQ-TREE/pll/models.c \
			$(srcdir)/IQ-TREE/pll/newick.c \
			$(srcdir)/IQ-TREE/pll/newviewGenericSpecial.c \
			$(srcdir)/IQ-TREE/pll/genericParallelization.c \
			$(srcdir)/IQ-TREE/pll/optimizeModel.c \
			$(srcdir)/IQ-TREE/pll/parsePartition.c \
			$(srcdir)/IQ-TREE/pll/queue.c \
			$(srcdir)/IQ-TREE/pll/randomTree.c \
			$(srcdir)/IQ-TREE/pll/recom.c \
			$(srcdir)/IQ-TREE/pll/restartHashTable.c \
			$(srcdir)/IQ-TREE/pll/searchAlgo.c \
			$(srcdir)/IQ-TREE/pll/ssort.c \
			$(srcdir)/IQ-TREE/pll/stack.c \
			$(srcdir)/IQ-TREE/pll/topologies.c \
			$(srcdir)/IQ-TREE/pll/trash.c \
			$(srcdir)/IQ-TREE/pll/treeIO.c \
			$(srcdir)/IQ-TREE/pll/utils.c \
			$(srcdir)/IQ-TREE/lbfgsb/lbfgsb_new.cpp \
			$(srcdir)/IQ-TREE/utils/eigendecomposition.cpp \
			$(srcdir)/IQ-TREE/utils/gzstream.cpp \
			$(srcdir)/IQ-TREE/utils/optimization.cpp \
			$(srcdir)/IQ-TREE/utils/stoprule.cpp \
			$(srcdir)/IQ-TREE/utils/tools.cpp \
			$(srcdir)/IQ-TREE/utils/pllnni.cpp \
			$(srcdir)/IQ-TREE/utils/checkpoint.cpp \
			$(srcdir)/IQ-TREE/utils/MPIHelper.cpp \
			$(srcdir)/IQ-TREE/model/modelmarkov.cpp \
			$(srcdir)/IQ-TREE/model/modelbin.cpp \
			$(srcdir)/IQ-TREE/model/modeldna.cpp \
			$(srcdir)/IQ-TREE/model/modelfactory.cpp \
			$(srcdir)/IQ-TREE/model/modelprotein.cpp \
			$(srcdir)/IQ-TREE/model/modelset.cpp \
			$(srcdir)/IQ-TREE/model/modelsubst.cpp \
			$(srcdir)/IQ-TREE/model/partitionmodel.cpp \
			$(srcdir)/IQ-TREE/model/partitionmodelplen.cpp \
			$(srcdir)/IQ-TREE/model/rategamma.cpp \
			$(srcdir)/IQ-TREE/model/rategammainvar.cpp \
			$(srcdir)/IQ-TREE/model/rateheterogeneity.cpp \
			$(srcdir)/IQ-TREE/model/rateinvar.cpp \
			$(srcdir)/IQ-TREE/model/ratemeyerdiscrete.cpp \
			$(srcdir)/IQ-TREE/model/ratemeyerhaeseler.cpp \
			$(srcdir)/IQ-TREE/model/ratekategory.cpp \
			$(srcdir)/IQ-TREE/model/ratefree.cpp \
			$(srcdir)/IQ-TREE/model/ratefreeinvar.cpp \
			$(srcdir)/IQ-TREE/model/rateheterotachy.cpp \
			$(srcdir)/IQ-TREE/model/rateheterotachyinvar.cpp \
			$(srcdir)/IQ-TREE/model/modelcodon.cpp \
			$(srcdir)/IQ-TREE/model/modelmorphology.cpp \
			$(srcdir)/IQ-TREE/model/modelmixture.cpp \
			$(srcdir)/IQ-TREE/model/modelunrest.cpp \
			$(srcdir)/IQ-TREE/model/modelliemarkov.cpp \
			$(srcdir)/IQ-TREE/model/modelpomo.cpp \
			$(srcdir)/IQ-TREE/model/modelpomomixture.cpp \
			$(srcdir)/IQ-TREE/model/modelfactorymixlen.cpp \
			$(srcdir)/IQ-TREE/tree/constrainttree.cpp \
			$(srcdir)/IQ-TREE/tree/candidateset.cpp \
			$(srcdir)/IQ-TREE/tree/iqtree.cpp \
			$(srcdir)/IQ-TREE/tree/matree.cpp \
			$(srcdir)/IQ-TREE/tree/memslot.cpp \
			$(srcdir)/IQ-TREE/tree/mexttree.cpp \
			$(srcdir)/IQ-TREE/tree/mtree.cpp \
			$(srcdir)/IQ-TREE/tree/mtreeset.cpp \
			$(srcdir)/IQ-TREE/tree/ncbitree.cpp \
			$(srcdir)/IQ-TREE/tree/node.cpp \
			$(srcdir)/IQ-TREE/tree/phylonode.cpp \
			$(srcdir)/IQ-TREE/tree/phylonodemixlen.cpp \
			$(srcdir)/IQ-TREE/tree/phylosupertree.cpp \
			$(srcdir)/IQ-TREE/tree/phylosupertreeplen.cpp \
			$(srcdir)/IQ-TREE/tree/phylotree.cpp \
			$(srcdir)/IQ-TREE/tree/phylotreemixlen.cpp \
			$(srcdir)/IQ-TREE/tree/phylotreepars.cpp \
			$(srcdir)/IQ-TREE/tree/phylotreesse.cpp \
			$(srcdir)/IQ-TREE/tree/quartet.cpp \
			$(srcdir)/IQ-TREE/tree/supernode.cpp \
			$(srcdir)/IQ-TREE/tree/tinatree.cpp \
			$(srcdir)/IQ-TREE/tree/parstree.cpp \
			$(srcdir)/IQ-TREE/tree/phylokernelfma.cpp \
			$(srcdir)/IQ-TREE/tree/phylotreeavx.cpp \
			$(srcdir)/IQ-TREE/pll/avxLikelihood.c \
			$(srcdir)/IQ-TREE/nclextra/mpdablock.cpp \
			$(srcdir)/IQ-TREE/nclextra/msetsblock.cpp \
			$(srcdir)/IQ-TREE/nclextra/msplitsblock.cpp \
			$(srcdir)/IQ-TREE/nclextra/modelsblock.cpp \
			$(srcdir)/IQ-TREE/pda/circularnetwork.cpp \
			$(srcdir)/IQ-TREE/pda/greedy.cpp \
			$(srcdir)/IQ-TREE/pda/gurobiwrapper.cpp \
			$(srcdir)/IQ-TREE/pda/hashsplitset.cpp \
			$(srcdir)/IQ-TREE/pda/pdnetwork.cpp \
			$(srcdir)/IQ-TREE/pda/pruning.cpp \
			$(srcdir)/IQ-TREE/pda/split.cpp \
			$(srcdir)/IQ-TREE/pda/splitgraph.cpp \
			$(srcdir)/IQ-TREE/pda/splitset.cpp \
			$(srcdir)/IQ-TREE/pda/lpwrapper.c \
			$(srcdir)/IQ-TREE/pda/ecopd.cpp \
			$(srcdir)/IQ-TREE/pda/graph.cpp \
			$(srcdir)/IQ-TREE/pda/pdtree.cpp \
			$(srcdir)/IQ-TREE/pda/pdtreeset.cpp \
			$(srcdir)/IQ-TREE/pda/ecopdmtreeset.cpp \
			$(srcdir)/IQ-TREE/whtest/eigen.c \
			$(srcdir)/IQ-TREE/whtest/eigen_sym.c \
			$(srcdir)/IQ-TREE/whtest/random.c \
			$(srcdir)/IQ-TREE/whtest/weisslambda.c \
			$(srcdir)/IQ-TREE/whtest/weisslambda_sub.c \
			$(srcdir)/IQ-TREE/whtest/whtest.c \
			$(srcdir)/IQ-TREE/whtest/whtest_sub.c \
			$(srcdir)/IQ-TREE/whtest/whtest_wrapper.cpp \
			$(srcdir)/IQ-TREE/sprng/lcg64.c \
			$(srcdir)/IQ-TREE/sprng/makeseed.c \
			$(srcdir)/IQ-TREE/sprng/memory.c \
			$(srcdir)/IQ-TREE/sprng/store.c \
			$(srcdir)/IQ-TREE/sprng/primes-lcg64.c \
			$(srcdir)/IQ-TREE/sprng/checkid.c \
			$(srcdir)/IQ-TREE/vectorclass/instrset_detect.cpp \
			$(srcdir)/IQ-TREE/main/main.cpp \
			$(srcdir)/IQ-TREE/main/phyloanalysis.cpp \
			$(srcdir)/IQ-TREE/main/phylotesting.cpp \
			$(srcdir)/IQ-TREE/obsolete/parsmultistate.cpp \
			$(srcdir)/iqtree.cpp

CLEANFILES = astral.jar
clean-local:
	pushd $(srcdir)/IQ-TREE && make clean; git clean -xfd; popd
	pushd $(srcdir)/ASTRAL && git clean -xfd; popd
