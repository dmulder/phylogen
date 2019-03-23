python_PYTHON = phylogen.py

git_checkout:
	cd $(top_srcdir) && git submodule init && git submodule update && cd $(srcdir)

all: git_checkout astral.jar iqtree

astral.jar:
	pushd $(srcdir)/ASTRAL; ./make.sh; popd
	cp $(srcdir)/ASTRAL/astral*.jar $(srcdir)/astral.jar
	cp -r $(srcdir)/ASTRAL/lib $(srcdir)/

iqtree:
	pushd $(srcdir)/IQ-TREE && cmake CMakeLists.txt && make && popd
	cp $(srcdir)/IQ-TREE/iqtree $(srcdir)/iqtree

CLEANFILES = astral.jar iqtree
