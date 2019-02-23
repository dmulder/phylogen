python_PYTHON = phylogen.py

all: astral.jar iqtree

astral.jar:
	pushd $(srcdir)/ASTRAL; ./make.sh; popd
	cp $(srcdir)/ASTRAL/astral*.jar $(srcdir)/astral.jar
	cp -r $(srcdir)/ASTRAL/lib $(srcdir)/

iqtree:
	pushd $(srcdir)/IQ-TREE && cmake CMakeLists.txt && make && popd
	cp $(srcdir)/IQ-TREE/iqtree $(srcdir)/iqtree

CLEANFILES = astral.jar iqtree
