python_PYTHON = phylogen.py

all: astral.jar

astral.jar:
	pushd $(srcdir)/ASTRAL; ./make.sh; popd
	cp $(srcdir)/ASTRAL/astral*.jar $(srcdir)/astral.jar
	cp -r $(srcdir)/ASTRAL/lib $(srcdir)/

CLEANFILES = astral.jar
