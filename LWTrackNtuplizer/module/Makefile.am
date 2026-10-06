AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -isystem$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -L$(ROOTSYS)/lib

pkginclude_HEADERS = \
  lwtrackntuplizer.h

lib_LTLIBRARIES = \
  liblwtrackntuplizer.la

liblwtrackntuplizer_la_SOURCES = \
  lwtrackntuplizer.cc

liblwtrackntuplizer_la_LIBADD = \
  -lCore \
  -lMathCore \
  -lPhysics \
  -lRIO \
  -lTree \
  -lfun4all \
  -lg4detectors \
  -lglobalvertex \
  -lphool \
  -lSubsysReco \
  -ltrack \
  -ltrackbase_historic

BUILT_SOURCES = testexternals.cc

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.cc
testexternals_LDADD   = liblwtrackntuplizer.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

clean-local:
	rm -f $(BUILT_SOURCES)
