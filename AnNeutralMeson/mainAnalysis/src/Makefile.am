AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libanneutral.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include \
  -isystem$(OPT_SPHENIX)/include

pkginclude_HEADERS = \
  AnNeutralMeson.h \
  AnNeutralMeson_micro_dst.h \
  AnNeutralMeson_nano.h \
  ClusterSmallInfo.h \
  ClusterSmallInfoContainer.h

ROOTDICTS = \
  ClusterSmallInfo_Dict.cc \
  ClusterSmallInfoContainer_Dict.cc

pcmdir = $(libdir)
nobase_dist_pcm_DATA = \
  ClusterSmallInfo_Dict_rdict.pcm \
  ClusterSmallInfoContainer_Dict_rdict.pcm

libanneutral_la_SOURCES = \
  $(ROOTDICTS) \
  AnNeutralMeson.cc \
  AnNeutralMeson_micro_dst.cc \
  AnNeutralMeson_nano.cc \
  ClusterSmallInfo.cc \
  ClusterSmallInfoContainer.cc

libanneutral_la_LIBADD = \
  -lphool \
  -lmbd_io \
  -lSubsysReco \
  -lfun4all \
  -lg4dst \
  -ltriggeremulator \
  -lpileupfilter \
  -luspin

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

#just to get the dependency
%_Dict_rdict.pcm: %_Dict.cc ;

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libanneutral.la

testexternals.C:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

# Rule for generating table CINT dictionaries.
%_Dict.cc: %.h %LinkDef.h
	rootcint -f $@ @CINTDEFS@ -c $(DEFAULT_INCLUDES) $(AM_CPPFLAGS) $^

clean-local:
	rm -f $(BUILT_SOURCES)
