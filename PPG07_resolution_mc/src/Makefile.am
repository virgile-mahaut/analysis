AUTOMAKE_OPTIONS = foreign

lib_LTLIBRARIES = \
    libppg07mcana.la

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -isystem$(ROOTSYS)/include

pkginclude_HEADERS = \
  MC_calo_pythia_pi0.h

libppg07mcana_la_SOURCES = \
  MC_calo_pythia_pi0.cc

libppg07mcana_la_LIBADD = \
  -lcalo_io \
  -lphool \
  -lmbd_io \
  -lSubsysReco \
  -lfun4all \
  -lg4dst \
  -lg4detectors_io \
  -lphg4hit

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals

testexternals_SOURCES = testexternals.C
testexternals_LDADD = libppg07mcana.la

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
