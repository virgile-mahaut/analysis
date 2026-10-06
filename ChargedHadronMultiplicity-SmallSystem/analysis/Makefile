CFLAGS = -g3 -pthread -m64 -Wno-deprecated -I src/ 
CFLAGS += $(shell root-config --cflags --libs)
CFLAGS += -lEG

PROGRAMS = TrackletAna_Combinatoric
GITIGNORE = ../.gitignore
GITIGNORE_BEGIN = \# BEGIN AUTO-GENERATED BINARIES
GITIGNORE_END = \# END AUTO-GENERATED BINARIES

.PHONY: all print-binaries update-gitignore clean

all : update-gitignore $(PROGRAMS)

TrackletAna_Combinatoric : src/TrackletAna_Combinatoric.cxx src/Hit.h src/Tracklet.h src/Tracklet_combinatoric.h
	g++ ${CFLAGS} src/TrackletAna_Combinatoric.cxx src/Hit.h src/Tracklet.h src/Tracklet_combinatoric.h -o TrackletAna_Combinatoric

print-binaries :
	@printf 'analysis/%s\n' $(PROGRAMS)

update-gitignore :
	@awk 'BEGIN {inside = 0} \
		$$0 == "$(GITIGNORE_BEGIN)" {inside = 1; next} \
		$$0 == "$(GITIGNORE_END)" {inside = 0; next} \
		!inside {print}' $(GITIGNORE) > $(GITIGNORE).tmp
	@printf '%s\n' "$(GITIGNORE_BEGIN)" >> $(GITIGNORE).tmp
	@$(MAKE) --no-print-directory print-binaries >> $(GITIGNORE).tmp
	@printf '%s\n' "$(GITIGNORE_END)" >> $(GITIGNORE).tmp
	@mv $(GITIGNORE).tmp $(GITIGNORE)

clean :
	rm -f $(PROGRAMS)
