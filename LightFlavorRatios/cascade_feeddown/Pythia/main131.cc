#include <iostream>
#include <cstdlib>
#include "Pythia8/Pythia.h"

#ifndef HEPMC2
#include "Pythia8Plugins/HepMC3.h"
#else
#include "Pythia8Plugins/HepMC2.h"
#endif

using namespace Pythia8;

//==========================================================================

int main(int argc, char* argv[]) {

  // Defaults
  std::string outputFile = "main131.hepmc";
  int nEvents = 500000;
  int seed = 0; // 0 = time-based seed in Pythia

  // Parse command-line arguments
  if (argc > 1) {
    outputFile = argv[1];
  }
  if (argc > 2) {
    nEvents = std::atoi(argv[2]);
    if (nEvents <= 0) {
      std::cerr << "Error: number of events must be positive\n";
      return 1;
    }
  }
  if (argc > 3) {
    seed = std::atoi(argv[3]);
    if (seed < 0) {
      std::cerr << "Error: seed must be >= 0\n";
      return 1;
    }
  }
  if (argc > 4) {
    std::cerr << "Usage: " << argv[0]
              << " [output.hepmc] [nEvents] [seed]\n";
    return 1;
  }

  std::cout << "Writing HepMC output to: " << outputFile << "\n";
  std::cout << "Number of events: " << nEvents << "\n";
  std::cout << "Random seed: " << seed << "\n";

  // Interface for conversion from Pythia8::Event to HepMC event
  Pythia8ToHepMC toHepMC(outputFile);

  // Generator
  Pythia pythia;

  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));

  // Change the following line to use a different steering card if desired
  pythia.readFile("/sphenix/user/hjheng/sPHENIXRepo/analysis/LightFlavorRatios/macros/simulation/steeringCards/phpythia8_minBias_Detroit.cfg");

  // Initialize
  if (!pythia.init()) return 1;

  // Event loop
  for (int iEvent = 0; iEvent < nEvents; ++iEvent) {
    if (!pythia.next()) continue;
    toHepMC.writeNextEvent(pythia);
  }

  // Statistics
  pythia.stat();

  return 0;
}
