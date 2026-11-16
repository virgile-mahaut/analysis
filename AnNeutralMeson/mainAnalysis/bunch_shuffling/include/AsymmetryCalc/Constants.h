#ifndef __CONSTANTS_H__
#define __CONSTANTS_H__

#include <cstdlib>
#include <string>

namespace ASYM_CONSTANTS {
  static constexpr size_t nBeams = 2; // Yellow or blue beam
  const std::string beams[nBeams] = {"yellow", "blue"};
  static constexpr size_t nParticles = 2; // pi0 or eta
  const std::string particle[nParticles] = {"pi0", "eta"};
  static constexpr size_t nRegions = 2; // peak band or side_band invariant mass region
  const std::string regions[nRegions] = {"peak", "side"};
  static constexpr size_t nPtBins = 8; // 9 default
  static constexpr size_t nEtaBins = 8;
  static constexpr size_t nXfBins = 8;
  static constexpr size_t nDirections = 2;
  const std::string directions[nDirections] = {"forward", "backward"};
  static constexpr size_t nSpins = 2; // up or down spin
  const std::string spins[nSpins] = {"up", "down"};
  static constexpr size_t nBunches = 111; // Number of proton bunches at RHIC (Run 24)
  static constexpr size_t nPhiBins = 12;
  static constexpr size_t nConfigs = 4;
  const std::string configuration[nConfigs] = {"0mrad_mbd", "0mrad_photon", "15mrad_mbd", "15mrad_photon"};
}

#endif
