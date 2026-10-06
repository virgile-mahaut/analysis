# Standalone PYTHIA8 setup

## Obtain the PYTHIA8 package and compile the pythia generation code
- Go to the PYTHIA8 website https://www.pythia.org/ and download the stanalone tarball 
    - (Version 8.316 was used as of writing. If you use a different version, make sure to update the version number in the commands below accordingly)
- Unpack the tarball and configure to link the HepMC3 package/libraries (these instructions are from Cameron)
```Bash
tar -zxvf pythia8316.tgz
cd pythia8316
./configure  --with-hepmc3-include=/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.15/include --with-hepmc3-lib=/cvmfs/sphenix.sdcc.bnl.gov/alma9.2-gcc-14.2.0/release/release_new/new.15/lib64
make -j $(nproc)
```
- You can copy the `main131.cc` file to the `examples` directory and compile it. To compile, simply do
```Bash
make main131
```
- If successful, you should see an executable named `main131` in the `examples` directory. You can run it with
```Bash
./main131
```
or with arguments different from the default ones
```Bash
./main131 <output file name (default: main131.hepmc)> <number of events (default: 500000)> <random seed (default: 0)>
```

## Convert the hepmc output to a ROOT file 
- Compile the `Xi_fraction.cc` with g++ (make sure to link the ROOT and HepMC3 libraries correctly)
```Bash
g++ -O2 -std=c++17 Xi_fraction.cc -o Xi_fraction $(root-config --cflags --libs) $(HepMC3-config --cflags --libs)
```
- If successful, you should see an executable named `Xi_fraction`. You can run it with
```Bash
./Xi_fraction <input hepmc file name> <output ROOT file name>
```

## Submit condor jobs
- You can submit condor jobs to run the above two steps for a large number of events. You may need to change the `Initialdir` (to where you unpack and compile the pythia code), `outputDir`, `NeventsPerJob, and the number of jobs in the `Queue` for your needs
```
Universe           = vanilla
Executable         = runPythia.sh
getenv             = True
Initialdir         = /sphenix/user/hjheng/sPHENIXRepo/HF-analysis/simulation/pythia8316/examples
Output             = $(Initialdir)/condorLog/Pythia_ppMinBias_$(Process).out
Error              = $(Initialdir)/condorLog/Pythia_ppMinBias_$(Process).err
Log                = $(Initialdir)/condorLog/Pythia_ppMinBias_$(Process).log
PeriodicHold       = (NumJobStarts>=1 && JobStatus == 1 && !(ON_EVICT_CHECK_RequestMemory_REQUIREMENTS))
request_memory     = 2000MB
retry_request_memory_increase = RequestMemory + 2000
retry_request_memory_max = 10000MB
Priority           = 20
job_lease_duration = 3600
outputDir          = /sphenix/tg/tg01/hf/hjheng/HF-analysis/simulation/Pythia_ppMinBias
outputHepMCFile    = $(outputDir)/pythia_ppMinBias_$(Process).hepmc
NeventsPerJob      = 500000
outputROOTFile     = $(outputDir)/ppMinBias_Xi_fraction_$(Process).root
Arguments          = $(outputHepMCFile) $(NeventsPerJob) $(Process) $(outputROOTFile)
Queue 20
```
- Inspect `runPythia.sh` to make sure everthing is in order. When you are ready to submit jobs, simply do
```Bash
condor_submit submit_standalonePythia.job
```