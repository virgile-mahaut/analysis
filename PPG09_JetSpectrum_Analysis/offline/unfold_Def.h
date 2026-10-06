#include <iostream>
#include <TFile.h>
#include <TF1.h>

double calibptbins[] = {19, 24, 30, 36, 43, 51, 60, 70};
double truthptbins[] = {15, 19, 24, 30, 36, 43, 51, 60, 70, 81};

int n_underflowbin = 1;
int n_overflowbin = 1;

int calibnpt = sizeof(calibptbins) / sizeof(calibptbins[0]) - 1;
int truthnpt = sizeof(truthptbins) / sizeof(truthptbins[0]) - 1;

const int nrep = 1000;
