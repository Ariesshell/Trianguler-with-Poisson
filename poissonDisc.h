#pragma once
#include <vector>

typedef std::pair<float, float> P;
std::vector<P> sampling(double h, double w, double r, unsigned char* data, int numSamplesBeforeReject = 30);