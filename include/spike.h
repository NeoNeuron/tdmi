#ifndef _IFNET_SPIKE_H_
#define _IFNET_SPIKE_H_
#include <string>
#include <vector>
using namespace std;

//	Convert double spike train to binary sequence;
//	VECTOR<DOUBLE> spikes: original spike train;
//	VECTOR<BOOL> binary_spikes: binary sequence; true for having spike; false for no spike;
//	DOUBLE tmax: maximum time of time range;
//	DOUBLE dt: the size of time step that a binary value considered;
//	Return: none;
void Spike2Bool(vector<double> &spikes, vector<bool> &binary_spikes, double tmax, double dt);

// Truncate spikes within selected time range;
// VECTOR<DOUBLE> spikes: spike train;
// DOUBLE* range[2]: time range with unit millisecond;
// Return: None;
void Truncate(vector<double>& spikes, vector<double> &range);

#endif
