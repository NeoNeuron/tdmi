//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-01-30
//	Description: Define model of local field potential(LFP);
//***************
#ifndef _LFP_H_
#define _LFP_H_
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

struct neuron_type{
	int index;
	bool type;
};

// Calculate spatial dependent weight
// variables_map vm: variable map of parameters;
// vector<double>& spaital_weights: container of spatial weights 
void CalculateSpatialWeight(po::variables_map & vm, vector<double> & spatial_weights);

//	Local field potential model [version 0.11]
//	Description: point current source model sptial distribution;
//	string dir: prefix of data files
//	variables_map& vm: variable map of neuronal model parameter;
//	vector<double>& lfp: container of generated data;
//	vector<int>& neuron_list: list of indices of attributed neurons;
//	string LFP_type: type of LFP
//  vector<double>& spatial_weights: spatial weights of attributed neurons;
//  double* t_range: time period used in calculation, with unit ms, include the last point while not the first point
//	double sampling_dt: sampling time step;
//	Return: none;
void CalculateLFP(string dir, po::variables_map & vm, vector<double>& lfp, vector<int>& neuron_list, string LFP_type, vector<double>& spatial_weights, double* t_range, double sampling_dt);

#endif // _IFNET_LFP_H_
