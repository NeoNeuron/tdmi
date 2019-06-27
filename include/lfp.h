//***************
//	Copyright: Kyle Chen
//	Author: Kyle Chen
//	Date: 2018-01-30
//	Description: Define model of local field potential(LFP);
//***************
#ifndef _LFP_H_
#define _LFP_H_
#include "common_header.h"
#include <boost/program_options.hpp>
namespace po = boost::program_options;

using namespace std;

struct neuron_type{
	int index;
	bool type;
};

void CalculateSpatialWeight(po::variables_map & vm, vector<double> & spatial_weights);

//	Local field potential model [version 0.11]
//	Description: point current source model without sptial distribution;
//	DOUBLE* t_range: time period used in calculation, with unit ms, include the last point while not the first point
//	STRING current_file: total membrane current;
//	VECTOR<DOUBLE> lfp: local field potential data;
//	Return: none;
void CalculateLFP(string dir, po::variables_map & vm, vector<double>& lfp, vector<int>& neuron_list, string LFP_type, vector<double>& spatial_weights, double* t_range, double sampling_dt);

#endif // _IFNET_LFP_H_
