//***************
//  Copyright: Kyle Chen
//  Author: Kyle Chen
//  Date: 2018-01-30
//  Description: source file of lfp.h
//***************
#include "lfp.h"
#include "io.h"
#include "common_header.h"
using namespace std;

bool comp(const int x, const int y) {
  return x < y;
}

inline double L2(vector<double>& a, vector<double>& b) {
  return sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]));
}

void CalculateSpatialWeight(po::variables_map &vm, vector<double> &spatial_weights) {
  bool is_spatial = vm["network.isspatial"].as<bool>();
  spatial_weights.clear();
  int neuron_num = vm["network.size"].as<int>();
  spatial_weights.resize(neuron_num, 1.0);
  if (is_spatial) {
    // Calculate the spatial weights of selected neurons:
    vector<double> electrode_pos = {vm["electrode.x"].as<double>(), vm["electrode.y"].as<double>()};
    // read the coordinate file;
    vector<vector<double> > coordinates;
    Read2D(vm["prefix"].as<string>() + vm["network.file"].as<string>(), coordinates);
    double distance;
    int decay_order = vm["model.order"].as<int>();
    if (decay_order == 1) {
      for (int i = 0; i < neuron_num; i ++) {
        distance = L2(electrode_pos, coordinates[i]);
        spatial_weights[i] = 1.0 / distance;
      }
    } else if (decay_order == 2) {
      for (int i = 0; i < neuron_num; i ++) {
        distance = L2(electrode_pos, coordinates[i]);
        spatial_weights[i] = 1.0 / distance / distance;
      }
    } else {
      cout << "Not proper decay order\n";
    }
  }
  return;
}

void CalculateLFP(string dir, po::variables_map &vm, vector<double> 
    &lfp, vector<int>& neuron_list, string LFP_type, vector<double>
    &spatial_weights, double* t_range, double sampling_dt) {
  // preliminary parameters;
  double sampling_rate = 1 / sampling_dt; // Unit ms;

  // Preparing time series;
  size_t t_begin = t_range[0] * sampling_rate; // not included
  size_t t_end = t_range[1] * sampling_rate; // included
  size_t size_of_lfp = t_end - t_begin;
  lfp.clear();
  lfp.resize(size_of_lfp);

  // Load neuron model;
  double g_m =    vm["neuron.g_m" ].as<double>();
  double V_rest = vm["neuron.er"  ].as<double>();
  double V_e =    vm["neuron.ee"  ].as<double>();
  double V_i =    vm["neuron.ei"  ].as<double>();
  int    v_id =   vm["neuron.v_id"  ].as<int>();
  int    ge_id =  vm["neuron.ge_id" ].as<int>();
  int    gi_id =  vm["neuron.gi_id" ].as<int>();

  ifstream V_file, GE_file, GI_file;
  size_t shape[2];
  V_file.open(dir + "/V.bin", ios::binary);
  V_file.read((char*)&shape, 2*sizeof(size_t));
  vector<double> buffer_vec;
  if (t_begin == 0) buffer_vec.resize(shape[1]);
  else buffer_vec.resize(t_begin*shape[1]);
  // Preparing diff_list:
  vector<int> diff_list(neuron_list.size() + 1);
  // Sort neuron list:
  if (neuron_list.size() == 1) {
    diff_list[0] = neuron_list[0];
    diff_list[1] = shape[1] - neuron_list[0] - 1;
  } else {
    sort(neuron_list.begin(), neuron_list.end(), comp);
    for (int i = 0; i < diff_list.size(); i ++) {
      if (i == 0) diff_list[i] = neuron_list[0];
      else if (i == diff_list.size() - 1) diff_list[i] = shape[1] - neuron_list[i - 1] - 1;
      else diff_list[i] = neuron_list[i] - neuron_list[i - 1] - 1;
    }
  }

  // Classify the LFP type and load other dym data file;
  if (LFP_type == "tot") {
    GE_file.open(dir + "/GE.bin", ios::binary);
    GI_file.open(dir + "/GI.bin", ios::binary);
    GE_file.read((char*)buffer_vec.data(), 2*sizeof(size_t));
    GI_file.read((char*)buffer_vec.data(), 2*sizeof(size_t));

    // For t = [0, t_begin];
    V_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
    GE_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
    GI_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
    // For t = (t_begin, t_end]
    double dym_val[3] = {0.0, 0.0, 0.0};
    double temp_lfp;
    for (size_t i = t_begin; i < t_end; i++) {
      temp_lfp = 0;
      if (neuron_list.size() == shape[1]) {
        for (int j = 0; j < neuron_list.size(); j ++) {
          V_file.read((char*)&dym_val[0], sizeof(double));
          GE_file.read((char*)&dym_val[1], sizeof(double));
          GI_file.read((char*)&dym_val[2], sizeof(double));
          temp_lfp += (- g_m * (dym_val[v_id] - V_rest) - dym_val[ge_id] * (dym_val[v_id] - V_e) - dym_val[gi_id] * (dym_val[v_id] - V_i)) * spatial_weights[j];
        }
        lfp[i - t_begin] = temp_lfp / neuron_list.size();
      } else {
        for (int j = 0; j < diff_list.size() - 1; j ++) {
          V_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
          GE_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
          GI_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
          V_file.read((char*)&dym_val[0], sizeof(double));
          GE_file.read((char*)&dym_val[1], sizeof(double));
          GI_file.read((char*)&dym_val[2], sizeof(double));
          temp_lfp += (- g_m * (dym_val[v_id] - V_rest) - dym_val[ge_id] * (dym_val[v_id] - V_e) - dym_val[gi_id] * (dym_val[v_id] - V_i)) * spatial_weights[j];
        }
        V_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
        GE_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
        GI_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
        lfp[i - t_begin] = temp_lfp / neuron_list.size();
      }
    }
    GE_file.close();
    GI_file.close();
  } else if (LFP_type == "lek") {
    // For t = [0, t_begin];
    V_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
    // For t = (t_begin, t_end]
    double dym_val[3] = {0.0, 0.0, 0.0};
    double temp_lfp;
    for (size_t i = t_begin; i < t_end; i++) {
      temp_lfp = 0;
      if (neuron_list.size() == shape[1]) {
        for (int j = 0; j < neuron_list.size(); j ++) {
          V_file.read((char*)&dym_val[0], sizeof(double));
          temp_lfp += (- g_m * (dym_val[v_id] - V_rest)) * spatial_weights[j];
        }
        lfp[i - t_begin] = temp_lfp / neuron_list.size();
      } else {
        for (int j = 0; j < diff_list.size() - 1; j ++) {
          V_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
          V_file.read((char*)&dym_val[0], sizeof(double));
          temp_lfp += (- g_m * (dym_val[v_id] - V_rest)) * spatial_weights[j];
        }
        V_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
        lfp[i - t_begin] = temp_lfp / neuron_list.size();
      }
    }
  } else if (LFP_type == "exi") {
    GE_file.open(dir + "/GE.bin", ios::binary);
    GE_file.read((char*)buffer_vec.data(), 2*sizeof(size_t));

    // For t = [0, t_begin];
    V_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
    GE_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
    // For t = (t_begin, t_end]
    double dym_val[3] = {0.0, 0.0, 0.0};
    double temp_lfp;
    for (size_t i = t_begin; i < t_end; i++) {
      temp_lfp = 0;
      if (neuron_list.size() == shape[1]) {
        for (int j = 0; j < neuron_list.size(); j ++) {
          V_file.read((char*)&dym_val[0], sizeof(double));
          GE_file.read((char*)&dym_val[1], sizeof(double));
          temp_lfp += (- dym_val[ge_id] * (dym_val[v_id] - V_e)) * spatial_weights[j];
        }
        lfp[i - t_begin] = temp_lfp / neuron_list.size();
      } else {
        for (int j = 0; j < diff_list.size() - 1; j ++) {
          V_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
          GE_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
          V_file.read((char*)&dym_val[0], sizeof(double));
          GE_file.read((char*)&dym_val[1], sizeof(double));
          temp_lfp += (- dym_val[ge_id] * (dym_val[v_id] - V_e)) * spatial_weights[j];
        }
        V_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
        GE_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
        lfp[i - t_begin] = temp_lfp / neuron_list.size();
      }
    }
    GE_file.close();
  } else if (LFP_type == "inh") {
    GI_file.open(dir + "/GI.bin", ios::binary);
    GI_file.read((char*)buffer_vec.data(), 2*sizeof(size_t));

    // For t = [0, t_begin];
    V_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
    GI_file.read((char*)buffer_vec.data(), shape[1]*t_begin*sizeof(double));
    // For t = (t_begin, t_end]
    double dym_val[3] = {0.0, 0.0, 0.0};
    double temp_lfp;
    for (size_t i = t_begin; i < t_end; i++) {
      temp_lfp = 0;
      if (neuron_list.size() == shape[1]) {
        for (int j = 0; j < neuron_list.size(); j ++) {
          V_file.read((char*)&dym_val[0], sizeof(double));
          GI_file.read((char*)&dym_val[2], sizeof(double));
          temp_lfp += (- dym_val[gi_id] * (dym_val[v_id] - V_i)) * spatial_weights[j];
        }
        lfp[i - t_begin] = temp_lfp / neuron_list.size();
      } else {
        for (int j = 0; j < diff_list.size() - 1; j ++) {
          V_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
          GI_file.read((char*)buffer_vec.data(), diff_list[j]*sizeof(double));
          V_file.read((char*)&dym_val[0], sizeof(double));
          GI_file.read((char*)&dym_val[2], sizeof(double));
          temp_lfp += (- dym_val[gi_id] * (dym_val[v_id] - V_i)) * spatial_weights[j];
        }
        V_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
        GI_file.read((char*)buffer_vec.data(), *(diff_list.end() - 1)*sizeof(double));
        lfp[i - t_begin] = temp_lfp / neuron_list.size();
      }
    }
    GI_file.close();
  } else throw runtime_error("ERROR: wrong LFP type (Accessible types: tot, lek, exi, inh.)");
  V_file.close();
}
