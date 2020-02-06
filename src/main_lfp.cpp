// ***************
//  Copyright: Kyle Chen
//  Author: Kyle Chen
//  Created: 2019-04-06
//  Description: main file for lfp.h and lfp.cpp
//***************
#include "io.h"
#include "lfp.h"
#include <random>
#include <chrono>
#include "cnpy.h"
using namespace std;

//  Function of calculating LFP with point current source model in 1-D loop network case;
int main(int argc, const char* argv[]) {
  // Config program options:
  bool verbose, shuffle_flag;
  po::options_description generic("All Options");
  generic.add_options()
    ("help,h", "show help message")
    ("verbose,v", po::bool_switch(&verbose), "enable verbose")
    ;
  po::options_description io_config;
  io_config.add_options()
    ("prefix", po::value<string>()->default_value("./"), "prefix of output files")
    ("config,c", po::value<string>(), "config file")
    // output file
    ("ofile,o", po::value<string>(), "[positional] : output LFP file, as numpy *.npy format")
    ("shuffle-flag,f", po::bool_switch(&shuffle_flag), "flag for random shuffle")
    ("shuffle-seed", po::value<int>(), "random seed for random shuffle, optional")
    ;
  po::positional_options_description pos_desc;
  pos_desc.add("ofile", 1);
  po::variables_map vm;
  //po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_desc).run(), vm);
  //po::notify(vm);

  po::options_description model_config("Configs");
  model_config.add_options()
    // [electrode]
    ("electrode.x", po::value<double>(), "X coordinates of electrode")
    ("electrode.y", po::value<double>(), "Y coordinates of electrode")
    // [model]
    ("model.order", po::value<int>(), "decay order: \n1 for r^-1 spatial dependence; \n2 for r^-2;")
    // [LFP]
    ("LFP.type", po::value<string>(), "type of neurons contributed in LFP")
    // [neuron]
    ("neuron.g_m", po::value<double>(), "leaky conductance")
    ("neuron.er" , po::value<double>(), "resting potential") 
    ("neuron.ee" , po::value<double>(), "exc reversal potential") 
    ("neuron.ei" , po::value<double>(), "inh reversal potential") 
    ("neuron.v_id" , po::value<int>(), "dynamical id of V")
    ("neuron.ge_id", po::value<int>(), "dynamical id of GE")
    ("neuron.gi_id", po::value<int>(), "dynamical id of GI")
    // [network]
    ("network.size", po::value<int>(), "number of neurons")
    ("network.isspatial", po::value<bool>(), "spatial structure: true for spatially weighted, otherwise not;")
    ("network.file", po::value<string>(), "filename of coordinates file")
    // [time]
    ("time.tmin", po::value<double>(), "lower bound of time range")
    ("time.tmax", po::value<double>(), "upper bound of time range")
    ("time.dt",   po::value<double>(), "size of time step")
    ;
  po::options_description cml_options, config_options;
  cml_options.add(generic).add(io_config).add(model_config);
  config_options.add(model_config);
  po::store(po::command_line_parser(argc, argv).options(cml_options).positional(pos_desc).run(), vm);
//  po::store(po::command_line_parser(argc, argv).options(config).run(), vm);
  po::notify(vm);
  if (vm.count("help")) {
    printf("Usage: cal-lfp [-c option_argument] config ofile\n");
    cout << generic << '\n';
    cout << io_config << '\n';
    cout << model_config << '\n';
    return 1;
  }
  // Loading config.ini:
  ifstream config_file;
  if (vm.count("config")) {
    config_file.open(vm["config"].as<string>().c_str());
    po::store(po::parse_config_file(config_file, cml_options), vm);
    po::notify(vm);
    config_file.close();
  }
  if (verbose) {
    cout << ">> Configs loaded.\n";
  }
  
  string dir = vm["prefix"].as<string>();
  string ofilename = vm["ofile"].as<string>();
  // Analyze listing series;
  int neuron_num = vm["network.size"].as<int>();
  vector<int> list(neuron_num);
  for (int i = 0; i < neuron_num; i ++) list[i] = i;
  if (verbose) {
    printf(">> %d connected neuron contribute to LFP\n", neuron_num);
  }
  fflush(stdout);

  // Calculate the spatial weights
  auto start = chrono::system_clock::now();
  vector<double> spatial_weights;
  CalculateSpatialWeight(vm, spatial_weights);
  auto finish = chrono::system_clock::now();
  chrono::duration<double> escaped_seconds;
  if (verbose) {
    escaped_seconds = finish - start;
    printf("[-] Spatial weight generation : %3.3e s\n", escaped_seconds.count());
  }

  //  Choose objective time range;
  double t_range[2]; // t_range[0] = t_min; t_range[1] = t_max;
  t_range[0] = vm["time.tmin"].as<double>();
  t_range[1] = vm["time.tmax"].as<double>();
  if (verbose) {
    printf(">> Time range = (%.2f, %.2f] ms\n", t_range[0], t_range[1]);
  }
  fflush(stdout);

  // Processing LFP data;
  start = chrono::system_clock::now();
  vector<double> lfp;
  string LFP_type = vm["LFP.type"].as<string>();
  CalculateLFP(dir, vm, lfp, list, LFP_type, spatial_weights, t_range, vm["time.dt"].as<double>());

  finish = chrono::system_clock::now();
  // shuffle flag;
  if (shuffle_flag) {
    mt19937 g;
    if (vm.count("shuffle-seed")) {
      g.seed(vm["shuffle-seed"].as<int>());
    } else {
      random_device rd;
      g.seed(rd());
    }
    shuffle(lfp.begin(), lfp.end(), g);
  }
  //  Output lfp:
  cnpy::npy_save(dir+ofilename, &lfp[0], {lfp.size()},"w");
  // counting time;
  if (verbose) {
    escaped_seconds = finish - start;
    printf("[-] LFP generation : %3.3f s\n", escaped_seconds.count());
  }
  return 0;
}
