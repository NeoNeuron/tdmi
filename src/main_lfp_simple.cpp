// ***************
//  Copyright: Kyle Chen
//  Author: Kyle Chen
//  Date: 2019-04-06
//  Description: main file for lfp.h and lfp.cpp
//***************
#include "io.h"
#include "lfp.h"
#include "cnpy.h"
using namespace std;

//  Function of calculating LFP with point current source model in 1-D loop network case;
int main(int argc, const char* argv[]) {
  clock_t start, finish;
  start = clock();
  // Config program options:
  bool verbose;
  po::options_description desc("All Options");
  desc.add_options()
    ("help,h", "show help message")
    ("verbose,v", po::bool_switch(&verbose), "enable verbose")
    ("prefix", po::value<string>()->default_value("./"), "prefix of output files")
    ("config,c", po::value<string>(), "[positional] : config file")
    ("ofile,o", po::value<string>(), "[positional] : output LFP file")
    ("index,I", po::value<vector<int> >()->multitoken(), "[positional] : index(indices) of target neurons")
    ;
  po::positional_options_description pos_desc;
  pos_desc.add("config", 1);
  pos_desc.add("ofile", 1);
  pos_desc.add("index", -1);
  po::variables_map vm;
  po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_desc).run(), vm);
  po::notify(vm);

  po::options_description config("Configs");
  config.add_options()
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
    ("network.file", po::value<string>(), "path of coordinates file")
    // [time]
    ("time.tmin", po::value<double>(), "lower bound of time range")
    ("time.tmax", po::value<double>(), "upper bound of time range")
    ("time.dt",   po::value<double>(), "size of time step")
    ;
  if (vm.count("help")) {
    printf("Usage: cal-lfp-simple [-c option_argument] config ofile index[id1, id2, ...]\n");
    cout << desc << '\n';
    cout << config << '\n';
    return 1;
  }
  // Loading config.ini:
  ifstream config_file;
  if (vm.count("config")) {
    config_file.open(vm["config"].as<string>().c_str());
  } else {
    cout << "ERROR : lack of config file\n";
    return -1;
  }
  po::store(po::parse_config_file(config_file, config), vm);
  po::notify(vm);
  if (verbose) {
    cout << ">> Configs loaded.\n";
  }
  
  string dir = vm["prefix"].as<string>();
  string ofilename = vm["ofile"].as<string>();
  // Analyze listing series;
  int neuron_num = vm["network.size"].as<int>();
  vector<int> list = vm["index"].as<vector<int> >();
  if (verbose) {
    printf(">> %d connected neuron contribute to LFP\n", (int)list.size());
  }
  fflush(stdout);

  // Define uniform spatial weights
  vector<double> spatial_weights(neuron_num, 1.0);

  //  Choose objective time range;
  double t_range[2]; // t_range[0] = t_min; t_range[1] = t_max;
  t_range[0] = vm["time.tmin"].as<double>();
  t_range[1] = vm["time.tmax"].as<double>();
  if (verbose) {
    printf(">> Time range = (%.2f, %.2f] ms\n", t_range[0], t_range[1]);
  }
  fflush(stdout);

  // Processing LFP data;
  vector<double> lfp;
  string LFP_type = vm["LFP.type"].as<string>();
  CalculateLFP(dir, vm, lfp, list, LFP_type, spatial_weights, t_range, vm["time.dt"].as<double>());

  //  Output lfp:
  cnpy::npy_save(dir+ofilename, &lfp[0], {lfp.size()},"w");
  finish = clock();
  // counting time;
  if (verbose) {
    printf("[-] LFP generation : %3.3f s\n", (finish - start)*1.0 / CLOCKS_PER_SEC);
  }
  return 0;
}
