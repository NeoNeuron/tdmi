#include "io.h"
#include "spike.h"
#include "common_header.h"
#include <boost/program_options.hpp>
using namespace std;
namespace po = boost::program_options;

int myrandom(int i) {return rand()%i;}

// arguments:
//
// argv[1] = path of the raster data file;
// argv[2] = path of output raster file;
// argv[3] = index of target neuron;
// argv[4] = time range of spikes;
// argv[5] = binning size of binary time series of spike train;
// argv[6] = shuffle flag;
//
int main(int argc, const char* argv[]) {
	clock_t start, finish;
	start = clock();
	// Config program options:
  bool shuffle_flag;
	po::options_description desc("All Options");
	desc.add_options()
		("help,h", "show help message")
		("prefix", po::value<string>()->default_value("./"), "prefix of output files")
		// input files
		("ifile,i", po::value<string>(), "[positional] : input data files")
		// output file
		("ofile,o", po::value<string>(), "[positional] : output TDMI file")
		// index of spike
		("index,I", po::value<int>(), "index of target spike train")
		// time range
		("trange,t", po::value< vector<double> >()->multitoken(), "range of time series")
		("dt", po::value<double>(), "size of time window")
		("shuffle-flag,f", po::bool_switch(&shuffle_flag), "flag for random shuffle")
		;
	po::positional_options_description pos_desc;
	pos_desc.add("ifile", 1);
	pos_desc.add("ofile", 1);
	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_desc).run(), vm);
	po::notify(vm);
	if (vm.count("help")) {
		printf("Usage: cal-spike [-c option_argument] ifile ofile\n");
		cout << desc << '\n';
		return 1;
	}
  // load data inputing arguments;
	string dir = vm["prefix"].as<string>();
	string ifilename = vm["ifile"].as<string>();
	string ofilename  = vm["ofile"].as<string>();
	int id = vm["index"].as<int>();
	vector<double> range = vm["trange"].as< vector<double> >();
	double dt = vm["dt"].as<double>();
  vector<double> spikes;
  Read1D(dir + ifilename, spikes, id, 0);

  // Truncate the spiking series;
  Truncate(spikes, range);
  // Convert double to binary;
  vector<bool> binary_spikes;
  double tmax = range[1] - range[0];
  Spike2Bool(spikes, binary_spikes, tmax, dt);

  // shuffle flag;
	if (shuffle_flag) {
		srand(unsigned(time(0)));
		random_shuffle(binary_spikes.begin(), binary_spikes.end(), myrandom);
	}
  // Output spike train;
  Print1DBin(dir + ofilename, binary_spikes, "trunc");
	finish = clock();
	// Time counting:
	cout << "[-] spike generation takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
