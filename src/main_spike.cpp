#include "io.h"
#include "spike.h"
#include "common_header.h"
#include <boost/program_options.hpp>
#include <xtensor/xarray.hpp>
#include <xtensor/xadapt.hpp>
#include <xtensor/xnpy.hpp>
using namespace std;
namespace po = boost::program_options;

int myrandom(int i) {return rand()%i;}

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
		("ofile,o", po::value<string>(), "[positional] : output TDMI file, as numpy *.npy format")
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
	// convert vector<bool> to xarray<bool>;
	vector<size_t> shape = {binary_spikes.size()};
	xt::xarray<bool, xt::layout_type::row_major> x_binary_spikes(shape);
	for (int i = 0; i < binary_spikes.size(); i ++) {
		if (binary_spikes[i]) x_binary_spikes(i) = true;
		else x_binary_spikes(i) = false;
	}
  // Output spike train;
	xt::dump_npy(dir + ofilename, x_binary_spikes);
  //Print1DBin(dir + ofilename, binary_spikes, "trunc");
	finish = clock();
	// Time counting:
	cout << "[-] spike generation takes " << (finish - start)*1.0 / CLOCKS_PER_SEC << "s" << endl;
  return 0;
}
