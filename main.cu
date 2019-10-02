#include "header.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string>


void FreeSimulationStruct(SimulationStruct* sim, int n_simulations);
int read_mua_mus(SimulationStruct** simulations, char* sim_input, char* tissue_input);
void DoOneSimulation(SimulationStruct* simulation, int index, char* output, bool do_replay, bool output_each_pathlength, bool do_output_average_pathlength, bool do_output_bin);
void show_usage(string name);


int main(int argc, char* argv[])
{
	SimulationStruct* simulations;
	int n_simulations;
	unsigned long long seed = (unsigned long long) time(NULL); // Default, use time(NULL) as seed

	if (argc < 4) {
		show_usage(argv[0]);
		return 0;
	}

	bool do_replay = false;
	bool output_each_pathlength = false;
	bool do_output_average_pathlength = false;
	bool do_output_bin = false;

	for (int i = 1; i < argc; i++) {
		if (string(argv[i]) == "-h") {
			show_usage(argv[0]);
			return 0;
		}
		else if (string(argv[i]) == "-R") {
			do_replay = true;
		}
		else if (string(argv[i]) == "-P") {
			output_each_pathlength = true;
		}
		else if (string(argv[i]) == "-A") {
			do_output_average_pathlength = true;
		}
		else if (string(argv[i]) == "-B") {
			do_output_bin = true;
		}
	}
	if (do_replay && !output_each_pathlength) {
		do_output_average_pathlength = true;
	}

	if (output_each_pathlength && !do_replay) {
		cout << "-P option only work with -R option!\n";
		return 0;
	}
	if (do_output_average_pathlength && !do_replay) {
		cout << "-A option only work with -R option!\n";
		return 0;
	}
	if (do_output_bin && !output_each_pathlength) {
		cout << "-B option only work with -P option!\n";
		return 0;
	}

	n_simulations = read_mua_mus(&simulations, argv[1], argv[2]); // read the input file

	// debug
	for (int s = 1; s <= simulations[0].num_detector; s++)
	{
		cout << simulations[0].detInfo[s].raduis << '\t';
	}
	cout << endl;
	system("PAUSE");

	if (n_simulations == 0)
	{
		printf("Something wrong with read_simulation_data!\n");
		return 1;
	}
	else
	{
		//printf("Successfully read data!\n");
	}

	clock_t time1, time2;

	// Start the clock
	time1 = clock();

	//perform all the simulations
	for (int i = 0; i < n_simulations; i++)
	{
		// Run a simulation
		printf("simulating %d\n", i);
		DoOneSimulation(&simulations[i], i, argv[3], do_replay, output_each_pathlength, do_output_average_pathlength, do_output_bin); //Wang modified
	}

	time2 = clock();
	printf("Simulation time: %.2f sec\n", (double)(time2 - time1) / CLOCKS_PER_SEC);

	FreeSimulationStruct(simulations, n_simulations);

	//system("PAUSE");
	return 0;
}

void show_usage(string name)
{
	size_t found = name.find_last_of("/\\");
	string exename = name.substr(found + 1);
	cerr << "Usage: " << exename << " sim_set.json input.txt output.txt <option(s)>\n"
		<< "Options:\n"
		<< "\t-h,--help\t Show this help message\n"
		<< "\t-R,\t\t Replay the detected photon after first simulation\n"
		// << "\t-W,\t\t Doing white Monte Carlo\n"
		<< "\t-P,\t\t Output the pathlength for each photon, otherwise output the calculated average pathlength\n"
		<< "\t-A,\t\t Calaulate and output the average pathlength\n"
		<< "\t-B,\t\t Output the pathlength file in binary format\n"
		<< endl;
}

void FreeSimulationStruct(SimulationStruct* sim, int n_simulations)
{
	for (int i = 0; i < n_simulations; i++) {
		free(sim[i].layers);
		free(sim[i].detInfo);
		free(sim[i].critical_arr);
	}
	free(sim);
}