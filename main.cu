#include "header.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <string>


void FreeSimulationStruct(SimulationStruct* sim, int n_simulations);
int read_mua_mus(SimulationStruct** simulations, char* sim_input, char* tissue_input);
void DoOneSimulation(SimulationStruct* simulation, int index, char* output, char* fiber1_position); //Wang modified
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

	n_simulations = read_mua_mus(&simulations, argv[1], argv[2]); // read the input file

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
		DoOneSimulation(&simulations[i], i, argv[3], argv[1]); //Wang modified
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