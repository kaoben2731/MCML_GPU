#include "header.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include "json.hpp"

using json = nlohmann::json;

void output_fiber(SimulationStruct* sim, float *data, char* output)
{
	ofstream myfile;
	myfile.open(output, ios::app);

	double scale1 = (double)0xFFFFFFFF * (double)sim->number_of_photons;
	if (NORMAL)
	{
		for (int i = 0; i < sim->num_detector; i++)
		{
			myfile << double(data[i] / scale1) << "\t";
		}
	}
	else
	{
		for (int i = 0; i < sim->num_detector; i++)
		{
			myfile << double(data[i] / scale1) << " ";
		}
	}
	myfile << endl;
	myfile.close();
}

int read_mua_mus(SimulationStruct** simulations, char* sim_input, char* tissue_input) //Wang modified
{
	// load the simulation parameters and form into json format
	ifstream inFile;
	inFile.open(sim_input);

	stringstream simStrStream;
	simStrStream << inFile.rdbuf();
	string simStr = simStrStream.str();

	json sim_input_struct = json::parse(simStr);

	// read and set the parameters from json file
	unsigned long long number_of_photons = (unsigned long long)sim_input_struct["number_photons"];
	const int n_simulations = sim_input_struct["number_simulation"];
	//const int n_simulations =1;

	double detector_reflectance = sim_input_struct["detector_reflectance"];
	if (detector_reflectance < 0 || detector_reflectance>1) {
		cout << "Detector reflectance: " << detector_reflectance << " is out of range !\n";
		return 0;
	}

	int n_layers = sim_input_struct["number_layers"];
	if (n_layers > PRESET_NUM_LAYER) {
		cout << "Number of layer is too large!\nThis program only allow number lower than " << PRESET_NUM_LAYER << " !\n";
		return 0;
	}

	float upper_n = sim_input_struct["upper_n"];
	float lower_n = sim_input_struct["lower_n"];

	int num_detector = sim_input_struct["probes"]["num_SDS"];
	if (num_detector > PRESET_NUM_DETECTOR) {
		cout << "Number of SDS is too large!\nThis program only allow number lower than " << PRESET_NUM_DETECTOR << " !\n";
		return 0;
	}
	if (num_detector != sim_input_struct["probes"]["detectors"].size()) {
		cout << "Number of SDS : "<<num_detector <<" not match the probe setting : "<< sim_input_struct["probes"]["detectors"].size() <<" !\n";
		return 0;
	}

	float lower_thickness = 10.0;

	float start_weight;

	// init the parameter array
	// mua for i-th simulation, j-th layer
	float ** thickness = new float*[n_simulations];
	float ** mua = new float*[n_simulations];
	float ** mus = new float*[n_simulations];
	float ** n = new float*[n_simulations];
	float ** g = new float*[n_simulations];
	for(int i=0;i<n_simulations;i++)
	{
		thickness[i] = new float[n_layers];
		mua[i] = new float[n_layers];
		mus[i] = new float[n_layers];
		n[i] = new float[n_layers];
		g[i] = new float[n_layers];
		for (int j = 0; j < n_layers; j++) {
			thickness[i][j] = 0;
			mua[i][j] = 0;
			mus[i][j] = 0;
			n[i][j] = 0;
			g[i][j] = 0;
		}
	}

	// read tissue parameter
	fstream myfile;
	myfile.open(tissue_input);  //Wang modified

	for (int i = 0; i < n_simulations; i++) {
		// for upper layers
		int j = 0; // layer index
		while ( j < n_layers - 1) {
			myfile >> thickness[i][j] >> mua[i][j] >> mus[i][j] >> n[i][j] >> g[i][j];
			j++;
		}
		// for the last layer
		myfile >> mua[i][j] >> mus[i][j] >> n[i][j] >> g[i][j];
		
		// check for input correctness
		if (mus[i][j] <= 0) {
			cout << "Mus for simulation " << i << " layer " << j << " Error with value " << mus[i][j] << " !\n";
			return 0;
		}
		if (n[i][j] <= 1) {
			cout << "n for simulation " << i << " layer " << j << " Error with value " << n[i][j] << " !\n";
			return 0;
		}
		if (g[i][j] > 1) {
			cout << "g for simulation " << i << " layer " << j << " Error with value " << g[i][j] << " !\n";
			return 0;
		}
	}
	myfile.close();

	// Allocate memory for the SimulationStruct array
	*simulations = (SimulationStruct*)malloc(sizeof(SimulationStruct)*n_simulations);
	if (*simulations == NULL) { perror("Failed to malloc simulations.\n"); return 0; }//{printf("Failed to malloc simulations.\n");return 0;}

	for (int i = 0; i < n_simulations; i++)
	{
		(*simulations)[i].number_of_photons = number_of_photons;
		(*simulations)[i].num_layers = n_layers;
		(*simulations)[i].num_detector = num_detector;
		(*simulations)[i].detector_reflectance = detector_reflectance;

		// read probe setting
		(*simulations)[i].detInfo = new DetectorInfoStruct[num_detector + 1];
		// for source
		(*simulations)[i].detInfo[0].NA = sim_input_struct["probes"]["source"]["NA"];
		(*simulations)[i].detInfo[0].raduis = sim_input_struct["probes"]["source"]["radius"];
		(*simulations)[i].detInfo[0].position = 0;
		(*simulations)[i].detInfo[0].angle = ANGLE*PI / 180;
		// for detector
		for (int d = 0; d < num_detector; d++) {
			(*simulations)[i].detInfo[d + 1].NA = sim_input_struct["probes"]["detectors"][d]["NA"];
			(*simulations)[i].detInfo[d + 1].raduis = sim_input_struct["probes"]["detectors"][d]["radius"];
			(*simulations)[i].detInfo[d + 1].position = sim_input_struct["probes"]["detectors"][d]["pos"];
			(*simulations)[i].detInfo[d + 1].angle = ANGLE*PI / 180;
		}
		// presetting critical angle for detectors
		(*simulations)->critical_arr = new float[num_detector + 1];
		for (int d = 1; d <= num_detector; d++) {
			(*simulations)->critical_arr[d] = asin((*simulations)[i].detInfo[d].NA / n_detector);
		}

		// Allocate memory for the layers (including one for the upper and one for the lower)
		(*simulations)[i].layers = (LayerStruct*)malloc(sizeof(LayerStruct)*(n_layers + 2));
		if ((*simulations)[i].layers == NULL) { perror("Failed to malloc layers.\n"); return 0; }//{printf("Failed to malloc simulations.\n");return 0;}

		// Set upper refractive index (medium)
		(*simulations)[i].layers[0].n = upper_n;

		float z_min_cumulate = 0;
		float z_max_cumulate = 0;

		// for setting layers
		for (int l = 1; l <= n_layers; l++) {
			z_max_cumulate += thickness[i][l-1];
			(*simulations)[i].layers[l].n = n[i][l - 1];
			(*simulations)[i].layers[l].mua = mua[i][l - 1];
			(*simulations)[i].layers[l].g = g[i][l - 1];
			(*simulations)[i].layers[l].z_min = z_min_cumulate;
			(*simulations)[i].layers[l].z_max = z_max_cumulate;
			(*simulations)[i].layers[l].mutr = 1.0f / (mua[i][l - 1] + mus[i][l - 1]);
			z_min_cumulate += thickness[i][l - 1];
			//cout << "layer " << l << ", mua= " << (*simulations)[i].layers[l].mua << ", z_min=" << (*simulations)[i].layers[l].z_min << ", z_max=" << (*simulations)[i].layers[l].z_max << endl;
			//system("pause");
		}
		// set the depth of the lower layer
		if ((*simulations)[i].layers[n_layers].z_max < lower_thickness) {
			(*simulations)[i].layers[n_layers].z_max = lower_thickness;
		}

		// Set lower refractive index (medium)
		(*simulations)[i].layers[n_layers + 1].n = lower_n;

		//calculate start_weight
		double n1 = (*simulations)[i].layers[0].n;
		double n2 = (*simulations)[i].layers[1].n;
		double r = (n1 - n2) / (n1 + n2);
		r = r*r;
		start_weight = (unsigned int)((double)0xffffffff * (1 - r));
		//printf("Start weight=%u\n",start_weight);
		(*simulations)[i].start_weight = start_weight;
	}

	// free the memory
	for (int i = 0; i<n_simulations; i++)
	{
		delete[] thickness[i];
		delete[] mua[i];
		delete[] mus[i];
		delete[] n[i];
		delete[] g[i];
	}
	delete[] thickness;
	delete[] mua;
	delete[] mus;
	delete[] n;
	delete[] g;

	return n_simulations;
}

void generate_filename(char *filename, char* prefix,int SDS_number)
{
	// make the output file name
	//char prefix[100] = "pathlength_SDS_";
	char postfix[10] = ".txt";
	string SDS_num = to_string(SDS_number);
	int i = 0;
	for (; prefix[i] != '\0'; i++) {
		filename[i] = prefix[i];
	}
	for (int j = 0; j < SDS_num.length(); j++) {
		filename[i] = SDS_num[j];
		i++;
	}
	for (int j = 0; postfix[j] != '\0'; j++) {
		filename[i] = postfix[j];
		i++;
	}
	filename[i] = '\0';
}

// SDS_to_output: exact number of SDS to output, 1 for SDS1 , 0 for output all SDS
void output_SDS_pathlength(SimulationStruct* simulation, float ***pathlength_weight_arr, int *temp_SDS_detect_num, int SDS_to_output)
{
	/*if (SDS_to_output != 0) {
		cout << "SDS_to_output= " << SDS_to_output << endl;
	}
	else {
		cout << "SDS_to_output= all" << endl;
	}*/

	int start_SDS_index, end_SDS_index;
	if (SDS_to_output==0){
		start_SDS_index = 0;
		end_SDS_index = simulation->num_detector;
	}
	else {
		start_SDS_index = SDS_to_output - 1;
		end_SDS_index = SDS_to_output;
	}

	for (int s = start_SDS_index; s < end_SDS_index; s++) {
		char output[100];
		generate_filename(output, "pathlength_SDS_", s + 1);
		//cout << "output to: " << output << ", add " << temp_SDS_detect_num[s] << " photons" << endl;

		// output the pathlength
		ofstream myfile;
		myfile.open(output, ios::app);
		for (int i = 0; i < temp_SDS_detect_num[s]; i++) {
			for (int j = 0; j <= simulation->num_layers + 1; j++) {
				myfile << pathlength_weight_arr[s][i][j] << '\t';
			}
			myfile << endl;
		}
		myfile.close();
		temp_SDS_detect_num[s] = 0;
	}
}

void output_sim_summary(SimulationStruct* simulation, SummaryStruct sumStruc, bool do_replay)
{
	double scale1 = (double)0xFFFFFFFF * (double)sumStruc.number_of_photons;

	ofstream myfile;
	myfile.open("summary.json", ios::app);
	myfile << "{\n";
	myfile << "\"sim_time\": " << (double)(sumStruc.time2 - sumStruc.time1) / CLOCKS_PER_SEC << ",\n";
	if (do_replay) {
		myfile << "\"replay_time\": " << (double)(sumStruc.time3 - sumStruc.time2) / CLOCKS_PER_SEC << ",\n";
	}
	myfile << "\"each_photon_weight\": " << scale1 << ",\n";
	myfile << "\"SDS_detected_number\": [";
	for (int d = 0; d < simulation->num_detector-1; d++) {
		myfile << " " << sumStruc.total_SDS_detect_num[d] << " ,";
	}
	myfile << " " << sumStruc.total_SDS_detect_num[simulation->num_detector - 1] << " ]\n";
	myfile << "}";
	myfile.close();
}

void output_A_rz(SimulationStruct* sim, unsigned long long *data)
{

	double scale1 = (double)0xFFFFFFFF * (double)sim->number_of_photons;
	double scale2; // scale for different r and z
	
	for (int s = 0; s < sim->num_detector; s++) {
		char output[100];
		generate_filename(output, "A_rz_SDS_", s + 1);

		ofstream myfile;
		myfile.open(output, ios::app);

		for (int z = 0; z < record_nz; z++)
		{
			for (int r = 0; r < record_nr; r++)
			{
				// divided by the small grid volume
				scale2 = scale1 * 2 * PI*(r + 0.5)*record_dr*record_dr*record_dz;
				myfile << double(data[s*record_nr*record_nz + z*record_nr + r] / scale2) << "\t";

				// NOT divided by the small grid volume
				//myfile << double(data[z*record_nr + r] / scale1) << "\t";
			}
			myfile << endl;
		}
		myfile.close();
	}
}

void output_A0_z(SimulationStruct* sim, unsigned long long *data)
{
	double scale1 = (double)0xFFFFFFFF * (double)sim->number_of_photons;
	double scale2 = scale1 * 2 * PI*(0 + 0.5)*record_dr*record_dr*record_dz; // scale for different r and z

	for (int s = 0; s < sim->num_detector; s++) {
		char output[100];
		generate_filename(output, "A0_z_SDS_", s + 1);

		ofstream myfile;
		myfile.open(output, ios::app);

		for (int z = 0; z < record_nz; z++)
		{
			myfile << double(data[s*record_nz + z] / scale2) << endl;
		}
		myfile.close();
	}
}