#include "header.h"
#include <string>

void output_fiber(SimulationStruct* sim, float *data, char* output)
{
	ofstream myfile;
	myfile.open(output, ios::app);

	double scale1 = (double)0xFFFFFFFF * (double)sim->number_of_photons;
	if (NORMAL)
	{
		for (int i = 0; i < NUM_OF_DETECTOR; i++)
		{
			myfile << double(data[i] / scale1) << "\t";
		}
	}
	else
	{
		for (int i = 0; i < NUM_OF_DETECTOR; i++)
		{
			myfile << double(data[i] / scale1) << " ";
		}
	}
	myfile << endl;
	myfile.close();
}


int read_mua_mus(SimulationStruct** simulations, char* input) //Wang modified
{
	// parameters to be modified
	unsigned long long number_of_photons = NUMBER_PHOTONS;
	const int n_simulations = NUMBER_SIMULATION;

	int n_layers = NUM_LAYER;                                   //Zhan modified - 5 layers(scalp+skull+CSF+gray matter+white matter); // Wang modified - 4 layers(epi+dermis+sub fat+muscle); double layer, default value = 2
	float medium_n = 1.457;								// float medium_n = 1.33;   // refractive index of medium // YU-modified
														//float medium_n = 1.457; //Wang-modified
	float lower_thickness = 10.0;						// YU-modified
														//float tissue_n = 1.60;                            // refractive index of tissue
														//float g_factor = 0.84;                            // anisotropic					//YU-modified

	float start_weight;
	//float upper_thickness; //YU-modified


	// read the file 
	fstream myfile;
	myfile.open(input);  //Wang modified

	float thickness_1[n_simulations], mua_1[n_simulations], mus_1[n_simulations], n_1[n_simulations], g_1[n_simulations];
	float thickness_2[n_simulations], mua_2[n_simulations], mus_2[n_simulations], n_2[n_simulations], g_2[n_simulations];
	float thickness_3[n_simulations], mua_3[n_simulations], mus_3[n_simulations], n_3[n_simulations], g_3[n_simulations];
	float thickness_4[n_simulations], mua_4[n_simulations], mus_4[n_simulations], n_4[n_simulations], g_4[n_simulations];
	float thickness_5[n_simulations], mua_5[n_simulations], mus_5[n_simulations], n_5[n_simulations], g_5[n_simulations];
	float thickness_6[n_simulations], mua_6[n_simulations], mus_6[n_simulations], n_6[n_simulations], g_6[n_simulations];

	for (int i = 0; i < n_simulations; i++)
		//myfile >> thickness_1[i] >> mua_1[i] >> mus_1[i] >> n_1[i] >> g_1[i] >> thickness_2[i] >> mua_2[i] >> mus_2[i] >> n_2[i] >> g_2[i] >> thickness_3[i] >> mua_3[i] >> mus_3[i] >> n_3[i] >> g_3[i] >> thickness_4[i] >> mua_4[i] >> mus_4[i] >> n_4[i] >> g_4[i] >> thickness_5[i] >> mua_5[i] >> mus_5[i] >> n_5[i] >> g_5[i] >> mua_6[i] >> mus_6[i] >> n_6[i] >> g_6[i];
		myfile >> thickness_1[i] >> mua_1[i] >> mus_1[i] >> n_1[i] >> g_1[i] >> thickness_2[i] >> mua_2[i] >> mus_2[i] >> n_2[i] >> g_2[i] >> thickness_3[i] >> mua_3[i] >> mus_3[i] >> n_3[i] >> g_3[i] >> thickness_4[i] >> mua_4[i] >> mus_4[i] >> n_4[i] >> g_4[i] >> mua_5[i] >> mus_5[i] >> n_5[i] >> g_5[i];
	myfile.close();

	/*fstream myfile;
	myfile.open ("mua_mus.txt");
	float up_mua[n_simulations],up_mus[n_simulations],down_mua[n_simulations],down_mus[n_simulations];
	myfile >> upper_thickness;
	for(int i = 0; i < n_simulations; i++)
	myfile >> up_mua[i] >> up_mus[i] >> down_mua[i] >> down_mus[i];
	myfile.close();*/

	/*fstream file;
	file.open("index.txt");
	float wavelength[n_simulations], tissue_n[n_simulations], medium_n[n_simulations];
	for(int i = 0; i < n_simulations; i++)
	file >> wavelength[i] >> tissue_n[i] >> medium_n[i];

	file.close();*/

	// Allocate memory for the SimulationStruct array
	*simulations = (SimulationStruct*)malloc(sizeof(SimulationStruct)*n_simulations);
	if (*simulations == NULL) { perror("Failed to malloc simulations.\n"); return 0; }//{printf("Failed to malloc simulations.\n");return 0;}

	for (int i = 0; i < n_simulations; i++)
	{
		(*simulations)[i].number_of_photons = number_of_photons;
		(*simulations)[i].n_layers = n_layers;

		// Allocate memory for the layers (including one for the upper and one for the lower)
		(*simulations)[i].layers = (LayerStruct*)malloc(sizeof(LayerStruct)*(n_layers + 2));
		if ((*simulations)[i].layers == NULL) { perror("Failed to malloc layers.\n"); return 0; }//{printf("Failed to malloc simulations.\n");return 0;}

																								 // Set upper refractive index (medium)
		(*simulations)[i].layers[0].n = medium_n;	//(*simulations)[i].layers[0].n = medium_n[i]; //YU-modified

													// Set the parameters of tissue (1st layer)
		(*simulations)[i].layers[1].n = n_1[i];
		(*simulations)[i].layers[1].mua = mua_1[i];
		(*simulations)[i].layers[1].g = g_1[i];			//(*simulations)[i].layers[1].g     = g_factor; //YU-modified 
		(*simulations)[i].layers[1].z_min = 0;
		(*simulations)[i].layers[1].z_max = thickness_1[i];	//(*simulations)[i].layers[1].z_max = upper_thickness; //YU-modified
		(*simulations)[i].layers[1].mutr = 1.0f / (mua_1[i] + mus_1[i]);

		// Set the parameters of tissue (2nd layer)  //Wang modified
		(*simulations)[i].layers[2].n = n_2[i];
		(*simulations)[i].layers[2].mua = mua_2[i];
		(*simulations)[i].layers[2].g = g_2[i];			//(*simulations)[i].layers[2].g     = g_factor; //YU-modified
		(*simulations)[i].layers[2].z_min = thickness_1[i];		//(*simulations)[i].layers[2].z_min = upper_thickness; //YU-modified
		(*simulations)[i].layers[2].z_max = thickness_1[i] + thickness_2[i];		//(*simulations)[i].layers[2].z_max = 1.0;            // set as infinity
		(*simulations)[i].layers[2].mutr = 1.0f / (mua_1[i] + mus_2[i]);

		// Set the parameters of tissue (3rd layer)  //Wang modified
		(*simulations)[i].layers[3].n = n_3[i];
		(*simulations)[i].layers[3].mua = mua_3[i];
		(*simulations)[i].layers[3].g = g_3[i];
		(*simulations)[i].layers[3].z_min = thickness_1[i] + thickness_2[i];
		(*simulations)[i].layers[3].z_max = thickness_1[i] + thickness_2[i] + thickness_3[i];
		(*simulations)[i].layers[3].mutr = 1.0f / (mua_3[i] + mus_3[i]);

		// Set the parameters of tissue (4th layer)  //Zhan modified
		(*simulations)[i].layers[4].n = n_4[i];
		(*simulations)[i].layers[4].mua = mua_4[i];
		(*simulations)[i].layers[4].g = g_4[i];
		(*simulations)[i].layers[4].z_min = thickness_1[i] + thickness_2[i] + thickness_3[i];
		(*simulations)[i].layers[4].z_max = thickness_1[i] + thickness_2[i] + thickness_3[i] + thickness_4[i];
		(*simulations)[i].layers[4].mutr = 1.0f / (mua_4[i] + mus_4[i]);

		// Set the parameters of tissue (5th layer)  //Zhan modified
		(*simulations)[i].layers[5].n = n_5[i];
		(*simulations)[i].layers[5].mua = mua_5[i];
		(*simulations)[i].layers[5].g = g_5[i];
		(*simulations)[i].layers[5].z_min = thickness_1[i] + thickness_2[i] + thickness_3[i] + thickness_4[i];
		//(*simulations)[i].layers[5].z_max = thickness_1[i] + thickness_2[i] + thickness_3[i] + thickness_4[i] + thickness_5[i];
		(*simulations)[i].layers[5].z_max = lower_thickness;
		(*simulations)[i].layers[5].mutr = 1.0f / (mua_5[i] + mus_5[i]);

		/*
		// Set the parameters of tissue (6th layer)  //Zhan modified
		(*simulations)[i].layers[6].n = n_6[i];
		(*simulations)[i].layers[6].mua = mua_6[i];
		(*simulations)[i].layers[6].g = g_6[i];
		(*simulations)[i].layers[6].z_min = thickness_1[i] + thickness_2[i] + thickness_3[i] + thickness_4[i] + thickness_5[i];
		(*simulations)[i].layers[6].z_max = lower_thickness;
		(*simulations)[i].layers[6].mutr = 1.0f / (mua_6[i] + mus_6[i]);
		*/

		// Set lower refractive index (medium)
		(*simulations)[i].layers[n_layers + 1].n = medium_n;		//(*simulations)[i].layers[n_layers+1].n = medium_n[i]; //YU-modified

																	//calculate start_weight
		double n1 = (*simulations)[i].layers[0].n;
		double n2 = (*simulations)[i].layers[1].n;
		double r = (n1 - n2) / (n1 + n2);
		r = r*r;
		start_weight = (unsigned int)((double)0xffffffff * (1 - r));
		//start_weight = 1-r;  
		//printf("Start weight=%u\n",start_weight);
		(*simulations)[i].start_weight = start_weight;
	}
	return n_simulations;
}

void generate_filename(char *filename, int SDS_number)
{
	// make the output file name
	char output[100] = "pathlength_SDS_";
	char postfix[10] = ".txt";
	string SDS_num = to_string(SDS_number);
	int i = 0;
	for (; output[i] != '\0'; i++) {
		filename[i] = output[i];
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
void output_SDS_pathlength(float ***pathlength_weight_arr, int *temp_SDS_detect_num, int SDS_to_output)
{
	if (SDS_to_output != 0) {
		cout << "SDS_to_output= " << SDS_to_output << endl;
	}
	else {
		cout << "SDS_to_output= all" << endl;
	}

	int start_SDS_index, end_SDS_index;
	if (SDS_to_output==0){
		start_SDS_index = 0;
		end_SDS_index = NUM_OF_DETECTOR;
	}
	else {
		start_SDS_index = SDS_to_output - 1;
		end_SDS_index = SDS_to_output;
	}

	for (int s = start_SDS_index; s < end_SDS_index; s++) {
		char output[100];
		generate_filename(output, s + 1);
		//cout << "output to: " << output << ", add " << temp_SDS_detect_num[s] << " photons" << endl;

		// output the pathlength
		ofstream myfile;
		myfile.open(output, ios::app);
		for (int i = 0; i < temp_SDS_detect_num[s]; i++) {
			for (int j = 0; j <= NUM_LAYER + 1; j++) {
				myfile << pathlength_weight_arr[s][j][i] << '\t';
			}
			myfile << endl;
		}
		myfile.close();
		temp_SDS_detect_num[s] = 0;
	}
}

void output_sim_summary(SimulationStruct* sim, int *total_SDS_detect_num)
{
	ofstream myfile;
	myfile.open("sumary.txt", ios::app);

	double scale1 = (double)0xFFFFFFFF * (double)sim->number_of_photons;
	
	myfile << "Photon weight should be divided by " << scale1 << endl << "Detected photon number:" << endl;
	for (int i = 0; i < NUM_OF_DETECTOR; i++) {
		myfile << "SDS " << i << ":\t" << total_SDS_detect_num[i] << endl;
	}
	myfile.close();
}

void output_A_rz(SimulationStruct* sim, unsigned long long *data, char* output)
{
	ofstream myfile;
	myfile.open(output, ios::app);

	double scale1 = (double)0xFFFFFFFF * (double)sim->number_of_photons;
	double scale2; // scale for different r and z

	for (int z = 0; z < record_nz; z++)
	{
		for (int r = 0; r < record_nr; r++)
		{
			// divided by the small grid volume
			scale2 = scale1 * 2 * PI*(r + 0.5)*record_dr*record_dr*record_dz;
			myfile << double(data[z*record_nr + r] / scale2) << "\t";

			// NOT divided by the small grid volume
			//myfile << double(data[z*record_nr + r] / scale1) << "\t";
		}
		myfile << endl;
	}

	myfile.close();
}

void output_A0_z(SimulationStruct* sim, unsigned long long *data, char* output)
{
	ofstream myfile;
	myfile.open(output, ios::app);

	double scale1 = (double)0xFFFFFFFF * (double)sim->number_of_photons;
	double scale2 = scale1 * 2 * PI*(0 + 0.5)*record_dr*record_dr*record_dz; // scale for different r and z

	for (int z = 0; z < record_nz; z++)
	{
		myfile << double(data[z] / scale2) << endl;
	}

	myfile.close();
}