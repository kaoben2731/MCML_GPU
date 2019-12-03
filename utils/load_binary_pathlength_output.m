% Load the MCML output binary file
% Benjamin Kao, 20191203

%% inputs
% sim_set_file: the file path of MCML input sim_setup.json file
% summary_file: the file path of MCML output summary.json file
% SDS_index: which SDS is this pathlength file
% pathlength_file: the file path of MCML output pathlength.bin file

%% outputs
% output: the content of the pathlength file

function output=load_binary_pathlength_output(sim_set_file,summary_file,SDS_index,pathlength_file)
sim_set=jsondecode(fileread(sim_set_file));
sim_sum=jsondecode(fileread(summary_file));

fid=fopen(pathlength_file);
output=fread(fid,[sim_set.number_layers+2 sim_sum.SDS_detected_number(SDS_index)],'float');
fclose all;
output=output';
end