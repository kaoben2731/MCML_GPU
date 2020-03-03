%{
Load the MCML output A_rz and A0_z binary file

Benjamin Kao
Last update: 2020/02/23
%}

function output=fun_load_A_arr(input_dir,SDS_num)

%% param
n_dr=1000;
n_dz=500;

%% main

%% load A_rz
fid=fopen(fullfile(input_dir,['A_rz_SDS_' num2str(SDS_num) '.bin']));
A_rz=fread(fid,[n_dr,n_dz],'float');
fclose all;

%% load A0_z
fid=fopen(fullfile(input_dir,['A0_z_SDS_' num2str(SDS_num) '.bin']));
A0_z=fread(fid,[n_dz,1],'float');
fclose all;

output=A_rz';
output(:,1)=output(:,1)+A0_z;

end
