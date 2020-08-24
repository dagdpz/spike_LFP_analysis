function [pop keys] = ph_load_population(folder,filename)
% folder='W:\Projects\Pulv_eye_gaze_position\ephys\ephys_analysis_v3_July2016_coordinates';
% filename='population_mem_';

allfiles=dir([folder filesep filename '*']);
allfiles=allfiles(~[allfiles.isdir]);
n_population=0;
for f=1:numel(allfiles)
   load([folder filesep allfiles(f).name])
   %temporary for noise correlation
   [spike_field.session]=deal(f);
   pop(n_population+1:n_population+numel(spike_field))=spike_field;
   n_population=numel(pop);
end

if isempty(allfiles)   
pop=struct(); 
end
end