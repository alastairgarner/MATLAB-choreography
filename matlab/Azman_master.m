%% Azman_master
% Author - Alastair Garner, alastairgarner@outlook.com
clear all; clc;

parameterFile = 'AG_config.yaml'; 
paramFile = dir(fullfile('.','**',[parameterFile,'*']));
cd(paramFile.folder)
params = yaml.ReadYaml(fullfile(paramFile.folder,paramFile.name));
addpath(params.directories.code)

%%
% initialise_folders

%%
if ~isempty(params.directories.raw_data)
    fix_rawdata(params.directories.raw_data,params.directories.choreography_inputs);
end
%%

run_choreography(params.directories.choreography_inputs,params.directories.choreography_results,params);
% run_choreography_select(params.directories.choreography_inputs,params.directories.choreography_results,params);
%%

% extract_outline(pathinp,pathres)

%%

% finalise_data()



