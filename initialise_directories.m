%% initialise_directories

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% August 2019; Last revision: 

function initialise_directories()

addpath(genpath('./code'))

parameterFile = 'default_config.yaml'; 
paramFile = dir(fullfile('.','**',[parameterFile,'*']));
params = yaml.ReadYaml(fullfile(paramFile.folder,paramFile.name));

% initialise_folders
temp = rmfield(params.directories,'master');
fields = struct2cell(temp);
for ii = 1:numel(fields)
    if ~isdir(fields{ii})
        mkdir(fields{ii})
    end
end
