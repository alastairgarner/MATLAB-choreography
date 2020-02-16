%% initialise_directories

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision: 


function initialise_directories()

%
parameterFile = 'default_config.yaml'; 
paramFile = dir(fullfile('.','**',[parameterFile,'*']));
params = yaml.ReadYaml(fullfile(paramFile.folder,paramFile.name));
cd(params.directories.master)
addpath(params.directories.code)

% initialise_folders
temp = rmfield(params.directories,'master');
fields = struct2cell(temp);
for ii = 1:numel(fields)
    if ~isdir(fields{ii})
        mkdir(fields{ii})
    end
end
