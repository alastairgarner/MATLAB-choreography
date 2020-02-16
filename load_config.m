%% load_config.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = load_config(filename)

if nargin == 0
    filename = 'default_config.yaml';
end

parameterFile = filename;
[~,parameterFile,~] = fileparts(parameterFile);
paramFile = dir(fullfile('.','**',[parameterFile,'*']));

if isempty(paramFile)
    fprintf('config file not found \n')
    return
end

params = yaml.ReadYaml(fullfile(paramFile.folder,paramFile.name));
cd(params.directories.master)
addpath(params.directories.code)

end
