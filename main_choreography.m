%% main_choreography

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% August 2019; Last revision: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('code'))

configfile = 'AZ_config';

params = load_config(configfile);

% Fix and copy across raw data to chore_input
if ~isempty(params.directories.mwt)
    fix_rawdata(params);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_choreography(params);





