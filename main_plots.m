%% main_plots

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load parameters

addpath(genpath('code'))

configfile = 'AZ_config';

params = load_config(configfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run plot function
generate_plots(params,'FigureType', 'pdf',...
    'PlotOrder', 'default',...
    'FixedDuration', 30,...
    'RidgeLimits', [0 180],...
    'Buffer',15) % set custom x limits for Ridgeline plots here
