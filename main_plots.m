%% main_plots

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% August 2019; Last revision: 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MANUAL OPTIONS

% General options
opt.general.figureFormat = 'pdf'; % Format you want to save figures to
opt.general.Buffer = 15; % Ignore tracking for the first X seconds for each object
opt.general.AverageWindow = 30; % Window (seconds) over which to calculate average values
opt.general.FilterWindow = [nan nan]; % only include data that lies between these time points (seconds).
                                        % [nan nan] uses all available data
opt.general.SmoothenWindow = 21; % Number of frames over which to smoothen data [-10 +10]
opt.general.FilterFunc = @(Start,End) (Start <= 60 & End >= 120);

% 'Pause' figures
opt.pause.Threshold = 0.01; % Minimum speed threshold to determine 'pause' (mm/s)
opt.pause.MinDuration = 0.5; % Minimum duration for 'pause' (seconds)
opt.pause.MinInterval = 0.5; % Merge 'pause' events within this interval

% 'Ridge' figures
opt.ridge.XLimits = [0 180]; % X limits for ridgeline plots

% Boxplot options
opt.box.PlotOrder = 'default';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load parameters

addpath(genpath('code'))

configfile = 'AZ_config';

params = load_config(configfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run plot function
generate_plots(params,opt)

% generate_plots(params,opt,...
%     'FigureType', 'pdf',...
%     'PlotOrder', 'default',...
%     'FixedDuration', 30,...
%     'RidgeLimits', [0 180],...
%     'Buffer',15) % set custom x limits for Ridgeline plots here

fig = figure(1, 'WindowStyle', 'Docked')


