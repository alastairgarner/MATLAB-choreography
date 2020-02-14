%% initialise_folders

if ~exist('datadir','var')
    mpath = mfilename('fullpath');
%     sp = split(mpath,{'/','\'} );
%     datadir = fullfile(sp{1:end-1});    
    datadir = fileparts(mpath);
end

pathfunc = fullfile(datadir,'matlab');
pathraw = fullfile(datadir,'raw_input');
pathinp = fullfile(datadir,'choreography_input');
pathres = fullfile(datadir,'choreography_results');
pathfigs = fullfile(datadir,'figures');
pathfin = fullfile(datadir,'data_final');
pathchore = fullfile(datadir,'packages','choreography');
% pathtemp = fullfile(pwd,'work_mem');

addpath(pathfunc);
addpath(genpath(pathfunc));