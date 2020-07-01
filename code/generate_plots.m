%% generate_plots.m

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision: 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function generate_plots(params,varargin)
% TO Do
%  - param file input
%  - save type (svg, pdf, m)
%  - custom time frame
%  - plot order (default, custom)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debug settings
% varargin = {'FigureType','pdf','PlotOrder','default','FixedDuration', 30, 'RidgeLimits', [],'Buffer',15};

% Parse Function Inputs

defaultWindow = 21;
defaultFilterFunc = @(Start,End) (Start <= 60 & End >= 120);
defaultPlotOrder = 'default';
% defaultConfig = 'default_config';
defaultFixedDuration = nan;

defaultFigType = 'svg';
validFigType = {'svg','pdf','fig'};
checkFigType = @(x) any(validatestring(x,validFigType));

defaultFilterWindow = [nan nan];
checkFilterWindow = @(x) isnumeric(x) & numel(x)==2;

p = inputParser;
addOptional(p,'PlotOrder',defaultPlotOrder);
addOptional(p,'FilterFunction',defaultFilterFunc);
% addOptional(p,'ConfigFile',defaultConfig);
addOptional(p,'FilterWindow',defaultFilterWindow,checkFilterWindow);
addOptional(p,'SmoothWindowWidth',defaultWindow);
addOptional(p,'FigureType',defaultFigType,checkFigType);
addOptional(p,'FixedDuration',defaultFixedDuration);
addOptional(p,'RidgeLimits',[]);
addOptional(p,'Buffer',5);

parse(p,varargin{:})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set local parameters

save_type = p.Results.FigureType;
time_window = p.Results.FilterWindow;
plot_order = p.Results.PlotOrder;
% config_file = p.Results.ConfigFile;

opt.smooth_window = p.Results.SmoothWindowWidth;
opt.filter_function = p.Results.FilterFunction;
opt.tStart = time_window(1);
opt.tEnd = time_window(2);
opt.tDur = p.Results.FixedDuration;
opt.tBuffer = p.Results.Buffer;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load params, data and filter data

data = choreHandler(params);
data = data.group_by();

for ii = 1:numel(data)
    data(ii) = data(ii).load_data_choreography;
end

data = data.apply_time_filter(opt.tStart,opt.tEnd,opt.tDur,opt.tBuffer);
data = data.process_choreography(opt);

data = data.format_title();
fig_dir = data.init_figure_directory();

if strcmpi(plot_order,'custom')
    genos = [data.genotype_display];
    genos = [string([1:numel(genos)]);genos];
    fprintf('\n\n\n\n\n\n\n\n\n\n\n');
    fprintf('%s - %s \n',genos{:});
    fprintf('\n')
    res = input('Set the genotype order (list of numbers):','s');
    res = strsplit(res,'\s*[ ,]\s*','DelimiterType','RegularExpression');

    params.group_order = double(string(res));
    if any(isnan(plot_order))
        params.group_order = [1:numel(data)];
    end
else
    params.group_order = [1:numel(data)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boxplots

metrics = {"path_distance_smooth",...
    "path_distance_crude",...
    "mean_speed_smooth",...
    "mean_speed_crude",...
    "vect_distance",...
    "path_distance_byarea",...
    "mean_curve_smooth"};

labels = {"Path length (mm) - Smoothened",...
    "Path length (mm) - Raw",...
    "Average speed (mm/s) - Smooth",...
    "Average speed (mm/s) - Raw",...
    "Distance from Origin (mm)",...
    "Path length (mm) - By area",...
    "Average curve (degrees)"};

for ii = 1:length(metrics)
    [printArray,bh,sh] = data.plot_boxplot(metrics{ii},labels{ii},params);
    fig_name = fullfile(fig_dir,metrics{ii});
    save_figure(gcf,fig_name,save_type);
    save_data_csv(printArray,fig_name);
    close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Path Figures

% gridsize = 3;
% for ii = 1:numel(data)
%     fh = data(ii).plot_paths(gridsize);
%     fig_name = fullfile(fig_dir,strcat(data(ii).driver,"@",data(ii).effector,'@path'));
%     save_figure(gcf,fig_name,save_type);
%     close
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Path Figures - Multi-page (all)
for ii = 1:numel(data)
    data(ii).plot_paths_full(opt,1)
    fig_name = fullfile(fig_dir,strcat(data(ii).driver,"@",data(ii).effector,'@full_path.pdf'));
    d = dir('./.temp/*.pdf');
    input_files = strcat({d.folder},filesep,{d.name});
    if isfile(fig_name)
        delete(fig_name);
    end
    % MATLAB's own libraries mess with ghostscript
    setenv('LD_LIBRARY_PATH','') % clears MATALB-based PATH variables
    append_pdfs(fig_name,input_files{:});
    delete(input_files{:})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Path Figures - Multi-page (all)
for ii = 1:numel(data)
    data(ii).plot_paths_full(opt,0)
    fig_name = fullfile(fig_dir,strcat(data(ii).driver,"@",data(ii).effector,'@path.pdf'));
    d = dir('./.temp/*.pdf');
    input_files = strcat({d.folder},filesep,{d.name});
    if isfile(fig_name)
        delete(fig_name);
    end
    % MATLAB's own libraries mess with ghostscript
    setenv('LD_LIBRARY_PATH','') % clears MATALB-based PATH variables
    append_pdfs(fig_name,input_files{:});
    delete(input_files{:})
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Ridgeline Figures

for ii = 1:numel(data)
    data(ii).plot_ridgeline()
    if ~isempty(p.Results.RidgeLimits)
        xlim([p.Results.RidgeLimits])
    end
    geno = strcat(data(ii).driver,"@",data(ii).effector,"@ridge_speed");
    fig_name = fullfile(fig_dir,geno);
    save_figure(gcf,fig_name,save_type);
    close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


