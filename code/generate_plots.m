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
function generate_plots(params,opt,varargin)
% TO Do
%  - param file input
%  - save type (svg, pdf, m)
%  - custom time frame
%  - plot order (default, custom)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Debug settings
% varargin = {'FigureType','pdf','PlotOrder','default','FixedDuration', 30, 'RidgeLimits', [],'Buffer',15};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parse Function Inputs

% defaultWindow = 21;
% defaultFilterFunc = @(Start,End) (Start <= 60 & End >= 120);
% defaultPlotOrder = 'default';
% % defaultConfig = 'default_config';
% defaultFixedDuration = nan;
% 
% defaultFigType = 'svg';
% validFigType = {'svg','pdf','fig'};
% checkFigType = @(x) any(validatestring(x,validFigType));
% 
% defaultFilterWindow = [nan nan];
% checkFilterWindow = @(x) isnumeric(x) & numel(x)==2;
% 
% p = inputParser;
% addOptional(p,'PlotOrder',defaultPlotOrder);
% addOptional(p,'FilterFunction',defaultFilterFunc);
% % addOptional(p,'ConfigFile',defaultConfig);
% addOptional(p,'FilterWindow',defaultFilterWindow,checkFilterWindow);
% addOptional(p,'SmoothWindowWidth',defaultWindow);
% addOptional(p,'FigureType',defaultFigType,checkFigType);
% addOptional(p,'FixedDuration',defaultFixedDuration);
% addOptional(p,'RidgeLimits',[]);
% addOptional(p,'Buffer',5);
% 
% parse(p,varargin{:})

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set local parameters

% save_type = p.Results.FigureType;
% time_window = p.Results.FilterWindow;
% plot_order = p.Results.PlotOrder;

save_type = opt.general.figureFormat;
time_window = opt.general.FilterWindow;
plot_order = opt.box.PlotOrder;

opt.smooth_window = opt.general.SmoothenWindow;
opt.filter_function = opt.general.FilterFunc;
opt.tStart = time_window(1);
opt.tEnd = time_window(2);
opt.tDur = opt.general.AverageWindow;
opt.tBuffer = opt.general.Buffer;

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

ridgeLimits = opt.ridge.XLimits;
for ii = 1:numel(data)
    data(ii).plot_ridgeline()
    if ~isempty(ridgeLimits)
        xlim([ridgeLimits])
    end
    geno = strcat(data(ii).driver,"@",data(ii).effector,"@ridge_speed");
    fig_name = fullfile(fig_dir,geno);
    save_figure(gcf,fig_name,save_type);
    close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% 

n = 10;
cmap = cbrewer('qual', 'Paired', n);
cmap = repmat(cmap,2,1);

LineStyles = {'-','--',':'};
hold on
for ii = 1:numel(data)
    ixStyle = ceil(ii/n);
    
    temp = data(ii).data_choreography;
    temp = temp([temp.animal_filter]);
    
    vals = [temp.speedPerFrame];
    
    threshold = [0:100] ./ 100 .* 0.5;
    value = mean(vals > threshold',2);
    line(threshold,value,'Color',cmap(ii,:),'LineWidth', 2,'LineStyle', LineStyles{ixStyle});
end
hold off

legend({data.genotype_display},'Location','eastoutside')
xlim([0,0.3])
xlabel('Speed threshold (mm/s)')
ylabel('Proportion of time above threshold')
pbaspect([1,1,1])

geno = strcat("speed_distribution");
fig_name = fullfile(fig_dir,geno);
save_figure(gcf,fig_name,save_type);
close


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Pause Analyses

lim = vertcat(data.data_choreography);
lim = ceil(max([lim.elapstime]));

printArray = [];
for ii = 1:numel(data)
    [time,id] = data(ii).below_threshold(opt,'speedPerFrame');
    filt = [data(ii).data_choreography.animal_filter];
    objNumber = [data(ii).data_choreography(filt).aniID];
    [tMin, tMax] = cellfun(@(x) bounds(x), {data(ii).data_choreography(filt).elapstime});
    idTot = [1:numel(tMin)] .* [1;1];
    
    printArray = [printArray, [ones(1,numel(id(1,:))).*ii; objNumber(id(1,:)); time] ];
    
    %%% Generate Ethogram
    LineWidth = 3;
    
    hold on
    plot([tMin; tMax],idTot,...
        'Color',[1 1 1].*0.7,...
        'LineWidth',LineWidth);
    
    plot(time,id,...
        'Color','red',...
        'LineWidth',LineWidth);
    hold off
    
    ylim([0 max(idTot(1,:))+1])
    xlim([0 lim])
    pbaspect([1,1,1])
    set(gca,'YTick',[])
    xlabel('Time (s)')
    title(data(ii).genotype_display)
    
    geno = strcat("pause_etho@",data(ii).driver,"@",data(ii).effector);
    fig_name = fullfile(fig_dir,geno);
    save_figure(gcf,fig_name,save_type);
    close
    
    %%% Generate proportion timeseries
    bins = [1:0.05:max(tMax)];
        
    paused = [time(1,:)' < bins] - [time(2,:)' < bins];
    counts = accumarray(id(1,:)',1);
    pausedCells = mat2cell(paused,counts,size(paused,2));
    pausedCells = cellfun(@(x) sum(x,1), pausedCells, 'UniformOutput', false);
    paused = vertcat(pausedCells{:});
    
    value = mean(paused,1);
    err = std(paused,1);
    
    errX = [bins, fliplr(bins)];
    errY = [value, fliplr(value)] + [err, -fliplr(err)];
    
    hold on
    p = patch(errX,errY, [1 1 1].*0.75,...
        'EdgeColor', 'none');
    l = line(bins,value, 'Color','black', 'LineWidth', 2);
    hold off
    
    pbaspect([2,1,1])
    ylim([0 1])
    xlim([0 lim])
    ylabel('Pause proportion')
    xlabel('Time (s)')
    title(data(ii).genotype_display)
    
    geno = strcat("pause_prop@",data(ii).driver,"@",data(ii).effector);
    fig_name = fullfile(fig_dir,geno);
    save_figure(gcf,fig_name,save_type);
    close 
end

locs = find([1 diff(printArray(1,:))]);
labs = repelem({""}, 1, size(printArray,2));
labs(locs) = {data.genotype_display};

pA = [labs; num2cell(printArray)];
fid = fopen(fullfile(fig_dir,'pause_data.csv'),'w');
fprintf(fid,'%s,%s,%s,%s,%s\n',string({'genotype','group','id','start','stop'}));
fprintf(fid,'%s,%d,%d,%f,%f\n',pA{:});
fclose(fid);
