%% additional_plots

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load parameters

addpath(genpath('code'))

configfile = 'AZ_config';

params = load_config(configfile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varargin = {'FigureType','pdf','PlotOrder','default','FixedDuration', 30, 'RidgeLimits', []};

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
addOptional(p,'RidgeLimits',[])

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load params, data and filter data

data = choreHandler(params);
data = data.group_by();

for ii = 1:numel(data)
    data(ii) = data(ii).load_data_choreography;
end

data = data.apply_time_filter(opt.tStart,opt.tEnd,opt.tDur);
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
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Optional Plots
if ~isdir(fullfile(fig_dir,"other"))
    mkdir(fullfile(fig_dir,"other"))
end

%%% Number of objects
for ii = 1:numel(data)
    temp = data(ii).data_choreography;
    temp = temp([temp.animal_filter]);
    
    tmax = max([temp.elapstime]);
    et = {temp.elapstime};
    bins = [0:0.25:ceil(tmax)];
    
    binned = cellfun(@(x) unique(discretize(x,bins)), et, 'UniformOutput', false);
    starts = cellfun(@min,binned);
    
    num_objects = accumarray([binned{:}]',1,[numel(bins),1])';
    cum_objects = cumsum(accumarray(starts',1,[numel(bins),1]))';
    
    line(bins,num_objects,...
        'Color','black',...
        'LineWidth',3);
    
    ylabel("Instantaneous Object #")
    
    yyaxis right
    line(bins,cum_objects,...
        'Color',[0.8500 0.3250 0.0980],...
        'LineWidth',3);
    ax = gca;
    ax.YAxis(2).Color = [0.8500 0.3250 0.0980];
    ylabel("Cumulative Object #")
    
    xlabel("Time (s)")
    
    pbaspect([2,1,1])
    xlim([0 tmax])
        
    geno = strcat("tracked@",data(ii).driver,"@",data(ii).effector);
    fig_name = fullfile(fig_dir,'other',geno);
    save_figure(gcf,fig_name,save_type);
    close
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Distribution of 'curve'
tmax = max([temp.elapstime]);
et = {temp.elapstime};
bins = [0:2.5:60];

hold on
for ii = 1:numel(data)
    temp = data(ii).data_choreography;
    vals = [temp.curve_smooth];
    
    % calculate kernel density estimation for the violin
    [density, value] = ksdensity(vals, 'bandwidth', []);
    density = density(value >= min(vals) & value <= max(vals));
    value = value(value >= min(vals) & value <= max(vals));
    value(1) = min(vals);
    value(end) = max(vals);

    % all data is identical
    if min(vals) == max(vals)
        density = 1;
    end
    width = 0.4/max(density);
    pos = ii;
    
    % Boxplot point
    quartiles = quantile(vals,[0.25 0.75]);
    ave = median(vals);
    
    scal = 0.075;
    boxY = [quartiles fliplr(quartiles)];
    boxX = pos + [-1 -1 1 1].*scal;
    
    %%% Plots
    ViolinPlot =  ... % plot color will be overwritten later
        patch([pos+density*width pos-density(end:-1:1)*width], ...
             [value value(end:-1:1)],...
             [1,0,0],'EdgeColor','none',...
             'FaceColor',[1 1 1].*0.5,'FaceAlpha',.6,...
             'DisplayName','violin');
    
    patch(boxX,boxY,[0 0 150]/255,...
        'EdgeColor','none')
    
    line(pos + [-1 1].*scal, [1 1].*ave,...
        'lineWidth', 2,...
        'Color', 'white');
end
hold off

pbaspect([2,1,1])

set(gca,...
    'XTick',[1:numel(data)],...
    'XTickLabels',{data.genotype_display})
ax = gca;
ax.XAxis.TickLabelRotation = 45;
ylabel('Curve (deg.)')

fig_name = fullfile(fig_dir,'other','violin_curve');
save_figure(gcf,fig_name,save_type);
close
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Risse plot
hold on
for ii = 1:numel(data)
    temp = data(ii).data_choreography;
    
    et = {temp.elapstime};
    curve = {temp.curve};
    curve_smooth = cellfun(@(x) movmean(x,21),curve, 'UniformOutput', false);
    
    tmax = max([et{:}]);    
    threshold = [5:0.5:40]';
    
    func = @(x,y) sum([diff(x > threshold,[],2)] > 0,2) ./ range(y);
    bend_rate = cellfun(func, curve_smooth, et, 'UniformOutput', false);
    bend_rate = mean([bend_rate{:}]');
    
    line(threshold,bend_rate,'Color','black') 
end
hold off

pbaspect([2,1,1])
xlim([5 30])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Curve ridgeline
threshold = 10;

for ii = 1:numel(data)
    temp = data(ii).data_choreography;
    temp = temp([temp.animal_filter]);
    lens = cellfun(@numel, {temp.curve});
    [~,I] = sort(lens);
    temp = temp(I);
    
    et = {temp.elapstime};
    curve = {temp.curve};
    curve_smooth = cellfun(@(x) movmean(x,21),curve, 'UniformOutput', false);
    
    tmax = max([et{:}]);
    edges = [0:0.2:ceil(tmax)];
    
    bins = cellfun(@(x) edges(edges>min(x) & edges<max(x)), et, 'UniformOutput', false);
    curve_interp = cellfun(@(x,y,z) interp1(x,y,z), et, curve_smooth, bins, 'UniformOutput', false);
    
    curve_interp_adj = cellfun(@(x) x./8, curve_interp, 'UniformOutput', false);
    
    len_obj = length(curve_smooth);
    hold on
    for jj = 1:len_obj
        Y = curve_interp_adj{jj};
        X = bins{jj};
        y = (-jj)+[Y,0,0];
        x = [X,max(X),min(X)];
        
        line(X,(-jj)+Y,'Color','white', 'LineWidth', 0.5) % background white line
        patch(x,y,[1 1 1].*0.5,'EdgeColor','none','LineWidth',0.05); % grey patch
        
        f = curve_interp{jj} > threshold;
        f2 = f | [f(2:end) 0] | [0 f(1:end-1)];
        
        Y(~f) = 0;
%         bins{jj}(~f) = 0;
        
        y = Y(f2);
        x = X(f2);
        y = (-jj)+[y,0,0];
        x = [x,max(x),min(x)];
        try
            patch(x,y,[0.8 0 0],'EdgeColor','none');
        end
    end
    hold off
    pbaspect([3,1,1])
    
    ylim([-len_obj,0]);
    xlim([0 max(cellfun(@max,et))]);
    set(gca,'YTickLabels',[],'YTick',[]);
    
    geno = strcat("ridge_curve@",data(ii).driver,"@",data(ii).effector);
    fig_name = fullfile(fig_dir,'other',geno);
    save_figure(gcf,fig_name,save_type);
    close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
hold on
for ii = 1:numel(data)
    temp = data(ii).data_choreography;
    
    et = {temp.elapstime};
    curve = {temp.curve};
    curve_smooth = cellfun(@(x) movmean(x,21),curve, 'UniformOutput', false);
    
    tmax = max([et{:}]);
    edges = [0:0.25:ceil(tmax)];
    
    bins = cellfun(@(x) edges(edges>min(x) & edges<max(x)), et, 'UniformOutput', false);
    curve_interp = cellfun(@(x,y,z) interp1(x,y,z), et, curve_smooth, bins, 'UniformOutput', false)
    
    curve_interp = [curve_interp{:}];
    [~,I] = ismember([bins{:}], edges);
    curve_mean = accumarray(I', curve_interp',[numel(edges),1],@mean);
    curve_std = accumarray(I', curve_interp',[numel(edges),1],@std);
    
    patch_x = [edges, fliplr(edges)];
    patch_y = [curve_mean;curve_mean]' + [curve_std',-curve_std'];

    hold on
    cellfun(@(x,y) line(x,y,'Color',[1 1 1]*0.7), et, curve_smooth)
    patch('X',patch_x,'Y',patch_y, 'FaceColor', 'red',...
        'EdgeColor', 'none',...
        'FaceAlpha',0.2);
    line(edges, curve_mean, 'Color','red', 'LineWidth', 2.5)
    hold off
    
    pbaspect([2,1,1])
    xlim([0 tmax])
    
    ylabel('Curve (deg.)')
    xlabel('Time (s)')
    title(data(ii).genotype_display)
    
    geno = strcat("curve@",data(ii).driver,"@",data(ii).effector);
    fig_name = fullfile(fig_dir,'other',geno);
    save_figure(gcf,fig_name,save_type);
    close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:numel(data)
    temp = data(ii).data_choreography;
    f = [temp.animal_filter];
    
    et = {temp(f).elapstime};
    et = unique([et{:}]);
    ar = {temp(f).area};
    ar_norm = cellfun(@(x) x./mean(x), ar, 'UniformOutput', false);
    
    c = num2cell([1:numel(ar_norm)]./numel(ar_norm).*0.7);
    
    hold on
    cellfun(@(y,c) plot(y,'Color',[1 1 1].*c), ar_norm, c)
    hold off
    
    ticks = get(gca,'XTick');
    set(gca,'XTickLabel',et(ticks(2:end)));
    
    ylim([0 2])
    pbaspect([2,1,1])
    
    ylabel('Normalised Area')
    xlabel('Time (s)')
    title(data(ii).genotype_display)
    
    geno = strcat("area@",data(ii).driver,"@",data(ii).effector);
    fig_name = fullfile(fig_dir,'other',geno);
    save_figure(gcf,fig_name,save_type);
    close
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Pause Analyses
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% hold on
% for ii = 1:numel(data)
%     temp = data(ii).data_choreography;
%     
%     j = 4
%     
%     et = temp(j).elapstime;
%     x = temp(j).x;
%     y = temp(j).y
%     xDiff = [0 diff(x)];
%     yDiff = [0 diff(y)];
%     
%     dist = temp(j).distance;
%     
%     ar = mean(temp(j).area);
%     r = sqrt(ar/pi);
%     d = 2*r;
%     
%     dist_prev = sqrt(...
%                     movsum(xDiff,[7 7]).^2 ...
%                     + movsum(yDiff,[7 7]).^2 ...
%                     );
%                 
% %     dist_prev = sqrt(...
% %                     movsum(xDiff,[15 1]).^2 ...
% %                     + movsum(yDiff,[15 1]).^2 ...
% %                     );
% %     dist_prev = dist_prev./d
%     
%     hold on
%     plot(et,dist_prev)
%     plot(et,movmean(dist_prev,21))
%     hold off
%     pbaspect([4,1,1])
%     
%     
%     
%     distByCirc = dist./d
%     
%     cumsum(dist)/d
%     
%     th = 0:pi/50:2*pi;
%     xunit = r * cos(th) + x(1);
%     yunit = r * sin(th) + y(1);
%     
%     frame = 21;
%     
%     hold on
%     plot(x,y)
%     plot(smooth(x,frame,'sgolay',2),smooth(y,frame,'sgolay',2))
%     plot(xunit, yunit, 'Color', 'red');
%     hold off
%     
%     axis equal
% 

