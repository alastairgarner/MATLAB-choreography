%% Azman_plots


%%
cd("/home/alastair/OneDrive/projects/vanMeyel_pipeline")
clear all; clc;

params.group_order = [6,7,11,8,2,1,5,3,4,9,10];

params.smooth_window = 21;
params.filter_function = @(tStart,tEnd,tDur) (tStart <= 60 & tEnd >= 120);
params.tStart = 60;
params.tEnd = 120;
params.tDur = 60;

addpath('./matlab')
data = choreHandler('AZ_config.yaml');
data = data.group_by();

for ii = 1:numel(data)
    data(ii) = data(ii).load_data_choreography;
end

data = data.apply_time_filter(params.tStart,params.tEnd,params.tDur);
data = data.process_choreography(params);

data = data.format_title();
fig_dir = data.init_figure_directory();


%% Boxplots

metrics = {"path_distance_smooth",...
    "path_distance_crude",...
    "mean_speed_smooth",...
    "mean_speed_crude",...
    "vect_distance",...
    "path_distance_byarea"};

labels = {"Path length (mm) - Smoothened",...
    "Path length (mm) - Raw",...
    "Average speed (mm/s) - Smooth",...
    "Average speed (mm/s) - Raw",...
    "Distance from Origin (mm)",...
    "Path length (mm) - By area"};

for ii = 1:length(metrics)
    [bh,sh] = data.plot_boxplot(metrics{ii},params);
    ax = gca;
    ax.YLabel.String = labels{ii};
    fig_name = fullfile(fig_dir,strcat(metrics{ii},".pdf"));
    print(fig_name,'-dpdf','-painters','-fillpage');
    close
end

%% Path figure

gridsize = 3;
for ii = 1:numel(data)
    fh = data(ii).plot_paths(gridsize);
    geno = strcat(data(ii).driver,"@",data(ii).effector,"@fig_path.pdf");
    fig_name = fullfile(fig_dir,geno);
    print(fig_name,'-dpdf','-painters','-fillpage');
    close
end

%% Ridgeline Plot

for ii = 1:numel(data)
    if ~any([data(ii).data_choreography.animal_filter])
        continue
    end
    filter = find([data(ii).data_choreography(:).animal_filter]);
    speed = {data(ii).data_choreography(filter).speed};
    et = {data(ii).data_choreography(filter).elapstime};
    filt = {data(ii).data_choreography(filter).time_filter};
    speeds = cellfun(@(x,y) x(y),speed,filt,'UniformOutput', false);
    et = cellfun(@(x,y) x(y),et,filt,'UniformOutput', false);
    
    speeds_norm = cellfun(@(x) movmean(x,21)./0.25, speeds, 'UniformOutput', false);
    
    len_obj = length(speeds_norm);
    cmap = flipud(viridis(len_obj));
%     cmap = cmocean('ice',len_obj)
    hold on
    for jj = 1:len_obj
        y = (-jj)+[speeds_norm{jj},0,0];
%         y = [speeds_norm{jj}+(jj-1),(jj-1)+zeros(1,len_vect)];
        x = [et{jj},max(et{jj}),min(et{jj})];
%         x = [et{jj},fliplr(et{jj})];
        patch(x,y,'red','EdgeColor','white','LineWidth',0.1,'FaceColor',cmap(jj,:));
    end
    hold off
    set(gcf,'PaperOrientation','landscape');
    pbaspect([3,1,1])
    
    ylim([-len_obj,0]);
    set(gca,'YTickLabels',[],'YTick',[]);
    
    geno = strcat(data(ii).driver,"@",data(ii).effector,"@fig_speeds.pdf");
    fig_name = fullfile(fig_dir,geno);
    print(fig_name,'-dpdf','-painters','-fillpage');
    close
end

%% Eliminate jitter

scaler = .5;
n = 1;
filter = find([data(1).data_choreography(:).animal_filter]);
for ii = 1:10
    temp = data(1).data_choreography(filter(ii));
    [pathlength,idx] = choreHandler.get_distance_by_area(temp,scaler);

    plot(temp.x,temp.y,'k.',temp.x_smooth,temp.y_smooth,'b.');
    viscircles([temp.x_smooth(idx);temp.y_smooth(idx)]',repelem(sqrt(scaler*mean(temp.area)/pi),length(idx))');
    pbaspect([1,1,1])
    ylim(mean(ylim)+[-4 4])
    xlim(mean(xlim)+[-4 4])
    file_name = strcat("by_area_demo_",num2str(filter(ii)));
    print(file_name,'-dpdf','-painters','-fillpage');
    close
    n = n+1;
end

scaler = .3;
n = 1;
filter = find([data(1).data_choreography(:).animal_filter]);
for ii = 1:10
    temp = data(1).data_choreography(filter(ii));
    [pathlength,idx] = choreHandler.get_distance_by_area(temp,scaler);

    plot(temp.x,temp.y,'k.',temp.x_smooth,temp.y_smooth,'b.',temp.x_smooth(idx),temp.y_smooth(idx),'r','LineWidth',1.5);
    pbaspect([1,1,1])
    ylim(mean(ylim)+[-4 4])
    xlim(mean(xlim)+[-4 4])
    file_name = strcat("by_area_line_",num2str(filter(ii)));
    print(file_name,'-dpdf','-painters','-fillpage');
    close
    n = n+1;
end

scal = fliplr([0:.05:1]);
stem(scal,sqrt(scal/pi))

A = 1;
B = 0.3183./A;
pi*A*B

r = 10;
A = sqrt(3*(r^2));
B = A/3;
[r^2,A*B]

[.865*r,A/2]

r^2 = A*B



%%

cmocean

for ii = 1:numel(data)
    filter = find([data(ii).data_choreography(:).animal_filter]);
    speed = {data(ii).data_choreography(filter).speed};
    et = {data(ii).data_choreography(filter).elapstime};
    filt = {data(ii).data_choreography(filter).time_filter};
    speeds = cellfun(@(x,y) x(y),speed,filt,'UniformOutput', false);
    et = cellfun(@(x,y) x(y),et,filt,'UniformOutput', false);
    
    speeds_norm = cellfun(@(x) movmean(x,21), speeds, 'UniformOutput', false);
    
    hold on
    for jj = 1:10
        x = speeds_norm{jj};
        y = fft(x);     
        m = abs(y);
        plot(m)
    end
    hold off
    xlim([0 100])
    
    
    winsize = 54;
    f = linspace(0,3*pi,winsize)
    sinwv = (sin(f)+1)/2

    for jj = 1:5
        x = speeds_norm{jj}./.3;
        [c,lags] = xcorr(x,sinwv);
        figure()
        hold on
        plot(lags,c)
        plot(x)
        xlim([0 2500])
    end
    hold off
    
    corrvect = [];
    for jj = winsize:length(x)
        C = corrcoef(sinwv,x(jj-winsize+1:jj));
        corr(sinwv', x(jj-winsize+1:jj)')
        corrvect = [corrvect C(2,1)];
    end
    plot([zeros(1,round(winsize/2)),corrvect],'b')
    hold on 
    plot(x,'r')
    hold off

    for t = windowSize+1:N
        C = corrcoef(dataMatrix(t-windowSize:t, :));
        idx = setdiff(1:M, [indexColumn]);
        correlationTS(t, :) = C(indexColumn, idx);
    end
end

%% Violin plots

genos = [data.genotype_display];
metric_to_plot = "path_distance_byarea";

dists = arrayfun(@(x) [x.data_choreography([x.data_choreography.animal_filter]).(metric_to_plot)], data,'UniformOutput', false);
groups = arrayfun(@(x) [x.data_choreography([x.data_choreography.animal_filter]).group], data,'UniformOutput', false);
dists = [dists{:}];
groups = [groups{:}];
geno_n = accumarray(groups',1);

go = [6,7,11,8,2,1,5,3,4,9,10];

figure()
boxplot(dists,groups,"GroupOrder",string(go))


%%

for ii = 1:100
    temp = data(1).data_choreography(ii)
    try
        x = temp.x_smooth(1);
        y = temp.y_smooth(1);
        area = mean(temp.area)
        plot(temp.x,temp.y,'k.',temp.x_smooth,temp.y_smooth,'r.')
        viscircles([x,y],sqrt(area/pi))
        pause(0.5)
    catch
    end
end


%%
cond{1} = @(x,y) sum(x > 30);
cond{2} = @(x,y) sum(x > 60);
cond{3} = @(x,y) sum(y(:,1) <= 60 & y(:,2) >= 120);
cond{4} = @(x,y) sum(y(:,1) >= 30 & y(:,1) <= 90 & x > 60);

x = time_tracked;
y = startfin;
for ii = 1:length(cond)
    out{ii} = cellfun(cond{ii}, x, y);
end
plot(vertcat(out{:})','.', 'MarkerSize', 20)

%%

temp = data(1).data_choreography
vertcat()
filters = cellfun(@(x) x.data_choreography.)
data(1).data_choreography.path_distance
startfin = arrayfun(@(x) vertcat(x.data_choreography.startfin), data,'UniformOutput', false);

arrayfun(@(x,y) sum(y(:,1) <= 60 & y(:,2) >= 120), , 'UniformOutput', false)

ii = 1;
temp = data(1).data_choreography


plot(temp.elapstime,temp.x,'k',...
    temp.elapstime,temp.x_smooth,'r')

figure()
plot(temp.x,temp.y,'k.',...
    temp.x_smooth,temp.y_smooth,'r.')
pbaspect([1,1,1])


ii = 1;
temp = data(1).data_choreography(ii)

[mi,ma] = bounds(temp.elapstime);


%%
% plot speed timeseries for each animal
% plot path (colour by speed)

initialise_folders
datadir;

fils = dir(fullfile(pathres,'**'));
fils = fils(~[fils.isdir]);
filt = contains({fils.name},'speed');
fils = fils(filt);

parts = arrayfun(@(x) split(x.name,["@","."]), fils, 'UniformOutput', false);
tstamps = cellfun(@(x) x{1}, parts, 'UniformOutput', false);
genos = cellfun(@(x) x{2}, parts, 'UniformOutput', false);

[genosU ia ic] = unique(genos);

delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%f%f%[^\n\r]';


for ii = 1:length(genos)
    
    datmat = [];
    filename = fullfile(fils(ii).folder,fils(ii).name);
    titlename = format_title(filename);
    sp = split(fils(ii).name,["@","."]);
    tstamp = sp{1};
    eff = sp{3};
    driv = sp{2};
    feat = sp{end-1};
    outdir = fullfile(pathfigs,eff,driv,tstamp);
    outname = [feat,'_fig'];
    outfull = fullfile(outdir,outname);
    
    if ~isdir(outdir)
        mkdir(outdir);
    end
    
    fileID = fopen(filename,'r');
    dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
    fclose(fileID);

    %datmat = vertcat(datmat,[dA{2:4}]);
    datmat = [dA{2:4}];
    
    mint = min(datmat(:,2));
    maxt = max(datmat(:,2));
    
    [uniani aa ac] = unique(datmat(:,1));
    length(uniani);
    
    n = 1;
    for jj = 1:length(uniani)
        filt = ac == jj;
        datshort = datmat(filt,:);
        
        subplot(6,1,n);
        p = plot(datshort(:,2),datshort(:,3));
        xlim([mint,maxt])
        ylim([0,1.5])
        text(5,1,string(uniani(jj)))
        
        mult = jj/6;
        if floor(mult)==mult
            n = 1;
            print(strcat(outfull,'_',string(ceil(mult))),'-dpdf','-painters','-fillpage')
            close;
            pr = 0;
        else
            pr = 1;
            n = n+1;
        end
    end
    
    if pr == 1
        print(strcat(outfull,'_',string(ceil(mult))),'-dpdf','-painters','-fillpage');
        close
    end
    
end


%%


fils = dir(fullfile(pathres,'**'));
fils = fils(~[fils.isdir]);
filt = contains({fils.name},{'.x.','.y.'});
fils = fils(filt);

parts = arrayfun(@(x) split(x.name,["@","."]), fils, 'UniformOutput', false);
tstamps = cellfun(@(x) x{1}, parts, 'UniformOutput', false);
genos = cellfun(@(x) x{2}, parts, 'UniformOutput', false);

[genosU ia ic] = unique(genos);

delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%f%f%[^\n\r]';

for ii = 1:length(genosU)
    
    filt = ic == ii;
    filt2 = contains({fils.name},{'.x.','.y.'});
    tfiles = fils(filt'&filt2);
    
    sp = split(tfiles(1).name,["@","."]);
    tstamp = sp{1};
    eff = sp{3};
    driv = sp{2};
    feat = sp{end-1};
    outdir = fullfile(pathfigs,eff,driv,tstamp);
    outname = ['path_fig'];
    outfull = fullfile(outdir,outname);
    
    if ~isdir(outdir)
        mkdir(outdir);
    end
    
    datmat = [];
    for jj = 1:2
        filename = fullfile(tfiles(jj).folder,tfiles(jj).name);

        fileID = fopen(filename,'r');
        dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
        fclose(fileID);

        %datmat = vertcat(datmat,[dA{2:4}]);
        if jj == 1
            datmat = [dA{2:4}];
        else
            datmat = [datmat dA{4}];
        end
    end
        
    mint = min(datmat(:,2));
    maxt = max(datmat(:,2));
    
    [uniani aa ac] = unique(datmat(:,1));
    anis = length(uniani);
    
    n = 1;
    for jj = 1:anis
        f1 = ac == jj;
        x = datmat(f1,4);
        y = datmat(f1,3);
        x = x-x(1);
        y = y-y(1);
        
        subplot(4,3,n)
        p = plot(x,y);
        xlim([-14 14]);
        ylim([-14 14]);
        text(-12,12,string(uniani(jj)),'FontSize',6)
        pbaspect([1 1 1]);
        
        mult = jj/12;
        if floor(mult)==mult
            n = 1;
            print(strcat(outfull,'_',string(ceil(mult))),'-dpdf','-painters','-fillpage');
            close;
            pr = 0;
        else
            pr = 1;
            n = n+1;
        end
        
    end
    if pr == 1
        print(strcat(outfull,'_',string(ceil(mult))),'-dpdf','-painters','-fillpage');
        close
    end
end
    
   
    

