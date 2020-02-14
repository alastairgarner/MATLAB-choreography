%% individual_figures

clear all
clc

initialise_folders

fils = dir(fullfile(pathres,'**'));
fils = fils(~[fils.isdir]);

file_details = get_file_details(fils);

[datesU,ia,ic] = unique([file_details.date]);

[indx,tf] = listdlg('ListString',datesU);

filter = any(ic == indx,2);
fils = fils(filter);

file_details = get_file_details(fils);
genotypes = [file_details.full_genotype];

%% 
[genotypesU ia ic] = unique(genotypes);

delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%f%f%[^\n\r]';

datStruct = [];
for ii = 1:length(genotypesU)
    
    filt = ic == ii;
    tfiles = fils(filt);
    
    f_deets = get_file_details(tfiles(1));
    outdir = fullfile(pathfigs,f_deets.effector,f_deets.driver, strcat(f_deets.date,'_',f_deets.time));
    outname = ['path_fig'];
    outfull = fullfile(outdir,outname);
    
    if ~isdir(outdir)
        mkdir(outdir);
    end
    
    % LOAD DATA
    dStruct = load_choreography_data(tfiles);
    [dStruct.complete]
%     [cellfun(@min, [dStruct.time]);cellfun(@max, [dStruct.time])]

    dStruct = dStruct([dStruct.complete]);
    
    lengths = cellfun(@length, {dStruct.animal});
    diffs = [1 diff([dStruct.animal]) ~= 0];
    unique_animals = cumsum(diffs);
    unique_animals = mat2cell(unique_animals,[1],lengths);
    [dStruct.unique_numb] = unique_animals{:};
    
    dStruct = make_animalStruct(dStruct);
    
    %%
    unique_ani = [dStruct.unique_numb];
    time_elapsed = [dStruct.time];
    obj_speed = [dStruct.speed];

    unique_ani = unique([unique_ani{:}]);
    time_min = cellfun(@min, time_elapsed);
    time_max = cellfun(@max, time_elapsed);
    ave_speed = cellfun(@mean, obj_speed);
    obj_speed_mean = cellfun(@(x) movmean(x,15), obj_speed, 'UniformOutput', false);

    timediff = cellfun(@(x) [0 diff(x)], time_elapsed, 'UniformOutput', false);
    dist = cellfun(@(x,y) x.*y, obj_speed, timediff, 'UniformOutput', false);
    dist_cumu = cellfun(@cumsum, dist, 'UniformOutput', false);
    dist_tot = cellfun(@max, dist_cumu );


    datStruct(ii).genotype = f_deets(1).full_genotype; 
    datStruct(ii).filename = tfiles(1).name;
    datStruct(ii).animal = unique_ani;
    datStruct(ii).tstart = time_min;
    datStruct(ii).tstop = time_max;
    datStruct(ii).tduration = time_max - time_min;
    datStruct(ii).meanspeed = ave_speed;
    datStruct(ii).movmeanspeed = obj_speed_mean;
    datStruct(ii).dist_cumulative = dist_cumu;
    datStruct(ii).dist_total = dist_tot;
    datStruct(ii).time_elapsed = [dStruct.time];
    datStruct(ii).x = [dStruct.x];
    datStruct(ii).y = [dStruct.y];
end


%%

shortStruct = datStruct;

filterType = 2;

maxspeed = 0.4;
for ii = 1:length(shortStruct)
    
    switch filterType
        case 1
            dur = shortStruct(ii).tduration;
            f = find(dur > 170);
            gridsize = 4;
            distspace = 20;
            maxlim = distspace*(gridsize+1);
        case 2
            windo = [60 120];
            f = shortStruct(ii).tstart < windo(1) & shortStruct(ii).tstop > windo(2);
            f = find(f);
            gridsize = 3;
            distspace = 6;
            maxlim = distspace*(gridsize+1);
    end
    
    if length(f) > 9
        f = randsample(f,9);
    end
    nlines = length(f);
    cols = ceil(nlines/3);
    rows = ceil(nlines/cols);
    
    idx = [repmat([1:cols]',rows,1), repelem(flipud([1:rows]'),cols,1)];
    
    cd = [uint8(jet(64)*255) uint8(ones(64,1))].';
    hold on
    n = 1;
    for jj = f
        if filterType == 2
            et = shortStruct(ii).time_elapsed{jj};
            filt = et > windo(1) & et < windo(2);
        else
            filt = repelem(true,1, length( shortStruct(ii).time_elapsed{jj} ));
        end
        x = shortStruct(ii).x{jj}(filt);
        y = shortStruct(ii).y{jj}(filt);
        xy = [x;y]';
%         xy = xy - xy(1,:);
%         et = shortStruct(ii).et{jj};
%         et = et - et(1);
%         etround = round(et)+1;
        xy = xy-mean(xy)+([distspace distspace].*idx(n,:));
        
        bins = linspace(0,maxspeed,64);
        speeds = shortStruct(ii).movmeanspeed{jj}(filt);
        Y = discretize(speeds,bins);
        Y(isnan(Y)) = 64;
        
%         x = [xy(1:end-1,1),xy(2:end,1)];
%         y = [xy(1:end-1,2),xy(2:end,2)];
%         p = line(x',y','LineWidth',2);
%         colcell = num2cell(cd(1:3,Y(2:end))',2);
%         [p.Color] = colcell{:};
%         
%         p = plot(xy(:,1),xy(:,2));
%         pause(0.01)
%         set(p.Edge,'ColorBinding','flat','ColorData',cd(1:4,Y))

        p = scatter(xy(:,1),xy(:,2),5,'.');
        p.CData = cd(1:3,Y)';       
        n = n+1;
    end
    hold off
    
    axis equal
    colormap(jet)
    xlim([0 maxlim])
    ylim([0 maxlim])
    pbaspect([cols,rows,1])
    
    titl = format_title(shortStruct(ii).filename);
%     titl = strrep(shortStruct(ii).genotype,'@','>hEAAT\wedge');
%     titl = strrep(titl,'None','');
    title(titl)
    c = colorbar;
    caxis([0 maxspeed])
    
    set(gca,'XTick',[0 maxlim],'YTick',[0 maxlim]);
    xlabel('Distance (mm)')    
    
    outname = fullfile(pathfigs,strcat(shortStruct(ii).genotype,'_path3_fig'));
    
    
    print(outname,'-dpdf','-painters','-fillpage');
    
    close

end

%%
shortStruct = datStruct;
genotypes = [shortStruct.genotype]

% genoOrder = [3,7,4,5,6];
% genoOrder = [4,8,3,11,7,5,6,9,10];
% genoOrder = [2,5,9,4,12,8,6,7,10,11];
% shortStruct = shortStruct(genoOrder);

filterType = 2;

speeds = []; dists = []; groups = []; lengths = [];
for ii = 1:length(shortStruct)
    switch filterType
        case 1
            dur = shortStruct(ii).tduration;
            f = find(dur > 170);
        case 2
            windo = [60 120];
            f = shortStruct(ii).tstart < windo(1) & shortStruct(ii).tstop > windo(2);
            f = find(f);
    end
    speeds = [speeds shortStruct(ii).meanspeed(f)];
    dists = [dists shortStruct(ii).dist_total(f)];
    groups = [groups repelem(ii,1,length(f))];
    lengths = [lengths length(f)];
end
lengths

genos = format_title({shortStruct.filename});

%%

figure(1)
h = boxplot(speeds,groups,'Colors',[.3 .3 .3],'Symbol','rx','Labels',genos,'LabelOrientation','inline');
hold on
p = plotSpread(speeds,'distributionIdx',groups,'distributionColors',[100/255 100/255 255/255]);
hold off
pbaspect([1,1,1])
ylim([0 0.25])
ylabel('Average speed (mm/s)')

set(h,'LineWidth',2);
obs = findobj(gca,'Marker','.');
set(obs,'MarkerSize',12);
set(gca,'XTickLabel',genos,'XTickLabelRotation',45,'LineWidth',2);

text([1:length(genos)],repelem(0.015,1,length(genos)),num2str(lengths'),'HorizontalAlignment','center');

outname = fullfile(pathfigs,['allgenos_speed_fig']);
print(outname,'-dpdf','-painters','-fillpage');
close


figure(1)
h = boxplot(dists,groups,'Colors',[.3 .3 .3],'Symbol','rx','Labels',genos,'LabelOrientation','inline');
hold on
p = plotSpread(dists,'distributionIdx',groups,'distributionColors',[100/255 100/255 255/255]);
hold off
pbaspect([1,1,1])
ylim([0 40])
ylabel('Distance moved (mm)')

set(h,'LineWidth',2);
obs = findobj(gca,'Marker','.');
set(obs,'MarkerSize',12);
set(gca,'XTickLabel',genos,'XTickLabelRotation',45,'LineWidth',2);

text([1:length(genos)],repelem(2,1,length(genos)),num2str(lengths'),'HorizontalAlignment','center');

outname = fullfile(pathfigs,['allgenos_distance_fig']);
print(outname,'-dpdf','-painters','-fillpage');
close


%% Stats - One-way ANOVA, Tukey-Kramer Correction

metrictotest = speeds;
fname = ['speed_ANOVA_Tukey.csv'];

variableNames = {'Genotype_1','Genotype_2','Lower_Confidence_Interval','Estimate','Upper_Confidence_Interval','P_value'};

[~,~,stats] = anova1(metrictotest,groups,'off');
[c,m,h,gnames] = multcompare(stats);

T = array2table([genos(c(:,1))', genos(c(:,2))', num2cell(c(:,3:6))]);
T.Properties.VariableNames = variableNames;

writetable(T,fname,'Delimiter',',','QuoteStrings',true)
type(fname)

%% Stats - Kruskal-Wallis, Dunn Correction

metrictotest = speeds;
fname = ['speed_Kruskal_Dunn.csv'];

variableNames = {'Genotype_1','Genotype_2','Lower_Confidence_Interval','Estimate','Upper_Confidence_Interval','P_value'};

[p,tbl,stats] = kruskalwallis(metrictotest,groups,'off');
[c,m,h,gnames] = multcompare(stats,'CType','dunn-sidak');

T = array2table([genos(c(:,1))', genos(c(:,2))', num2cell(c(:,3:6))]);
T.Properties.VariableNames = variableNames;

writetable(T,fname,'Delimiter',',','QuoteStrings',true)
type(fname)