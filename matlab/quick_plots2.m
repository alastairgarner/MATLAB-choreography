%% quick_plots2

clear all
clc

initialise_folders

date_filt = [20190328 20190615];

fils = dir(fullfile(pathres,'\*\*\*\*'));
fils = fils(~[fils.isdir]);

% dStruct = load_choreography_data(fils,[1;2])
% filt = contains({fils.name},{'.x.','.y.'});
% fils = fils(filt);

expr = '\d\d\d\d\d\d\d\d';
tstamps = regexp([fils.name],expr,'match');
tstamps = str2num(vertcat(tstamps{:}));
date_filt = tstamps > date_filt(1) & tstamps < date_filt(2);
fils = fils(date_filt);

parts = arrayfun(@(x) split(x.name,["@","."]), fils, 'UniformOutput', false);
tstamps = cellfun(@(x) x{1}, parts, 'UniformOutput', false);
genos = cellfun(@(x) [x{2} '@' x{3}], parts, 'UniformOutput', false);

[genosU ia ic] = unique(genos);

delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%f%f%[^\n\r]';

for ii = 1:length(genosU)
    
    filt = ic == ii;
%     filt2 = contains({fils.name},{'.x.','.y.'});
    filt2 = contains({fils.name},{'.speed.'});
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
    
    datCat = {};
    tstamps = {};
    for jj = 1:length(tfiles)
        filename = fullfile(tfiles(jj).folder,tfiles(jj).name);
        fileID = fopen(filename,'r');
        dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
        fclose(fileID);
        datCat = [datCat dA];        
    end
    
    datCat = reshape(datCat,length(dA),[],1)';
    datmat = cell2mat(datCat(:,2:4));
    uninum = cumsum([1; diff(datmat(:,1))]~=0);
    datmat = [uninum datmat];
    
    mint = min(datmat(:,3));
    maxt = max(datmat(:,3));
    
%     nums = datmat(:,1);
%     diffs = [1 [abs(diff(nums))>0]' ]';
%     diffs(1) = 1;
%     nums = cumsum(diffs)';
    
    [uniani aa ac] = unique(datmat(:,1));
    
    time_min = accumarray(ac,datmat(:,3),[],@min);
    time_max = accumarray(ac,datmat(:,3),[],@max);
    ave_speed = accumarray(ac,datmat(:,4),[],@mean);
    
    sps = accumarray(ac,1);
    speeds = mat2cell(datmat(:,4),sps,1)';
    speeds = cellfun(@(x) movmean(x,15), speeds, 'UniformOutput', false);
    
    timediff = [0 diff(datmat(:,3))']';
    f = abs(timediff) > 0.15;
    timediff(f) = 0;
    dist = timediff.*datmat(:,4);
    
    total_dist = accumarray(ac,dist);
    
    datTab = [uniani time_min time_max ave_speed total_dist];
    
    
    datStruct(ii).genotype = strcat(driv,'@',eff); 
    datStruct(ii).filename = tfiles(1).name;
    datStruct(ii).animal = uniani;
    datStruct(ii).tstart = time_min';
    datStruct(ii).tstop = time_max';
    datStruct(ii).tduration = time_max' - time_min';
    datStruct(ii).meanspeed = ave_speed';
    datStruct(ii).movmeanspeed = speeds;
    datStruct(ii).totdist = total_dist';
    
    %%
    
    f = find( ic == ii & contains({fils.name},'.x.')' );
    tfiles = fils(f);
    
    datCat = {};
    for jj = 1:length(tfiles)
        filename = fullfile(tfiles(jj).folder,tfiles(jj).name);
        fileID = fopen(filename,'r');
        dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
        fclose(fileID);
        datCat = [datCat dA];
    end
    datCat = reshape(datCat,length(dA),[],1)';
    datmat = cell2mat(datCat(:,2:4));
    uninum = cumsum([1; diff(datmat(:,1))]~=0);
    datmat = [uninum datmat];
    
    xs = datmat;
    
    f = find( ic == ii & contains({fils.name},'.y.')' );
    tfiles = fils(f);
    
    datCat = {};
    for jj = 1:length(tfiles)
        filename = fullfile(tfiles(jj).folder,tfiles(jj).name);
        fileID = fopen(filename,'r');
        dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
        fclose(fileID);
        datCat = [datCat dA];
    end
    datCat = reshape(datCat,length(dA),[],1)';
    datmat = cell2mat(datCat(:,2:4));
    uninum = cumsum([1; diff(datmat(:,1))]~=0);
    datmat = [uninum datmat];
    
    datmat = [xs datmat(:,end)];
    
    [uniani aa ac] = unique(datmat(:,1));
    
    sps = accumarray(ac,1);
    datStruct(ii).xy = mat2cell(datmat(:,4:5),sps,2)';
    datStruct(ii).et = mat2cell(datmat(:,3),sps,1)';
    
end


%%

% filt = contains({datStruct.genotype},{'Repo','None'});
% shortStruct = datStruct(filt);
shortStruct = datStruct;

maxspeed = 0.4;
for ii = 1:length(shortStruct)
    
    dur = shortStruct(ii).tduration;
    f = find(dur > 170);
    if length(f) > 16
        f = randsample(f,16);
    end
    nlines = length(f);
    cols = ceil(nlines/4);
    rows = ceil(nlines/cols);
    
    idx = [repmat([1:cols]',rows,1), repelem(flipud([1:rows]'),cols,1)];
    
    cd = [uint8(jet(64)*255) uint8(ones(64,1))].';
    hold on
    n = 1;
    for jj = f
        xy = shortStruct(ii).xy{jj};
%         xy = xy - xy(1,:);
%         et = shortStruct(ii).et{jj};
%         et = et - et(1);
%         etround = round(et)+1;
        xy = xy-mean(xy)+([20 20].*idx(n,:));
        
        bins = linspace(0,maxspeed,64);
        speeds = shortStruct(ii).movmeanspeed{jj};
        Y = discretize(speeds,bins);
        Y(isnan(Y)) = 64;
        
%         p = plot(xy(:,1),xy(:,2));
%         pause(0.01)
%         set(p.Edge,'ColorBinding','flat','ColorData',cd(1:4,Y))

        p = scatter(xy(:,1),xy(:,2),2,'.');
        p.CData = cd(1:3,Y)';       
        n = n+1;
    end
    hold off
    
    axis equal
    colormap(jet)
    xlim([0 100])
    ylim([0 100])
    pbaspect([1,1,1])
    
    titl = format_title(shortStruct(ii).filename);
%     titl = strrep(shortStruct(ii).genotype,'@','>hEAAT\wedge');
%     titl = strrep(titl,'None','');
    title(titl)
    c = colorbar;
    caxis([0 maxspeed])
    
    
    outname = fullfile(pathfigs,[shortStruct(ii).genotype,'_path2_fig']);
    
    
    print(outname,'-dpdf','-painters','-fillpage');
    
    close

end


%%

speeds = []; dists = []; groups = [];
for ii = 1:length(shortStruct)
    durs = shortStruct(ii).tduration;
    f = durs > 170;
    speeds = [speeds shortStruct(ii).meanspeed(f)];
    dists = [dists shortStruct(ii).totdist(f)];
    groups = [groups repelem(ii,1,sum(f))];
end
% 
% genos = {shortStruct.genotype};
% genos = cellfun(@(x) strrep(x,'None',''), genos, 'UniformOutput', false);
% genos = cellfun(@(x) strrep(x,'@','>hEAAT\wedge'), genos, 'UniformOutput', false);

genos = format_title({shortStruct.filename});


figure(1)
boxplot(speeds,groups,'Colors',[0 0 0],'Symbol','rx','Labels',genos,'LabelOrientation','inline')
hold on
plotSpread(speeds,'distributionIdx',groups,'distributionColors',[.7 .7 1])
hold off
pbaspect([1,1,1])
ylim([0 0.25])
ylabel('Average speed (mm/s)')

ax = gca;
ax.XTickLabel = genos;
ax.XTickLabelRotation = 45;

outname = fullfile(pathfigs,['allgenos_speed_fig']);
print(outname,'-dpdf','-painters','-fillpage');
close


figure(1)
boxplot(dists,groups,'Colors',[0 0 0],'Symbol','rx','Labels',genos,'LabelOrientation','inline')
hold on
plotSpread(dists,'distributionIdx',groups,'distributionColors',[.7 .7 1])
hold off
pbaspect([1,1,1])
ylim([0 40])
ylabel('Distance moved (mm)')

ax = gca;
ax.XTickLabel = genos;
ax.XTickLabelRotation = 45;

outname = fullfile(pathfigs,['allgenos_distance_fig']);
print(outname,'-dpdf','-painters','-fillpage');
close

