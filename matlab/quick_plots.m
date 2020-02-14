%% quick_plots

clear all
clc

%% Load and set up

initialise_folders
datadir;

fils = dir(fullfile(pathres,'\*\*\*\*'));
fils = fils(~[fils.isdir]);


parts = arrayfun(@(x) split(x.name,["@","."]), fils, 'UniformOutput', false);
tstamps = cellfun(@(x) x{1}, parts, 'UniformOutput', false);
genos = cellfun(@(x) strcat(x{2},'@',x{3}), parts, 'UniformOutput', false);

[genosU ia ic] = unique(genos);

delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%f%f%[^\n\r]';

datstruct = struct();
for ii = 1:length(genosU)
    
    f = find( ic == ii & contains({fils.name},'.speed')' );
    fils_geno = fils(f);
    
    sp = split(fils_geno(1).name,["@","."]);
    tstamp = sp{1};
    eff = sp{3};
    driv = sp{2};
    feat = sp{end-1};
    outdir = fullfile(pathfigs,eff,driv);
    outname = [feat,'_average_fig'];
    outfull = fullfile(outdir,outname);
    
    if ~isdir(outdir)
        mkdir(outdir);
    end
    
    datCat = {};
    tstamps = {};
    for jj = 1:length(fils_geno)
        filename = fullfile(fils_geno(jj).folder,fils_geno(jj).name);
        fileID = fopen(filename,'r');
        dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
        fclose(fileID);
        datCat = [datCat dA];        
    end
    datCat = reshape(datCat,length(dA),[],1)';
    datmat = cell2mat(datCat(:,2:4));
    
    mint = min(datmat(:,2));
    maxt = max(datmat(:,2));
    
    nums = datmat(:,1);
    diffs = [1 [abs(diff(nums))>0]' ]';
    diffs(1) = 1;
    nums = cumsum(diffs)';
    
    [uniani aa ac] = unique(nums);
    length(uniani);
    
    time_min = accumarray(ac,datmat(:,2),[],@min);
    time_max = accumarray(ac,datmat(:,2),[],@max);
    ave_speed = accumarray(ac,datmat(:,3),[],@mean);
    
    timediff = [0 diff(datmat(:,2))']';
    f = abs(timediff) > 0.15;
    timediff(f) = 0;
    dist = timediff.*datmat(:,3);
    
    total_dist = accumarray(ac,dist);
    
    datTab = [uniani' time_min time_max ave_speed total_dist];
    
    
    datStruct(ii).genotype = strcat(driv,'@',eff); 
    datStruct(ii).animal = uniani;
    datStruct(ii).tstart = time_min';
    datStruct(ii).tstop = time_max';
    datStruct(ii).meanspeed = ave_speed';
    datStruct(ii).totdist = total_dist';
    
    %% 
    
    f = find( ic == ii & contains({fils.name},'.x.')' );
    fils_geno = fils(f);
    
    datCat = {};
    for jj = 1:length(fils_geno)
        filename = fullfile(fils_geno(jj).folder,fils_geno(jj).name);
        fileID = fopen(filename,'r');
        dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
        fclose(fileID);
        datCat = [datCat dA];
    end
    datCat = reshape(datCat,length(dA),[],1)';
%     mult = cellfun(@length, datCat(:,1));
    
    xs = cell2mat(datCat(:,[2:4]));
    
    f = find( ic == ii & contains({fils.name},'.y.')' );
    fils_geno = fils(f);
    
    datCat = {};
    for jj = 1:length(fils_geno)
        filename = fullfile(fils_geno(jj).folder,fils_geno(jj).name);
        fileID = fopen(filename,'r');
        dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
        fclose(fileID);
        datCat = [datCat dA];
    end
    datCat = reshape(datCat,length(dA),[],1)';
%     mult = cellfun(@length, datCat(:,1));

    ys = cell2mat(datCat(:,[2:4]));
    
%     unique(xs(:,1)) unique(ys(:,1))
    
    [~,ia,ib] = intersect(xs(:,1:2),ys(:,1:2),'rows');
    datmat = [xs(ia,:), ys(ib,3)];
    
    nums = datmat(:,1);
    diffs = [1 [abs(diff(nums))>0]' ]';
    diffs(1) = 1;
    nums = cumsum(diffs)';
    
    [uniani aa ac] = unique(nums);
    
    sps = accumarray(ac,1);
    datStruct(ii).xy = mat2cell(datmat(:,3:4),sps,2)';
    datStruct(ii).et = mat2cell(datmat(:,2),sps,1)';

end

%% Select data

filt = contains({datStruct.genotype},{'Repo','None'});
shortStruct = datStruct(filt);

%% Plot
    
for ii = 1:length(shortStruct)
    
    dur = cellfun(@max, shortStruct(ii).et) - cellfun(@min, shortStruct(ii).et);
    f = find(dur > 60);
    if length(f) > 7
        f = randsample(f,7);
    end
    
    cd = [uint8(jet(181)*255) uint8(ones(181,1))].';
    hold on
    for jj = f
        xy = shortStruct(ii).xy{jj};
        xy = xy - xy(1,:);
        et = shortStruct(ii).et{jj};
        et = et - et(1);
        etround = round(et)+1;
        
        p = scatter(xy(:,1),xy(:,2),'.');
        p.CData = cd(1:3,etround)';       
        
    end
    hold off
    
    colormap(jet)
    xlim([-15 15])
    ylim([-15 15])
    pbaspect([1,1,1])
    title(shortStruct(ii).genotype)
    c = colorbar;
    caxis([0 180])
    
    
    outname = fullfile(pathfigs,[shortStruct(ii).genotype,'_path2_fig']);
    
    
    print(outname,'-dpdf','-painters','-fillpage');
    
    close

end
    

%% Plot 2

speeds = {};
for ii = 1:length(shortStruct)
    durs = shortStruct(ii).tstop - shortStruct(ii).tstart;
    f = durs > 60;
    speeds{ii} = shortStruct(ii).meanspeed(f);
    dists{ii} = shortStruct(ii).totdist(f);
end

genos = {shortStruct.genotype};

figure(1)
plotSpread(speeds,'showMM',4)
ylim([0 0.2])
pbaspect([1 1 1])
ylabel('Larva speed (mm/s)')
ax = gca;
ax.XTickLabel = genos;
ax.XTickLabelRotation = 45;

outname = fullfile(pathfigs,['allgenos_speed_fig']);
print(outname,'-dpdf','-painters','-fillpage');
close


figure(2)
plotSpread(dists,'showMM',4)
ylim([0 35])
pbaspect([1 1 1])
ylabel('Distance moved (mm)')
ax = gca;
ax.XTickLabel = genos;
ax.XTickLabelRotation = 45;

outname = fullfile(pathfigs,['allgenos_distance_fig']);
print(outname,'-dpdf','-painters','-fillpage');
close


