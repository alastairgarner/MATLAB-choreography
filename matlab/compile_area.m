%% compile_area.m

%% 
search_term = fullfile(params.directories.choreography_results,'**','**','*area.dat')

d = dir(search_term);

[~,ix,~] = unique(fullfile({d.folder},{d.name}));
d = d(ix);

delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%f%f%s%[^\n\r]';

for ii = 1:numel(d)
    filepath = fullfile(d(ii).folder,d(ii).name);
    
    fid = fopen(filepath);
    datA = textscan(fid, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow, 'ReturnOnError', false);
    fclose(fid);
    
    ids = datA{2};
    ar = datA{4};
    et = datA{3};
    [uids,~,ids] = unique(ids,'stable');
    
    area_mean = accumarray(ids,ar,[],@mean);
    tstart = accumarray(ids,et,[],@min);
    tend = accumarray(ids,et,[],@max);
    
    tstamps = datA{1}(1:numel(uids));
    vals = num2cell([uids,area_mean,tstart,tend]);
    
    to_print = [tstamps,vals]';
    
    outpath = strrep(filepath,'.dat','Summary.txt');
    
    fid = fopen(outpath,'w');
    fprintf(fid,'%s %d %f %.3f %.3f\n',to_print{:});
    fclose(fid);
    
end
    
    
  


