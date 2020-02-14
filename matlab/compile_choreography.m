%% compile_choreography.m

%% 
search_term = fullfile(params.directories.choreography_results,'**','**','*.dat');

d = dir(search_term);

[~,ix,~] = unique(fullfile({d.folder},{d.name}));
d = d(ix);

delimiter = ' ';
startRow = 0;
formatSpec = '%s%s%s%s%s%[^\n\r]';

cols = {'curve', 'midline', 'speed', 'crabspeed', 'x', 'y', 'area'};
expr = '\d{8}_\d{6}';
tstamps = regexp({d.name},expr,'match','once');
[ux,~,ix] = unique(tstamps,'stable');

for ii = 1:numel(ux)
    fprintf(strcat(string(ii),'\n'))
    f = ix == ii;
    dd = d(f);
    print_mat = [];
    
    for jj = 1:numel(cols)
        seacrhterm = ['.',cols{jj},'.'];
        ff = contains({dd.name},seacrhterm);
        filepath = fullfile(dd(ff).folder,dd(ff).name);
        
        fid = fopen(filepath);
        datA = textscan(fid, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
        fclose(fid);
        
        if jj == 1
            print_mat = [print_mat datA([2,1,3,4])];
        elseif strcmp(cols{jj},'area')
            ids = double(string(datA{2}));
            [~,~,ids] = unique(ids,'stable');
            vals = double(string(datA{4}));
            counts = accumarray(ids,1);
            vals = accumarray(ids,vals,[],@mean);
            print_mat = [print_mat {cellstr(string(repelem(vals,counts,1)))}];
        else
            print_mat = [print_mat datA(4)];
        end
    end
    
    print_mat = [print_mat{:}]';
    
    [pat,fil,ext] = fileparts(filepath);
    out_dir = '/home/alastair/Downloads/For_Jane/For_Jane/compiled_chore/';
    sp = split(fil,'@100@');
    outpath = fullfile(out_dir,[sp{1},'@100_compiledChore.txt']);
    
    fid = fopen(outpath,'w');
    fprintf(fid,'%s %s %s %s %s %s %s %s %s %s\n',string({'id', 'timestamp', 'time'}),cols{:});
    fprintf(fid,'%s %s %s %s %s %s %s %s %s %s\n',print_mat{:});
    fclose(fid);
end
    
    
  


