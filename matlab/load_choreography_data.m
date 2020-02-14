%% load_choreography_data

function [dataStructure] = load_choreography_data(file_structure,geno_number)

% geno_number = [1; 3; 5];
% file_structure = tfiles
fils = file_structure;

if nargin > 1
    pick_geno = geno_number;
end

fils = fils(~[fils.isdir]);

filenames = string({fils.name});
expr = ['(?<date>^\d\d\d\d\d\d\d\d)[_]'...
    '(?<time>\d\d\d\d\d\d)[@]'...
    '(?<driver>[\w()]*)[@]'...
    '(?<effector>[\w()]*)[@]'...
    '(?<tracker>[t]\d*)[@]'...
    '[#](?<prot1>\w*)'...
    '[#](?<prot2>\w*)'...
    '[#](?<prot3>\w*)'...
    '[#](?<prot4>\w*)[@]'...
    '\d*[@]'...
    '[.](?<filetype>\w*)[.]'];
descript = regexp(filenames,expr,'names');
descript = vertcat(descript{:});
descript_list = unique([descript.filetype]);
descript_list = cellstr(descript_list);

genos = strcat([descript.driver],'@',[descript.effector])';
[genosU,ia,ic] = unique(genos);

if exist('pick_geno')
    if size(pick_geno,1) > 1
        pick_geno = pick_geno';
    end
    filter = any(ic == pick_geno,2);
    fils = fils(filter);
    descript = descript(filter);
end
    
delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%f%f%[^\n\r]';

tstamps = strcat([descript.date],'_',[descript.time]);
[tstampsU,ia,Ix] = unique(tstamps,'stable');

dStruct = struct();
for ii = 1:length(tstamps)
    idx = Ix(ii);
    fname = fullfile(fils(ii).folder,fils(ii).name);
    
    fileID = fopen(fname,'r');
    dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
    fclose(fileID);
    
    datmat = [dA{2:4}];
    
    if any(ia == ii)
        dStruct(idx).driver = descript(ii).driver;
        dStruct(idx).effector = descript(ii).effector;
        dStruct(idx).tracker = descript(ii).tracker;        
        dStruct(idx).date = descript(ii).date;
        dStruct(idx).timestamp = descript(ii).time;
        dStruct(idx).animal = datmat(:,1)';
        dStruct(idx).time = datmat(:,2)';
        dStruct(idx).complete = false;
    end
    
    dStruct(idx).(descript(ii).filetype) = datmat(:,3)';
    
    if strcmp(descript(ii).filetype,'outline')
        fileID = fopen(fname,'r');
        dA = textscan(fileID, '%s%s%[^\n\r]', 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false,'MultipleDelimsAsOne',1);
        fclose(fileID);
        dStruct(idx).(descript(ii).filetype) = dA{3};
    end
    
end

for ii = 1:length(dStruct)
    filt = cellfun(@(x) ~isempty(dStruct(ii).(x)), descript_list );
    if all(filt)
        dStruct(ii).complete = true;
    end
end

dataStructure = dStruct;
