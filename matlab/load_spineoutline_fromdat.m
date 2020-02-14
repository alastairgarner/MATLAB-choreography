%% load_spineoutline_fromdat

function structure = load_spineoutline_fromdat(timestamp_directory)

ftype = {'.spine','.outline'};
d = dir(timestamp_directory);
f = contains({d.name},ftype);
d = d(f);
formatSpec = '%s%f%f%[^\n\r]';

filename = fullfile(d(1).folder, d(1).name);
fileID = fopen(filename,'r');
dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow, 'ReturnOnError', false);
fclose(fileID);

o = [dA{2},dA{3}];
outlines = cellfun(@str2num, dA{4}, 'UniformOutput', false);

filename = fullfile(d(2).folder, d(2).name);
fileID = fopen(filename,'r');
dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow, 'ReturnOnError', false);
fclose(fileID);

s = [dA{2},dA{3}];
spines = cellfun(@str2num, dA{4}, 'UniformOutput', false);

[C,ia,ib] = intersect(o,s,'rows');
s = s(ib,:);
spines = spines(ib);
clear dA

structure = cell2struct([num2cell(s),spines,outlines],...
    {'aninum','et','spine','outline'},2);