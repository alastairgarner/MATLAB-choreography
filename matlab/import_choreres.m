%% import_choreres

function [et,aninum,feat] = import_choreres(fileName)

% fileName = 'C:\Users\alyga\OneDrive - McGill University\r\vonMeyel_pipeline\choreography_results\fastlarva@fastlarva\#n#n#n#n\20181116_133755\20181116_133755@fastlarva@fastlarva@t95@#n#n#n#n@100@.x.dat'

sp = split(fileName,'.');
feat = sp{end-1};

delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%f%f%s%[^\n\r]';

fileID = fopen(fileName,'r');
datA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow, 'ReturnOnError', false);
fclose(fileID);

aninum = datA{2};
et = datA{3};
feat = datA{4};

