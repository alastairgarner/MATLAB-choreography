%% finalise_data
% Author - Alastair Garner, alastairgarner@outlook.com


function finalise_data(InputDirectory,OutputDirectory)

currentDir = pwd;
% InputDirectory = pathres;
% OutputDirectory = pathfin;

%%

expr = '\d\d\d\d\d\d\d\d_\d\d\d\d\d\d';

fileType = '*.dat';
inpFiles = dir(horzcat(InputDirectory,'\*\*\*\',fileType));
inpStamps = regexp([inpFiles.name],expr,'match');

outType = '*.outline.dat';
outFiles = dir(horzcat(OutputDirectory,'\*\*\*\',outType));

if length(outFiles) ~= 0
    outStamps = regexp([outFiles.name],expr,'match');
    [Lia Locb] = ismember(inpStamps,outStamps);
    inpFiles = inpFiles(~Lia);

    if length(inpFiles) == 0
        fprintf('\n All files have already been run \n')
        return
    end
    inpStamps = regexp([inpFiles.name],expr,'match');
end

%% 

delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%s%s%s%s%s%s%s%s%s%d%d%d%s%[^\n\r]';

[uStamps ia ic] = unique(inpStamps);

for ii = 1:length(uStamps)
    f = ic == ii;
    iFiles = inpFiles(f);
    f2 = find(~contains({iFiles(:).name},'outline'));
    
    dat = [];
    nam = {};
    n = 1;
    for jj = f2
        sp = split(iFiles(jj).name,'.');
        fname = fullfile(iFiles(jj).folder,iFiles(jj).name);
        [et,aninum,feat] = import_choreres(fname);
        if n == 1
            dat = [et,aninum,feat];
            nam(n) = sp(end-1);
        else
            dat = [dat,feat];
            nam(n) = sp(end-1);
        end
        
        n = n+1;
    end

%     import_choreres()
    
    
    
end





