%% extract_outline
% Author - Alastair Garner, alastairgarner@outlook.com


function extract_outline(InputDirectory,OutputDirectory)


currentDir = pwd;
% InputDirectory = pathinp;
% OutputDirectory = pathres;

%%

expr = '\d\d\d\d\d\d\d\d_\d\d\d\d\d\d';

fileType = '*.blobs';
inpFiles = dir(fullfile(InputDirectory,'**',fileType));
inpStamps = regexp([inpFiles.name],expr,'match');

outType = '*.outline.dat';
outFiles = dir(fullfile(OutputDirectory,'**',outType));

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

%% Loop

delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%s%s%s%s%s%s%s%s%s%d%d%d%s%[^\n\r]';

[uStamps ia ic] = unique(inpStamps);

for ii = 1:length(uStamps)
    
    %% import blobs files
    f = ic == ii;
    flist = inpFiles(f);
    
    dA = {};
    for jj = 1:length(flist)
        filename = fullfile(flist(jj).folder,flist(jj).name);
        fileID = fopen(filename,'r');
        datA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'HeaderLines', startRow, 'ReturnOnError', false);
        fclose(fileID);
        dA = vertcat(dA,datA);
    end
    
    %% reshape data
    cntr = vertcat(dA{:,1});
    ref = find(contains(cntr,'%'));
    idx = setdiff(1:length(cntr),ref);
    
    et = vertcat(dA{:,2});
    anis = et(ref);
    reps = diff([ref' length(cntr)+1]'-1);
    nums = repelem(anis,reps);
    anis = sort(anis);
    
    pts = [vertcat(dA{:,12}),vertcat(dA{:,13}),vertcat(dA{:,14})];
    outlines = vertcat(dA{:,15});
    
    nums = nums(idx);
    et = et(idx);
    pts = pts(idx,:);
    outlines = outlines(idx);
    
    clear datA dA
    
    %% make output filename
    
    sp = split(flist(1).name,'@');
    outputFold = fullfile(OutputDirectory,[sp{2},'@', sp{3}],sp{[5,1]});
    sp = split(flist(1).name,'.blob');
    outName = fullfile(outputFold,[sp{1}(1:end-7),'.outline.dat']);
    
    d = dir(fullfile(outputFold,'*.outline.dat'));
    try
        delete(fullfile(d(1).folder,d(1).name));
    end
    
    %% save new file 
    % to choreography results
    
    fSpec = '%4.0f %3.4f %3.0f %3.0f %3.0f %s \n';
    nums = reshape(nums,length(nums),1);
    outarr = [num2cell(nums), num2cell(et),num2cell(pts),outlines]';

    fileID = fopen(outName,'w');
    fprintf(fileID,fSpec,outarr{:});
    fclose(fileID);
    
    fprintf(['outlines for ',uStamps{ii} ,' generated \n'])
    
end



