%% fix_rawdata

% Author - Alastair Garner, alastairgarner@outlook.com

function fix_rawdata(params)

InputDirectory = params.directories.mwt;
OutputDirectory = params.directories.choreography_input;

currentDir = pwd;

fileType = '*.summary';
inpFiles = dir(fullfile(InputDirectory,'**',fileType));
outFiles = dir(fullfile(OutputDirectory,'**',fileType));

%% check for new files

if ~numel(inpFiles)
    return
end

if length(outFiles) ~= 0
    expr = '\d\d\d\d\d\d\d\d_\d\d\d\d\d\d';
    inpStamps = regexp([inpFiles.folder],expr,'match');
    outStamps = regexp([outFiles.name],expr,'match');
    [Lia Locb] = ismember(inpStamps,outStamps);
    inpFiles = inpFiles(~Lia);
end
    
%% copy and rename

iserror = [];
for ii = 1:length(inpFiles)
    spName = split(inpFiles(ii).name,'@');
    sp = split(inpFiles(ii).folder,{'\','/'});
    expression = '\d\d\d\d\d\d\d\d_\d\d\d\d\d\d';
    isdate = regexp(sp{end},expression);
    if isdate == 1
        folderLoc = inpFiles(ii).folder;
        newLoc = fullfile(OutputDirectory,sp{end-1:end});
        status = rmdir(newLoc, 's');
        status = mkdir(newLoc);
        status = copyfile(folderLoc,newLoc);
        
        cd(newLoc)
        localFiles = dir();
        localFiles = localFiles(~[localFiles.isdir]);
        
        prefix = arrayfun(@(x) regexp(x.name,expression),localFiles,'UniformOutput', false);
        if all(cellfun(@isempty,prefix))
            arrayfun(@(x) movefile(x.name,[sp{end},'@',x.name]), localFiles);
        end
        
        cd(currentDir)
        fprintf(['files for ' sp{end}, ' copied \n']);
    else
        fprintf(['copying files for ' sp{end} ' failed \n']);
    end
end

