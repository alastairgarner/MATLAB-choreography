%% Azman_Choreography

dataDir = uigetdir('Select RAW Data folder');
currentDir = cd;

%% Get raw files

fileType = '*.summary';
chorFiles = dir(horzcat(dataDir,'\',rawDir,'\*\*\',fileType));

%% Preprocess

iserror = [];
for i = 1:length(chorFiles)
    spName = split(chorFiles(i).name,'@');
    sp = split(chorFiles(i).folder,'\');
    expression = '\d\d\d\d\d\d\d\d_\d\d\d\d\d\d';
    isdate = regexp(sp{end},expression);
    if isdate == 1
        folderLoc = chorFiles(i).folder;
        newLoc = fullfile(dataDir,inpDir,sp{end-1:end});
        status = rmdir(newLoc, 's');
        status = mkdir(newLoc);
        status = copyfile(folderLoc,newLoc);
        cd(newLoc)
        localFiles = dir();
        arrayfun(@(x) movefile(x.name,[sp{end},'@',x.name]), localFiles(3:end));
        cd(currentDir)
    else
        iserror = [iserror isdate];
    end
end

%% Choreography command line

cmdLine{1} = 'java -Xincgc -Xms8000m -Xmx8000m -jar C:\Users\alyga\Desktop\Chore.jar -t 5 -s 0.1 -p 0.0528 -S --nanless -o Dts1234 -O speed -N all --target ';
cmdLine{2} = 'java -Xincgc -Xms8000m -Xmx8000m -jar C:\Users\alyga\Desktop\Chore.jar -t 5 -s 0.1 -p 0.0528 -S --nanless -o Dtx1234 -O x -N all --target ';
cmdLine{3} = 'java -Xincgc -Xms8000m -Xmx8000m -jar C:\Users\alyga\Desktop\Chore.jar -t 5 -s 0.1 -p 0.0528 -S --nanless -o Dty1234 -O y -N all --target ';
outFold = fullfile(dataDir,resDir);

%% LoopChore

% Just generate speed files

fileType = '*.summary';
chorFiles = dir(horzcat(dataDir,'\',inpDir,'\*\*\',fileType));
folderList = struct();
folderList = {chorFiles.folder}';

delimiter = ' ';
startRow = 0;
formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';

for i = 2:length(folderList)
    sp = split(chorFiles(i).name,'@');
    inputFold = chorFiles(i).folder;
    outputFold = fullfile(outFold,[sp{2},'@', sp{3}],sp{[5,1]});
    
    fullLine = {};
    for jj = 1:length(cmdLine)
        fullLine{jj} = [cmdLine{jj}, outputFold, ' ', inputFold];
    end
    
    if isdir(outputFold)
        rmdir(outputFold,'s');
        mkdir(outputFold);
    else
        mkdir(outputFold);
    end    
    
    for jj = 1:length(cmdLine)
        status = system(fullLine{jj});
    end
    cd(outputFold);
    aniFiles = dir();
    fileNames = {aniFiles.name}';
    fileNames = fileNames(~[aniFiles.isdir]);
    fileTypes = cellfun(@(x) split(x,'.'), fileNames, 'UniformOutput', false);
    fileTypes = cellfun(@(x) x{2}, fileTypes, 'Uniformoutput', false);
    [fileTypesU Ua Uc] = unique(fileTypes);
    
    for jj = 1:length(fileTypesU)
        if length(fileNames) > 2
            datCat = {};
            for j = find(Uc == jj)'
                fileID = fopen(fileNames{j},'r');
                dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
                fclose(fileID);

                datCol = string(repmat(sp{1},length(dA{1}),1));
                dA = [{datCol} dA];
                datCat = vertcat(datCat,dA);
            end

            datCat = cellfun(@(col) horzcat(col{:}), num2cell(datCat, 2), 'UniformOutput', false);
            datCat = vertcat(datCat{:});

            spOut = split(fileNames{Ua(jj)},'.');
            outName = [spOut{1},'.',spOut{2},'.',spOut{4}];
            fileID = fopen(outName,'w');
            [nrows,ncols] = size(datCat);

            for row = 1:nrows
                fprintf(fileID,'%s %s %s %s %s %s %s %s %s\n',datCat(row,:));
            end

            fclose(fileID);

            delete(fileNames{Uc == jj});
            
        end
    end
    cd(currentDir)
end

%%



