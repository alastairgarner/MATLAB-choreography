%% run_choreography
% Author - Alastair Garner, alastairgarner@outlook.com

function run_choreography(InputDirectory,OutputDirectory,params)

narginchk(1,4);


%% parameters

% InputDirectory = pathinp;
% OutputDirectory = pathres;

pval = params.choreography_config.pval;
choreLoc = params.choreography_config.choreLoc;
features = params.choreography_config.features;
M = params.choreography_config.M;
t = params.choreography_config.t;
s = params.choreography_config.s;
currentDir = pwd;

%% list of features

featName = {'x','y','speed','midline','curve','crabspeed','area'};
featCode = {'x','y','s','m','c','r','e'};

[C ia ic] = intersect(featName,features,'stable');
[Lia Locb] = ismember(features,featName);


%% construct commands

if isunix
    java_call = 'java -jar';
elseif ismac
    java_call = 'java -jar';
elseif ispc
%                 java_call = 'java -Xincgc -Xms8000m -Xmx8000m -jar';
    java_call = 'java -jar';
else
    fprintf('\n Operating System not supported')
    return
end

cmd = {};
n = 1;
for ii = Locb
    if contains(featCode(ii),{'s','x','y','e'})
        cmd{n} = [java_call ' "' choreLoc '" -t ' num2str(t) ' -s ' num2str(s) ' -M ' num2str(M) ' -q --quiet -p ' num2str(pval) ' -S --nanless -o Dt' featCode{ii} '1234 -O ' featName{ii} ' -N all --target '];
    elseif contains(featCode(ii),{'r','c','m'})
        cmd{n} = [java_call ' "' choreLoc '" -t ' num2str(t) ' -s ' num2str(s) ' -M ' num2str(M) ' -q --quiet -p ' num2str(pval) ' --plugin Reoutline::exp --plugin Respine::0.23::tapered=0.28,1,2 --plugin SpinesForward::rebias --minimum-biased 3mm -S --nanless -o Dt' featCode{ii} '1234 -O ' featName{ii} ' -N all --target '];
    elseif contains(featCode(ii),{'S'})
        cmd{n} = [java_call ' "' choreLoc '" -t ' num2str(t) ' -s ' num2str(s) ' -M ' num2str(M) ' -q --quiet -p ' num2str(pval) ' --plugin Reoutline::exp --plugin Respine::0.23::tapered=0.28,1,2 --plugin SpinesForward::rebias --minimum-biased 3mm -S --nanless --plugin Extract::spine --target '];
    end
    n = n+1;
end

if ~params.choreography_config.runQuiet
    cmd = strrep(cmd,'-q --quiet ','');
end

%% check for files already run

expr = '\d\d\d\d\d\d\d\d_\d\d\d\d\d\d';
outFold = OutputDirectory;

inpType = '*.summary';
inpFiles = dir(fullfile(InputDirectory,'**',inpType));
inpStamps = regexp([inpFiles.name],expr,'match');

outType = '*.dat';
outFiles = dir(fullfile(OutputDirectory,'**',outType));

if length(outFiles) ~= 0

    outStamps = regexp([outFiles.name],expr,'match');
    
    [Lia, ~] = ismember(inpStamps,outStamps);
    [~, Locb] = ismember(outStamps,inpStamps);
    
    [C,~,ic] = unique(outStamps,'stable');
    outlength = accumarray(ic,1);
    redo = outlength < length(features) | outlength > length(features)+2;
    [~,redo] = ismember(inpStamps,C(redo));
    
    inpFiles = inpFiles(~Lia|redo);

    if length(inpFiles) == 0
        fprintf('\n All files have already been run \n')
        return
    end
    inpStamps = regexp([inpFiles.name],expr,'match');

end
%% run choreography

delimiter = ' ';
startRow = 0;
formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';

folderList = struct();
folderList = {inpFiles.folder}';


for ii = 1:length(folderList)
    sp = split(inpFiles(ii).name,'@');
    inputFold = inpFiles(ii).folder;
    outputFold = fullfile(outFold,[sp{2},'@', sp{3}],sp{[5,1]});
    
    fullLine = {};
    for jj = 1:length(cmd)
        fullLine{jj} = [cmd{jj},'"',outputFold,'" "', inputFold,'"'];
    end
    
    if isdir(outputFold)
        rmdir(outputFold,'s');
        mkdir(outputFold);
    else
        mkdir(outputFold);
    end    
    
    for jj = 1:length(cmd)
        status = system(fullLine{jj});
    end
    cd(outputFold);
    aniFiles = dir();
    fileNames = {aniFiles.name}';
    fileNames = fileNames(~[aniFiles.isdir]);
    fileTypes = cellfun(@(x) split(x,'.'), fileNames, 'UniformOutput', false);
    fileTypes = cellfun(@(x) x{2}, fileTypes, 'Uniformoutput', false);
    [fileTypesU Ua Uc] = unique(fileTypes);
    
    fprintf('\n')
    
    for jj = 1:length(fileTypesU)
        if length(fileNames) > 2
            
            fprintf(['        ...merging "' fileTypesU{jj} '" data... \n\n']);
            
            datCat = {};
            for j = find(Uc == jj)'
%                 fileID = fopen(fileNames{j},'r');
%                 dA = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'HeaderLines', startRow, 'ReturnOnError', false);
%                 fclose(fileID);
%                 datCol = string(repmat(sp{1},length(dA{1}),1));
%                 dA = [{datCol} dA];
%                 datCat = vertcat(datCat,dA);

                dA = fileread(fileNames{j});
                datRow = strsplit(dA,'\n');
                datCat = [datCat datRow(1:end-1)];
            end
            
            datCat = cellfun(@(x) [sp{1} ' ' x], datCat, 'UniformOutput', false);
            
            spOut = split(fileNames{Ua(jj)},'.');
            outName = [spOut{1},'.',spOut{2},'.',spOut{4}];
            fileID = fopen(outName,'w');
            fprintf(fileID, '%s\n', datCat{:});
            fclose(fileID);

            delete(fileNames{Uc == jj});
            
            fprintf(['\b\b success \n']);
            
        end
    end
    Fs = dir('*.dat');
    
    yaml.WriteYaml("config.yaml",params);
    
    cd(currentDir)
    
    fprintf(['Choreography complete for timestamp ' sp{1} '\n']);
end


