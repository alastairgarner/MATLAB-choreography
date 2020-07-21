%% run_choreography

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision: 


function run_choreography(params)
%Some text here
%
% Some more text
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(1,4);

%% parameters
InputDirectory = params.directories.choreography_input;
OutputDirectory = params.directories.choreography_output;

pval = params.choreography_config.pval;
choreLoc = params.choreography_config.choreLoc;
features = params.choreography_config.features;
M = params.choreography_config.M;
t = params.choreography_config.t;
s = params.choreography_config.s;
currentDir = pwd;

%% list of features

featName = {'x','y','speed','midline','curve','crabspeed',...
    'area','morpwidth','kink','cast','bias','dir','spine'};
featCode = {'x','y','s','m','c','r',...
    'e','M','k','c','b','d','S'};

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
    elseif contains(featCode(ii),{'r','c','m','M','k','b','d'})
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

if numel(inpFiles) ~= numel(inpStamps)
    fprintf('\n Bad file found \n')
    return
end

outType = '*.dat';
outFiles = dir(fullfile(OutputDirectory,'**',outType));

%%%%%%%%%%%%%%%%%%%%%%
if length(outFiles) ~= 0
    outStamps = regexp([outFiles.name],'\d\d\d\d\d\d\d\d_\d\d\d\d\d\d','match');

    feats = params.choreography_config.features;
    inpComp = [repelem(inpStamps,numel(feats));repmat(feats,1,numel(inpStamps))]';
    outComp = regexp([outFiles.name],'(\d\d\d\d\d\d\d\d_\d\d\d\d\d\d)[a-zA-Z_0-9@#]*[.](\w+)[.]','tokens');
    outComp = regexp([outFiles.name],'(\d{8}_\d{6})[a-zA-Z_0-9@#()]*[.](\w+)[.]','tokens');
    outComp = reshape([outComp{:}],2,[])';

    f = ~ismember(string(inpComp),string(outComp),'rows');
    inpComp = unique(inpComp(f),'stable');
    
    f = ismember(inpStamps,inpComp);
    inpFiles = inpFiles(f);
    
    if length(inpFiles) == 0
        fprintf('\n All files have already been run \n')
        return
    end
    inpStamps = regexp({inpFiles.name},expr,'match','once');
end

%%%%%%%%%%%%%%%%%%%%%%

% if length(outFiles) ~= 0
% 
%     outStamps = regexp([outFiles.name],expr,'match');
%     
%     [Lia, ~] = ismember(inpStamps,outStamps);
%     [~, Locb] = ismember(outStamps,inpStamps);
%     
%     [C,~,ic] = unique(outStamps,'stable');
%     outlength = accumarray(ic,1);
%     redo = outlength < length(features) | outlength > length(features)+2;
%     [~,redo] = ismember(inpStamps,C(redo));
%     
%     inpFiles = inpFiles(~Lia|redo);
% 
%     if length(inpFiles) == 0
%         fprintf('\n All files have already been run \n')
%         return
%     end
%     inpStamps = regexp([inpFiles.name],expr,'match');
% 
% end
%%

% filterByDriver = false;
% if filterByDriver
%     fid = fopen('screen_genotypes.txt');
%     myGenotypes = textscan(fid,'%s\n');
%     myGenotypes = string([myGenotypes{:}]);
%     fclose(fid);
%     
%     expr = '\d{8}_\d{6}@([a-zA-Z_0-9]+)@*';
%     driver = regexp({inpFiles.name},expr,'tokens','once');
%     driver = string(driver)';
% 
%     genoFilter = ismember(driver,myGenotypes);
%     inpFiles = inpFiles(genoFilter);
% end

%% run choreography

% 20170928_111031

delimiter = ' ';
startRow = 0;
formatSpec = '%f%f%f%f%f%f%f%[^\n\r]';

folderList = struct();
folderList = {inpFiles.folder}';

for ii = 1:length(folderList)
    sp = split(inpFiles(ii).name,'@');
    
    if ~startsWith(inpFiles(ii).folder,'/media')
        isExternal = false;
        inputFold = inpFiles(ii).folder;
        outputFold = fullfile(outFold,[sp{2},'@', sp{3}],sp{[5,1]});
    else
        isExternal = true;
        if isdir('./_temp')
            rmdir('./_temp','s');
        end
        inputFold = './_temp/chore_input';
        mkdir(inputFold);
        copyfile(inpFiles(ii).folder,inputFold);
        outputFold = fullfile('./_temp',[sp{2},'@', sp{3}],sp{[5,1]});
        copyDest = fullfile(outFold,[sp{2},'@', sp{3}],sp{[5,1]});
    end
    
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fileNamesAdj = regexprep(fileNames,'\d{5}[.]','');
    fileTypes = regexp(fileNamesAdj,'[.]','split','once');
    fileTypes = cellfun(@(x) x(2),fileTypes);
%     fileTypes = cellfun(@(x) x(ismember(x,featName)), fileTypes);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [fileTypesU Ua Uc] = unique(fileTypes);
    
    fprintf('\n')
    for jj = 1:length(fileTypesU)
        if length(fileNames) > 2
            
            fprintf(['        ...merging "' fileTypesU{jj} '" data... \n\n']);
            
            datCat = {};
            for kk = find(Uc == jj)'
                dA = fileread(fileNames{kk});
                datRow = strsplit(dA,'\n');
                datCat = [datCat datRow(1:end-1)];
            end
            
            datCat = cellfun(@(x) [sp{1} ' ' x], datCat, 'UniformOutput', false);
            
            spOut = split(fileNames{Ua(jj)},'.');
            outName = strcat(spOut{1},'.',fileTypesU(jj));
            fileID = fopen(string(outName),'w');
            fprintf(fileID, '%s\n', datCat{:});
            fclose(fileID);

            delete(fileNames{Uc == jj});
            
            fprintf(['\b\b success \n']);
            
        end
    end
    Fs = dir('*.dat');
    
    yaml.WriteYaml("config.yaml",params);
    
    cd(currentDir)
    
    if isExternal
        if isdir(copyDest)
            rmdir(copyDest,'s');
        end
        mkdir(copyDest);
        copyfile(outputFold,copyDest)
    end
    
    fprintf(['Choreography complete for timestamp ' sp{1} '\n']);
end

try
    rmdir('./_temp','s');
end