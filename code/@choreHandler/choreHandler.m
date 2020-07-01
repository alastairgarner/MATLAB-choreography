%% choreHandler

% Other m-files required: 
% Subfunctions: 
% MAT-files required: 
% See also:

% Author: Alastair Garner
% email: alastairgarner@outlook.com
% Website: https://github.com/alastairgarner/
% August 2019; Last revision: 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef choreHandler
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    %%
    properties
        filepaths = "";
        date = "";
        time = "";
        driver = "";
        effector = "";
        protocol1 = "";
        protocol2 = "";
        protocol3 = "";
        protocol4 = "";
        rig = "";
        datatypes = "";
        group = [];
        plot_order = [];
        iscontrol = false;
        filetypes = "";
        genotype_display = "";
        
        data_choreography
        
        
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS - NORMAL
    methods
        %% CONSTRUCTOR
        function obj = choreHandler(config_file)
            if nargin == 0
                params = choreHandler.load_parameters();
            else
                if isstr(config_file)
                    params = choreHandler.load_parameters(config_file);
                else
                    params = config_file;
                end
            end

            files = choreHandler.get_filelist(params);
            file_details = choreHandler.get_file_details(files);
            [C,ia,ic] = unique([file_details.full_timestamp],'stable');
            
            objects = {};
            for ii = 1:length(C)
                idx = find(ic == ii);
                filepaths = string(fullfile({files(idx).folder},{files(idx).name}))';
                res = regexp(filepaths,'[.]\d{5}[.]');
                if any(~cellfun(@isempty, res))
                    continue
                end
                obj.filepaths = filepaths;
                obj.date = file_details(idx(1)).date;
                obj.time = file_details(idx(1)).time;
                obj.driver = file_details(idx(1)).driver;
                obj.effector = file_details(idx(1)).effector;
                obj.protocol1 = file_details(idx(1)).prot1;
                obj.protocol2 = file_details(idx(1)).prot2;
                obj.protocol3 = file_details(idx(1)).prot3;
                obj.protocol4 = file_details(idx(1)).prot4;
                obj.rig = file_details(idx(1)).tracker;
                obj.filetypes = [file_details(idx).filetype];
                objects{ii} = obj;
            end
            obj = vertcat(objects{:});
        end
        
        %%
        function obj = load_data_choreography(obj)
            fprintf('\n')
            delimiter = ' ';
            startRow = 0;
            formatSpec = '%s%f%f%f%[^\n\r]';
            
            filenames = obj.filepaths;
            timestamps = regexp(filenames,'(\d{8}[_]\d{6})','tokens','once');
            [uni_ts,~,idx_ts] = unique([timestamps{:}]);

            clear ani_cell
            for jj = 1:length(uni_ts)
                fprintf('\t ...loading choreography data for %s \n', uni_ts(jj))
                
                filenames_ts = filenames(idx_ts == jj);
                filetypes = regexp(filenames_ts,'[.](\w+)[.]','tokens','once');
                filt = startsWith([filetypes{:}],["x","y"]);
                if sum(filt)<2
                    fprintf('\t x or y data missing\n')
                    continue
                end
                clear ani_struct
                for ii = 1:length(filenames_ts)
                    filetype = regexp(filenames_ts(ii),'[.](\w+)[.]','tokens','once');
                    if filetype == "outline"
                        continue
                    end

                    fileID = fopen(filenames_ts(ii));
                    dataArray = textscan(fileID, formatSpec,...
                        'Delimiter', delimiter,...
                        'MultipleDelimsAsOne', true,...
                        'HeaderLines' ,startRow,...
                        'ReturnOnError', false,...
                        'CollectOutput',true);
                    fclose(fileID);
                    mat = dataArray{2};

                    [C,~,ic] = unique(mat(:,1));
                    if ~exist('ani_struct')
                        C = num2cell(C);
                        [ani_struct(1:length(C)).aniID] = C{:};
                        [ani_struct(:).timestamp] = dataArray{1}{1:length(C)};
                        elem_count = accumarray(ic,1);
                        et = mat2cell(mat(:,2)',1,elem_count);
                        [ani_struct(:).elapstime] = et{:};
                        C = [C{:}];
                    end
                    existing_ID = [ani_struct(:).aniID];
                    idx = find(ismember(existing_ID,C));
                    
                    existing_counts = arrayfun(@(x) length(x.elapstime), ani_struct);
                    blank_cells = cellfun(@(x) nan(1,x), num2cell(existing_counts), 'UniformOutput', false);

                    elem_count = accumarray(ic,1);
                    beh_data = mat2cell(mat(:,3)',1,elem_count);
                    blank_cells(idx) = beh_data;
                    
                    [ani_struct(:).(filetype)] = blank_cells{:};
                end
                ani_cell(jj) = {ani_struct};
            end
            fields = cellfun(@(x) string(strjoin(sort(fieldnames(x)))), ani_cell)';
            filt = max(sum(fields == fields')) == sum(fields == fields');
            ani_cell = ani_cell(filt);
            data_chore = [ani_cell{:}]';
            uniID = num2cell(1:length(data_chore));
            [data_chore(:).uniID] = uniID{:};
            obj.data_choreography = data_chore;
        end
        
        %%
        function obj = process_choreography(obj,params)
            for ii = 1:length(obj)
                grp = obj(ii).group;
                temp = obj(ii).data_choreography;

                dist_func = @(x,y) sqrt(sum([[0;0],diff([x;y],[],2)].^2));
                vect_func = @(x,y) sqrt(sum(diff([x([1,end]);y([1,end])],[],2).^2));
                
                curve_smooth = arrayfun(@(f) movmean(f.curve(f.time_filter),params.smooth_window), temp,'UniformOutput',false);
                mean_curve_smooth = cellfun(@mean,curve_smooth,'UniformOutput', false);
                
                x_crude = arrayfun(@(f) f.x(f.time_filter), temp,'UniformOutput',false);
                y_crude = arrayfun(@(f) f.y(f.time_filter), temp,'UniformOutput',false);
                distance_crude = cellfun(dist_func, x_crude, y_crude, 'UniformOutput', false);
                x_smooth = arrayfun(@(f) movmean(f.x(f.time_filter),params.smooth_window), temp,'UniformOutput',false);
                y_smooth = arrayfun(@(f) movmean(f.y(f.time_filter),params.smooth_window), temp,'UniformOutput',false);
                distance_smooth = cellfun(dist_func, x_smooth, y_smooth, 'UniformOutput', false);
                
                path_distance_crude = num2cell(cellfun(@sum, distance_crude));
                path_distance_smooth = num2cell(cellfun(@sum, distance_smooth));
                vect_distance = cellfun(vect_func, x_smooth, y_smooth, 'UniformOutput', false,'ErrorHandler',@(x,y,z) nan);
                
                mean_speed_crude = num2cell(cellfun(@(x,y) x./range(y), path_distance_crude, {temp.elapstime}'));
                mean_speed_smooth = num2cell(cellfun(@(x,y) x./range(y), path_distance_smooth, {temp.elapstime}'));
                
                [mi,ma] = arrayfun(@(x) bounds(x.elapstime), temp);
                startfin = num2cell([mi,ma],2);

                [temp.group] = deal(grp);
                [temp.curve_smooth] = curve_smooth{:};
                [temp.mean_curve_smooth] = mean_curve_smooth{:};
                [temp.x_smooth] = x_smooth{:};
                [temp.y_smooth] = y_smooth{:};
                [temp.distance] = distance_smooth{:};
                [temp.path_distance_smooth] = path_distance_smooth{:};
                [temp.path_distance_crude] = path_distance_crude{:};
                [temp.vect_distance] = vect_distance{:};
                [temp.startfin] = startfin{:};
                [temp.mean_speed_smooth] = mean_speed_smooth{:}; 
                [temp.mean_speed_crude] = mean_speed_crude{:};
                
                path_distance_byarea = num2cell(arrayfun(@(x) choreHandler.get_distance_by_area(x,.5), temp));
                [temp.path_distance_byarea] = path_distance_byarea{:};
                
                obj(ii).data_choreography = temp;
            end
        end
        
        %% Sort By Genotype
        function obj = group_by(obj,group)
            if nargin == 0
                groupby = "genotype";
            end
            [obj(:).iscontrol] = deal(false);
            
            genos = strcat(obj.get_full_genotype,'@',obj.get_full_protocol);
            [C,ia,ic] = unique(genos,'stable');
            groups = num2cell(ic,2);

            [obj(:).group] = groups{:};
            [obj(:).plot_order] = groups{:};
            
            groups = [groups{:}];
            for ii = unique(groups)
                idx = find(groups == ii);
                data_grouped(ii) = obj(idx(1));
                data_grouped(ii).filepaths = vertcat(obj(idx).filepaths);
                data_grouped(ii).date = [obj(idx).date];
                data_grouped(ii).time = [obj(idx).time];
                data_grouped(ii).data_choreography = vertcat(obj(idx).data_choreography);
            end
            obj = data_grouped;
        end
        
        %% Apply Filter (by time)
        function obj = apply_time_filter(obj,tStart,tEnd,tDur,tBuffer)
            if isnan(tStart)
                tStart = floor(min(arrayfun(@(x) min([x.data_choreography(:).elapstime]), obj)));
            end
            if isnan(tEnd)
                tEnd = ceil(max(arrayfun(@(x) max([x.data_choreography(:).elapstime]), obj)));
            end
            
            for ii = 1:length(obj)
                temp = obj(ii).data_choreography;
                
                et = {temp.elapstime};
                time_filt = cellfun(@(x) x >= tStart & x <= tEnd, et, 'UniformOutput', false);
                
                %%%
                durs = cellfun(@(x,y) max([nan,range(x(y))]), et, time_filt);
                if isnan(tDur)
                    ani_filt = repelem({true}, numel(et));
                else
                    ani_filt = num2cell(durs > tDur + tBuffer);
                    time_filt = cellfun(@(x,y) y & x >= [min(x(y))+tBuffer] & x <= [min(x(y))+tDur+tBuffer], et, time_filt,...
                        'UniformOutput', false, 'ErrorHandler', @(x,y,z) repelem(false,1,length(y)));
%                     time_filt = cellfun(@(x,y) y & x <= [min(x(y))+tDur], et, time_filt,...
%                         'UniformOutput', false, 'ErrorHandler', @(x,y,z) repelem(false,1,length(y)));
                end
                                                
                %%%
%                 filt = cellfun(@(x,y) ~isempty(x(y)), et, time_filt);
% 
%                 time_filt(filt) = cellfun(@(x,y) y & x <= [min(x(y))+tDur], et(filt), time_filt(filt), 'UniformOutput', false, 'ErrorHandler', @(x,y,z) repelem(false,1,length(y)));
%                 filt2 = filt;
% 
%                 filt2(filt) = cellfun(@(x,y) range(x(y))>(tDur-1), et(filt), time_filt(filt));
%                 
%                 ani_filt = num2cell(filt2);
                %%%
                [temp(:).animal_filter] = ani_filt{:};
                [temp(:).time_filter] = time_filt{:};
                
                obj(ii).data_choreography = temp;
            end
        end
        
        %%
        function obj = format_title(obj)
            filename = arrayfun(@(x) x.filepaths(1), obj)';
            
            expr = '\d\d\d\d\d\d\d\d_\d\d\d\d\d\d[@]([\w()]*[@][\w()]*)';

            token = regexp(filename,expr,'tokens','once');
            names = string(token);

            names = strrep(names,'None','');

            expr = '([@][a-zA-Z0-9]*)[_](\w*)';
            exprrep = '$1^{$2}';
            namesfix = regexprep(names,expr,exprrep);

            expr = '(^[a-zA-Z0-9]*)([_])(\w*[@])';
            exprrep = '$1; $3';
            namesfix = regexprep(namesfix,expr,exprrep);

            namesfix = strrep(namesfix,'@','>');
            namesfix = strrep(namesfix,'_','-');

            namesfix = regexprep(namesfix,'SM2>','SM2-null>');
            expr = '(SM2)[-]([\w()-]*)([>])';
            exprrep = 'dEaat1^{$2}$3';
            namesfix = regexprep(namesfix,expr,exprrep);

            namesfix = strrep(namesfix,'-NEDA','::Venus');
            namesfix = strrep(namesfix,'SM2','dEaat1^{null}');
            namesfix = strrep(namesfix,'>S103A','>hEAAT1^{S103A}');
            namesfix = strrep(namesfix,'>WT','>hEAAT1^{WT}');
            namesfix = strrep(namesfix,'>M128R','>hEAAT1^{M128R}');

            titlename = num2cell(namesfix);

            [obj.genotype_display] = titlename{:};
        end
        
        %% get full genotype
        function full_genotype = get_full_genotype(obj)
            full_genotype = strcat([obj.driver],'@',[obj.effector]);
        end

        %% get full protocol
        function full_protocol = get_full_protocol(obj)
            full_protocol = strcat([obj.protocol1],...
                '#',[obj.protocol2],...
                '#',[obj.protocol3],...
                '#',[obj.protocol4]);
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
        
        %% Initialise figure directories
        function figure_dir = init_figure_directory(obj)
            filelist = vertcat(obj(:).filepaths);
            filelist = strrep(filelist,pwd,'.');
            n_timestamps = length([obj.time]);

            [dmin,dmax] = bounds(str2double(unique([obj(:).date])));
            full_protocol = strcat(obj(1).protocol1,"#",obj(1).protocol2,"#",obj(1).protocol3,"#",obj(1).protocol4,"@x",string(n_timestamps));
            figure_dir = fullfile(".","figures",strcat(string(dmin), "_", string(dmax),"@",full_protocol));
            if ~isdir(figure_dir)
                mkdir(figure_dir)
            end

            fileID = fopen(fullfile(figure_dir,'experiment_list.txt'),'w');
            fprintf(fileID,'%s\n',filelist);
            fclose(fileID);
        end
        
        %% Boxplots
        function [printArray,box_handle,spread_handle] = plot_boxplot(obj,metric_to_plot,y_label,params)
            genos = [obj.genotype_display];

            dists = arrayfun(@(x) [x.data_choreography([x.data_choreography.animal_filter]).(metric_to_plot)], obj,'UniformOutput', false);
            groups = arrayfun(@(x) [x.data_choreography([x.data_choreography.animal_filter]).group], obj,'UniformOutput', false);
            objectID = arrayfun(@(x) [x.data_choreography([x.data_choreography.animal_filter]).aniID], obj,'UniformOutput', false);
            dists = [dists{:}];
            groups = [groups{:}];
            objectID = [objectID{:}];
            [unique_groups,~,alt_groups] = unique(groups,'stable');
            genos = genos(unique_groups);
            if max(groups) > length(unique_groups)
                groups = alt_groups';
                genos = genos(unique_groups);
            end
            geno_n = accumarray(groups',1);
            
            if length(unique(genos)) == length(params.group_order)
                [B,I] = sort(params.group_order);
                groups = I(groups);
                genos = genos(params.group_order);
                geno_n = geno_n(params.group_order);
            end
            
            figure()
            h = boxplot(dists,groups,'Colors',[.3 .3 .3],'Symbol','rx','Labels',genos,'LabelOrientation','inline','symbol','','GroupOrder',string([1:length(genos)]));
            hold on
            p = plotSpread(dists,'distributionIdx',groups,'distributionColors',[100/255 100/255 255/255]);
            
            hold off
            pbaspect([1,1,1])
            ylabel(metric_to_plot)
            
            xlim = get(gca,"Xlim");
            xpos = ([1:length(genos)]/range(xlim)) +min(xlim);
            text(xpos,repelem(0.03,1,length(genos)),num2str(geno_n),'Units','normalized','HorizontalAlignment','center');
            set(h,'LineWidth',2);
            obs = findobj(gca,'Marker','.');
            set(obs,'MarkerSize',12);
            set(gca,'XTickLabel',genos,'XTickLabelRotation',45,'LineWidth',2);
            ax = gca;
            ax.YLabel.String = y_label;
            
            % generate table to print
            genotypes = obj.get_full_genotype();
            printArray = [cellstr(genotypes(groups));num2cell(groups);num2cell(objectID);num2cell(dists)]';
            
            box_handle = h;
            spread_handle = p;
        end
        
        %%
        function fig_handle = plot_paths(obj, gridsize)
            temp = obj.data_choreography;
            if ~any([temp.animal_filter])
                fig_handle = figure();
                return
            end
            
            if ~exist("gridsize")
                gridsize = 3;
            end
            
            maxspeed = 0.4;
            distspace = 6;
            maxlim = distspace*(gridsize+1);

            f = find([temp.animal_filter]);
            if length(f) > 9
                f = randsample(f,9);
            end
            nlines = length(f);
            cols = ceil(nlines/3);
            rows = ceil(nlines/cols);

            idx = [repmat([1:cols]',rows,1), repelem(flipud([1:rows]'),cols,1)];

            cmap = magma;
            cd = [uint8(viridis(64)*255) uint8(ones(64,1))].';
            figure()
            hold on
            n = 1;
            for jj = f
                et = temp(jj).elapstime(temp(jj).time_filter);
                x = temp(jj).x_smooth;
                y = temp(jj).y_smooth;

                xy = [x;y]';
                xy = xy-mean(xy)+([distspace distspace].*idx(n,:));

                bins = linspace(0,maxspeed,64);
                speeds = temp(jj).distance./[0 diff(et)];
                Y = discretize(speeds,bins);
                Y(isnan(Y)) = 64;

                p = scatter(xy(:,1),xy(:,2),20,'.');
                p.CData = cd(1:3,Y)';       
                fig_handle(n) = p;
                n = n+1;
            end
            hold off

            axis equal
            colormap(viridis)
            xlim([0 maxlim])
            ylim([0 maxlim])
            pbaspect([cols,rows,1])

            titl = obj.genotype_display;
            title(titl)
            c = colorbar;
            caxis([0 maxspeed])

            set(gca,'XTick',[0 maxlim],'YTick',[0 maxlim]);
            xlabel('Distance (mm)')
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%
        function plot_paths_full(obj,params,FullDuration)
            temp = obj.data_choreography;
            if ~any([temp.animal_filter])
                return
            end

            if isdir('./.temp')
                rmdir('./.temp','s')    
            end
            mkdir('./.temp')
            
            object_num = find([temp.animal_filter]);
            [~,timestamp_idx,~] = unique({temp.timestamp});
            maxspeed = 0.4;

            cols = 4;
            rows = 5;

            gridsize = 8;
            ylims = [0 rows].*gridsize;
            xlims = [0 cols].*gridsize;
            indices = [repmat([1:cols],1,rows);repelem(fliplr([1:rows]),cols)]';
            centres = [[1:max([cols,rows])].*8];
            line_xy = centres(indices) - gridsize/2;
            text_xy = centres(indices) - [gridsize,0];

            fig = figure('Units','inches','Visible','off');
            set(fig,'Position',[0,0,8.5,11]);
            ax = axes('Units', 'inches');
            set(ax,'InnerPosition',[0.25,0.25,8,10])

            cmap = magma;
            cd = [uint8(viridis(64)*255) uint8(ones(64,1))].';

            n = 1;
            ix = 1;
            gridelems = [rows*cols];
            for ii = 1:ceil(numel(object_num)/gridelems)
                objects_remain = numel(object_num) - [ii-1]*gridelems;
                objects_remain = min([objects_remain,gridelems]);
                hold on
                for jj = 1:objects_remain
                    ob = object_num(ix);
                    if ~FullDuration
                        et = temp(ob).elapstime(temp(ob).time_filter);
                        x = temp(ob).x_smooth;
                        y = temp(ob).y_smooth;
                        speeds = temp(ob).distance./[0 diff(et)];
                        id = temp(ob).aniID;
                    else
                        et = temp(ob).elapstime;
                        x = movmean(temp(ob).x,params.smooth_window);
                        y = movmean(temp(ob).y,params.smooth_window);
                        id = temp(ob).aniID;
                        dist = sqrt(sum([[0;0],diff([x;y],[],2)].^2));
                        speeds = dist./[0 diff(et)];
                    end
                    etR = [ceil(et(1)*10):floor(et(end)*10)]./10;
                    x = interp1(et,x,etR);
                    y = interp1(et,y,etR);
                    speeds = interp1(et,speeds,etR);
                    et = etR;
                    
                    xy = [x;y]';
                    xy = xy-mean(xy)+line_xy(n,:);

                    bins = linspace(0,maxspeed,64);
                    Y = discretize(speeds,bins);
                    Y(isnan(Y)) = 64;

                    p = scatter(xy(:,1),xy(:,2),20,'.');
                    p.CData = cd(1:3,Y)';
                    fig_handle(n) = p;
                    
                    if ismember(ix,timestamp_idx)
                       txt = strcat(string(id), ' - ', temp(ob).timestamp);
                    else
                        txt = string(id);
                    end
                    text(text_xy(n,1),text_xy(n,2),txt,...
                        'HorizontalAlignment','left',...
                        'VerticalAlignment','top',...
                        'Interpreter','none');

                    n = n+1;
                    ix = ix+1;
                end
                %%% Add scale
                offset = 2;
                len = 2;
                line_x = [0,0;0,1]*len + offset;
                line_y = [0,1;0,0]*len + offset;
                lab = sprintf('%.0fmm',len);
                
                line(line_x',line_y','Color','black', 'LineWidth',2)
                text(mean(line_x(1,:)),mean(line_y(1,:)),lab,...
                    'HorizontalAlignment','right',...
                    'VerticalAlignment', 'middle');
                text(mean(line_x(2,:)),mean(line_y(2,:)),lab,...
                    'HorizontalAlignment','center',...
                    'VerticalAlignment', 'top');
                hold off
                
                %%% Save current figure
                colormap(viridis)
                xlim(xlims)
                ylim(ylims)

                titl = obj.genotype_display;
                title(titl)

                set(gca,'XTick',[0 xlims(2)],'YTick',[0 ylims(2)]);
                ax.XAxis.Visible = 'off';
                ax.YAxis.Visible = 'off';

                filename = ['./.temp/tmp_',sprintf('%04.0f',ii)];
                print(gcf,filename,'-dpdf','-painters');

                delete(fig_handle)
                delete(findobj(gcf,'Type','Text'))
                delete(findobj(gcf,'Type','Line'))
                n = 1;
            end
            delete(ax)
            
            c = colorbar();
            caxis([0 maxspeed]);
            c.Label.String = 'Speed (mm/s)';
            pbaspect([1,1,1]);
            filename = ['./.temp/tmp_',sprintf('%04.0f',ii+1)];
            print(gcf,filename,'-dpdf','-painters');
            
            close all
        end
        
        %%
        function plot_ridgeline(obj)
            if ~any([obj.data_choreography.animal_filter])
                return
            end

            filter = find([obj.data_choreography(:).animal_filter]);
            speed = {obj.data_choreography(filter).speed};
            et = {obj.data_choreography(filter).elapstime};
            %%% Filter by filter window
%             filt = {obj.data_choreography(filter).time_filter};
%             speeds = cellfun(@(x,y) x(y),speed,filt,'UniformOutput', false);
%             et = cellfun(@(x,y) x(y),et,filt,'UniformOutput', false);
            %%% Get full time tracked
            speeds = speed;
            %%%
            speeds_norm = cellfun(@(x) movmean(x,21).*5, speeds, 'UniformOutput', false);
            
            len_obj = length(speeds_norm);
            cmap = flipud(viridis(len_obj));
            hold on
            for jj = 1:len_obj
                y = (-jj)+[speeds_norm{jj},0,0];
                x = [et{jj},max(et{jj}),min(et{jj})];
%                 patch(x,y,'red','EdgeColor','white','LineWidth',0.1,'FaceColor',cmap(jj,:));
            patch(x,y,'red','EdgeColor','white','LineWidth',0.05,'FaceColor',cmap(jj,:));
%                 patch(x,y,'red','EdgeColor','none','FaceColor',cmap(jj,:));
            end
            hold off
            set(gcf,'PaperOrientation','landscape');
            pbaspect([3,1,1])

            ylim([-len_obj,0]);
            xlim([0 max(cellfun(@max,et))]);
            set(gca,'YTickLabels',[],'YTick',[]);
        end
        
    end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% METHODS - Static
    methods (Static)
        function params = load_parameters(config_file)
            if nargin == 0
                parameterFile = 'default_config.yaml'; 
                fprintf("\n No config file specified - loading deafault config\n")
            else
                parameterFile = config_file;
            end
            
            paramFile = dir(fullfile('.','**',[parameterFile,'*']));
            params = yaml.ReadYaml(fullfile(paramFile.folder,paramFile.name));
            try
                cd(params.directories.master)
            catch
                error('Non-existant directory specified - please update config.yaml file')
            end
            addpath(params.directories.code)
        end
        
        %%
        function file_selection = get_filelist(params)
            files = dir(fullfile(params.directories.choreography_output,'**'));
            files = files(~[files.isdir]);
            [fs,idx] = choreHandler.get_file_details(files);
            files = files(idx);
            
            flist = [fs.date; fs.full_genotype]';
            [filelist,I] = sortrows(flist,[1,2],["descend","ascend"]);
            files = files(I);
            [uni_files,ia,ic] = unique(filelist,'rows');
            files_prompt = strcat(uni_files(:,1),'@',uni_files(:,2));

            [indx,tf] = listdlg('ListString',files_prompt);
            filt = ismember(ic,indx);
            filelist = filelist(filt,:);
            file_selection = files(filt);
       end
        %%
        function [file_struct,index_goodfiles] = get_file_details(filelist)
            filenames = string({filelist.name});
            expr = ['(?<date>^\d\d\d\d\d\d\d\d)[_]'...
                '(?<time>\d\d\d\d\d\d)[@]'...
                '(?<driver>[\w()]*)[@]'...
                '(?<effector>[\w()]*)[@]'...
                '(?<tracker>[t]\d*)[@#]+'...
                '(?<prot1>\w*)'...
                '[#](?<prot2>\w*)'...
                '[#](?<prot3>\w*)'...
                '[#](?<prot4>\w*)[@]'...
                '\d*[@]'...
                '[.](?<filetype>\w*)[.]'];

            file_struct = regexp(filenames,expr,'names','forcecelloutput');
            index_goodfiles = ~cellfun(@isempty, file_struct);
            file_struct = vertcat(file_struct{:});

            fullgeno = strcat([file_struct.driver],'@',[file_struct.effector]);
            fulltimestamp = strcat([file_struct.date],'_',[file_struct.time]);
            for ii = 1:length(file_struct)
                file_struct(ii).full_genotype = fullgeno(ii);
                file_struct(ii).full_timestamp = fulltimestamp(ii);
            end
        end
        
        %%
        function [distance_by_area,index] = get_distance_by_area(blob,scaler)
            if ~blob.animal_filter
                distance_by_area = nan;
                index = 1;
                return
            end
            
            x = [blob.x_smooth];
            y = [blob.y_smooth];
            area = mean(blob.area);
            rad = sqrt(scaler*area/pi);

            tf = true;
            n=1;
            idx = 1;
            I = 1;
            while tf
                xy = [x(I:end); y(I:end)] - [x(I); y(I)];
                dist = sqrt(sum(xy.^2,1));
                first_idx = find(dist > (rad),1);
                if isempty(first_idx)
                    tf = false;
                    idx(n+1) = length(x);
                    continue
                end
                I = I - 1 + first_idx;
                n = n+1;
                idx(n) = I;
            end
            index = idx;
            distance_by_area = sum(sqrt(sum(diff([x(idx); y(idx)],[],2).^2)));
        end
        
    end
end