%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% ClickSnr %%%%%%%%%%%%%%%%
%%%%%%%%% Data Preprocessing - EEG %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paths + Preprocessing parameters

%%%%%%%%%%%%%%%%%%
% Initial setups %
%%%%%%%%%%%%%%%%%%

    clear; clc
    
    
    % Set main paths
    Paths.computer_username     = 'emmanuelcoulon'; %'florarosenzweig/Dropbox/Mon Mac (mac-114-318.local)';
    Paths.projectName           = 'ClickSnrFF';
    
    Paths.matlab                = fullfile ('/Users',Paths.computer_username,'Documents/MATLAB');
    Paths.toolboxes             = fullfile (Paths.matlab, 'Toolboxes');
    Paths.PTB                   = fullfile ('/Users',Paths.computer_username,'Documents/PTB');
    Paths.source                = fullfile (Paths.PTB, Paths.projectName, 'output/source');
    Paths.lw                    = fullfile (Paths.PTB, Paths.projectName, 'output/lw');
    Paths.ProjectPath           = fullfile (Paths.matlab, 'PROJECTS', Paths.projectName);
    Paths.EEGPreprocessedStruct = fullfile (Paths.ProjectPath, 'Preprocessing/EEG_preprocessed.mat');
    Paths.ElecLabels            = fullfile (Paths.ProjectPath, 'Preprocessing/electrodelabels.csv');
    Paths.analysis              = fullfile (Paths.ProjectPath,'Analysis');
    Paths.analysisLw            = fullfile (Paths.analysis,'LW');
    Paths.figures               = fullfile (Paths.ProjectPath,'Figures');

    % Add the toolboxes to Matlab path
    if exist(Paths.toolboxes,'dir')
        addpath(genpath(Paths.toolboxes))
        LW_init(); 
    else
        warning('provide valid path to the toolboxes folder');
    end
    
    % Load EEGPreprocessedStruct
    if isfile (Paths.EEGPreprocessedStruct)
        load(Paths.EEGPreprocessedStruct);
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocessing parameters %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Conditions
    Cfg.condNames           = {'ClickSnrFF' };
    condNames               = Cfg.condNames;
    
    %Butterworth filter
    Cfg.filter.filterType   = 'highpass';
    Cfg.filter.cutOff       = 0.1;
    Cfg.filter.Order        = 4;
    
    Cfg.trigCode            = {'5'}; % Look at the continuous view in Letswave
    
    % Rereference
    Cfg.Rereference.refElec   = {'Mast1','Mast2'};
    
    %BadEpochRejection
    Cfg.badEpoThr          = 100; % +/- ??V

    %ICA
    Cfg.ICA.numbIC         = 60; % # of IC
    
    %LwSave
    Cfg.lwSave             = 1; 
    
    
    
    
    
%% Load new data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get the source and lw output folders % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('-----------------------------\n get a participant''s directory (eg.: sub-001)\n')

    % Select the subject
    Paths.Sub.dir          = uigetdir (Paths.source);
    Paths.Sub.dataDir      = fullfile (Paths.Sub.dir,'ses-001/eeg');
    [~,Log.subject]        = fileparts(Paths.Sub.dir);
    
    fprintf(['--> ', Log.subject, ' selected'])
    
    % Select the lw output folder
    if ~exist (Paths.lw, 'dir');     mkdir(Paths.lw); end
    Paths.Sub.out         = fullfile (Paths.lw, Log.subject);
    MsgBox                = {'Would you like to save the preprocessing in this specific directory',...
                            ['Current Directrory : ',Paths.Sub.out]};
    outputPath            = questdlg (MsgBox, 'Analysis Directory','Yes','No','Yes');
    
    switch outputPath
        case 'No'
            waitfor (msgbox ('Please indicate where you would like to save the lw analysis'))
            Paths.Sub.out = uigetdir (Paths.Sub.out);
    end;    clear outputPath MsgBox
    
    % If the lw doesn't exist, create it
    if ~exist (Paths.Sub.out,'dir');    mkdir(Paths.Sub.out); end
    cd (Paths.Sub.out)
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% select and load the bdf file % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    % Select the right bdf file
    fprintf('\n\n-----------------------------\n get bdf file\n')
    [Paths.Sub.bdfName, Paths.Sub.bdfPath] = uigetfile({'*.bdf'}, 'pick a file', Paths.Sub.dataDir);
    
    % Load bdf
    FLW_import_data.get_lwdata('filename',Paths.Sub.bdfName,'pathname',Paths.Sub.bdfPath,'is_save',1);
    
    % Load the imported bdf in the workspace
    [~,Paths.Sub.bdfName,] = fileparts(Paths.Sub.bdfName); % remove the bdf extension from the name
    option                 = struct('filename',fullfile (Paths.Sub.out, [Paths.Sub.bdfName, '.lw6']));
    lwdata                 = FLW_load.get_lwdata(option);
    
    fprintf('--> bdf file loaded \n')
    

%% Preprocessing

fprintf('\n\n-----------------------------\n Start of the preprocessing\n')

Paths.preprocessinglw       = fullfile(Paths.lw,Log.subject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% High-Pass Butterworth filter %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    option = struct ('filter_type', Cfg.filter.filterType,...
                     'low_cutoff',Cfg.filter.cutOff,...
                     'filter_order',Cfg.filter.Order,...
                     'suffix','butt','is_save',Cfg.lwSave);
    lwdata = FLW_butterworth_filter.get_lwdata(lwdata,option);  

    fprintf ('--> Signal filtered \n')
    
    
%%%%%%%%%%%%%%%%%%%%%%%%   
% New electrode labels %
%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Open the electorelabels.csv from the preprocessing folder
    Paths.ElecLabels        = fullfile (Paths.ProjectPath, 'Preprocessing/electrodelabels.csv');
    option                  = delimitedTextImportOptions ("NumVariables", 1);
    option.DataLines        = [1, Inf]; % Get all the rows containing the electrode labels
    Cfg.elecLabels          = readtable (Paths.ElecLabels, option); 

    % Get the name associated to the old channels 
    for iChan = 1:size (Cfg.elecLabels,1) % Only run through the number of electrodes from the csv file (not the EXG3,EXG4,etc.)
        Cfg.oldChan {iChan} = lwdata.header.chanlocs(iChan).labels;
    end

    option = struct('old_channel',{Cfg.oldChan},...
                    'new_channel',{table2cell(Cfg.elecLabels)},...
                    'suffix','chanlabels','is_save',Cfg.lwSave);
    lwdata = FLW_electrode_labels.get_lwdata(lwdata,option);    

    fprintf('--> New electrode labels applied \n')
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Removal of unused channels %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    option = struct('type','channel','items',{table2cell(Cfg.elecLabels)}, ...
                    'suffix','sel_chan','is_save',Cfg.lwSave);
    lwdata = FLW_selection.get_lwdata(lwdata,option);

    fprintf('--> Unused channels removed \n')
    
  
%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Channels interpolation %
%%%%%%%%%%%%%%%%%%%%%%%%%%

    Cfg.Interpolation.quest1 = questdlg ('Would you like to interpolate channels?', ...
                                                   'Channel Interpolation','Yes','No','No');

    switch Cfg.Interpolation.quest1

        case 'No'
            disp('No channel interpolated')
            Log.Interpolation.InterpElec  = 'None';
            Log.Interpolation.ClosestElec = 'None';

        case 'Yes'
            letswave7
            pause(60) % let the user 60 seconds to do the interpolation
            Cfg.Interpolation.interpFinished = 'No';

            while strcmp(Cfg.Interpolation.interpFinished,'No') == 1

                Cfg.Interpolation.quest2 = questdlg('Have you finished the channel interpolation ?', ...
                                                              'Channel Interpolation','Yes','No','No');
                switch Cfg.Interpolation.quest2

                    case 'No'
                        pause(30)

                    case 'Yes'
                        % Find the latest file in the LW directory
                        lwDir = dir ('*.lw6');

                        for i = 1: size(lwDir,1)
                            time(i,1) = lwDir(i).datenum;            
                        end

                        [~,timeIdx] = max(time);

                        % if the latest file begins with 'chan_interp', load it
                        if strncmp (string (lwDir(timeIdx).name),'chan_interp', 11) == 1
                            option = struct('filename',fullfile(Paths.Sub.out,lwDir(timeIdx).name));
                            lwdata = FLW_load.get_lwdata(option);

                        else
                            warning ('The file that was created last is not a chan_interp file')                           
                            Log.Interpolation.InterpElec  = 'None';
                            Log.Interpolation.ClosestElec = 'None';
                        end

                        % Complete a dialogbox to indicate the parameters of the interpolation (further used in the ReadMe.txt)
                        chan_interp = inputdlg({'Enter the channel(s) that was/were interpolated', ...
                                                'Enter the number of closest electrodes selected'},'Interpolation');

                        Log.Interpolation.InterpElec  = chan_interp {1};
                        Log.Interpolation.ClosestElec = chan_interp {2};

                        Cfg.Interpolation.interpFinished = 'Yes'; % stops the while loop

                        fprintf('--> Channel(s) correctly interpolated \n')                
                end
            end  
    end

    close all force 
    clear chan_interp time timeIdx lwDir    


%%%%%%%%%%%%%%%%%%%%%%%
% Epochs Segmentation %
%%%%%%%%%%%%%%%%%%%%%%%

    % Get sequence durations
    Cfg.SequenceDir     = dir([Paths.Sub.dataDir,'/*.tsv']); 
    Cfg.Sequence        = struct2dataset (bids.util.tsvread(...
                                    fullfile(Cfg.SequenceDir.folder, Cfg.SequenceDir.name))); 
    Cfg.Sequence (any(isnan(Cfg.Sequence.onset),2),:) = []; % Delete NaNs

    % Check that all Sequence.duration are identicals
    for iSeq = 1:length(Cfg.Sequence)-1
        if ~isequal(Cfg.Sequence.duration(iSeq),Cfg.Sequence.duration(iSeq+1))
            warning('All sequence duration are not identical')
            error('Epoch segmentation cannot be correctly computed')
        end 
    end

    option      = struct('event_labels',Cfg.trigCode, ...
                         'x_start',0, ...
                         'x_end',Cfg.Sequence.duration(1), ...
                         'x_duration',Cfg.Sequence.duration(1), ...
                         'suffix','ep','is_save',Cfg.lwSave);
   
    lwdata   = FLW_segmentation.get_lwdata(lwdata,option);
    
    fprintf('--> Epochs segmented \n')
    
    
%%%%%%%%%%%%%%%  
% Rereference %
%%%%%%%%%%%%%%%

    Cfg.Rereference.list   = {'Mastoids','Mast Right', 'Mast Left', 'All'};
    Cfg.Rereference.answer = listdlg('PromptString',{'Which electrode to select as reference?'},...
                                     'SelectionMode','single','ListString',Cfg.Rereference.list);                                

    if strcmp(Cfg.Rereference.answer,'Mastoids')
        Cfg.Rereference.refElec   = {'Mast1','Mast2'};
        
    elseif strcmp(Cfg.Rereference.answer,'Mast Right')  
        Cfg.Rereference.refElec   = {'Mast2'};
        
    elseif strcmp(Cfg.Rereference.answer,'Mast Left')   
        Cfg.Rereference.refElec   = {'Mast1'};
        
    else
        Cfg.Rereference.refElec   = table2cell(Cfg.elecLabels);
    end
    
    
    Cfg.Rereference.applyElec = table2cell(Cfg.elecLabels);
    
    option      = struct('reference_list',{Cfg.Rereference.refElec}, ...
                         'apply_list',{Cfg.Rereference.applyElec}, ...
                         'suffix','reref','is_save',Cfg.lwSave);
    lwdata      = FLW_rereference.get_lwdata(lwdata,option);
       
    fprintf('--> Signal rereferenced \n')  
    

%%%%%%%%%%%%%%%%%% 
% Run Merged ICA %
%%%%%%%%%%%%%%%%%%

    % Find latest reref file and load it
    lwDir = dir ('reref*.lw6');
    
    if size (lwDir,1) == 1
         
        clear option
        option.filename{1} = fullfile (lwDir.folder,lwDir.name);
        lwdataset          = FLW_load.get_lwdataset(option);  
    else
        error('Please seperate the two analysis in different folders. There shouldn''t be two reref files in the folder...');
    end
    
    option      = struct('ICA_mode',2, ...
                         'algorithm',1, ...
                         'num_ICs',Cfg.ICA.numbIC, ...
                         'suffix','ica_merged','is_save',Cfg.lwSave);
    lwdataset   = FLW_compute_ICA_merged.get_lwdataset(lwdataset,option);

    fprintf('--> %i ICs calculated \n',Cfg.ICA.numbIC)



%%%%%%%%%%%%% 
% Apply ICA %
%%%%%%%%%%%%%

    fprintf('\n\n-----------------------------\n In Letswave7, manually apply the ICA matrices and delete the artefact components \n')

    letswave7
    pause(60)
    Cfg.ICAFinished = 'No';

    while strcmp(Cfg.ICAFinished,'No') ==1 

        Cfg.ICAQuest = questdlg ('Have you applied the ICA?','ICA','Yes','No','No');

        switch Cfg.ICAQuest
            case 'No'
                pause (30)
                
            case 'Yes'
                Cfg.ICAFile       = dir('sp_filter*.lw6'); % Get the name and path of the sp_filter files
                Cfg.ICAFinished   = 'Yes';
        end
    end

    
    close all force % closes Letswave7
    
    
%%%%%%%%%%%  
% Average %
%%%%%%%%%%%

    Cfg.ICADir = dir ([Paths.preprocessinglw,'/sp_filter*.lw6']);

    for iCond = 1: size(condNames,2)
        
        % Load sp_filter data
        option  = struct('filename',fullfile(Paths.preprocessinglw,Cfg.ICADir(iCond).name));
        lwdata  = FLW_load.get_lwdata(option);
        
        % Average the epochs
        option = struct('operation','average','suffix','avg','is_save',Cfg.lwSave);
        lwdata = FLW_average_epochs.get_lwdata(lwdata,option);
        
        % Put the lwdata in the EEG_preprocessed structure
        Log.subject = strrep(Log.subject,'-',''); 
        EEG_preprocessed.(condNames{iCond}).(Log.subject) = lwdata;
        
    end
    
    
    fprintf('--> Signal averaged \n') 
    
    
%% Write ReadMe file to detail the current analysis

    Log.ReadMe.commentQuest = questdlg ('Would you like to add a specific comment about this analysis?', ...
                                        'Add a comment','Yes','No','No');

    switch Log.ReadMe.commentQuest
        case 'No'
            Log.ReadMe.comment = {' '};
        case 'Yes'
            Log.ReadMe.comment = inputdlg ('Please enter your comment here','ReadMe comment');
    end

    ReadMe = fopen('ReadMe.txt','w');
    
    fprintf(ReadMe,['\nSummary of the Preprocessing \n---------------------------\n\n', ...
                    'Subject ID : ', Log.subject, '\n\n', ...
                    'Date of analysis : ', datestr(now),'\n\n', ...
                    'Comment : ', Log.ReadMe.comment{1},'\n\n\n', ...
                    'Butterworth filter  \n------------------ \n \t filter type : ',Cfg.filter.filterType  , '\n \t low cutoff : ',num2str(Cfg.filter.cutOff) , '\n\t filter order : ',num2str(Cfg.filter.Order),'\n\n', ...
                    'Removal of unused channels \n--------------------------\n\n', ...
                    'Channel interpolation   \n----------------------\n \t interpolated channel(s) : ', Log.Interpolation.InterpElec, '\n \t number of closest electrodes selected : ',Log.Interpolation.ClosestElec,'\n\n', ...
                    'Independent component analysis \n------------------------------ \n \t number of ICs : ', num2str(Cfg.ICA.numbIC),...
                    ]);
    
    fclose(ReadMe);    
    
         

%% Save the new EEG_prep

cd (fullfile(Paths.ProjectPath,'preprocessing'))
save('EEG_preprocessed.mat','EEG_preprocessed','Paths','Cfg', '-v7.3')

disp('EEG_preprocessed structure saved')

%% To do 

% - Regarder dans les composants et il y a souvent le heartbeat. On peut
% souvent le voir dans l'acc??l??rom??tre et ??a permet de facilement le
% retrouver dans les IC. En g??n??ral, c'est vers 1Hz et ??a a la forme d'une
% vague qui traverse tout le cr??ne. 

% - Bien mesurer avec l'acc??l??rom??tre

% - Quand on fait l'ICA, on peut faire de base plus de composants, surtout
% que nos participants sont bien calmes et puis on vire environ 5% des ICs.

% - Voir si on trouverait pas un moyen de d??finir l'ordre du preprocessing
% au d??but du script. Si l'on fait ??a, enlever les squeeze pour garder les
% dimensions n??cessaires

