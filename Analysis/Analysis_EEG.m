%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% ClickSnr %%%%%%%%%%%%%%%%
%%%%%%%%%%% Data Analysis - EEG %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Paths + Preprocessing parameters

%%%%%%%%%%%%%%%%%%
% Initial setups %
%%%%%%%%%%%%%%%%%%

    clear ; clc ; close

    % Load preprocessed data
    Paths.computer_username     = 'emmanuelcoulon';
    Paths.projectName           = 'ClickSnrFF';
    load(['/Users/',Paths.computer_username,'/Documents/MATLAB/PROJECTS/',Paths.projectName,'/Preprocessing/EEG_preprocessed.mat'])  

    % Add the toolboxes to Matlab path
    if exist(Paths.toolboxes,'dir')
        addpath(genpath(Paths.toolboxes))
        LW_init(); 
    else
        warning('provide valid path to the toolboxes folder');
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%
% Analysis parameters %
%%%%%%%%%%%%%%%%%%%%%%%

    Cfg.subjects            = fieldnames(EEG_preprocessed); % Subject names
    subjects                = Cfg.subjects;                 
    Cfg.fs                  = 512;                          % Sampling frequency
    Cfg.percentSign         = 0.05;                         % p-value
    Cfg.sign                = norminv(1 - Cfg.percentSign); % Transformation of the p-value in zscore
    Cfg.frex                = 1/(0.2*12) * [1:12];          % All meter-(un)related frequencies
    Cfg.whichMeterRel       = [3,6,9,12];                   % Meter-related frequencies
    Cfg.whichMeterUnrel     = [1,2,4,5,7,8,10,11];          % Meter-unrelated frequencies
    Cfg.frontalPool         = {'F1','Fz','F2','FC1','FCz','FC2','C1','Cz','C2'};
    Cfg.frontalPoolStr      = join (string (Cfg.frontalPool));
    Cfg.bl_snr.lowBin       = 2;
    Cfg.bl_snr.highBin      = 5;
    Cfg.lowPass.filterType  = 'lowpass';
    Cfg.lowPass.cutoff      = 30;
    Cfg.lowPass.filterOrder = 4;
    Cfg.chunkSegm.onset     = 0;
    Cfg.chunkSegm.duration  = 2.4;
    Cfg.chunkSegm.interval  = 2.4;
    Cfg.LwAnalysisSave      = 0;
    
    
    % Add the "-" in the "sub-00x" name (useful to navigate through the BIDS directories)
    for iSubjects = 1:length(subjects)
        Cfg.subDir(iSubjects) = insertAfter(subjects(iSubjects),"sub","-");
    end

    % Let the user decide whether he would like to launch the analysis for the fft, the ERPS or both
    Question.AnalysisType = questdlg('Which analysis would you like to process?',...
                                     'Analysis','Fft','ERP','Both','Both');



%% Spectral analysis


if or(strcmp(Question.AnalysisType,'Fft') == 1, ...
      strcmp(Question.AnalysisType,'Both') == 1) 
    
%%%%%%%    
% FFT %
%%%%%%%

    
    for iSubjects = 1:length(Cfg.subDir)

        % Extraction of the data
        lwdata = EEG_preprocessed.(subjects{iSubjects});

        % Save the fft_analysis in the lw folder of each participant
        if ~exist (Paths.analysisLw, 'dir');     mkdir(Paths.analysisLw); end
        cd(Paths.analysisLw)
        

        % FFT
        option  = struct('output','amplitude',...
                         'half_spectrum',1,...
                         'suffix','fft','is_save',Cfg.LwAnalysisSave);
        lwdata  = FLW_FFT.get_lwdata(lwdata,option);
        

        % Baseline correction
        option  = struct('xstart',Cfg.bl_snr.lowBin,...
                         'xend',Cfg.bl_snr.highBin,...
                         'num_extreme',0,...
                         'operation','subtract',...
                         'suffix','bl_snr','is_save',Cfg.LwAnalysisSave);
        lwdata  = FLW_baseline_SNR.get_lwdata(lwdata,option);
        

        % Put the lwdata in the EEG_fft structure
        EEG_fft.(subjects{iSubjects}).blfft      = lwdata;

        % Remove empty dimensions in the lw.data structure
        EEG_fft.(subjects{iSubjects}).blfft.data = squeeze(EEG_fft.(subjects{iSubjects}).blfft.data);
        
    end
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Means and pools of electrodes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
       
    for iSubjects = 1:length(Cfg.subDir)

       % Reload lwdata
       lwdata                                   = EEG_fft.(subjects{iSubjects}).blfft;
        
       % Get the number of electrodes
       Cfg.(subjects{iSubjects}).nElec          = length(EEG_preprocessed.(subjects{iSubjects}).header.chanlocs);
       nElec                                    = Cfg.(subjects{iSubjects}).nElec;

       % Get the mean fft across all electrodes
       lwdata.data(nElec+1,:)                   = mean (EEG_fft.(subjects{iSubjects}).blfft.data,1);
       lwdata.header.chanlocs(nElec+1).labels   = 'MeanOfAllElectrodes';  
           
       % Get the mean fft of the prefrontal electrodes pool
       Cfg.elecLabels                           = {lwdata.header.chanlocs(1:end-1).labels};             % Get electrode labels (without the added 'mean')
       Cfg.frontalPoolIdx                       = find (ismember(Cfg.elecLabels, Cfg.frontalPool));     % Find the index corresponding to these electrodes
       lwdata.data (nElec+2,:)                  = mean (lwdata.data(Cfg.frontalPoolIdx,:),1); 
       lwdata.header.chanlocs(nElec+2).labels   = Cfg.frontalPoolStr;
       
       % Save new version of lwdata
       EEG_fft.(subjects{iSubjects}).blfft      = lwdata;
       tempLwdata (:,:,iSubjects)               = lwdata.data;
       
       % Update the Cfg.elecLabels with the new pools
       for iChan = 1:size (EEG_fft.(subjects{iSubjects}).blfft.header.chanlocs,2)    
            Cfg.elecLabels {iChan} = lwdata.header.chanlocs(iChan).labels;
       end
       
    end
        
    % Get the mean of all participants
    EEG_fft.all.blfft.data = mean(tempLwdata,3);
    subjects{end+1}        = 'all'; 

    clear tempLwdata
    
    
    
    
    
%%%%%%%%%%%  
% Zscores %
%%%%%%%%%%%

EEG_zscores=[];
   
    for iSubjects = 1:length(subjects) 

        if iSubjects <= length(Cfg.subDir)

            % Create the frequency array
            datasize        = EEG_fft.(subjects{iSubjects}).blfft.header.datasize(6);   % Get the number of frequency datapoints
            Cfg.freqRes     = EEG_fft.(subjects{iSubjects}).blfft.header.xstep;         % Frequency resolution
            Cfg.freq        = [0 : datasize-1] * Cfg.freqRes;                           % Frequency array

            % Find the frequency index corresponding to the 12 frequencies of interest
            Cfg.frexidx     = dsearchn(Cfg.freq', Cfg.frex');   

        end

        for iChan = 1:nElec +2

            % Extraction of the amps for the 12 frex of interest
            EEG_zscores.(subjects{iSubjects}).rawfftAmps(iChan,:) ...
                        = EEG_fft.(subjects{iSubjects}).blfft.data(iChan,Cfg.frexidx);

            % Zscore
            EEG_zscores.(subjects{iSubjects}).zscores(iChan,:) ...
                        = zscore(EEG_zscores.(subjects{iSubjects}).rawfftAmps(iChan,:));

            % Get meter-(un)related frequencies
            EEG_zscores.(subjects{iSubjects}).zscoresMetRel(iChan,:) ...
                        = EEG_zscores.(subjects{iSubjects}).zscores(iChan,Cfg.whichMeterRel) ;
            EEG_zscores.(subjects{iSubjects}).zscoresMetUnrel(iChan,:) ...
                        = EEG_zscores.(subjects{iSubjects}).zscores(iChan,Cfg.whichMeterUnrel) ;

            % Get the mean of the meter-(un)related frequencies
            EEG_zscores.(subjects{iSubjects}).meanZscoresMetRel(iChan,:) ...
                        = mean (EEG_zscores.(subjects{iSubjects}).zscoresMetRel(iChan,:),2);
            EEG_zscores.(subjects{iSubjects}).meanZscoresMetUnrel(iChan,:) ...
                        = mean (EEG_zscores.(subjects{iSubjects}).zscoresMetUnrel(iChan,:),2) ;


            %%% Repeat the zscore procedure without taking the frequency associated with the unitary event into account

            % Zscore (witout 12)
            EEG_zscores.(subjects{iSubjects}).zscoresNo12(iChan,:) ...
                        = zscore(EEG_zscores.(subjects{iSubjects}).rawfftAmps(iChan,1:end-1));

            % Get meter-(un)related frequencies (witout 12)
            EEG_zscores.(subjects{iSubjects}).zscoresNo12MetRel(iChan,:) ...
                        = EEG_zscores.(subjects{iSubjects}).zscoresNo12(iChan,Cfg.whichMeterRel(1:end-1)) ;
            EEG_zscores.(subjects{iSubjects}).zscoresNo12MetUnrel(iChan,:) ...
                        = EEG_zscores.(subjects{iSubjects}).zscoresNo12(iChan,Cfg.whichMeterUnrel) ;

            % Get the mean of the meter-(un)related frequencies (witout 12)
            EEG_zscores.(subjects{iSubjects}).meanZscoresNo12MetRel(iChan,:) ...
                        = mean (EEG_zscores.(subjects{iSubjects}).zscoresNo12MetRel(iChan,:),2);
            EEG_zscores.(subjects{iSubjects}).meanZscoresNo12MetUnrel(iChan,:) ...
                        = mean (EEG_zscores.(subjects{iSubjects}).zscoresNo12MetUnrel(iChan,:),2) ;
        end  
    end  
end
    

%% Temporal analysis

if or(strcmp(Question.AnalysisType,'ERP') == 1,...
      strcmp(Question.AnalysisType,'Both') == 1)
    
%%%%%%%%    
% ERPs %
%%%%%%%%


    for iSubjects = 1:length(Cfg.subDir) % Not take the last subject which reprensents the mean of all participants

        % Extraction of the data
        lwdata = EEG_preprocessed.(subjects{iSubjects});

        % Save the ERP_analysis in the lw folder of each participant
        cd(Paths.analysisLw)

        % Low-pass Butterworth 
        option  = struct('filter_type',Cfg.lowPass.filterType,...
                         'high_cutoff',Cfg.lowPass.cutoff,...
                         'filter_order',Cfg.lowPass.filterOrder,...
                         'suffix','butt','is_save',Cfg.LwAnalysisSave);
        lwdata  = FLW_butterworth_filter.get_lwdata(lwdata,option);

        % Segmentation in succesive chunks
        option  = struct('chunk_onset',Cfg.chunkSegm.onset,...
                         'chunk_duration',Cfg.chunkSegm.duration,...
                         'chunk_interval',Cfg.chunkSegm.interval,...
                         'suffix','chunk','is_save',Cfg.LwAnalysisSave);
        lwdata  = FLW_segmentation_chunk.get_lwdata(lwdata,option);

        % Average
        option  = struct('operation','average','suffix','avg','is_save',Cfg.LwAnalysisSave);
        lwdata  = FLW_average_epochs.get_lwdata(lwdata,option);

        % Put the lwdata in the EEG_fft structure + remove empty dimensions
        EEG_ERP.(subjects{iSubjects})       = lwdata;
        EEG_ERP.(subjects{iSubjects}).data  = squeeze(EEG_ERP.(subjects{iSubjects}).data);

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Means and pools of electrodes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


    for iSubjects = 1:length(Cfg.subDir)

       % Reload lwdata
       lwdata = EEG_ERP.(subjects{iSubjects});
        
       % Get the number of electrodes
       Cfg.(subjects{iSubjects}).nElec          = length(EEG_preprocessed.(subjects{iSubjects}).header.chanlocs);
       nElec                                    = Cfg.(subjects{iSubjects}).nElec;

       % Get the mean ERP across all electrodes
       lwdata.data(nElec+1,:)                   = mean(lwdata.data,1);
       lwdata.header.chanlocs(nElec+1).labels   = 'MeanOfAllElectrodes';  

       % Get the mean fft of the prefrontal electrodes pool     
       Cfg.elecLabels                           = {lwdata.header.chanlocs(1:end-1).labels};             % Get electrode labels (without the added 'mean')
       Cfg.frontalPoolIdx                       = find (ismember(Cfg.elecLabels, Cfg.frontalPool));     % Find the index corresponding to these electrodes
       lwdata.data (nElec+2,:)                  = mean (lwdata.data(Cfg.frontalPoolIdx,:),1);
       lwdata.header.chanlocs(nElec+2).labels   = Cfg.frontalPoolStr;
       
       % Save new version of lwdata
       EEG_ERP.(subjects{iSubjects})            = lwdata;
       tempLwdata (:,:,iSubjects)               = lwdata.data;
       
       % Update the Cfg.elecLabels with the new pools
       for iChan = 1:size (EEG_fft.(subjects{iSubjects}).blfft.header.chanlocs,2)    
            Cfg.elecLabels {iChan} = lwdata.header.chanlocs(iChan).labels;
       end

    end

    % Get the mean of all participants
    EEG_ERP.all.data = mean(tempLwdata,3);


end; clear tempLwdata

%% Save

cd(Paths.analysis)

save('EEG_analysis.mat','Cfg','EEG_fft','EEG_zscores','EEG_ERP','Paths','-v7.3');

disp('EEG_analysis.mat saved')


