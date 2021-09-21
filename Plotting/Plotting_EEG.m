%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% PitchChange %%%%%%%%%%%%%%%
%%%%%%%%%%% Data Plotting - EEG %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Paths + Plotting parameters

%%%%%%%%%%%%%%%%%%
% Initial setups %
%%%%%%%%%%%%%%%%%%

    clear ; clc ; close
    
    global Cfg EEG_fft EEG_zscores EEG_ERP Paths

   
    Paths.computer_username     = 'emmanuelcoulon';
    Paths.projectName           = 'ClickSnrFF';
    
    addpath(genpath(['/Users/',Paths.computer_username,'/Documents/MATLAB/PROJECTS/',Paths.projectName,'/Plotting']))

%%%%%%%%%%%%%%%%%%%%%%%
% Select data to plot %
%%%%%%%%%%%%%%%%%%%%%%%  

  % Type of analysis
    plotIdx = questdlg ('Which data would you like to plot?','Plotting',...
                              'Spectral analysis','Temporal analysis','Both','Both');
                          
  % Load the corresponding data
    if or (strcmp(plotIdx,'Spectral analysis') == 1,...
           strcmp(plotIdx,'Both') == 1)

           disp(' --> Loading spectral analysis') 
           load(['/Users/',Paths.computer_username,'/Documents/MATLAB/PROJECTS/',Paths.projectName,'/Analysis/SpectralAnalysis.mat'])        
    end  

    if or (strcmp(plotIdx,'Temporal analysis') == 1,...
           strcmp(plotIdx,'Both') == 1)

           disp('--> Loading temporal analysis')
           load(['/Users/',Paths.computer_username,'/Documents/MATLAB/PROJECTS/',Paths.projectName,'/Analysis/TemporalAnalysis.mat'])     
    end
                          
  % Subject(s)
    if or (strcmp(plotIdx,'Spectral analysis') == 1, ...
           strcmp(plotIdx,'Both') == 1)
       
           subList         = Cfg.subjects;        
           subIdx          = listdlg('PromptString',{'Please select the subject(s) you',...
                                  'would like to plot'},'ListString',subList,'SelectionMode','multiple'); 
    end
    
  % Electrode(s)
    ElecIdx   = listdlg('PromptString',{'Please select the electrode(s) you would like to plot'},...
                        'ListString',Cfg.elecLabels,'SelectionMode','single');
    
   
%%%%%%%%%%%%%%%%%%%%%%%
% Plotting parameters %
%%%%%%%%%%%%%%%%%%%%%%%  

    Cfg.figure.fontsize         = 16 ;
    Cfg.figure.linewidth        = 1.5;
    
    subjects                    = Cfg.subjects; 
    fontsize                    = Cfg.figure.fontsize;
    linewidth                   = Cfg.figure.linewidth;
    
    condNames                   = Cfg.condNames;
    
    % Define figure size depending on the number of conditions
    if size(condNames,2) == 1
        Cfg.figPos = [1 1 400 250];
    elseif size(condNames,2) == 4
        Cfg.figPos = [1 1 1200 650];   
    end


%% Spectral Analysis

if or (strcmp(plotIdx,'Spectral analysis') == 1,...
      strcmp(plotIdx,'Both') == 1)
  
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % figure and main title %
    %%%%%%%%%%%%%%%%%%%%%%%%%

    figFFT = figure('position',Cfg.figPos,'Color',[1 1 1]);

    if length(subIdx) == 1
        Cfg.figure.subTitle = subList{subIdx};
    elseif length(subIdx) == length(subList)
        Cfg.figure.subTitle = 'All subjects';
    else
        Cfg.figure.subTitle = ['Subjects ', num2str(subIdx)];
    end

    annotation ('textbox', [0.1,0.9,0.8,0.08], ...
                'String',strjoin(['EEG - ', Cfg.figure.subTitle, ' - Electrode(s): ', Cfg.elecLabels{ElecIdx}]),...
                'fontsize',fontsize+3, 'EdgeColor','none', 'FontWeight','bold', 'HorizontalAlignment','center'); 

    %%%%%%%%%%%%%%%%%%%%%%%          
    % Plot all conditions %
    %%%%%%%%%%%%%%%%%%%%%%%
    
    for iCond = 1: size(condNames,2) 

      % Zscores
        zscoresPlot (iCond,subIdx,ElecIdx)
      % FFT
        fftPlot (iCond,subIdx,ElecIdx)     
    end
end




%% Temporal Analysis



