%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% ClickSnr %%%%%%%%%%%%%%%%
%%%%%%%%%%% Data Plotting - EEG %%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Paths + Plotting parameters

%%%%%%%%%%%%%%%%%%
% Initial setups %
%%%%%%%%%%%%%%%%%%

    clear ; clc ; close
    
    global Cfg EEG_fft EEG_zscores EEG_ERP Paths

    % Load data
    Paths.computer_username     = 'emmanuelcoulon';
    Paths.projectName           = 'ClickSnrFF';
    
    load(['/Users/',Paths.computer_username,'/Documents/MATLAB/PROJECTS/',Paths.projectName,'/Analysis/EEG_analysis.mat'])  
 
    
%%%%%%%%%%%%%%%%%%%%%%%
% Plotting parameters %
%%%%%%%%%%%%%%%%%%%%%%%  

    Cfg.figure.fontsize         = 16 ;
    Cfg.figure.linewidth        = 1.5;
    
    subjects                    = Cfg.subjects; 
    fontsize                    = Cfg.figure.fontsize;
    linewidth                    = Cfg.figure.linewidth;
    
    % Type of analysis to plot
    plotIdx = questdlg ('Which data would you like to plot?','Plotting',...
                              'Spectral analysis','Temporal analysis','Both','Both');
                          
    % Subject to plot
    if or (strcmp(plotIdx,'Spectral analysis') == 1, strcmp(plotIdx,'Both') == 1)
        subList         = Cfg.subjects;
        % subList{end+1}  = 'All subjects';          
        subIdx          = listdlg('PromptString',{'Please select the subject(s) you',...
                                  'would like to plot'},'ListString',subList,'SelectionMode','multiple');           
    end
    
    % Electrode(s) to plot
    ElecIdx   = listdlg('PromptString',{'Please select the electrode(s) you would like to plot'},...
                        'ListString',Cfg.elecLabels,'SelectionMode','single');
    
 
%% Extraction of the appropriate zscores from the EEG_zscores structure

    for iSubjects = 1:length(subjects)

            % Meter-related frequencies
            zscores(1,iSubjects)        = EEG_zscores.(subjects{iSubjects}).meanZscoresMetRel(ElecIdx);
            % Meter-unrelated frequencies
            zscores(2,iSubjects)        = EEG_zscores.(subjects{iSubjects}).meanZscoresMetUnrel(ElecIdx); 

            % Meter-related frequencies (without 12th frequency)
            zscoresNo12(1,iSubjects)    = EEG_zscores.(subjects{iSubjects}).meanZscoresNo12MetRel(ElecIdx);
            % Meter-unrelated frequencies (without 12th frequency)
            zscoresNo12(2,iSubjects)    = EEG_zscores.(subjects{iSubjects}).meanZscoresNo12MetUnrel(ElecIdx);        
    end


%% Spectral analysis plots

if or(strcmp(plotIdx,'Spectral analysis')==1 , strcmp(plotIdx,'Both')==1)
    
    % Create the figure
    figFFT = figure('position',[1 1 1200 650],'Color',[1 1 1]);
    
    % Main title
    annotation ('textbox',[0.1,0.9,0.8,0.08],'String',['EEG - ', subList{subIdx}],...
                'fontsize',fontsize+4,'EdgeColor','none','FontWeight','bold','HorizontalAlignment','center'); 
   
     %%%%%%%%%%%
     % Zscores % 
     %%%%%%%%%%%
     
        axes('Position',[0.1 0.1 0.35 0.8])

        % Zscores with the 12th frequency
        plot ([1.1,1.9],squeeze(zscores(:,subIdx)),'LineWidth',linewidth,'Marker','o'); hold on
        
        % Zscores without the 12th frequency
        for iSubjects = subIdx
           scatter ([0.95,2.05],zscoresNo12(:,iSubjects),'LineWidth',linewidth,'Marker','o','SizeData',50,'jitter','on','jitterAmount',0.04)
        end
        
        % Plot layout
        box off
        set(gca,'LineWidth', linewidth,'fontsize',fontsize,'xlim',[0.75 2.25],'xcolor','none','ylim',[-1.2 1.2],'ytick',[-1 0 1])
        line(xlim(), [0,0], 'LineWidth', linewidth, 'Color', 'k');
        legend(subjects(1:end),'Position',[0.303846108967363,0.059951923076923,0.136923076923077,0.116875]) 
    
        
    %%%%%%%    
    % Fft %
    %%%%%%% 
    
        axes('Position',[0.5 0.1 0.35 0.8])
        
        if length(subIdx) > 1 
            data2plot = EEG_fft.all.blfft.data(ElecIdx,:);
        else
            data2plot = EEG_fft.(subjects{subIdx}).blfft.data(ElecIdx,:);
        end

        stem(Cfg.freq,data2plot,'marker','none','LineWidth',linewidth,'Color','k')
        hold on 
        stem(Cfg.frex(Cfg.whichMeterRel),data2plot(Cfg.frexidx(Cfg.whichMeterRel)),'r','marker','none','LineWidth',1.75)
        stem(Cfg.frex(Cfg.whichMeterUnrel),data2plot(Cfg.frexidx(Cfg.whichMeterUnrel)),'b','marker','none','LineWidth',1.75)
        
        % Plot layout
        box off
        set (gca,'Tickdir', 'out','fontsize',fontsize,'LineWidth',linewidth,'xlim',[0 7],'ylim',[-0.1 0.7]);
        xlabel ({'Frequency (in Hz)'},'fontsize',fontsize-1,'Position',[3.5,-0.240])
        ylabel ({'µV'},'fontsize',fontsize,'Position',[-1 0.3 -1])  
end

%% Save figure

cd(Paths.figure)
set(gcf,'PaperPositionMode','auto')
print('test.jpg','-djpeg','-r800')
                
%% To Do
% Trouver une solution pour les couleurs du scatter et ne pas utiliser du
% rouge et bleu comme celui de la fft
