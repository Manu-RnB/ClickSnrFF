


%% Intro

clear ; clc
% close all force % closes Letswave7 if previously open

    Paths.computer_username     = 'emmanuelcoulon';
    Paths.projectName           = 'ClickSnrFF';
    
    load(['/Users/',Paths.computer_username,'/Documents/MATLAB/PROJECTS/',Paths.projectName,'/Analysis/EEG_fft.mat'])  
    load(['/Users/',Paths.computer_username,'/Documents/MATLAB/PROJECTS/',Paths.projectName,'/Analysis/EEG_zscores.mat'])
    load(['/Users/',Paths.computer_username,'/Documents/MATLAB/PROJECTS/',Paths.projectName,'/Analysis/EEG_ERP.mat'])  
    load(['/Users/',Paths.computer_username,'/Documents/MATLAB/PROJECTS/',Paths.projectName,'/Analysis/EEG_parameters.mat'])  
    Paths.Figures = fullfile(Paths.ProjectPath,'Figures');

    
    
    
%%




% Choose the type of plot (Spectral or temporal analysis)
Data2plot=questdlg('Which data would you like to plot?','Plotting','Spectral analysis','Temporal analysis','Spectral and Temporal analysis','Spectral and Temporal analysis');

% Select the subject to plot (only for the spectral analysis)
if or( strcmp(Data2plot,'Spectral analysis')==1 , strcmp(Data2plot,'Spectral and Temporal analysis')==1)

    subList = subjects;
    subList{end} = 'All subjects'; % Rename the mean of all participants (just so that it looks better)
    [subIdx] = listdlg('PromptString',{'Please select the subject(s) you',...
                                      'would like to plot'},'ListString',subList,'SelectionMode','single');
else
    subIdx=length(subjects);
end
                              
% Select the channel(s) to plot
chanIdx = 68;



    fontsize = 16;
    linewidth = 1.5;
    


    
%% Extraction of the appropriate zscores from the EEG_zscores structure

% Extract the mean zscore at the meter-(un)related frequencies (for the specified subject at the specified channel)
for iSubjects=1:length(subjects)
    

        % Meter-related frequencies
        zscores(1,iSubjects)=EEG_zscores.(subjects{iSubjects}).mean_zscores_metRel(chanIdx);
        % Meter-unrelated frequencies
        zscores(2,iSubjects)=EEG_zscores.(subjects{iSubjects}).mean_zscores_metUnrel(chanIdx); 

        % Meter-related frequencies (without 12th frequency)
        zscoresNo12(1,iSubjects)=EEG_zscores.(subjects{iSubjects}).mean_zscoresNo12_metRel(chanIdx);
        % Meter-unrelated frequencies (without 12th frequency)
        zscoresNo12(2,iSubjects)=EEG_zscores.(subjects{iSubjects}).mean_zscoresNo12_metUnrel(chanIdx); 
        
end


%% fft

if or( strcmp(Data2plot,'Spectral analysis')==1 , strcmp(Data2plot,'Spectral and Temporal analysis')==1)
    
    % Create the figure
    figFFT=figure('position',[1 1 1200 650],'Color',[1 1 1]);
    
    % Main title
    annotation('textbox', [0.1,0.9,0.8,0.08],'String',['EEG - ', subList{subIdx}],'fontsize',fontsize+4,'EdgeColor','none','FontWeight','bold','HorizontalAlignment','center'); 
   
    axes('Position',[0.1 0.1 0.35 0.8])
    
    
        % With the 12th frequency
        plot ([1.1,1.9],squeeze(zscores(:,:)),'LineWidth',linewidth,'Marker','o'); hold on
        
        legend(subjects(1:end),'Position',[0.813846108967362,0.796875,0.136923076923077,0.116875])
    
        % Without the 12th frequency
        for iSubjects=1:subIdx-1
           scatter ([0.95,2.05],zscoresNo12(:,iSubjects),'LineWidth',linewidth,'Marker','o','SizeData',50,'jitter','on','jitterAmount',0.04)
        end
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Plot layout %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    box off
    set(gca,'LineWidth', linewidth,'fontsize',fontsize,'xlim',[0.75 2.25],'xcolor','none','ylim',[-1.2 1.2],'ytick',[-1 1])
    line(xlim(), [0,0], 'LineWidth', linewidth, 'Color', 'k');
    
    
    axes('Position',[0.5 0.1 0.35 0.8])
    
    
        data2plot=EEG_fft.(subjects{subIdx}).blfft.data(chanIdx,:);

        stem(freq,data2plot,'marker','none','LineWidth',linewidth,'Color','k')
        hold on 
        stem(frex(whichMeterRel),data2plot(frexidx(whichMeterRel)),'r','marker','none','LineWidth',1.75)
        stem(frex(whichMeterUnrel),data2plot(frexidx(whichMeterUnrel)),'b','marker','none','LineWidth',1.75)
        
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%% Plot layout %%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        box off
        set (gca,'Tickdir', 'out','fontsize',fontsize,'LineWidth',linewidth,'xlim',[0 7],'ylim',[-0.1 0.7]);


            xlabel ({'Frequency (in Hz)'},'fontsize',fontsize-1,'Position',[3.5,-0.240])

            ylabel ({'µV'},'fontsize',fontsize,'Position',[-1 0.3 -1])

    
end


                set(gcf,'PaperPositionMode','auto')
                print('test.jpg','-djpeg','-r800')