function fftPlot(iCond,subIdx,ElecIdx)

%%%%%%%%%%%%%%%%%%%%%%%
% Plotting parameters %
%%%%%%%%%%%%%%%%%%%%%%%

 % Define global variables 
   global Cfg EEG_fft EEG_zscores EEG_ERP Paths 
   
   
 % Define plot positions depending on the number of conditions
   if size(Cfg.condNames,2) == 1      
       AxesPos = [0.6 0.1 0.35 0.8];
              
   elseif size (Cfg.condNames,2) == 4      
        axesPos{1}  = [0.32 0.62 0.18 0.24];
        axesPos{2}  = [0.32 0.25 0.18 0.24];
        axesPos{3}  = [0.55 0.62 0.18 0.24];
        axesPos{4}  = [0.55 0.25 0.18 0.24];       
   end   
   

%%%%%%%%%%%%%%%%%%%% 
% Get Data to Plot % 
%%%%%%%%%%%%%%%%%%%%   
   
   if size(subIdx,2) == 1 % One participant
        data2plot = EEG_fft.(Cfg.condNames{iCond}).(Cfg.subjects{subIdx}).blfft.data(ElecIdx,:);
   
   else % Multiple participants
       
        % Extract the data2plot from the structure 
        for iSub = 1:length(subIdx)
            data2plot_temp (iSub,:) = EEG_fft.(Cfg.condNames{iCond}).(Cfg.subjects{subIdx(iSub)}).blfft.data(ElecIdx,:);  
        end
        
        % Mean fft of the selected participants
        data2plot = mean(data2plot_temp,1);
   end
   
%%%%%%%%%%%%%%%%
% Plot the fft %
%%%%%%%%%%%%%%%%

    axes('Position',axesPos{iCond});
    
    stem(Cfg.freq,data2plot,'marker','none','LineWidth',Cfg.figure.linewidth-0.35,'Color',[0.15 0.15 0.15])
    hold on 
    
    % Meter-related frequencies in red
    stem(Cfg.frex(Cfg.whichMeterRel),data2plot(Cfg.frexidx(Cfg.whichMeterRel)),...
            'Color',[238,44,44]/255,'marker','none','LineWidth',Cfg.figure.linewidth+0.6)
    
    % Meter-unrelated frequencies in blue    
    stem(Cfg.frex(Cfg.whichMeterUnrel),data2plot(Cfg.frexidx(Cfg.whichMeterUnrel)),...
            'Color',[24,116,205]/255,'marker','none','LineWidth',Cfg.figure.linewidth+0.6)
        
        
    % Tomas
    % Red [214, 52, 24]/255
    % Blue [45, 114, 224]/255

    % Ceren
    % firebrick2 238,44,44
    % Dodgerblue3 24,116,205 

%%%%%%%%%%%%%%%
% Plot layout %
%%%%%%%%%%%%%%%
    
    box off
    set (gca,'Tickdir', 'out','fontsize',Cfg.figure.fontsize,...
         'LineWidth',Cfg.figure.linewidth,'xlim',[0 7],'ylim',[-0.1 0.7]);
         
    % Xlabels
    if isequal(iCond,2) || isequal(iCond,4)
        xlabel ({'Frequency (in Hz)'},'fontsize',Cfg.figure.fontsize-1,'Position',[3.5,-0.240])
    end

    % Ylabels
    if isequal(iCond,1) || isequal(iCond,2)
        ylabel ({'µV'},'fontsize',Cfg.figure.fontsize,'Position',[-1 0.3 -1])
    end
   
   
   
end

