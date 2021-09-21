function zscoresPlot(iCond,subIdx,ElecIdx)

%%%%%%%%%%%%%%%%%%%%%%%
% Plotting parameters %
%%%%%%%%%%%%%%%%%%%%%%%

 % Define global variables 
   global Cfg EEG_fft EEG_zscores EEG_ERP Paths 

   
 % Define plot positions depending on the number of conditions
   if size(Cfg.condNames,2) == 1   
       disp('1 condition')
       axesPos = [0.1 0.1 0.35 0.8];
              
   elseif size (Cfg.condNames,2) == 4      
        axesPos{1}  = [0.08 0.62 0.18 0.26];
        axesPos{2}  = [0.08 0.25 0.18 0.26];
        axesPos{3}  = [0.8 0.62 0.18 0.26];
        axesPos{4}  = [0.8 0.25 0.18 0.26];       
   end
   
 % Define colors
   Colors_Palette1{1}    = [0.9765    0.2549    0.2667];
   Colors_Palette1{2}    = [0.1529    0.6078    0.6314];
   Colors_Palette1{3}    = [0.2627    0.6667    0.5451];
   Colors_Palette1{4}    = [0.9529    0.4471    0.1725];
   Colors_Palette1{5}    = [0.9765    0.7804    0.3098];
   Colors_Palette1{6}    = [0.3412    0.4588    0.5647];
   Colors_Palette1{7}    = [0.9725    0.5882    0.1176];
   Colors_Palette1{8}    = [0.5647    0.7451    0.4275];
   Colors_Palette1{9}    = [0.9765    0.5176    0.2902];
   Colors_Palette1{10}   = [0.3020    0.5647    0.5569];
   
   Colors_Palette2{1}    = [239, 71, 111]./255;
   Colors_Palette2{2}    = [255, 209, 102]./255;
   Colors_Palette2{3}    = [6, 214, 160]./255;
   Colors_Palette2{4}    = [7, 59, 76]./255;
   Colors_Palette2{5}    = [46, 31, 39]./255;
   Colors_Palette2{6}    = [94, 76, 90]./255;
   Colors_Palette2{7}    = [133, 77, 39]./255;
   Colors_Palette2{8}    = [2, 195, 154]./255;
   Colors_Palette2{9}    = [16, 73, 17]./255;
   Colors_Palette2{10}   = [16, 43, 63]./255;
   
%%%%%%%%%%%%%%%%%%%% 
% Get Data to Plot % 
%%%%%%%%%%%%%%%%%%%%

  % Meter-(un)related zscores with(out) the frequency of the unitary event
    for iSub = 1:length(subIdx)
        
        zscores(:,iSub)     = [EEG_zscores.(Cfg.condNames{iCond}).(Cfg.subjects{subIdx(iSub)}).meanZscoresMetRel(ElecIdx), ...
                               EEG_zscores.(Cfg.condNames{iCond}).(Cfg.subjects{subIdx(iSub)}).meanZscoresMetUnrel(ElecIdx)];
        zscoresNo12(:,iSub) = [EEG_zscores.(Cfg.condNames{iCond}).(Cfg.subjects{subIdx(iSub)}).meanZscoresNo12MetRel(ElecIdx), ...
                               EEG_zscores.(Cfg.condNames{iCond}).(Cfg.subjects{subIdx(iSub)}).meanZscoresNo12MetUnrel(ElecIdx)];    
    end

%%%%%%%%%%%%%%%%%%%%
% Plot the zscores %
%%%%%%%%%%%%%%%%%%%%

    if size(Cfg.condNames,2) == 1  
        axes('Position',axesPos);
    else 
        axes('Position',axesPos{iCond});
    end
    
%     plot([1:2], zscores, 'LineWidth',Cfg.figure.linewidth,'Marker','o','Color',Colors{1});
%     hold on
    
    if length(subIdx) == 1
%         plot([1:2], zscores, 'LineWidth',Cfg.figure.linewidth,'Marker','o','Color',Colors_Palette2{1});
%         hold on
%         plot([1:2],zscoresNo12, 'LineWidth',Cfg.figure.linewidth,'Marker','o','Color',Colors_Palette2{2}) 
        plot([1:2], zscores, 'LineWidth',Cfg.figure.linewidth,'Marker','o','Color',[0.65 0.65 0.65]);
        hold on
        plot([1:2],zscoresNo12, 'LineWidth',Cfg.figure.linewidth,'Marker','o','Color',[0.6 0.6 0.6]) 
        
    else 
        for iSub = 1:length(subIdx)
            plot([1:2], zscores(:,iSub), 'LineWidth',Cfg.figure.linewidth,'Marker','o','Color',[0.65 0.65 0.65]);
            hold on
            scatter ([0.9,2.1],zscoresNo12(:,iSub),...
                'LineWidth',Cfg.figure.linewidth,'Marker','o','SizeData',50,'jitter','on','jitterAmount',0.04,'MarkerEdgeColor',[0.6 0.6 0.6])     
        end  
    end
    
      
%%%%%%%%%%%%%%%
% Plot layout %
%%%%%%%%%%%%%%%

    box off
    set(gca,'LineWidth', Cfg.figure.linewidth,'fontsize',Cfg.figure.fontsize,'xlim',[0.75 2.25],'xcolor','none','ylim',[-1.2 1.2],'ytick',[-1 1])
    line(xlim(), [0,0], 'LineWidth', Cfg.figure.linewidth, 'Color', 'k');
    title(Cfg.condNames{iCond},'fontsize',Cfg.figure.fontsize)
    
  % Xlabel 
    if isequal(iCond,1) || isequal(iCond,2) || size(Cfg.condNames,2) == 1 
        ylabel ({'Amplitude';'(Zscore)'},'fontsize',Cfg.figure.fontsize,'Position',[0.73 0 -1])
    end

  % Ylabel
    if size (Cfg.condNames,2) == 4 && isequal(iCond,2)
        annotation('textbox', [0.0783,0.1646,0.0894,0.087],'String',{'Meter-related';'frequencies'},'EdgeColor','none','fontsize',Cfg.figure.fontsize-2,'Color',[0.9333    0.1725    0.1725])
        annotation('textbox', [0.1783,0.1646,0.0894,0.087],'String',{'Meter-unrelated';'frequencies'},'EdgeColor','none','fontsize',Cfg.figure.fontsize-2,'Color',[0.0941    0.4549    0.8039])   
    
    elseif size (Cfg.condNames,2) == 4 && isequal(iCond,4)
        annotation('textbox', [0.7975,0.1646,0.0894,0.087],'String',{'Meter-related';'frequencies'},'EdgeColor','none','fontsize',Cfg.figure.fontsize-2,'Color',[0.9333    0.1725    0.1725])
        annotation('textbox', [0.8991,0.1646,0.0894,0.087],'String',{'Meter-unrelated';'frequencies'},'EdgeColor','none','fontsize',Cfg.figure.fontsize-2,'Color',[0.0941    0.4549    0.8039])
    
    elseif size (Cfg.condNames,2) == 1
        % Redefine the correct position when this scenario first occurs
        annotation('textbox', [0.7975,0.1646,0.0894,0.087],'String',{'Meter-related';'frequencies'},'EdgeColor','none','fontsize',Cfg.figure.fontsize-2,'Color',[0.9333    0.1725    0.1725])
        annotation('textbox', [0.8991,0.1646,0.0894,0.087],'String',{'Meter-unrelated';'frequencies'},'EdgeColor','none','fontsize',Cfg.figure.fontsize-2,'Color',[0.0941    0.4549    0.8039])
    end
    
    
  % Legend
    if size (Cfg.subjects(subIdx),1) == 1
        legend ('with freq 12','without freq 12','Position',[0.172,0.775,0.107,0.084],'fontsize',Cfg.figure.fontsize-3)
    else 
        % legend(Cfg.subjects(subIdx),'Position',[0.080,0.002,0.192,0.143])   
    end
  
  
   
end

