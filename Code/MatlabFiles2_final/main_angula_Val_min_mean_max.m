%updted 12_9_2017
clc
clear all
close all
fclose all
opengl software
% close all
%Input data
%set(gca, 'DefaultFontSize', 12);
M=[1];
CR=[100];
% S=[4 8 16 32];%number of part in each dir
S=[8];%number of part in each dir
dteta=5;
statType=2;%min-max(1)   std(2)
SVEtype='S';% V=voronoi or S=square
dataRootLoc='../OUTPUT/Square/';
% dataT2run=[1 2 3];
dataT2run=[2];
OutFodler='../PLOTs';

%plot format
legfs = 18;
axisfs = 24;
titlefs = 22;
LINEWIDTH=2;
%data types
dataType_read{1}='(TStrn_Angle).OTS';
dataType_read{2}='(SStrn_Angle).OSS';
dataType_read{3}='(ModeMixity_Angle).OMM';

dataType_write{1}='(TStrn_Angle)';
dataType_write{2}='(SStrn_Angle)';
dataType_write{3}='(ModeMixity_Angle)';

dataType_plot{1}='\tilde{s}_n/\sigma^{TH}';
dataType_plot{2}='\tilde{s}_t/\sigma^{TH}';
dataType_plot{3}='\beta';

legendName{1}={'$$\textrm{Min}$$';'$$\textrm{Mean}$$';'$$\textrm{Max}$$'};
legendName{2}={'$$\textrm{- S.Div}$$';'$$\textrm{Mean}$$';'$$\textrm{+ S.Div}$$'};
%------------------------------------------------------------------
if(length(dataT2run)>length(dataType_read))
    fprintf('Error in wanted data to plot');
    pauss;
end
if exist(OutFodler,'dir')~=7
    mkdir(OutFodler);
end
for im=1:length(M)
    for icr=1:length(CR)
        for is=1:length(S)
            if strcmp(SVEtype,'S')
                id0=strcat(dataRootLoc,'M',num2str(M(im)),'CR',...
                    num2str(CR(icr)),'S',num2str(S(is)));
            elseif strcmp(SVEtype,'V')
                id0=strcat(dataRootLoc,'M',num2str(M(im)),'CR',...
                    num2str(CR(icr)),'V',num2str(S(is)));
            else
                error('SVE type not implemented\n');
            end
            nTeta=180/dteta+1;
            Nsve=S(is)*S(is);
            for idt_r=1:length(dataT2run)
                val=zeros(nTeta,2,Nsve);
                nonExsistSVE=[];
                for isve=1:Nsve
                    [is isve]
                    dt_r=dataT2run(idt_r);
                    id=strcat(id0,'/SVE',num2str(isve),dataType_read{dt_r});
                    if exist(id,'file')==2
                        file = fopen(id);
                        line=fgets(file);%RVE	 SVE	 dteta	 Nteta
                        line=fgets(file);
                        temp=sscanf(line,'%d %d %f %d');
                        idRve=temp(1);
                        idSve=temp(2);
                        if nTeta~=temp(4) || dteta~=temp(3)
                            error('Error in number of teta\n');
                        end
                        line=fgets(file);%Angle	 dataType
                        
                        for iteta=1:nTeta
                            line=fgets(file);
                            val(iteta,:,isve)=sscanf(line,'%f %f');%i=itete  j=ifld  k=isve
                        end
                        fclose(file);
                    else
                       nonExsistSVE=[nonExsistSVE;isve];
                    end
                end
                val(:,:,nonExsistSVE)=[];
                lowerBand=zeros(nTeta,1);
                upperBand=zeros(nTeta,1);
                meanBand=zeros(nTeta,1);
                for it=1:nTeta
                    temp=val(it,2,:);
                    Min=min(temp);
                    Max=max(temp);
                    Mean=mean(temp);
                    Std=std(temp);
                    if statType==1
                        lowerBand(it)=Min;
                        upperBand(it)=Max;
                        meanBand(it)=Mean;
                    elseif statType==2
                        lowerBand(it)=Mean-Std;
                        upperBand(it)=Mean+Std;
                        meanBand(it)=Mean;
                    end
                end
                fclose all
                close all
                fgm=figure('Visible','on');
                %                 figure('Visible','off');
                plot(val(:,1,1),lowerBand,'Color','k','LineStyle','--', 'Linewidth', LINEWIDTH);
                hold on
                plot(val(:,1,1),meanBand,'Color','k','LineStyle','-', 'Linewidth', LINEWIDTH);
                plot(val(:,1,1),upperBand,'Color','k','LineStyle','--', 'Linewidth', LINEWIDTH);
                inBetween = [lowerBand', fliplr(upperBand')];
                x2 = [val(:,1,1)', fliplr(val(:,1,1)')];
                fillhandle=fill(x2, inBetween,'k');
                transparency=.2;
                set(fillhandle,'EdgeColor','k','FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color
                
                axis([0 inf .1 .55])
                set(gca,'XTick',[0:30:180]);
                set(gca,'YTick',[.1:.15:.55]);
                xAx = get(gca, 'XAxis');%to change axis tick size
                set(xAx, 'FontSize', 12);
                yAx = get(gca, 'YAxis');%to change axis tick size
                set(yAx, 'FontSize', 12);
                set(gca,'TickLabelInterprete','latex')
                
                xlabel('$$ \theta $$', 'Interpreter', 'latex', 'FontSize', axisfs, 'VerticalAlignment','Top');
                ylabel(['$$ ',dataType_plot{dt_r}, ' $$'], 'Interpreter', 'latex', 'FontSize', axisfs, 'VerticalAlignment','Bottom');
%                                 ylabel('$$ \textrm{Value} $$', 'Interpreter', 'latex', 'FontSize', axisfs, 'VerticalAlignment','Bottom');
                lgd=legend(legendName{statType},'Interpreter', 'latex','Location','best');
%                 title(lgd,['$$',dataType_plot{dt_r},'$$'],'Interpreter', 'latex');
                legend('boxoff');
                lgd.FontSize = legfs;
                
                
                if strcmp(SVEtype,'S')
                    id=strcat(OutFodler,'/M',num2str(M(im)),'CR',...
                        num2str(CR(icr)),'S',num2str(S(is)),'STATtyp',...
                        num2str(statType),'_AllSVE',dataType_write{dt_r},'.png');
                    id2=strcat(OutFodler,'/M',num2str(M(im)),'CR',...
                        num2str(CR(icr)),'S',num2str(S(is)),'STATtyp',...
                        num2str(statType),'_AllSVE',dataType_write{dt_r},'.fig');
                elseif strcmp(SVEtype,'V')
                    id=strcat(OutFodler,'/M',num2str(M(im)),'CR',...
                        num2str(CR(icr)),'V',num2str(S(is)),'STATtyp',...
                        num2str(statType),'_AllSVE',dataType_write{dt_r},'.png');
                    id2=strcat(OutFodler,'/M',num2str(M(im)),'CR',...
                        num2str(CR(icr)),'V',num2str(S(is)),'STATtyp',...
                        num2str(statType),'_AllSVE',dataType_write{dt_r},'.fig');
                else
                    error('SVE type is not implemented!\n');
                end
                
                print('-dpng',id);
                savefig(gcf, id2);
                %                 fprintf('RVE=%d\t SVE=%d\t data=%d\n',idRve,idSve,dt_r);
                %                 hold off
            end
            
        end
    end
end
