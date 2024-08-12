%3_16_2018 developed by B. Bahmani

function [  ] = plot_pdfXfield_diffSVE_constContRatio(...
    DataAll,ioFrmt,contrastRatio )
%plot collected data for an specific contrast ratio (collect over RVE types regard to each SVE size)
%plot PDF format
%field:tens/shear/betta

numData=length(DataAll);
SVElength_all=zeros(numData,1);


%sellect different cases for SVE size
for i=1:numData
    SVElength_all(i)=DataAll{i}.SVElength;
end
SVElength_all=unique(SVElength_all);
% SVElength_all=sort(SVElength_all,'descend');
corrRVE=cell(1,length(SVElength_all));%corresponding RVE for each of SVElength


for iSVE=1:length(SVElength_all)%loop over each distinct SVE size
    targ_tens=[];
    targ_shear=[];
    targ_betta=[];
    
    for iData=1:numData%loop over all the data to select those with the considering same SVE size
        
        if abs(DataAll{iData}.SVElength-SVElength_all(iSVE))<1e-14 ...
                && abs(DataAll{iData}.contrastRatio-contrastRatio)<1e-14
            
            corrRVE{iSVE}=[corrRVE{iSVE};iData];%save the RVE id for the specific SVE size
            if DataAll{iData}.fieldActivity(1)==1
                targ_tens=[targ_tens;DataAll{iData}.TensileStrength_sveXstat];%collecting all relevent data
            end
            if DataAll{iData}.fieldActivity(2)==1
                targ_shear=[targ_shear;DataAll{iData}.ShearStrength_sveXstat];%collecting all relevent data
            end
            if DataAll{iData}.fieldActivity(3)==1
                targ_betta=[targ_betta;DataAll{iData}.BettaStrength_sveXstat];%collecting all relevent data
            end
        end
        
        
    end
    
    %now we have all the required data for this specific SVE size
    
    numStatVar=length(ioFrmt.fieldType_stat);%num of statistic field (min/max/mean/std)
    
    for istat=1:numStatVar%loop over each stat field
        if ~isempty(targ_tens)
            [pdf_val, x_val] = ksdensity(targ_tens(:,1+istat), 'function', 'PDF');
            idataFiled=1;
            figNum_tens=(idataFiled-1)*numStatVar+istat;
            figure(figNum_tens);
            hold on;
            plot(x_val, pdf_val, 'linewidth', 1);
            xlabelText = ['$$ ', ioFrmt.StatNotation_plot{istat},...
                ioFrmt.fieldType_plot{1}, ') $$'];
            xlabel(xlabelText, 'Interpreter', 'latex', 'FontSize', ioFrmt.axisfs, 'VerticalAlignment','Top');
            ylabel('PDF');
            %%'SVE Size 25 (M1-M4)'
            %fprintf('RVE=%d\t data=%d\t stat=%d\n',idRve,dt_r,istat);
        end
        if ~isempty(targ_shear)
            [pdf_val, x_val] = ksdensity(targ_shear(:,1+istat), 'function', 'PDF');
            idataFiled=2;
            figNum_tens=(idataFiled-1)*numStatVar+istat;
            figure(figNum_tens);
            hold on;
            plot(x_val, pdf_val, 'linewidth', 1);
            xlabelText = ['$$ ', ioFrmt.StatNotation_plot{istat},...
                ioFrmt.fieldType_plot{2}, ') $$'];
            xlabel(xlabelText, 'Interpreter', 'latex', 'FontSize', ioFrmt.axisfs, 'VerticalAlignment','Top');
            ylabel('PDF');
            %%'SVE Size 25 (M1-M4)'
            %fprintf('RVE=%d\t data=%d\t stat=%d\n',idRve,dt_r,istat);
        end
        
        if ~isempty(targ_betta)
            [pdf_val, x_val] = ksdensity(targ_shear(:,1+istat), 'function', 'PDF');
            idataFiled=3;
            figNum_tens=(idataFiled-1)*numStatVar+istat;
            figure(figNum_tens);
            hold on;
            plot(x_val, pdf_val, 'linewidth', 1);
            xlabelText = ['$$ ', ioFrmt.StatNotation_plot{istat},...
                ioFrmt.fieldType_plot{3}, ') $$'];
            xlabel(xlabelText, 'Interpreter', 'latex', 'FontSize', ioFrmt.axisfs, 'VerticalAlignment','Top');
            ylabel('PDF');
            %%'SVE Size 25 (M1-M4)'
            %fprintf('RVE=%d\t data=%d\t stat=%d\n',idRve,dt_r,istat);
        end
        
    end
    
end


%isert legend to plots

legendName=cell(length(SVElength_all),1);
for iSVElist=1:length(SVElength_all)
    collRVE=[];
    for i=1:length(corrRVE{iSVElist})
        collRVE=strcat(collRVE,num2str(DataAll{corrRVE{iSVElist}(i)}.RVEnum),',');
    end
    legendName{iSVElist}=strcat('SVE size ',num2str(SVElength_all(iSVElist)),' microStruct:',collRVE);
end

for istat=1:numStatVar
    if ~isempty(targ_tens)
        idataFiled=1;
        figNum_tens=(idataFiled-1)*numStatVar+istat;
        figure(figNum_tens);
        legend(legendName,'Location','bestoutside');
        %         title(TestTypeT{itype});
        id=strcat(ioFrmt.OutFodler,'/RVEMixed_',ioFrmt.fieldType_stat{istat},...
            ioFrmt.fieldType_write_pdfXstat{1});
        print('-dpng',id);
    end
    
    if ~isempty(targ_shear)
        idataFiled=2;
        figNum_tens=(idataFiled-1)*numStatVar+istat;
        figure(figNum_tens);
        legend(legendName,'Location','bestoutside');
        %         title(TestTypeT{itype});
        id=strcat(ioFrmt.OutFodler,'/RVEMixed_',ioFrmt.fieldType_stat{istat},...
            ioFrmt.fieldType_write_pdfXstat{2});
        print('-dpng',id);
    end
    
    if ~isempty(targ_betta)
        idataFiled=2;
        figNum_tens=(idataFiled-1)*numStatVar+istat;
        figure(figNum_tens);
        legend(legendName,'Location','bestoutside');
        %         title(TestTypeT{itype});
        id=strcat(ioFrmt.OutFodler,'/RVEMixed_',ioFrmt.fieldType_stat{istat},...
            ioFrmt.fieldType_write_pdfXstat{3});
        print('-dpng',id);
    end
    
end


end

