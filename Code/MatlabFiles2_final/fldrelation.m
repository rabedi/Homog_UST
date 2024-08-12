function [  ] = fldrelation(fldId,normFactor,field_X,field_Y)
global  DataAll
num_case=length(DataAll);
for icase=1:num_case
    [numSve,E,Sdiv,COV,COR]=comput_COV_COR(fldId,icase);
    printData(numSve,E,Sdiv,COV,COR,icase,fldId);
    plotField(field_X,field_Y,icase);
end

end

function [numSve,E,Sdiv,COV,COR]=comput_COV_COR(fldId,icase)
global  DataAll obj_Fields_Name
num_fld=length(fldId);
E=zeros(num_fld,1);
Sdiv=zeros(num_fld,1);
P=zeros(num_fld);
COV=zeros(num_fld);
COR=zeros(num_fld);
for ifld=1:num_fld
    numSve=length((DataAll{icase}.data_sveXfield(:,fldId(ifld))));
    E(ifld)=mean(DataAll{icase}.data_sveXfield(:,fldId(ifld)));%/obj_Fields_Name.normalization{field_normalization(ifld)}{1};
    Sdiv(ifld)=std(DataAll{icase}.data_sveXfield(:,fldId(ifld)));%/obj_Fields_Name.normalization{field_normalization(ifld)}{1};
    for jfld=1:num_fld
        temp=sum(DataAll{icase}.data_sveXfield(:,fldId(ifld)).*...
            DataAll{icase}.data_sveXfield(:,fldId(jfld)));
        P(ifld,jfld)=temp;
        COV(ifld,jfld)=temp/numSve+E(ifld)*E(jfld);
    end
end
for ifld=1:num_fld
    for jfld=1:num_fld
        COR(ifld,jfld)=COV(ifld,jfld)/sqrt(COV(ifld,ifld))/...
            sqrt(COV(jfld,jfld));
    end
end
end

function printData(numSve,E,Sdiv,COV,COR,icase,fldId)
global DataAll
name=strcat('../OutputPlot/COR_COV_RVE',num2str(DataAll{icase}.RVE_lx),'X',...
    num2str(DataAll{icase}.RVE_ly),...
    'SVE',num2str(DataAll{icase}.SVE_lx),'X',num2str(DataAll{icase}.SVE_ly),'BC_',...
    DataAll{icase}.BC_type,'.txt');
file=fopen(name,'w');
fprintf(file,'#number of Data Point\n');
fprintf(file,'%d\n',numSve);
fprintf(file,'#number of fields\n');
fprintf(file,'%d\n',length(E));
fprintf(file,'#field names\n');
for i=1:length(fldId)
    temp=DataAll{icase}.fieldsID_to_plot(fldId(i));
    temp=temp{1};
    temp=temp{1};
    temp=temp{1};
    fprintf(file,'%s\t',temp);
end
fprintf(file,'\n');

fprintf(file,'#Mean vector\n');
for i=1:length(E)
    fprintf(file,'%f\n',E(i));
end
fprintf(file,'#Standard Div vector\n');
for i=1:length(Sdiv)
    fprintf(file,'%f\n',Sdiv(i));
end
fprintf(file,'#Cov mat\n');
for i=1:size(COV,1)
    for j=1:size(COV,2)
        fprintf(file,'%f\t',COV(i,j));
    end
    fprintf(file,'\n');
end
fprintf(file,'#Cor mat\n');
for i=1:size(COR,1)
    for j=1:size(COR,2)
        fprintf(file,'%f\t',COR(i,j));
    end
    fprintf(file,'\n');
end

end

function plotField(field_X,field_Y,icase)
global DataAll obj_Fields_Name
pppp = plt_plot_plotProperties;
for iy=1:length(field_Y)
    for ix=1:length(field_X)
        ppdp = getLineProperties_fldrelation();
        xDat=DataAll{icase}.data_sveXfield(:,field_X(ix));
        yDat=DataAll{icase}.data_sveXfield(:,field_Y(iy));
        plt_plotData_plotXYbasedOnDataSpec(ppdp, xDat, yDat);
        
        
        
        
        
        % as mentioned above, we want to use latex interpreter
isLatex = 1;
% this is for eps figures which we don't want -> set it to 0
isPSFrag = 0;
        idfld_glob=DataAll{icase}.fieldsID_to_plot{field_X(ix)}(2);
        idfld_glob=idfld_glob{1};
        tempFldName=obj_Fields_Name.fieldLib{idfld_glob}(3);
        xlab = ['$$',tempFldName{1},'$$'];
        idfld_glob=DataAll{icase}.fieldsID_to_plot{field_Y(iy)}(2);
        idfld_glob=idfld_glob{1};
        tempFldName=obj_Fields_Name.fieldLib{idfld_glob}(3);
        ylab = ['$$',tempFldName{1},'$$'];
        
        % for axis limits nan is NOT setting the limit for min and/or max values ->
        % using a real number DOES set the limit
        xMinIn = nan;
        xMaxIn = nan;
        
        % if you want to set x lim min or max values uncomment
        %        ONE
        %or
        %       BOTH
        %lines below
        
        %xMinIn = 0.2;
        xMaxIn = nan;
        
        yMinIn = nan;
        yMaxIn = nan;
        
        pppp.setAxesLimitsLabelsFontSize(gca, xlab, ylab, isPSFrag, xMinIn, xMaxIn, yMinIn, yMaxIn, isLatex)';
        % legend
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % legend location
        %pppp.legend.location = 'bestoutside';
% %         titleLgd='';
% %         leg='';
% %         pppp.setLegend(gca, leg, isLatex,titleLgd);
        
        
        % printing ...
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fldx=DataAll{icase}.fieldsID_to_plot{field_X(ix)}(1);
        fldx=fldx{1}{1};
        fldy=DataAll{icase}.fieldsID_to_plot{field_Y(iy)}(1);
        fldy=fldy{1}{1};
        outPlotName=strcat(fldx,'_',fldy,'_RVE',num2str(DataAll{icase}.RVE_lx),...
            'X',num2str(DataAll{icase}.RVE_ly),'SVE',num2str(DataAll{icase}.SVE_lx),...
            'X',num2str(DataAll{icase}.SVE_ly),'_BC_',DataAll{icase}.BC_type);
        locTemp=strcat('../OutputPlot/','FldRel_',outPlotName,'.png');
        print('-dpng', locTemp);
        locTemp=strcat('../OutputPlot/','FldRel_',outPlotName,'.fig');
        savefig(gcf, locTemp);
        
    end
end

end

function ppdp = getLineProperties_fldrelation()
global lineDataBase;

ppdp = plt_plotDataProp;
ppdp.val_lineStyle='none';
ppdp.val_marker = '.';%lineDataBase.markerStyleAllTb(1,1);
% % ppdp.val_lineStyle = lineDataBase.lineStyleTb{ibc};
% % ppdp.val_lineColor = lineDataBase.colorNameClrTb{isz, 2};
% % if nargin>2
% %     ppdp.val_marker = lineDataBase.markerStyleAllTb{ifld,2};
% % end

end
