function [  ] = plotGrpoup_Val_teta( fieldName_teta,field_normalization, bc_type,...
    sve_lx,sve_ly,sveId,phaseId )
%plot a group of data on a same figure
%fieldName_teta:special names of the fields
%field_normalization:shows a factor correspond to each field id to
%normalized
%sveId:is a vector of SVE ids to plot simultanuasly
%phaseId:is a vector
%phaseId: is a vector 0,1,2,3 should start from 0!!!!
%%%fieldName_teta={'LT0LF';'LT1LF';'betta'};
global obj_Fields_Name DataAll
sveId=sort(sveId);
figure;
cntr=1;
%pppp stores legend, xaxes, yaxes default properties, see the uses below
pppp = plt_plot_plotProperties;
[diffCases,legendTag,outPlotName]=getDiffCasesIdentifier(fieldName_teta,...
    bc_type,sve_lx,sve_ly,sveId,phaseId );

for ibc=1:length(bc_type)
    for isz=1:length(sve_lx)
        dataId=getDataId(DataAll,bc_type{ibc},sve_lx(isz),sve_ly(isz));
        if dataId==-1
            continue;
        end
        baseFolder=DataAll{dataId}.RVE_folder_name;
        dataTeta=readData_teta(baseFolder,sveId,phaseId,fieldName_teta);
        
        for icase=1:length(dataTeta)
            if strcmp(dataTeta(icase).fieldName,'LT0LF') || ...
                    strcmp(dataTeta(icase).fieldName,'LT1LF')
                y=dataTeta(icase).fldVal/obj_Fields_Name.normalization{field_normalization(1)}{1};
            elseif strcmp(dataTeta(icase).fieldName,'betta')%betta is already normalized!(do not have any unit)
                y=dataTeta(icase).fldVal;
            else
                error('not implemented case!');
            end
            ppdp = getLineProperties_PDFbased(icase);%SAMPLE_getLineProperties
            %         if length(field_id)>1
            %             ppdp = getLineProperties_PDFbased(ibc,isz,ifld);%SAMPLE_getLineProperties
            %         else
            %             ppdp = getLineProperties_PDFbased(ibc,isz);%SAMPLE_getLineProperties
            %         end
            plt_plotData_plotXYbasedOnDataSpec(ppdp, dataTeta(icase).tetaArr, y);
            hold on;
            % %         idfld_glob=DataAll{dataId}.fieldsID_to_plot{field_id(ifld)}(2);
            % %         idfld_glob=idfld_glob{1};
            % %         tempFldName=obj_Fields_Name.fieldLib(idfld_glob);
            % %         tempFldName=tempFldName{1};
            % %         tempFldName=tempFldName(3);
            % %         leg{cntr}=getLegendTitle(dataTeta(icase).fieldName, bc_type{ibc},...
            % %             sve_lx(isz),sve_ly(isz),diffCases);
            leg{cntr}=getLegendTitle(dataTeta(icase).fieldName, bc_type{ibc},...
                sve_lx(isz),sve_ly(isz),dataTeta(icase).sveId,...
                dataTeta(icase).phase,diffCases);
            
            cntr = cntr + 1;
        end
    end
end

% as mentioned above, we want to use latex interpreter
isLatex = 1;
% this is for eps figures which we don't want -> set it to 0
isPSFrag = 0;

% setting x, y labels and limits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[appendAxisDataLabel2TitleX, appendAxisDataLabel2TitleY] =
if length(fieldName_teta)==1
% %     idfld_glob=DataAll{1}.fieldsID_to_plot{field_id}(2);
% %     idfld_glob=idfld_glob{1};
% %     tempFldName=obj_Fields_Name.fieldLib{idfld_glob}(3);
    if strcmp(fieldName_teta,'LT0LF')
        temp='(\tilde{s}_n)';
    elseif strcmp(fieldName_teta,'LT1LF')
        temp='(\tilde{s}_t)';
    else
        error('imposible case!\n');
    end
    ylab = ['$$',temp,' ',obj_Fields_Name.normalization{field_normalization(1)}{2},'$$'];
    xlab = '$$ \theta $^{\circ}$ $$';
else
    ylab = ['$$','\textrm{Value}',' ',obj_Fields_Name.normalization{field_normalization(1)}{2},'$$'];
    xlab = '$$ \theta ^{\circ} $$';
end


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
titleLgd=legendTag;
pppp.setLegend(gca, leg, isLatex,titleLgd);


% printing ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
locTemp=strcat('../OutputPlot/','Teta_',outPlotName,'.png');
print('-dpng', locTemp);
locTemp=strcat('../OutputPlot/','Teta_',outPlotName,'.fig');
savefig(gcf, locTemp);




hold off;
close all;
end

function dataId=getDataId(DataAll,bc_type,sve_lx,sve_ly)

dataId=-1;
tag=false;
i=1;
while tag==false && i<=length(DataAll)
    if strcmp(DataAll{i}.BC_type,bc_type) && abs(DataAll{i}.SVE_lx-sve_lx)<1e-5...
            && abs(DataAll{i}.SVE_ly-sve_ly)<1e-5
        dataId=i;
        tag=true;
    end
    i=i+1;
end

end

function ppdp = getLineProperties_PDFbased(icase)
global lineDataBase;

ppdp = plt_plotDataProp;
ppdp.val_lineStyle = lineDataBase.lineStyleTb{icase};
% % ppdp.val_lineColor = lineDataBase.colorNameClrTb{isz, 2};
% % if nargin>2
% %     ppdp.val_marker = lineDataBase.markerStyleAllTb{ifld,2};
% % end

end

function titleLgd=getLegendTitle(tempFldName, bc_type,sve_lx,sve_ly,sveId,...
    phaseId,diffCases )

%diffCases: is an array with 0 and 1 two show which of the legend variables
%are uniqe(0) or have different cases(1)

if diffCases(1)==1 && diffCases(2)==0 && diffCases(3)==0 && diffCases(4)==0 && diffCases(5)==0 && diffCases(6)==0
    if strcmp(tempFldName,'LT0LF')
        temp='(\tilde{s}_n)';
    elseif strcmp(tempFldName,'LT1LF')
        temp='(\tilde{s}_t)';
    else
        error('imposible case!\n');
    end
    titleLgd = ['$$',temp,'$$'];
    
elseif diffCases(1)==0 && diffCases(2)==1 && diffCases(3)==0 && diffCases(4)==0 && diffCases(5)==0 && diffCases(6)==0
    
    titleLgd = ['$$',bc_type,'$$'];
    
elseif diffCases(1)==0 && diffCases(2)==0 && diffCases(3)==1 && diffCases(4)==1 && diffCases(5)==0 && diffCases(6)==0
    
    titleLgd = ['$$',num2str(sve_lx),'\times',num2str(sve_ly),'$$'];
    
elseif diffCases(1)==0 && diffCases(2)==0 && diffCases(3)==0 && diffCases(4)==0 && diffCases(5)==1 && diffCases(6)==0
    
    titleLgd = ['$$',num2str(sveId),'$$'];
    
elseif diffCases(1)==0 && diffCases(2)==0 && diffCases(3)==0 && diffCases(4)==0 && diffCases(5)==0 && diffCases(6)==1
    
    titleLgd = ['$$',num2str(phaseId),'$$'];
    
else
    titleLgd = ['$$\textrm{Field:}',tempFldName,...
        ',','\textrm{BC:}',bc_type,',','\textrm{SVE Size:}',...
        num2str(sve_lx),'\times',num2str(sve_ly),',',num2str(sveId),...
        ',',num2str(phaseId),' $$'];
end


end

function [diffCases,legendTag,outPlotName]=getDiffCasesIdentifier(fieldName,...
    bc_type,sve_lx,sve_ly, sveId,phaseId)
%global DataAll obj_Fields_Name
diffCases=zeros(6,1);
legendTag='';
count=0;
tag=0;
outPlotName='';
dataId=1;
if length(fieldName)~=1
    diffCases(1)=1;
    count=count+1;
    tag='$$\textrm{Field}$$';
end
for ifld=1:length(fieldName)
    % %     idfld_glob=DataAll{dataId}.fieldsID_to_plot{fieldName(ifld)}(2);
    % %     idfld_glob=idfld_glob{1};
    % %     tempFldName=obj_Fields_Name.fieldLib(idfld_glob);
    % %     tempFldName=tempFldName{1};
    % %     tempFldName=tempFldName(2);
    if ifld==1
        if strcmp(fieldName(ifld),'LT0LF')
            temp='Nfstrength';
        elseif strcmp(fieldName(ifld),'LT1LF')
            temp='Sfstrength';
        else
            error('imposible case!\n');
        end
        outPlotName=strcat(outPlotName,temp);
    else
        if strcmp(fieldName(ifld),'LT0LF')
            temp='Nfstrength';
        elseif strcmp(fieldName(ifld),'LT1LF')
            temp='Sfstrength';
        else
            error('imposible case!\n');
        end
        outPlotName=strcat(outPlotName,'&',temp);
    end
    
end

if length(bc_type)~=1
    diffCases(2)=1;
    tag='$$\textrm{BC}$$';
    count=count+1;
end
for ibc=1:length(bc_type)
    if ibc==1
        outPlotName=strcat(outPlotName,'_',bc_type{ibc});
    else
        outPlotName=strcat(outPlotName,'&',bc_type{ibc});
    end
    
end

if length(sve_lx)~=1
    diffCases(3)=1;
    count=count+1;
    tag='$$\textrm{SVE Size}$$';
    
end

for ilx=1:length(sve_lx)
    if ilx==1
        outPlotName=strcat(outPlotName,'_','SVEsz',num2str(sve_lx(ilx)));
    else
        outPlotName=strcat(outPlotName,'&',num2str(sve_lx(ilx)));
    end
    
end

if length(sve_ly)~=1
    diffCases(4)=1;
    %     count=count+1;
end



if length(sveId)~=1
    diffCases(5)=1;
    count=count+1;
    tag='$$\textrm{SVE Id}$$';
    
end

for isve=1:length(sveId)
    if isve==1
        outPlotName=strcat(outPlotName,'_','SVEid',num2str(sveId(isve)));
    else
        outPlotName=strcat(outPlotName,'&',num2str(sveId(isve)));
    end
end



if length(phaseId)~=1
    diffCases(6)=1;
    count=count+1;
    tag='$$\textrm{Phase}$$';
    
end

for iph=1:length(phaseId)
    if iph==1
        outPlotName=strcat(outPlotName,'_PhId',num2str(phaseId(iph)));
    else
        outPlotName=strcat(outPlotName,'&',num2str(phaseId(iphs)));
    end
end



if count==1
    legendTag=tag;
end

end

function dataTeta=readData_teta(baseFolder,sveId,phaseId,fieldName)

%we assume sveId is sorted before!
cnt=0;


for ip=1:length(phaseId)
    for ifl=1:length(fieldName)
        id=strcat(baseFolder,'dataVSteta/','PH',num2str(phaseId),...
            fieldName{ifl},'.txt');%location of the file
        if exist(id,'file')~=2
            error('Error: Data File Does not exist!');
        end
        file=fopen(id);
        fldVal_minor=[];
        cnt_sv=0;
        while ~feof(file) && cnt_sv<length(sveId)
            temp=textscan(file,'%*s %d %*s %d %*s %d %*s %f');
            isve=temp{2}+1;
            nr=temp{3};
            temp=fgets(file);
            nc=3;
            if ismember(isve,sveId)
                [temp,file]=readBlock(file,nc,nr);
                cnt=cnt+1;
                cnt_sv=cnt_sv+1;
                dataTeta(cnt).sveId=isve;
                dataTeta(cnt).phase=phaseId(ip);
                dataTeta(cnt).fieldName=fieldName(ifl);
                dataTeta(cnt).tetaArr=temp(:,1);
                dataTeta(cnt).fldVal=temp(:,2);
            else
                [~,file]=readBlock(file,nc,nr);
            end
            
        end
        fclose(file);
    end
end

end
function [data,file]=readBlock(file,nc,nr)
data=zeros(nr,nc);
for i=1:nr
    line=fgets(file);
    data(i,:)=sscanf(line,'%f %f %d');
end
end
