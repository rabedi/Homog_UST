function [  ] = plotGrpoup_SVEsize_Val_sizeEff(statPlotFlag,CR)
%for using to plot size effects
%plot a group of data on a same figure
%field_id:id numbers of fields to plot
%field_normalization:shows a factor correspond to each field id to
%normalized
%statPlotFlag:it shows which of stats you want to plot on same figure
%statPlotFlag=[minFlag meanFlaf maxFlag] with 0 and 1 numbers

%for this kind of plot we should find average of the intrested field(s)
%over all sves with same size (constructed the RVE)
global obj_Fields_Name DataAll lineDataBase
[field_id,field_normalization, bc_type,sve_lx,caseInEachSVEsz,cases]=readConf_sizeEff();

%pppp stores legend, xaxes, yaxes default properties, see the uses below
pppp = plt_plot_plotProperties;
[diffCases,legendTag,outPlotName]=getDiffCasesIdentifier(field_id, bc_type,...
    statPlotFlag );

sve_lx=sort(sve_lx);
dataX=sve_lx;
if size(dataX,1)>size(dataX,2)
    dataX=dataX';
end
[fig_mm,lgdObj_mm,colorId,ppdp_mean,ppdp_range]=legendPlot(field_id,bc_type,1);
[fig_std,lgdObj_std,~,~,~]=legendPlot(field_id,bc_type,2);

for ifld=1:length(field_id)
    
    for ibc=1:length(bc_type)
        data_mean=[];
        data_range_mm=[];
        data_range_std=[];
        %getting data for all sve sizes
        
        for isz=1:length(sve_lx)
            y=[];
            microID=getMicroID(caseInEachSVEsz(isz),cases,sve_lx(isz),...
                bc_type(ibc),CR);
            
            
            for iM=1:length(microID)
                dataId=getDataId(DataAll,bc_type{ibc},sve_lx(isz),...
                    sve_lx(isz),CR,microID(iM));
                if dataId==-1
                    continue;
                end
                y=[y;DataAll{dataId}.data_sveXfield(:,field_id(ifld))/...
                    obj_Fields_Name.normalization{field_normalization(ifld)}{1}];
            end
            
            data_mean(isz)=mean(y);
            data_std(isz)=std(y);
            
            data_range_mm(isz,1)=min(y);
            data_range_mm(isz,2)=max(y);
            
            data_range_std(isz,1)=data_mean(isz)-data_std(isz);
            data_range_std(isz,2)=data_mean(isz)+data_std(isz);
        end
        
        %hold on;
        ppdp_range.val_marker=lineDataBase.markerStyleAllTb{13,2};
        ppdp_range.val_lineColor=lineDataBase.colorNameClrTb{colorId(ibc),2};
        ppdp_range.val_markerEdgeColor=ppdp_range.val_lineColor;
        ppdp_range.val_markerSize=24;
        
        ppdp_mean.val_marker=lineDataBase.markerStyleAllTb{13,2};
        ppdp_mean.val_lineColor=lineDataBase.colorNameClrTb{colorId(ibc),2};
        ppdp_mean.val_markerEdgeColor=ppdp_mean.val_lineColor;
        ppdp_mean.val_markerSize=24;
        
        
        figure(fig_mm);
        plt_plotData_plotXYbasedOnDataSpec(ppdp_range, dataX,  data_range_mm(:,1));
        figure(fig_mm);
        plt_plotData_plotXYbasedOnDataSpec(ppdp_range, dataX,  data_range_mm(:,2));
        figure(fig_mm);
        plt_plotData_plotXYbasedOnDataSpec(ppdp_mean, dataX,  data_mean);
        x2 = [dataX, fliplr(dataX)];
        inBetween = [data_range_mm(:,1)', fliplr(data_range_mm(:,2)')];
        figure(fig_mm);
        fillhandle=fill(x2, inBetween,lineDataBase.colorNameClrTb{colorId(ibc),2});
        transparency=ppdp_mean.val_transparency/5;
        figure(fig_mm);
        set(fillhandle,'EdgeColor',lineDataBase.colorNameClrTb{colorId(ibc),2},...
            'FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color
        
        %         axis([dataX(1) dataX(end) 0.5 3.5])
        %         set(gca,'XTick',dataX)
        
        
        
        figure(fig_std);
        plt_plotData_plotXYbasedOnDataSpec(ppdp_range, dataX,  data_range_std(:,1));
        figure(fig_std);
        plt_plotData_plotXYbasedOnDataSpec(ppdp_range, dataX,  data_range_std(:,2));
        figure(fig_std);
        plt_plotData_plotXYbasedOnDataSpec(ppdp_mean, dataX,  data_mean);
        x2 = [dataX, fliplr(dataX)];
        inBetween = [data_range_std(:,1)', fliplr(data_range_std(:,2)')];
        figure(fig_std);
        fillhandle=fill(x2, inBetween,lineDataBase.colorNameClrTb{colorId(ibc),2});
        transparency=ppdp_mean.val_transparency/5;
        figure(fig_std);
        set(fillhandle,'EdgeColor',lineDataBase.colorNameClrTb{colorId(ibc),2},...
            'FaceAlpha',transparency,'EdgeAlpha',transparency);%set edge color
        
        %                 axis([dataX(1) dataX(end) 0.5 3.5])
        %         set(gca,'XTick',dataX)
        
        
    end
end





% as mentioned above, we want to use latex interpreter
isLatex = 1;
% this is for eps figures which we don't want -> set it to 0
isPSFrag = 0;

% setting x, y labels and limits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[appendAxisDataLabel2TitleX, appendAxisDataLabel2TitleY] =
if length(field_id)==1
    idfld_glob=DataAll{1}.fieldsID_to_plot{field_id}(2);
    idfld_glob=idfld_glob{1};
    tempFldName=obj_Fields_Name.fieldLib{idfld_glob}(3);
    %     xlab = '$$ \textrm{SVE Size[mm]} $$';
    xlab = '$$ \textrm{SVE Size}(\delta) $$';
    ylab = ['$$',tempFldName{1},' ',obj_Fields_Name.normalization{field_normalization(ifld)}{2},'$$'];
    
else
    %     xlab = '$$ \textrm{SVE Size[mm]} $$';
    xlab = '$$ \textrm{SVE Size}(\delta) $$';
    ylab = '$$ Value $$';
end


% for axis limits nan is NOT setting the limit for min and/or max values ->
% using a real number DOES set the limit
% if you want to set x lim min or max values uncomment
%        ONE
%or
%       BOTH
%lines below

xMinIn = dataX(1) ;
xMaxIn = dataX(end);


% yMinIn = 0.5;
% yMaxIn = 3.5;
% seqY=[yMinIn:.5:yMaxIn];

yMinIn = nan;
yMaxIn = nan;
seqY=[];

seqX=dataX;
%seqY=[yMinIn:.5:yMaxIn];


show1OverSizeTick = 1;
figure(fig_mm);
pppp.setAxesLimitsLabelsFontSize(gca, xlab, ylab, isPSFrag, xMinIn, xMaxIn,...
    yMinIn, yMaxIn, isLatex,seqX,seqY, show1OverSizeTick);


figure(fig_std);
pppp.setAxesLimitsLabelsFontSize(gca, xlab, ylab, isPSFrag, xMinIn, xMaxIn,...
    yMinIn, yMaxIn, isLatex,seqX,seqY, show1OverSizeTick);


% legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend location
%pppp.legend.location = 'bestoutside';
% % % titleLgd=legendTag;
% % % pppp.setLegend(gca, leg, isLatex,titleLgd);
lgdObj_mm.Location='best';
lgdObj_std.Location='best';


% printing ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
locTemp=strcat('../PLOTs/','SzEff_mm_',outPlotName,'.png');
figure(fig_mm);
print('-dpng', locTemp{1});
locTemp=strcat('../PLOTs/','SzEff_mm_',outPlotName,'.fig');
figure(fig_mm);
savefig(gcf, locTemp{1});


locTemp=strcat('../PLOTs/','SzEff_std_',outPlotName,'.png');
figure(fig_std);
print('-dpng', locTemp{1});
locTemp=strcat('../PLOTs/','SzEff_std_',outPlotName,'.fig');
figure(fig_std);
savefig(gcf, locTemp{1});



hold off;
close all;
end


function [field_id,field_normalization, bc_type,sve_lx,...
    caseInEachSVEsz,cases]=readConf_sizeEff()
fname='Config_SizeEff.txt';
if exist(fname, 'file') ~= 2
    error('config for size eff does not exist\n');
end
file=fopen(fname);

buf=fgets(file);%#fld id(y-axis)
field_id=sscanf(fgets(file),'%d');

buf=fgets(file);%#fld_normalizationId
field_normalization=sscanf(fgets(file),'%d');

buf=fgets(file);%#num of SVE sizes(x-axis)
numSVEsz=sscanf(fgets(file),'%d');



buf=fgets(file);%#sizes
setUp='';
for i=1:numSVEsz
    setUp=strcat(setUp,'%f');
end
sve_lx=sscanf(fgets(file),setUp);
sve_ly=sve_lx;

buf=fgets(file);%#dataNumber of each size
setUp='';
for i=1:numSVEsz
    setUp=strcat(setUp,'%d');
end
caseInEachSVEsz=sscanf(fgets(file),setUp);

buf=fgets(file);%#numBC(as different categories in same plot)
numBCtype=sscanf(fgets(file),'%d');
buf=fgets(file);%#numBC(as different categories in same plot)
for i=1:numBCtype
    bc_type{i}=fscanf(file,'%s',1);
end


buf=fgets(file);%%\n
buf=fgets(file);%#num data
numData=sscanf(fgets(file),'%d');

buf=fgets(file);%#micStruc		CR			SVEsz		BC	s
for i=1:numData
    cases(i).microStruct=fscanf(file,'%d',1);
    cases(i).CR=fscanf(file,'%f',1);
    cases(i).SVEsz=fscanf(file,'%f',1);
    cases(i).BC=fscanf(file,'%s',1);
end



end

function dataId=getDataId(DataAll,bc_type,sve_lx,sve_ly,CR,M)

dataId=-1;
tag=false;
i=1;
while tag==false && i<=length(DataAll)
    if strcmp(DataAll{i}.BC_type,bc_type) && abs(DataAll{i}.SVE_lx-sve_lx)<1e-5...
            && abs(DataAll{i}.SVE_ly-sve_ly)<1e-5 && abs(DataAll{i}.CR-CR)<1e-5 ...
            && DataAll{i}.RVE_MicroStruc_type==M
        dataId=i;
        tag=true;
    end
    i=i+1;
end

end

function ppdp = getLineProperties_PDFbased(ibc,isz,ifld)
global lineDataBase;

ppdp = plt_plotDataProp;
ppdp.val_lineStyle = lineDataBase.lineStyleTb{ibc};
ppdp.val_lineColor = lineDataBase.colorNameClrTb{isz, 2};
if nargin>2
    ppdp.val_marker = lineDataBase.markerStyleAllTb{ifld,2};
end

end

function titleLgd=getLegendTitle(dataId, fldId_Glob,bc_type,stat_name,...
    diffCases)
%diffCases: is an array with 0 and 1 two show which of the legend variables
%are uniqe(0) or have different cases(1)

global DataAll obj_Fields_Name
idfld_glob=DataAll{dataId}.fieldsID_to_plot{fldId_Glob}(2);
idfld_glob=idfld_glob{1};
tempFldName=obj_Fields_Name.fieldLib(idfld_glob);
tempFldName=tempFldName{1};
tempFldName=tempFldName(3);
tempFldName=tempFldName{1};
if diffCases(1)==1 && diffCases(2)==0 && diffCases(3)==0
    
    titleLgd = ['$$',tempFldName,'$$'];
    
elseif diffCases(1)==0 && diffCases(2)==1 && diffCases(3)==0
    
    titleLgd = ['$$',bc_type,'$$'];
    
elseif diffCases(1)==0 && diffCases(2)==0 && diffCases(3)==1
    titleLgd = ['$$',stat_name,'$$'];
else
    titleLgd = ['$$\textrm{Field:}',tempFldName,...
        ',','\textrm{BC:}',bc_type,',','\textrm{stat:}',stat_name,' $$'];
end


end

function [diffCases,legendTag,outPlotName]=getDiffCasesIdentifier(field_id,...
    bc_type, statPlotFlag)
global DataAll obj_Fields_Name
diffCases=zeros(3,1);
legendTag='';
count=0;
tag=0;
outPlotName='';
dataId=1;
if length(field_id)~=1
    diffCases(1)=1;
    count=count+1;
    tag='$$\textrm{Field}$$';
end
for ifld=1:length(field_id)
    idfld_glob=DataAll{dataId}.fieldsID_to_plot{field_id(ifld)}(2);
    idfld_glob=idfld_glob{1};
    tempFldName=obj_Fields_Name.fieldLib(idfld_glob);
    tempFldName=tempFldName{1};
    tempFldName=tempFldName(2);
    if ifld==1
        outPlotName=strcat(outPlotName,tempFldName);
    else
        outPlotName=strcat(outPlotName,'&',tempFldName);
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

if sum(statPlotFlag)>1
    diffCases(3)=1;
    count=count+1;
    tag='$$\textrm{Statistic}$$';
    
end

% % % for istat=1:length(statPlotFlag)
% % %     if statPlotFlag(istat)==0
% % %         continue;
% % %     end
% % %     name=obj_Fields_Name.statVarables{istat}{3};
% % %     if istat==1
% % %         outPlotName=strcat(outPlotName,'_',name);
% % %     else
% % %         outPlotName=strcat(outPlotName,'&',name);
% % %     end
% % %
% % % end


if count==1
    legendTag=tag;
end

end



function [legendTitle]=plotAndLegend(dataX,dataY,dataId,ifld,field_ids,ibc,bc_name,...
    isz,istat,stat_name, diffCases)

if length(field_ids)>1
    ppdp = getLineProperties_PDFbased(ibc,isz,ifld);%SAMPLE_getLineProperties
else
    ppdp = getLineProperties_PDFbased(ibc,istat);%SAMPLE_getLineProperties
end
plt_plotData_plotXYbasedOnDataSpec(ppdp, dataX, dataY);
legendTitle=getLegendTitle(dataId, field_ids(ifld),bc_name,stat_name,...
    diffCases);

hold on

end


function [fig,lgdObj,colorId,ppdp_mean,ppdp_range]=legendPlot(fld_id_local,...
    cases_name,tag)
%tag:show the special case to consider (1) working with min and max  (2)
%working with mean+sdt and mean-std
global lineDataBase DataAll obj_Fields_Name
fig=figure;

cnt=0;

dataId=1;
idfld_glob=DataAll{dataId}.fieldsID_to_plot{fld_id_local}(2);
idfld_glob=idfld_glob{1};
fld_name=obj_Fields_Name.fieldLib(idfld_glob);
fld_name=fld_name{1};
fld_name=fld_name(3);
fld_name=fld_name{1};

ppdp_mean = plt_plotDataProp;
colorId=2;
if length(cases_name)==1
    ppdp_mean.val_lineColor=lineDataBase.colorNameClrTb{colorId,2};
end
plt_plotData_plotXYbasedOnDataSpec(ppdp_mean, nan,nan);
cnt=cnt+1;
% leg{cnt}=['$$','\textrm{Mean}\ ',fld_name,'$$'];
leg{cnt}=['$$','\textrm{Mean}\ ','$$'];
hold on;

ppdp_range=ppdp_mean;
ppdp_range.val_lineStyle=lineDataBase.lineStyleTb{2};%'---'
plt_plotData_plotXYbasedOnDataSpec(ppdp_range, nan,nan);
cnt=cnt+1;
if(tag==1)
%     leg{cnt}=['$$',fld_name,'\ \textrm{Range}','$$'];
leg{cnt}=['$$','\ \textrm{Range}','$$'];
elseif(tag==2)
%     leg{cnt}=['$$',fld_name,'\ +\setminus - \textrm{S.Div Range}','$$'];
leg{cnt}=['$$','\ +\setminus - \textrm{S.Div Range}','$$'];
else
    error('the tag is not implemented')
end

ppdp_case=ppdp_range;
ppdp_case.val_lineStyle=lineDataBase.lineStyleTb{1};%'-'
ppdp_case.val_lineWidth=ppdp_range.val_lineWidth-1;
if length(cases_name)~=1
    for i=1:length(cases_name)
        cnt=cnt+1;
        ppdp_case.val_lineColor=lineDataBase.colorNameClrTb{1+i,2};
        plt_plotData_plotXYbasedOnDataSpec(ppdp_case, nan,nan);
        leg{cnt}=['$$','\textrm{',cases_name{i},'}','$$'];
        colorId=[colorId,colorId(i)+1];
    end
end


pppp = plt_plot_plotProperties;
titleLgd='';
isLatex=1;
lgdObj=pppp.setLegend(gca, leg, isLatex,titleLgd);

end

function Mid=getMicroID(Nm,cases,sveSz,BC,CR)
cnt=0;
for ic=1:length(cases)
    if abs(cases(ic).CR-CR)<1e-5 && abs(cases(ic).SVEsz-sveSz)<1e-5 ...
            && strcmp(cases(ic).BC,BC)
        cnt=cnt+1;
        Mid(cnt)=cases(ic).microStruct;
    end
end

if Nm~=cnt
    error('Error\n');
end

end


