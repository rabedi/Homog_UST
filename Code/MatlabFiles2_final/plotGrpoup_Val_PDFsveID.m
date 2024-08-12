function [  ] = plotGrpoup_Val_PDFsveID(CR)
%plot a group of data on a same figure
%field_id:id numbers of fields to plot
%field_normalization:shows a factor correspond to each field id to
%normalized
global obj_Fields_Name DataAll
[field_id,field_normalization, bc_type,sve_lx,caseInEachSVEsz,cases]=...
    readConf_fldPDF();
figure;
cntr=1;
%pppp stores legend, xaxes, yaxes default properties, see the uses below
pppp = plt_plot_plotProperties;
[diffCases,legendTag,outPlotName]=getDiffCasesIdentifier(field_id, bc_type,sve_lx,sve_lx );
sve_lx=sort(sve_lx);
dataX=sve_lx;
if size(dataX,1)>size(dataX,2)
    dataX=dataX';
end

for ifld=1:length(field_id)
    for ibc=1:length(bc_type)
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

            [pdf_val, val] = ksdensity(y, 'function', 'PDF');
            if length(field_id)>1
                ppdp = getLineProperties_PDFbased(ibc,isz,ifld);%SAMPLE_getLineProperties
            else
                ppdp = getLineProperties_PDFbased(ibc,isz);%SAMPLE_getLineProperties
            end
            
            plt_plotData_plotXYbasedOnDataSpec(ppdp, val, pdf_val);
            hold on;
            idfld_glob=DataAll{dataId}.fieldsID_to_plot{field_id(ifld)}(2);
            idfld_glob=idfld_glob{1};
            tempFldName=obj_Fields_Name.fieldLib(idfld_glob);
            tempFldName=tempFldName{1};
            tempFldName=tempFldName(3);
            leg{cntr}=getLegendTitle(tempFldName{1}, bc_type{ibc},...
                sve_lx(isz),sve_lx(isz),diffCases);
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
if length(field_id)==1
    idfld_glob=DataAll{1}.fieldsID_to_plot{field_id}(2);
    idfld_glob=idfld_glob{1};
    tempFldName=obj_Fields_Name.fieldLib{idfld_glob}(3);
    xlab = ['$$',tempFldName{1},' ',obj_Fields_Name.normalization{field_normalization(ifld)}{2},'$$'];
    ylab = '$$ \textrm{PDF} $$';
else
    xlab = '$$ Value $$';
    ylab = '$$ textrm{PDF} $$';
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

pppp.setAxesLimitsLabelsFontSize(gca, xlab, ylab, isPSFrag, xMinIn, xMaxIn, yMinIn, yMaxIn, isLatex,[],[]);



% legend
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% legend location
%pppp.legend.location = 'bestoutside';
titleLgd=legendTag;
pppp.setLegend(gca, leg, isLatex,titleLgd);


% printing ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
locTemp=strcat('../PLOTs/','PDF_',outPlotName,'.png');
print('-dpng', locTemp{1});
locTemp=strcat('../PLOTs/','PDF_',outPlotName,'.fig');
savefig(gcf, locTemp{1});




hold off;
close all;
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

function titleLgd=getLegendTitle(tempFldName, bc_type,sve_lx,sve_ly,diffCases )
%diffCases: is an array with 0 and 1 two show which of the legend variables
%are uniqe(0) or have different cases(1)

if diffCases(1)==1 && diffCases(2)==0 && diffCases(3)==0 && diffCases(4)==0
    
    titleLgd = ['$$',tempFldName,'$$'];
    
elseif diffCases(1)==0 && diffCases(2)==1 && diffCases(3)==0 && diffCases(4)==0
    
    titleLgd = ['$$',bc_type,'$$'];
    
elseif diffCases(1)==0 && diffCases(2)==0 && diffCases(3)==1 && diffCases(4)==1
    
    titleLgd = ['$$',num2str(sve_lx),'\times',num2str(sve_ly),'$$'];
    
else
    titleLgd = ['$$\textrm{Field:}',tempFldName,...
        ',','\textrm{BC:}',bc_type,',','\textrm{SVE Size:}',...
        num2str(sve_lx),'\times',num2str(sve_ly),' $$'];
end


end

function [diffCases,legendTag,outPlotName]=getDiffCasesIdentifier(field_id, bc_type,sve_lx,sve_ly )
global DataAll obj_Fields_Name
diffCases=zeros(4,1);
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

if length(sve_lx)~=1
    diffCases(3)=1;
    count=count+1;
    tag='$$\textrm{SVE Size}$$';
    
end

for ilx=1:length(sve_lx)
    if ilx==1
        outPlotName=strcat(outPlotName,'_',num2str(sve_lx(ilx)));
    else
        outPlotName=strcat(outPlotName,'&',num2str(sve_lx(ilx)));
    end
    
end

if length(sve_ly)~=1
    diffCases(4)=1;
    %     count=count+1;
end

if count==1
    legendTag=tag;
end

end


function [field_id,field_normalization, bc_type,sve_lx,...
    caseInEachSVEsz,cases]=readConf_fldPDF()
fname='Config_fldPDF.txt';
if exist(fname, 'file') ~= 2
    error('config for size eff does not exist\n');
end
file=fopen(fname);


buf=fgets(file);%#numFlds
numFld=sscanf(fgets(file),'%d');

buf=fgets(file);%#fld id(y-axis)
setUp='';
for i=1:numFld
    setUp=strcat(setUp,'%d');
end
field_id=sscanf(fgets(file),setUp);


buf=fgets(file);%#fld_normalizationId
setUp='';
for i=1:numFld
    setUp=strcat(setUp,'%d');
end
field_normalization=sscanf(fgets(file),setUp);


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