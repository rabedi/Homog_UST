function [  ] = plotGrpoup_SVEid_Val( field_id,field_normalization, bc_type,sve_lx,sve_ly )
%plot a group of data on a same figure
%field_id:id numbers of fields to plot
%field_normalization:shows a factor correspond to each field id to
%normalized
global obj_Fields_Name DataAll
figure;
cntr=1;
%pppp stores legend, xaxes, yaxes default properties, see the uses below
pppp = plt_plot_plotProperties;
[diffCases,legendTag]=getDiffCasesIdentifier(field_id, bc_type,sve_lx,sve_ly );
for ifld=1:length(field_id)
    for ibc=1:length(bc_type)
        for isz=1:length(sve_lx)
            dataId=getDataId(DataAll,bc_type{ibc},sve_lx(isz),sve_ly(isz));
            if dataId==-1
                continue;
            end
            y=DataAll{dataId}.data_sveXfield(:,field_id(ifld))/obj_Fields_Name.normalization{field_normalization(ifld)}{1};
            x=DataAll{dataId}.data_sveIndex(:,3);
            
            if length(field_id)>1
                ppdp = getLineProperties_PDFbased(ibc,isz,ifld);%SAMPLE_getLineProperties
            else
                ppdp = getLineProperties_PDFbased(ibc,isz);%SAMPLE_getLineProperties
            end
            
            plt_plotData_plotXYbasedOnDataSpec(ppdp, x, y);
            hold on;
            idfld_glob=DataAll{dataId}.fieldsID_to_plot{field_id(ifld)}(2);
            idfld_glob=idfld_glob{1};
            tempFldName=obj_Fields_Name.fieldLib(idfld_glob);
            tempFldName=tempFldName{1};
            tempFldName=tempFldName(3);
            leg{cntr}=getLegendTitle(tempFldName{1}, bc_type{ibc},...
                sve_lx(isz),sve_ly(isz),diffCases);
            
            % % % %
            % % % %             % adding legeng entry
            % % % %             if length(field_id)==1
            % % % %                 leg{cntr} = ['$$','BC:',bc_type{ibc},',','SVE Size:',...
            % % % %                     num2str(sve_lx(isz)),'\times',num2str(sve_ly(isz)),' $$'];
            % % % %             else
            % % % %                 idfld_glob=DataAll{dataId}.fieldsID_to_plot{field_id(ifld)}(2);
            % % % %                 idfld_glob=idfld_glob{1};
            % % % %                 tempFldName=obj_Fields_Name.fieldLib{idfld_glob}(3);
            % % % %                 leg{cntr} = ['$$','Field:',tempFldName{1},...
            % % % %                     ',','BC:',bc_type{ibc},',','SVE Size:',...
            % % % %                     num2str(sve_lx(isz)),'\times',num2str(sve_ly(isz)),' $$'];
            % % % %             end
            % % % %
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
    xlab = '$$ \textrm{SVE} $$';
    ylab = ['$$',tempFldName{1},' ',obj_Fields_Name.normalization{field_normalization(ifld)}{2},'$$'];
else
    error('fsfd');
    ylab = '$$ textrm{PDF} $$';
    ylab = '$$ Value $$';
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
print('-dpng', 'sample.png');
savefig(gcf, 'sample.fig');




hold off;
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
        ',','\textrm{BC:}',bc_type{ibc},',','\textrm{SVE Size:}',...
        num2str(sve_lx),'\times',num2str(sve_ly),' $$'];
end


end

function [diffCases,legendTag]=getDiffCasesIdentifier(field_id, bc_type,sve_lx,sve_ly )
diffCases=zeros(4,1);
legendTag='';
count=0;
tag=0;
if length(field_id)~=1
    diffCases(1)=1;
    count=count+1;    
    tag='$$\textrm{Field}$$';
end

if length(bc_type)~=1
    diffCases(2)=1;
    tag='$$\textrm{BC}$$';
    count=count+1;
end

if length(sve_lx)~=1
    diffCases(3)=1;
    count=count+1;
    tag='$$\textrm{SVE Size}$$';
end

if length(sve_ly)~=1
    diffCases(4)=1;
%     count=count+1;
end

if count==1
    legendTag=tag;
end

end
