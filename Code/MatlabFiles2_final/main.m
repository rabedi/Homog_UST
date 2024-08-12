%3_19_2018 Developed by B. Bahmani
clc
clear all
close all
fclose all
set(0,'DefaultLegendAutoUpdate','off');

%%% RA: to load new data set use loadData='', it will read new data and
%%% save it in saveDataName.mat file
loadData='dataOld.mat';%'data1.mat';
loadData = '';
%readDataLocation='config_plot.txt';
readDataLocation='config_plot_Square.txt';

saveDataName='dataAllSV';
global lineDataBase DataAll obj_Fields_Name;

%if loadData is NOT '', data is loaded from its mat file, otherwise it goes
%to else and reads the data by the instructions in the readDataLocation
%file
if ~isempty(loadData)
    lineDataBase = plt_plotDataPropTable;
    obj_Fields_Name=Fields_Name;
    DataAll=load(loadData);
    DataAll=DataAll.DataAll;
else
    readRawData(readDataLocation,saveDataName);
end

% only sizes that are here are processed from the data read
% same with bc_type, etc. (field_id - one of the "13" fields that now
% exists, field_normalization is given in Field_Names (1 no normaization,
% 2,3 are GPa and MPa)
sve_lx=[3.125 6.25 12.5 25];
sve_ly=sve_lx;
bc_type={'D'};
field_id=[1];% this is a local: so if one data set has 5 fields for plotting 2 in here means second data of that 6 fields not in global field (49 fields)
field_normalization=[1];% it reads data in "normalization" of Fields_Name class
% include min, mean, max, not include sdiv (see statVarables in
% Fields_Name)
statPlotFlag=[1 1 1 0];
CR=100;




%%%you should make configs for the following tasks

if 1
    %first check the config file ->"Config_SizeEff.txt"
    plotGrpoup_SVEsize_Val_sizeEff(statPlotFlag,CR);%plotting size effects
end
if 1
    %first check the config file ->"config_fldPDF.txt"
    plotGrpoup_Val_PDFsveID(CR );%PDF base plots
end
if 0
    %first check the config file ->"config_convPDF.txt"
    plotGrpoup_Val_PDFsveID_Converg(CR);
end

% angle plot ONE SVE (angle versus strength)--- for Dr. Soghrati
% fieldName_teta='Nfstrength_mean';
% field_normalization=1;
% bc_type='D';
% sve_lx=12.5;
% sveId=2;
% phaseId=1;
% plotGrpoup_Val_teta( fieldName_teta,field_normalization, bc_type,...
%     sve_lx,sve_ly,sveId,phaseId )



function readRawData(fileLocation,saveDataName)%'..//config//config.txt'

global lineDataBase DataAll obj_Fields_Name;
lineDataBase = plt_plotDataPropTable;
obj_Fields_Name=Fields_Name;


ioFrmt=ioFormat;
confFile=fopen(fileLocation);


line=fgets(confFile);%#nCases	nCol
line=fgets(confFile);%nCases	nCol
line=sscanf(line,'%d %d');
numCase=line(1);
numColumn=line(2);


line=fgets(confFile);%#number of field to plot -1:means all
line=fgets(confFile);%num
numField=sscanf(line,'%d');


line=fgets(confFile);%#fields to plot
for i=1:numField
    fieldToPlot{i}=fscanf(confFile,'%s',1);
end
line=fgets(confFile);%array of numbers
% % % listId='';
% % % for i=1:numField
% % %     listId=strcat(listId,'%d');
% % % end
% % % fieldID=sscanf(line,listId);%list of filed to plot

line=fgets(confFile);%#FileRoot					BC_type									MicroStruct_type			RVE_lx			RVE_ly			SVE_lx		SVE_ly			overlapPersent 			tensOpt		shearOpt	bettaOpt	tetaOpt		tetaIncr

DataAll=cell(numCase,1);
caseCount=0;
while caseCount<numCase
    temp=fscanf(confFile,'%s',1);%rootFile
    if temp(1)=='#'%if it is comented
        line=fgets(confFile);
        line=[];
        continue;
    else
        caseCount=caseCount+1;
    end
    strData{1}=temp;
    strData{2}=fscanf(confFile,'%s',1);%BC
    line=fgets(confFile);
    numData=sscanf(line,'%d %f %f %f %f %f %f %f');
    DataAll{caseCount}=dataRead;
    DataAll{caseCount}=DataAll{caseCount}.setConfig(strData,numData,fieldToPlot);
    DataAll{caseCount}=DataAll{caseCount}.readData();
end


if ~isempty(saveDataName)
    save(saveDataName,'DataAll');
end


end