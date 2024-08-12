%3_16_2018 developed by B. Bahmani

classdef dataRead
    %SVE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        RVE_folder_name;% an string
        BC_type;% an string
        RVE_MicroStruc_type;%a number
        RVE_lx;%a number
        RVE_ly;%a num
        SVE_lx;%a num
        SVE_ly;%a num
        overlapingPercent;
        tetaInc;%increment value for teta
        CR;%contrast Ratio
        numSVE=-1;
        
        fieldNameTot;%all the fields in the file (may want or not)
        num_Fields_plot;%number of fields that you want to plot
        fieldsID_to_plot;%% a one dimension cell like cell{i}={'name_i',id_num_i}
        
        data_sveIndex;%is a matrix that holds index of each sve (x_ind y_ind  glob_ind)
        data_sveXfield;% a matrix that rows indicates SVE index and columns are target fields
        
    end
    
    methods
        %%
        function obj=setConfig(obj,strData,numData,fieldToPlot)
            global obj_Fields_Name;
            %fieldLib: it is an object of the "Fields_Name" class
            obj.RVE_folder_name=strData{1};
            obj.BC_type=strData{2};
            obj.RVE_MicroStruc_type=numData(1);
            obj.RVE_lx=numData(2);
            obj.RVE_ly=numData(3);
            obj.SVE_lx=numData(4);
            obj.SVE_ly=numData(5);
            obj.overlapingPercent=numData(6);
            obj.tetaInc=numData(7);
            obj.CR=numData(8);
            
            
            
            obj.num_Fields_plot=length(fieldToPlot);
            for i=1:obj.num_Fields_plot
                name_str=fieldToPlot{i};
                nameId_num=obj_Fields_Name.findFieldId_num(fieldToPlot{i});% name && id_num
                obj.fieldsID_to_plot{i}={name_str,nameId_num};
                temp_id(i)=nameId_num;
            end
            [temp_id, temp_loc]=sort(temp_id);
            tempCell=obj.fieldsID_to_plot;
            for i=1:obj.num_Fields_plot
                obj.fieldsID_to_plot{i}=...
                    {tempCell{temp_loc(i)}(1),temp_id(i)};
                %                 obj.fieldsID_to_plot{i}(2)=temp_id(i);
            end
            
            
        end
        
        
        %%
        function obj=readData(obj)
            
            id=strcat(obj.RVE_folder_name,'(AllFld_SVE).OAllFld');%location of the file
            if exist(id,'file')~=2
                error('Error: Data File Does not exist!');
            end
            file=fopen(id);
            
            line=fgets(file);%#NumGridPoint_X
            line=fgets(file);
            line=fgets(file);%#<field Type>
            line=fgets(file);
            
            fscanf(file,'%s',1);%void
            fscanf(file,'%s',1);%tenFULL
            temp=fscanf(file,'%d %d');
            numTotField=temp(2);
            
            
            line=fgets(file);%# Nfstrength
            line=fgets(file);
            line=fgets(file);%# Sfstrength
            line=fgets(file);
            line=fgets(file);%# num sve
            num=fscanf(file,'%d');%# num sve
            obj.numSVE=num;
            obj.data_sveIndex=zeros(obj.numSVE,3);
            obj.data_sveXfield=zeros(obj.numSVE,obj.num_Fields_plot);
            for i=1:numTotField
                obj.fieldNameTot{i}=fscanf(file,'%s',1);
            end
            
            %double check of consitency of fields
            for i=1:obj.num_Fields_plot
                flsID_name=obj.fieldsID_to_plot{i}(1);
                flsID_name=flsID_name{1};
                fldID_num=obj.fieldsID_to_plot{i}(2);
                fldID_num=fldID_num{1};
                
                if ~(strcmp(flsID_name,obj.fieldNameTot{fldID_num}))
                    error('The target field does not exist in data source');
                end
            end
            
            
            line=fgets(file);% #Xind
            for isve=1:obj.numSVE
% %                 fprintf('************** SVE %d\n',isve);
%                 isve
% %                 if isve==88
% %                     baha=1;
% %                 end
                
                obj.data_sveIndex(isve,1)=fscanf(file,'%d',1);
                obj.data_sveIndex(isve,2)=fscanf(file,'%d',1);
                obj.data_sveIndex(isve,3)=fscanf(file,'%d',1);
                offset=3;
                count=0+offset;
                j=1;
                for ifld=1:numTotField-offset
% %                     ifld
% %                     if ifld==28
% %                         baha=2;
% %                     end
                    count=count+1;
                    fieldTargId_num=obj.fieldsID_to_plot{j}(2);
                    if count==fieldTargId_num{1} && j<=obj.num_Fields_plot
%                         obj.data_sveXfield(isve,j)=fscanf(file,'%f',1);
                           tempN=fscanf(file,'%f',1);
%                           tempN=fread(file,'char');
                          obj.data_sveXfield(isve,j)=tempN;
                        if j<obj.num_Fields_plot
                            j=j+1;
                        end
% %                         fprintf('fild id is %d, value is %f\n',count,tempN);
                    else
                        temp=fscanf(file,'%f',1);
% %                         fprintf('fild id is %d, value is %f\n',count,temp);
                    end
                    
                end
            end
            
            
        end

    end
    
end

