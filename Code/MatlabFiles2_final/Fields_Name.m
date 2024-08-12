classdef Fields_Name
    
    properties
        num_fields=13;
        %num name_on_file  name_on_plot
        fieldLib={...
              {1,'Xind',''}
              {2,'Yind',''}
              {3,'Glob_id',''}
              {4,'Nfstrength_mean','\tilde{S}_n'}%{4,'Nfstrength_mean','\mathrm{E}(\tilde{s}_n)'}%
              {5,'Sfstrength_mean','\tilde{S}_t'}%{5,'Sfstrength_mean','\mathrm{E}(\tilde{s}_t)'}
              {6,'betta_mean','\mathrm{E}(\beta)'}
              {7,'c00','$$\textrm{C00}$$'}
              {8,'c10','$$\textrm{C10}$$'}
              {9,'c11','$$\textrm{C11}$$'}
              {10,'c20','$$\textrm{C20}$$'}
              {11,'c21','$$\textrm{C21}$$'}
              {12,'c22','$$\textrm{C22}$$'}
              {13,'bulkModulus','\kappa'}
              }
          
          normalization={...
              {1,''}
              {1e6,'[\textrm{Mpa}]'}
              {1e9,'[\textrm{Gpa}]'}
              }
          
          statVarables={...
              {1,'min','\textrm{min}'}
              {2,'mean','\textrm{mean}'}
              {3,'max','\textrm{max}'}
              {4,'stdDiv','textrm{stdDiv}'}
              }
    end
    
    methods
        function id_num=findFieldId_num(obj,targField_str)
            %targField_str: it is an string you want to find its number in
            %the list
            id_num=-1;
            res=false;
            itr=0;
            while res==false && itr<obj.num_fields
                itr=itr+1;
                if strcmp(obj.fieldLib{itr}(2),targField_str)
                    id_num=itr;
                    res=true;
                end
            end
        end
    end
    
end

