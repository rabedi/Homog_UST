%3_16_2018 developed by B. Bahmani

classdef ioFormat
    %IOFORMAT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        OutFodler='../PLOTs';
        %plot format
        legfs = 18;
        axisfs = 24;
        titlefs = 22;
        fieldType_read_sveXstat={'(Tstat_SVE).OTSTAT','(Sstat_SVE).OSSTAT',...
            '(MMstat_SVE).OMMSTAT'};
        
        fieldType_write_sveXstat={'(Tstat_SVE).png','(Sstat_SVE).png',...
            '(MMstat_SVE).png'};
        
        fieldType_write_pdfXstat={'(Tstat_PDF).png','(Sstat_PDF).png',...
            '(MMstat_PDF).png'};
        
        fieldType_plot={'\tilde{s}_n/\sigma^{TH}','\tilde{s}_t/\sigma^{TH}',...
            '\beta'};
        
        StatNotation_plot={'\mathrm{min}(','\mathrm{max}(','\mathrm{E}(',...
            '\sigma('};
        
        fieldType_stat={'MIN','MAX','MEAN','STDDIV'};
    end
    
    methods
    end
    
end

