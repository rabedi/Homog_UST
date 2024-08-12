classdef plt_plot_axisProperties
    properties
        % for label = 'auto'
        % we may need to add log, ..., and various other operators before
        % and after auto label
        % for now, I permit reading of (preLabel and postLabel)
        % both are set to none in default meaning that they won't get a
        % value.
        % in future option auto for these will figure out if:
        %   1. label is in reverse order
        %   2. what operations are done on the data
        preLabel = 'none';
        postLabel = 'none';

        
        scale = 'linear'; % other option is 'log'
        minorTick = 'on'; % other option is 'off'
        dir = 'normal'; % other option is 'reverse'
        limType = 0;
        % options are:
        % 0:    inactive
        % 1:    relative (e.g. [0, 1] goes from actual minimum data to
        % maximum [-0.2 1.2] goes 20% over on each side
        % 2:    absolute input values are used
        limMin = 0;
        limMax = 1;
        % these two values are used in computing with options 1 and 2 above
        
        % for limType = 1 the following two parameters are needed
        limTypeRelOption = 1;  %  for limType = 1: option:
        %  1: simply computes xMin, xMax based on limMin, limMax and actual min max of data
        %  2: first does 1 ^ but then goes a bit
        %  beyond the limits to get roughly
        %  limTypeRelNum of subdivisions of round 10
        %  based subdivisions (refer to subdivide
        %  function
        limTypeRelNum = 10;    %  only active for limType = 1 & limTypeRelOption = 2
        %  it's the suggested number of subdivisions
        label = 'auto'; % options are:
        % 1. given text, e.g., 'x1' it will apear as actual label
        % 2. 'auto' the label is obtained from the names of the input fields
        % 3. 'none' no label will appear

    end
    methods
        % reading the object from a file:
        % format for reading is
        % [ FLG VAL ]
        % where FLG VAL's are as follows:
        % field:  FLG           VAL format
        % scale             "scale"       'linear'(0) or 'log'(1)
        % minorTick         "mTick"       'on' (0)     or 'off' (1)
        % dir               "dir"         'normal'(0)  or 'reverse' (1)
        % limType           "limT"        0 ('off', 'no')
        %                                 1 ('rel')
        %                                 2 ('abs')
        % limMin, limMax    "lim"         limMin limMax (two numbers
        % limTypeRelOption  "limTOpt"     1 or 2
        % limTypeRelNum     "limTNum"     (number)
        % label             "label"       text: options are 
                                          % 'auto' auto from input names
                                          % 'none' no label appears
                                          % 'actual text'
        % 1. given text, e.g., 'x1' it will apear as actual label
        % 2. 'auto' the label is obtained from the names of the input fields
        % 3. 'none' no label will appear
        
%         function obj = plt_plot_axisProperties(fid)
%             obj = read(obj, fid);
%         end
        function objout = read(obj, fid)
            objout = plt_plot_axisProperties;
            buf = fscanf(fid, '%s', 1);
            if (strcmp(buf, '[') == 0)
                fprintf(1, 'reading plt_plot_axisProperties: data set does not start with [, instead %s\n', buf);
                pause;
            end
            
            buf = fscanf(fid, '%s', 1);
            while (strcmp(buf, ']') == 0)
                if (strcmp(buf, 'scale') == 1)
                    sym = fscanf(fid, '%s', 1);
                    if ((strcmp(sym, 'linear') == 1) || (strcmp(sym, '0') == 1))
                        objout.scale = 'linear';
                    elseif ((strcmp(sym, 'log') == 1) || (strcmp(sym, '1') == 1))
                        objout.scale = 'log';
                    else
                        sym
                        fprintf(1, 'invalid sym for scale in plt_plot_axisProperties\n');
                        pause;
                    end
                    
                elseif (strcmp(buf, 'mTick') == 1)
                    sym = fscanf(fid, '%s', 1);
                    if ((strcmp(sym, 'on') == 1) || (strcmp(sym, '1') == 1))
                        objout.minorTick = 'on';
                    elseif ((strcmp(sym, 'off') == 1) || (strcmp(sym, '0') == 1))
                        objout.minorTick = 'off';
                    else
                        sym
                        fprintf(1, 'invalid sym for minorTick in plt_plot_axisProperties\n');
                        pause;
                    end
                    
                elseif (strcmp(buf, 'dir') == 1)
                    sym = fscanf(fid, '%s', 1);
                    if ((strcmp(sym, 'normal') == 1) || (strcmp(sym, '0') == 1))
                        objout.dir = 'normal';
                    elseif ((strcmp(sym, 'reverse') == 1) || (strcmp(sym, '1') == 1))
                        objout.dir = 'reverse';
                    else
                        sym
                        fprintf(1, 'invalid sym for dir in plt_plot_axisProperties\n');
                        pause;
                    end
                    
                elseif (strcmp(buf, 'limT') == 1)
                    sym = fscanf(fid, '%s', 1);
                    if ((strcmp(sym, 'off') == 1) || (strcmp(sym, '0') == 1))
                        objout.limType = 0;
                    elseif ((strcmp(sym, 'rel') == 1) || (strcmp(sym, '1') == 1))
                        objout.limType = 1;
                    elseif ((strcmp(sym, 'abs') == 1) || (strcmp(sym, '2') == 1))
                        objout.limType = 2;
                    else
                        sym
                        fprintf(1, 'invalid sym for limType in plt_plot_axisProperties\n');
                        pause;
                    end
                    
                elseif (strcmp(buf, 'lim') == 1)
                    objout.limMin = fscanf(fid, '%lg', 1);
                    objout.limMax = fscanf(fid, '%lg', 1);
                elseif (strcmp(buf, 'limTOpt') == 1)
                    objout.limTypeRelOption = fscanf(fid, '%d', 1);
                    if ((objout.limTypeRelOption ~= 1) && (objout.limTypeRelOption ~= 2))
                        limTypeRelOption = objout.limTypeRelOption
                        fprintf(1, 'invalid limTypeRelOption\n');
                        pause;
                    end
                elseif (strcmp(buf, 'limTNum') == 1)
                    objout.limTypeRelNum = fscanf(fid, '%d', 1);
                elseif (strcmp(buf, 'label') == 1)
                    str = fscanf(fid, '%s', 1);
                    objout.label = readStringWithoutSpace(str, '*', ' ');
                elseif (strcmp(buf, 'preLabel') == 1)
                    str = fscanf(fid, '%s', 1);
                    objout.preLabel = readStringWithoutSpace(str, '*', ' ');
                elseif (strcmp(buf, 'postLabel') == 1)
                    str = fscanf(fid, '%s', 1);
                    objout.postLabel = readStringWithoutSpace(str, '*', ' ');
                elseif (strcmp(buf, ']') == 0)
                    buf
                    fprintf(1, 'invalid format in reading plt_plot_axisProperties %s \n', buf);
                    pause;
                end
                buf = fscanf(fid, '%s', 1);
            end
        end
        
        function [lims, setlims] = getAxisLimits(obj, minLimitActualPlot, maxnLimitActualPlot)
            if (obj.limType == 0)
                lims = [minLimitActualPlot, maxnLimitActualPlot];
                if ((isnan(minLimitActualPlot)) && (isnan(maxnLimitActualPlot)))
                    setlims = 0;
                else
                    setlims = 1;
                end
                return;
            elseif (obj.limType == 1)
                expandFactor = [obj.limMin obj.limMax];
                [lims(1) lims(2)] = subdivide(minLimitActualPlot, maxnLimitActualPlot, obj.limTypeRelNum, expandFactor, obj.limTypeRelOption);
                setlims = 1;
                return;
            elseif (obj.limType == 2)
                lims = [obj.limMin obj.limMax];
                setlims = 1;
                return;
            else
                limtType = obj.limType
                fprintf(1, 'invalid limtType\n');
                pause;
            end
        end
        
        function [labelOut, setLabel, appendAxisDataLabel2Title] = getAxisLabel(obj, labelIn, axis, isPSFrag, isLatex)
            if nargin < 3
                isPSFrag = 0;
            end
            if nargin < 4
                isLatex = 0;
            end
            
            appendAxisDataLabel2Title = 1;
            if (strcmp(obj.label, 'none') == 1)
                labelOut = '';
                setLabel = 0;
                return;
            end
            setLabel = 1;
            if (isPSFrag == 1)
                labelOut = [axis, 'Label'];
                return;
            end
            if (strcmp(obj.label, 'auto') == 1)
                labelOut = [];
                if (strcmp(obj.preLabel, 'none') == 0)
                    if (isLatex == 0)
                        labelOut = [obj.preLabel, labelOut];
                    else
                        pn = obj.preLabel;
                        if (strcmp(pn, 'log[') == 1)
                            pn = '$$\log$$(';
                        end
                        labelOut = [pn, labelOut];
                    end
                end
                 labelOut = [labelOut, labelIn];
                if (strcmp(obj.postLabel, 'none') == 0)
                    if (isLatex == 0)
                         labelOut = [labelOut, obj.postLabel];
                    else
                        pn = obj.postLabel;
                        if (strcmp(pn, ']') == 1)
                            pn = ')';
                        end
                         labelOut = [labelOut, pn];
                    end
                end
                 appendAxisDataLabel2Title = 0;
            else
                 labelOut = obj.label;
            end
        end
    end
end