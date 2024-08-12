classdef plt_plot_legendProperties
    properties
        box = 0; % 0 means no box
        location = 'Best';
        legendFontSize = 20;
    end
    methods
        % reading the object from a file:
        % format for reading is
        % [ FLG VAL ]
        % where FLG VAL's are as follows:
        % field:  FLG           VAL format
        % box             "box"       'off'(0) or 'on'(1)
        % location        "loc"       text
        % legendFontSize  "lfs"       number
        
        function objout = read(obj, fid)
            objout = plt_plot_legendProperties;
            buf = fscanf(fid, '%s', 1);
            if (strcmp(buf, '[') == 0)
                fprintf(1, 'reading plt_plot_legendProperties: data set does not start with [, instead %s\n', buf);
                pause;
            end
            
            buf = fscanf(fid, '%s', 1);
            while (strcmp(buf, ']') == 0)
                if (strcmp(buf, 'box') == 1)
                    sym = fscanf(fid, '%s', 1);
                    if ((strcmp(sym, 'off') == 1) || (strcmp(sym, '0') == 1))
                        objout.box = 0;
                    elseif ((strcmp(sym, 'on') == 1) || (strcmp(sym, '1') == 1))
                        objout.box = 1;
                    else
                        sym
                        fprintf(1, 'invalid sym for box in plt_plot_legendProperties\n');
                        pause;
                    end
                    
                elseif (strcmp(buf, 'loc') == 1)
                    objout.location = fscanf(fid, '%s', 1);
                elseif (strcmp(buf, 'lfs') == 1)
                    objout.legendFontSize = fscanf(fid, '%d', 1);
                elseif (strcmp(buf, ']') == 0)
                    buf
                    fprintf(1, 'invalid format in reading plt_plot_axisProperties %s \n', buf);
                    pause;
                end
                buf = fscanf(fid, '%s', 1);
            end
        end
        
        function h = setLegend(obj, axes, labels, isLatex,titleLgd)
            if (isLatex == 0)
                h = legend(axes,'Location',obj.location, labels);
                if ~isempty(titleLgd)
                    title(h,titleLgd);
                end
            else
                h = legend(axes,'Location',obj.location, labels, 'Interpreter', 'latex');
                if ~isempty(titleLgd)
                    title(h,titleLgd, 'Interpreter', 'latex');
                end
            end
            set(h, 'FontSize', obj.legendFontSize);
            if (obj.box == 0)
                legend('boxoff');
            end
            
        end
    end
end