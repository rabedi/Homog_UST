classdef plt_plotDataProp
   properties
       val_lineStyle = '-';
       val_lineColor = [0 0 0];
       val_lineWidth = 2;
       val_marker = 'none';
       val_markerEdgeColor = 'k';
       val_markerFaceColor = 'none';
       val_markerSize = 8;
       val_transparency=1.0;%[0 1]
       
       flg_lineStyle = 0;
       flg_lineColor = 0;
       flg_lineWidth = 0;
       flg_marker = 0;
       flg_markerEdgeColor = 0;
       flg_markerFaceColor = 0;
       flg_markerSize = 0;
       flg_markerFullDataSet = 1; % uses full data set which starts from those that have interior
   end
end

% flags meaning:    (Default value is 0 except marker colors)
%    0   flags are completely inactive  (the final property of the plot is
%    not determined from this data set)
%    employed
%   -1   corresponding value is obtained from the table from the input
%   number (e.g. data set number)
%   -2  corresponding val_ is directly employed
%   n >= 1 => values are again given from the table (BUT from number (n)
%   not from data set number (for flg_markerEdgeColor and
%   flg_markerFaceColor refer to below)

% other special flags:
% for markers:
%       flg_markerEdgeColor
%       flg_markerFaceColor
%   -3 means that color again is got from dataNumber (same as flag == -1) 
%   BUT with the difference that 
%   dataNumber == 1 (MakerEdge(Face) gets color of line)
%   dataNumber == 2 (MakerEdge(Face) gets no color)
%   -4 means it always gets the color of line (regardless of dataNumber)