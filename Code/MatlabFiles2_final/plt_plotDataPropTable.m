classdef plt_plotDataPropTable
   properties
       lineWidthTb = {2, 3, 4, 5, 6};
       lineStyleTb = {'-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.', ':', '-', '--', '-.', ':'};
       colorBWTb = {[0 0 0], ...
                      [0.3 0.3 0.3], ...
                      [0.5 0.5 0.5], ...
                      [0.7 0.7 0.7]};
       colorNameClrTb = { 'black', [0	0	0];  
                            'green', [0 135/255 0];
                            'red', [1	0	0]; % comment out for Coh Trajectory Conv
                            'blue',	[0	0	1];
                            'magenta',	[1	0	1];
                            'brown', [0.5	0.25	0];
                            'dark_blue', [0	0	0.5];
                            'orange', [1	102/255	0];
                            'teal', [0	1	1];
                            'red2', [203/255 0   51/255];
                            'green_blue', [0	0.5	0.25];
                            'light_gray', [0.75	0.75	0.75];
                            'peach', [1	0.5	0.5];
                            'yellow', [1	1	0]; 
                            'arghavani', [0.5	0	0.25];  
                            'dark_gray2', [0.5	0.5	0.5];   
                            'purple', [0.5	0	1];
                            'olive', [0.5	0.5	0];
                            'blue2', [0.25	0.5	0.5];
                            'rosy_pink', [1 0.5	0.75];	   
                            'purple2', [0.5	0	0.5];
                            'black', [0	0	0]
                            };

     markerStyleAllTb = {'Circle',                              'o';	
                    'Square',                              's';	
                    'Right-pointing triangle',             '>';
                    'Upward-pointing triangle',            '^';
                    'Downward-pointing triangle',          'v';
                    'Left-pointing triangle',              '<';
                    'Diamond',                             'd';	
                    'Five-pointed star (pentagram)',       'p';
                    'Six-pointed star (hexagram)',         'h';
                    'Cross',                               'x';
                    'Plus sign',                           '+';	
                    'Asterisk',                            '*';
                    'Point',                               '.'
                    'none',                                'none';  };
                        
     markerStyleNoAllTb = {'none',                  'none';  
                    'Circle',                              'o';	
                    'Square',                              's';	
                    'Cross',                               'x';
                    'Plus sign',                           '+';	
                    'Diamond',                             'd';	
                    'Five-pointed star (pentagram)',       'p';
                    'Right-pointing triangle',             '>';
                    'Upward-pointing triangle',            '^';
                    'Six-pointed star (hexagram)',         'h';
                    'Downward-pointing triangle',          'v';
                    'Left-pointing triangle',              '<';
                    'Asterisk',                            '*';
                    'Point',                               '.'};
        markerSizeTb = {8;
                        12;
                        13;
                        15;
                        17;
                        };

%         markerSizeTb = {5;
%                         7;
%                         9;
%                         11;
%                         13;
%                         15;
%                         };
   end
end
