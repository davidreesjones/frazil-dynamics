% MIT License
% 
% Copyright (c) [2017] [David Rees Jones]
%     David Rees Jones [david.reesjones@earth.ox.ac.uk] 
% 	  University of Oxford 	  
% 	  Earth Sciences
% 	  South Parks Road
% 	  Oxford, OX1 3AN
% 	  United Kingdom
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.

close all
clear

% Growth law figure
figure_growth_law % Figure 1

% Mixed layer frazil dynamics
figure_bifurcation_example % Figure 2
figure_CSD_evolution % Figure 4
figure_growth_example % Figure 5
figure_initial_example % Figure 6
figure_CSD_steadystates % Figure 7

% Frazil explosion regime diagram
figure_regime_diagram % Figure 3
% N.B. Regime diagram Figure 3 produced from stored result files regime_data_X.mat
% These files can be recomputed by running calculate_regime_diagram.m

% ISW plume frazil dynamics
plume_sensitivity_experiments % Run calculations
figure_plume_summary % Figure 8
clean_plume_data % Removes large data files calculated by plume model
