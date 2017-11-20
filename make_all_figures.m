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

% NOTES
% 1. Figure numbers are consistent with revised Cryosphere manuscript
% 2. Figure 1 (a sketch added to manuscript in revision) is not plotted

% REFERENCE:
% Rees Jones, D. W. and Wells, A. J.:
% Frazil-ice growth rate and dynamics in mixed layers and sub-ice-shelf plumes,
% The Cryosphere Discuss.,
% https://doi.org/10.5194/tc-2017-155,
% Accepted for publication in The Cryosphere, 2017.
% (Please cite final revised manuscript when available.)

close all
clear

% Growth law figure
figure_growth_law % Figure 2

% Mixed layer frazil dynamics
figure_bifurcation_example % Figure 3
figure_CSD_evolution % Figure 5
figure_growth_example % Figure 6
figure_initial_example % Figure 7
figure_CSD_steadystates % Figure 8

% Frazil explosion regime diagram
figure_regime_diagram % Figure 3
% N.B. Regime diagram Figure 3 produced from stored result files regime_data_X.mat
% These files can be recomputed by running calculate_regime_diagram.m

% ISW plume frazil dynamics
plume_sensitivity_experiments % Run calculations
figure_plume_summary % Figure 9
clean_plume_data % Removes large data files calculated by plume model
