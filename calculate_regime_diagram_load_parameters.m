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

%% Helper script for calculate_regime_diagram.m
N=64; ti_mm=0.05; ri0=0.2/1000;
ri_mm=logspace(log10(5e-3),log10(10),N); Vri=ri_mm/1000;
n_init=heaviside(2*ri0-Vri).*Vri*log(Vri(end)/Vri(1))/(N-1)/(2*ri0);
Qext=1200;
eps_turb=5e-3;
nu=1.95e-6;
tmax=40*60;
T_init=0;
D=1;
tilde_nmax = 4e6;
Npts=40;
Npts_N=60;