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

%% Run to calculate results of frazil laden plume model
% Visualize by running 
% (1) figure_setup.m 
% (2) figure_plume_summary.m
% (3) figure_plume_evolution.m

%% load run with no frazil to use as initial condition for frazil experiments
load('plume_ref_nofrazil.mat') 
yC=zeros(1,4);
xC=420e3; % frazil nucleation position
xmax=520e3; % end of calculations

%% Specifiy discretization in crystal size space
%ri_mm=[0.01,0.05,0.15,0.3,0.4,0.5,0.6,0.8,1,2]; % Approximately like SJ04
N=1e3;
ri_mm=logspace(log10(5e-3),log10(500),N); 
delta_ri=cat(2,diff(ri_mm),ri_mm(end)-ri_mm(end-1));

%% initial frazil distribution
C_initial=4e-9; % initial concentration
ri_min=0.1; % smallest nuclei
ri_max=1.0; % largest nuclei
C_star=C_initial/(ri_max-ri_min);
C_init=C_star*delta_ri; C_init(ri_mm>ri_max | ri_mm<ri_min) = 0; % top hat

for I=1:4; yC(I) = interp1(Xref,Yref(:,I),xC); end
yC(5:4+numel(ri_mm))=yC(1)*C_init; % multiply by DU

%% Run experiments showing sensitivity to growth rate and nucleation
run_counter=0;
for opt=1:2
    for tilde_nmax=[0,500,4e6]
        run_counter=run_counter+1;
        disp(run_counter)
        [vdydt,Y,vx,vz,vU,vD,vT,vS,vCi,vCik,vdelT,vm,vfp,vfpk,vpp,vppk,vTf,vdrho] = isw_plume_frazil_v1_0(xC,xmax,yC,ri_mm,opt,tilde_nmax);
    end
end