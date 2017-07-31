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


%% Calculate dependence on D (mixed layer depth)
clear
calculate_regime_diagram_load_parameters
vN=logspace(2,9,Npts_N);
vD=logspace(-1,2,Npts);
[sN,sD]=meshgrid(vN,vD);
st=zeros(length(vD),length(vN));
sC=zeros(length(vD),length(vN));
sR=zeros(length(vD),length(vN));
step_counter=0;
for J=1:length(vN)
    for I=1:length(vD)
        step_counter = step_counter + 1;
        [~,~,~,~,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,sN(I,J)*n_init,eps_turb,Qext,sD(I,J),tilde_nmax,1);
        sC(I,J)=flag_crit;
        st(I,J)=tsp;
        sR(I,J)=rsp;
        disp(['D Loop: completed ',num2str(step_counter),' out of ',num2str(Npts*Npts_N),' calculations'])
    end
end
save('regime_data_D.mat')

%% Calculate dependence on E (turbulent intensity)
clear
calculate_regime_diagram_load_parameters
vN=logspace(4,9,Npts_N);
vE=logspace(-6,0,Npts);
[sN,sE]=meshgrid(vN,vE);
st=zeros(length(vE),length(vN));
sC=zeros(length(vE),length(vN));
sR=zeros(length(vE),length(vN));
step_counter=0;
for J=1:length(vN)
    for I=1:length(vE)
        step_counter = step_counter + 1;
        [~,~,~,~,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,sN(I,J)*n_init,sE(I,J),Qext,D,tilde_nmax,1);
        sC(I,J)=flag_crit;
        st(I,J)=tsp;
        sR(I,J)=rsp;
        disp(['E Loop: completed ',num2str(step_counter),' out of ',num2str(Npts*Npts_N),' calculations'])
    end
end
save('regime_data_E.mat')

%% Calculate dependence on Q (external heat flux)
clear
calculate_regime_diagram_load_parameters
vN=logspace(4,7,Npts_N);
vQ=logspace(2,5,Npts);
[sN,sQ]=meshgrid(vN,vQ);
st=zeros(length(vQ),length(vN));
sC=zeros(length(vQ),length(vN));
sR=zeros(length(vQ),length(vN));
step_counter=0;
for J=1:length(vN)
    for I=1:length(vQ)
        step_counter = step_counter + 1;
        [~,~,~,~,flag_crit,tsp,rsp] = mixed_layer_frazil_v1_0(tmax,ri_mm,ti_mm,T_init,sN(I,J)*n_init,eps_turb,sQ(I,J),D,tilde_nmax,1);
        sC(I,J)=flag_crit;
        st(I,J)=tsp;
        sR(I,J)=rsp;
        disp(['Q Loop: completed ',num2str(step_counter),' out of ',num2str(Npts*Npts_N),' calculations'])

    end
end
save('regime_data_Q.mat')
