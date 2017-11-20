**********************************************************************
*************************  frazil_dynamics ***************************
**********************************************************************

MIT License

Copyright (c) [2017] [David Rees Jones]
    David Rees Jones [david.reesjones@earth.ox.ac.uk] 
	  University of Oxford 	  
	  Earth Sciences
	  South Parks Road
	  Oxford, OX1 3AN
	  United Kingdom

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


NB: This distribution of frazil_dynamics was developed and tested
on MATLAB 2016b. Running these routines on different versions of
MATLAB may lead to compatibility issues.

REFERENCE:
Rees Jones, D. W. and Wells, A. J.:
Frazil-ice growth rate and dynamics in mixed layers and sub-ice-shelf plumes,
The Cryosphere Discuss.,
https://doi.org/10.5194/tc-2017-155,
Accepted for publication in The Cryosphere, 2017.
(Please cite final revised manuscript when available.)

CODE SUMMARY:
1. mixed_layer_frazil_v1_0.m runs a frazil dynamics calculation in a simple mixed layer
2. isw_plume_frazil_v1_0.m runs a frazil dynamics calculation in an isw plume.
3. calculate_regime_diagram.m calculates data used in the regime diagram for frazil explosions
    (warning: likely to take an hour or more)

INSTRUCTIONS:
1. Run make_all_figures.m to produce the figures from the Cryosphere paper.
2. Note that the run-time may extend to several minutes (especially figure 8).
3. Figure 3 is plotted from saved data. Run calculate_regime_diagram.m to recalculate the data
    (warning: likely to take an hour or more).
