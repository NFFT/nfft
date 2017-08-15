% Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts
%
% This program is free software; you can redistribute it and/or modify it under
% the terms of the GNU General Public License as published by the Free Software
% Foundation; either version 2 of the License, or (at your option) any later
% version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
% details.
%
% You should have received a copy of the GNU General Public License along with
% this program; if not, write to the Free Software Foundation, Inc., 51
% Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
N=1000;	M=1000;
delta=5.520e+02+i*4.000e+02;
error=[
8	1.459e-02	1.459e-02	1.459e-02	1.748e-02	1.748e-02	1.748e-02
12	9.688e-03	9.688e-03	9.688e-03	1.212e-02	1.212e-02	1.212e-02
16	5.186e-03	5.186e-03	5.186e-03	7.868e-03	7.868e-03	7.868e-03
20	5.493e-03	5.493e-03	5.493e-03	6.925e-03	6.925e-03	6.925e-03
24	2.641e-03	2.641e-03	2.641e-03	3.190e-03	3.190e-03	3.190e-03
28	7.300e-04	7.300e-04	7.300e-04	1.181e-03	1.181e-03	1.181e-03
32	3.378e-04	3.378e-04	3.378e-04	5.849e-04	5.849e-04	5.849e-04
36	1.425e-04	1.425e-04	1.426e-04	1.829e-04	1.829e-04	1.829e-04
40	6.382e-05	6.382e-05	6.381e-05	8.337e-05	8.337e-05	8.337e-05
44	1.630e-05	1.630e-05	1.630e-05	2.820e-05	2.820e-05	2.819e-05
48	6.291e-06	6.291e-06	6.293e-06	8.636e-06	8.636e-06	8.618e-06
52	2.105e-06	2.105e-06	2.110e-06	2.701e-06	2.701e-06	2.710e-06
56	6.697e-07	6.697e-07	6.758e-07	1.062e-06	1.062e-06	1.065e-06
60	2.041e-07	2.041e-07	2.170e-07	2.543e-07	2.543e-07	2.552e-07
64	2.205e-08	2.205e-08	5.198e-08	3.413e-08	3.413e-08	6.138e-08
68	3.523e-09	3.523e-09	3.502e-08	9.928e-09	9.928e-09	3.745e-08
72	1.056e-09	1.056e-09	5.175e-08	1.489e-09	1.489e-09	5.231e-08
76	1.781e-10	1.781e-10	2.979e-08	3.771e-10	3.771e-10	2.992e-08
80	3.189e-11	3.189e-11	3.221e-08	5.446e-11	5.446e-11	3.220e-08
84	2.943e-12	2.943e-12	3.287e-08	3.739e-12	3.739e-12	3.287e-08
88	7.051e-13	7.051e-13	3.262e-08	1.039e-12	1.039e-12	3.262e-08
92	8.021e-14	8.041e-14	3.834e-08	1.675e-13	1.676e-13	3.834e-08
96	6.945e-15	7.029e-15	4.199e-08	1.540e-14	1.532e-14	4.199e-08
100	4.382e-16	7.516e-16	3.566e-08	8.655e-16	1.063e-15	3.566e-08
104	7.821e-17	4.575e-16	3.811e-08	8.676e-17	4.271e-16	3.811e-08
108	7.166e-17	5.269e-16	4.142e-08	7.515e-17	5.257e-16	4.142e-08
112	3.767e-17	4.911e-16	2.658e-08	3.993e-17	4.923e-16	2.658e-08
116	3.286e-17	2.782e-16	3.143e-08	3.264e-17	2.790e-16	3.143e-08
120	4.229e-17	4.084e-16	3.057e-08	4.211e-17	4.076e-16	3.057e-08
124	5.182e-17	3.743e-16	3.312e-08	5.084e-17	3.716e-16	3.312e-08
128	4.748e-17	3.891e-16	3.163e-08	4.748e-17	3.882e-16	3.163e-08
];
