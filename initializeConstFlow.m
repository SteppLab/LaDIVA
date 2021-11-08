
% Initialize Constant FLow parameters from Lucero 2016
constFlow.PGO = 0 ; % Ap = 0 
constFlow.Ae = 3.0000e-04;     % Ae = 3 cm^2 = 3/10000 m^2
constFlow.As = 3.0000e-04;     % As = 3 cm^2 = 3/10000 m^2
constFlow.mu = 1.8600e-05;     % mu = 0.000186 g/(cm s) = 0.000186*(100/1000) kg/(m s)
constFlow.rho =  1.1400;       % rho = 0.00114 g/cm^3 = 0.00114*(1000000/1000) kg/m^3
constFlow.c = 350;             % c = 35000 cm/s = 350 m/s
constFlow.L = 0.0160;          % L = 1.6 cm = 1.6/100 m
constFlow.T = 3.5000e-05;      % T = 0.35 mm = 0.35/10000 m 
constFlow.solver = 'TITZE84M';


Qg = solveFlow(delta_p,Ag,constFlow);