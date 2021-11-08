%%
% BCMSimulate: Function for produce the one-step simulation of the body 
% cover model of the vocal folds.
%
% Structure: BCMSimulate(BCMobj, Ps, Pe)
% where
%
% BCMObj: is an object from BodyCoverModel (handle) class,
% Ps: is the subglottal pressure in the trachea in Pascals,
% Pe: is the supraglottal pressure in the epilarynx in Pascals.
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body‐cover 
%     model of the vocal folds,” J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function Simulate(BCMobj, Ps, Pe)
  if ~BCMobj.SimParamOK % Simulation parameters are missing
    error('The simulation parameter are missing! See ''setSimulationParameter'' function.')
  end
  
  % Sampling time
  Ts = BCMobj.Ts;
  
  % Past dynamic state variables
  xu = BCMobj.xData(1); % [m] upper mass displacement
  xl = BCMobj.xData(2); % [m] lower mass displacement
  xb = BCMobj.xData(3); % [m] body mass displacement
  vu = BCMobj.xData(4); % [m/s] upper mass velocity
  vl = BCMobj.xData(5); % [m/s] lower mass velocity
  vb = BCMobj.xData(6); % [m/s] body mass velocity
  
  % Computation of internal and external forces
  [Fku,Fkl,Fkb,Fkc] = BCMobj.ElasticForces; % Elastic forces  
  [Fdu,Fdl,Fdb] = BCMobj.DampingForces; % Damping forces 
  [Fu_col,Fl_col] = BCMobj.CollisionForces; % Collision forces
  [Feu, Fel] = BCMobj.AeroPressure2DrivingForces(Ps, Pe); % Externar driving forces
  
  % Total Force Equations
  Fu = Fku + Fdu - Fkc + Feu + Fu_col;         % Total Force in the Upper Mass [N] 
  Fl = Fkl + Fdl + Fkc + Fel + Fl_col;         % Total Force in the Lower Mass [N]
  Fb = Fkb + Fdb - (Fku + Fdu + Fkl + Fdl);  % Total Force in the Body Mass [N] (Typo in the original paper)
  
  % Acceleration
  au = Fu/BCMobj.mu;
  al = Fl/BCMobj.ml;
  ab = Fb/BCMobj.mb;

  % computation of present dynamic state variables
  % Truncated Taylor series method is applied for discretizing and solving
  % the body cover model of the vocal folds.
  xData = zeros(9,1);
  xData(1) = xu + Ts*vu + 0.5*Ts^2*au;      % Position of upper mass [m]
  xData(2) = xl + Ts*vl + 0.5*Ts^2*al;      % Position of lower mass [m]
  xData(3) = xb + Ts*vb + 0.5*Ts^2*ab;      % Position of body mass [m]
  xData(4) = vu + Ts*au;      % Velocity of upper mass [m s^{-1}]
  xData(5) = vl + Ts*al;      % Velocity of lower mass [m s^{-1}]
  xData(6) = vb + Ts*ab;      % Velocity of body mass [m s^{-1}]
  xData(7) = au; 
  xData(8) = al; 
  xData(9) = ab; 
  
  % Actuallize the dynamic state for BCMobj
  BCMobj.xData = xData;
  BCMobj.n_IterCont = BCMobj.n_IterCont + 1;
end