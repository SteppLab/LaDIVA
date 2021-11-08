function Qg = getFlow(Ag,Ps,constFlow) % Modified for source only transmission line model
%function Qg = solveFlow(Ps_plus,Pe_minus,Am,constFlow)
%   config = inputParser;

%   config.addRequired('Ps_plus',@isnumeric); % Incident sub-glottal pressure [Pa]
%   config.addRequired('Pe_minus',@isnumeric);  % Incident supra-glottal pressure [Pa]
%   config.addRequired('Am',@isnumeric); % Membranous glottal area [m^2]
%   config.addOptional('Ae',0,@isnumeric);  % Supra-glottal area [m^2]
%   config.addOptional('As',0,@isnumeric);  % Sub-glottal area [m^2]
%   config.addOptional('rho',1.14,@isnumeric);% air density [kg m^-3]
%   config.addOptional('mu',18.36922e-6,@isnumeric); % air viscosity [Pa s]
%   config.addOptional('c',350,@isnumeric); % sound velocity [m s^-1]
%   config.addOptional('L',1,@isnumeric); % length of the vocal folds [m]
%   config.addOptional('T',1,@isnumeric); % thickness of the vocal folds [m]
%   config.addOptional('PGO',0,@isnumeric); % Posterior glottal area [m^2]
%   config.addOptional('PGD',0,@isnumeric); % Posterior glottal displacement [m]
%   config.addOptional('solver','galindo2016',@ischar); % solver to be used
%   
%   config.parse(varargin{:});
%   
%   Am = config.Results.Am;
%   Ap = config.Results.PGO;
%   Ae = config.Results.Ae;
%   As = config.Results.As;  
%   mu = config.Results.mu;
%   rho = config.Results.rho;
%   c = config.Results.c;  
%   L = config.Results.L;
%   T = config.Results.T;  
%   Ps_plus = config.Results.Ps_plus;
%   Pe_minus = config.Results.Pe_minus;

  %Ap = constFlow.PGO;
  Ae = constFlow.Ae;
  As = constFlow.As;  
  mu = constFlow.mu;
  rho = constFlow.rho;
  c = constFlow.c;  
  L = constFlow.L;
  T = constFlow.T;
  delta_p = Ps;
  
  switch upper(constFlow.solver)
    case 'TITZE84M' %?
      
      k_t = (1-Ag/Ae).^2 + (Ag/Ae).^2; % HW k_t = 1 for simplicity?

      Qg = sign(delta_p).*Ag.*sqrt(abs(delta_p)./(0.5*rho*k_t)); % HW consider no load condition in Titze 1984 Equation 36
      
    case 'LUCERO2015'
      Ag = Am + Ap;
      A_star = (Ae^-1 + As^-1)^-1;

      k_t = (1-Ag/Ae)^2 + (Ag/Ae)^2;

      r_s = (As - Ag)/(As + Ag); % Titze and Worley - JASA - 2009 - Modeling source-filter interaction in belting and pitched operatic male singing
      r_e = (Ae - Ag)/(Ae + Ag); % Titze and Worley - JASA - 2009 - Modeling source-filter interaction in belting and pitched operatic male singing
      r_s = 1.0;
      r_e = 0.8;
      k_t = 1.0;

      delta_p = (1+r_s)*Ps_plus - (1+r_e)*Pe_minus;

      gamma_m = 12*mu*(L^2)*T;

      Qg =  (2*(Ag^3)*delta_p)/(((rho*c*(Ag^3))/A_star) + gamma_m + sqrt((((rho*c*(Ag^3))/A_star) + gamma_m)^2 + (2*k_t*rho*(Ag^4)*abs(delta_p)))); % Equation 19
      
    case 'GALINDO2016'
      Ag = Am + Ap;
      A_star = (Ae^-1 + As^-1)^-1;

      k_t = (1-Ag/Ae)^2 + (Ag/Ae)^2;

      r_s = (As - Ag)/(As + Ag); % Titze and Worley - JASA - 2009 - Modeling source-filter interaction in belting and pitched operatic male singing
      r_e = (Ae - Ag)/(Ae + Ag); % Titze and Worley - JASA - 2009 - Modeling source-filter interaction in belting and pitched operatic male singing
      r_s = 1.0;
      r_e = 1.0; 
      k_t = 1.0;

      delta_p = (1+r_s)*Ps_plus - (1+r_e)*Pe_minus; % 2 (Ps_plus  -Pe_minus)

      gamma_m = 12*mu*(L^2)*T; 
      gamma_p = 8*mu*T*pi; 

      Rm = gamma_m/(Am^3); 
      Rp = gamma_p/(Ap^2); 

      Rg = (Rm^-1 + Rp^-1)^-1;

      if (Ap <= (eps*1e5))
        Qg =  (2*(Ag^3)*delta_p)/(((rho*c*(Ag^3))/A_star) + gamma_m + sqrt((((rho*c*(Ag^3))/A_star) + gamma_m)^2 + (2*k_t*rho*(Ag^4)*abs(delta_p)))); 
      else
        Qg = sign(delta_p)*(c*Ag/k_t)*(-Ag*(1/A_star+Rg/(rho*c))+sqrt((Ag*(1/A_star+Rg/(rho*c)))^2 + (2*k_t/(rho*c^2))*abs(delta_p)));
      end
    case 'OFF'
      Qg = 0;
  end
end

