%%
% AeroPressure2DrivingForces: Function implementing the computation of the 
% driving forces impiging on the cover masses due to the aerodinamic
% pressure in the glottis. The method allows for applying the expression
% introduced in the original formulation of the body cover model [1], or
% the new rules introduced later on [2]. The switching between both solvers
% is controlled by means of the function 'setDrivingForceSolver'. 
%
% Structure: [Feu, Fel] = AeroPressure2DrivingForces(BCMobj, Ps, Pe, Ae, Ph)
% where
%
% BCMObj: is an object from BodyCoverModel (handle) class,
% Ps: is the subglottal pressure in the trachea in Pascals,
% Pe: is the supraglottal pressure in the epilarynx in Pascals,
% Ae (optional, default 5e-4 [m^2]): Epilarynx (first supraglottal) tube area in meters,
% Ph (optional, default 0.0 [Pa]): hydrostatic pressure involved in collitions of cover masses in Pascals,
% Feu: is the resulting elastic force in the upper mass,
% Fel: is the resulting elastic force in the lower mass.
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body‐cover 
%     model of the vocal folds,” J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
% [2] I. R. Titze, “Regulating glottal airflow in phonation: Application of
%     the maximum power transfer theorem to a low dimensional phonation 
%     model,” J. Acoust. Soc. Am., vol. 111, no. 1, pp. 367–376, Jan. 2002.
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function [Feu, Fel] = AeroPressure2DrivingForces(BCMobj, Ps, Pe, varargin)
  % Check input variables
  if (nargin>=3)&&(nargin<=5)&&isnumeric(Ps)&&isnumeric(Pe)
    Ae = 5e-4; % [m^2] Epilarynx tube area
    Ph = 0.0; % [Pa] Hydrostatic pressure    
    if (nargin>=4)
      if isnumeric(varargin{1})&&(varargin{1}>=0)
        Ae = varargin{1}; % [m^2] Epilarynx tube area
      else
        error('Epilarynx area Ae must be a real positive numer.')
      end
    end
    if (nargin==5)
      if isnumeric(varargin{2})&&(varargin{2}>=0)
        Ph = varargin{2}; % [Pa] Hydrostatic pressure
      else
        error('Hydrostatic pressure must be a real positive numer.')
      end
    end
  else
    error('Incorrect input arguments! See function structure.')
  end
  
  % Current areas for the upper and lower masses
  au = BCMobj.au; % Upper glottal area [m^2]
  al = BCMobj.al; % Lower glottal area [m^2]
  am = max([0, min([al au])]);
  
  % Driving forces produced during contact phase of the cover masses
  if (am <= 0)
    if ~BCMobj.UseUpdated_ContactRules % Original formulation introduced in [1]
      if (au <= 0) && (al <= 0)     % Both masses colliding
        Pu = 0;                        % Upper mass pressure [N m^{-2}]
        Pl = 0;                        % Lower mass pressure [N m^{-2}]
      elseif (au <= 0) && (al > 0)  % Upper mass colliding
        Pu = 0;                        % Upper mass pressure [N m^{-2}]
        Pl = Ps;                      % Lower mass pressure [N m^{-2}]
      elseif (al <= 0) && (au > 0)  % Lower mass colliding
        Pu = Pe;                      % Upper mass pressure [N m^{-2}]
        Pl = 0;                        % Lower mass pressure [N m^{-2}]
      end
    else % New expressions introduced in [2]
      if (au <= 0) && (al <= 0)     % Both masses colliding
        Pu = Ph;                        % Upper mass pressure [N m^{-2}]
        Pl = Ph;                        % Lower mass pressure [N m^{-2}]
      else
%         zc = abs(al)/(abs(au)+abs(al))*(BCMobj.Tl+BCMobj.Tu); % Check this expression!
        zc = abs(BCMobj.xData(2)/(BCMobj.xData(2)-BCMobj.xData(1)))*(BCMobj.Tl+BCMobj.Tu); % improved expression!
        if (au <= 0) && (al > 0)  % Upper mass colliding
          if zc<BCMobj.Tl
            Pl=1/BCMobj.Tl*(zc*Ps+(BCMobj.Tl-zc)*Ph);
            Pu=Ph;
          else
            Pl=Ps;
            Pu=1/BCMobj.Tu*((zc-BCMobj.Tl)*Ps+(BCMobj.Tl+BCMobj.Tu-zc)*Ph);
          end
%           Pu = Ph;                        % Upper mass pressure [N m^{-2}]
%           Pl = Ps;                      % Lower mass pressure [N m^{-2}]
        elseif (al <= 0) && (au > 0)  % Lower mass colliding
          if zc<BCMobj.Tl
            Pl=1/BCMobj.Tl*(zc*Ph + (BCMobj.Tl-zc)*Pe);
            Pu=Pe;
          else
            Pl=Ph;
            Pu=1/BCMobj.Tu*((zc-BCMobj.Tl)*Ph + (BCMobj.Tl+BCMobj.Tu-zc)*Pe);
          end
%           Pu = Pe;                      % Upper mass pressure [N m^{-2}]
%           Pl = Ph;                        % Lower mass pressure [N m^{-2}]
        end
      end
    end
  % Driving forces produced during open phase of the glottis
  else
    if ~BCMobj.UseUpdated_AeroDrivingForces % Original formulation introduced in [1]
      Pu = Pe;                     % Upper mass pressure [N m^{-2}]
      Pl = Ps-(Ps-Pe)*(am/al)^2;   % Lower mass pressure [N m^{-2}]
    else % New expressions introduced in [2]
      % New aerodynamic terms
      spr= 1.2; %separation point ratio
      ad = max(0,min(spr*al,au));
      an = al+BCMobj.Tl*(au-al)/(BCMobj.Tl+BCMobj.Tu);
      zd = min(BCMobj.Tl+BCMobj.Tu,max(0,(ad-al)*(BCMobj.Tl+BCMobj.Tu)/(au-al)));
      ke = 2*ad/Ae*(1-ad/Ae);
%       ke = 0;
      Pkd=(Ps-Pe)/(1-ke);
      if al>=au     %convergent
        Pl = Ps - Pkd*(au^2/(al*an));
        Pu = Ps - Pkd*(au/an);
      else          %divergent
        if zd>BCMobj.Tl  %separation point above nodal point
          Pl = Ps - Pkd*(ad^2/(al*an));
          Pu = Ps - Pkd/BCMobj.Tu*(BCMobj.Tl+BCMobj.Tu-zd + (ad/an)*(zd-BCMobj.Tl));
        else      %separation point below nodal point
          Pl = Ps - Pkd/BCMobj.Tl*(BCMobj.Tl-zd+(ad/al)*zd);
          Pu = Ps - Pkd;
        end
      end
    end
  end
  
  % Aerodinamic pressure to driving forces on the cover masses
  Feu = Pu*BCMobj.Lg*BCMobj.Tu;  % Upper pressure force [N]
  Fel = Pl*BCMobj.Lg*BCMobj.Tl;  % Lower pressure force [N]
end