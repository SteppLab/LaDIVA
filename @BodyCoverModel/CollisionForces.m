%%
% CollisionForces: Function implementing the computation of the 
% collision forces produced during closure phase in accordance to the body 
% cover model of the vocal folds.
%
% Structure: [Fu_col,Fl_col] = CollisionForces(BCMobj)
% where
%
% BCMObj: is an object from BodyCoverModel (handle) class,
% Fu_col: is the collision force in the upper mass,
% Fl_col: is the collision force in the lower mass.
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body‐cover 
%     model of the vocal folds,” J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function [Fu_col,Fl_col] = CollisionForces(BCMobj)
  % Mass displacements
  xu = BCMobj.xData(1);
  xl = BCMobj.xData(2);
  xu_col = BCMobj.xu_col;
  xl_col = BCMobj.xl_col;
  deltau = (xu-xu_col)*(xu<=xu_col);
  deltal = (xl-xl_col)*(xl<=xl_col);
  
  % Collision forces
  Fu_col = - BCMobj.hu_col*(deltau + BCMobj.NonLinMode*BCMobj.etau_col*deltau^3); % [N] Collision force in the upper mass
  Fl_col = - BCMobj.hl_col*(deltal + BCMobj.NonLinMode*BCMobj.etal_col*deltal^3); % [N] Collision force in the lower mass
end