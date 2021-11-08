%%
% DampingForces: Function implementing the computation of the (linear and
% nonlinear) elastic forces involved in the body cover model of the vocal
% folds.
%
% Structure: [Fdu,Fdl,Fdb] = DampingForces(BCMobj)
% where
%
% BCMObj: is an object from BodyCoverModel (handle) class,
% Fdu: is the resulting elastic force in the upper mass,
% Fdl: is the resulting elastic force in the lower mass,
% Fdb: is the resulting elastic force in the body mass.
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body‐cover 
%     model of the vocal folds,” J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function [Fdu,Fdl,Fdb] = DampingForces(BCMobj)
  % Mass velocities
  vu = BCMobj.xData(4);
  vl = BCMobj.xData(5);
  vb = BCMobj.xData(6);
% Damping Equations
      Fdu = -BCMobj.du*(vu-vb);  % Damping force of upper mass [N]
      Fdl = -BCMobj.dl*(vl-vb);  % Damping force of lower mass [N]
      Fdb = -BCMobj.db*(vb);     % Damping force of body mass [N]
end