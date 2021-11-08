%%
% ElasticForces: Function implementing the computation of the (linear and
% nonlinear) elastic forces involved in the body cover model of the vocal
% folds.
%
% Structure: [Fku,Fkl,Fkb,Fkc] = ElasticForces(BCMobj)
% where
%
% BCMObj: is an object from BodyCoverModel (handle) class,
% Fku: is the resulting elastic force in the upper mass,
% Fkl: is the resulting elastic force in the lower mass,
% Fkb: is the resulting elastic force in the body mass,
% Fkc: is the resulting coupling force between the upper and lower mases.
%
% References:
% [1] B. H. Story and I. R. Titze, “Voice simulation with a body‐cover 
%     model of the vocal folds,” J. Acoust. Soc. Am., vol. 97, no. 2, 
%     pp. 1249–1260, Feb. 1995.
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function [Fku,Fkl,Fkb,Fkc] = ElasticForces(BCMobj)
  % Mass displacements
  xu = BCMobj.xData(1);
  xl = BCMobj.xData(2);
  xb = BCMobj.xData(3);
  xu0 = BCMobj.xu0;
  xl0 = BCMobj.xl0;
  xb0 = BCMobj.xb0;
  
  % Elastic forces
  Fku = - BCMobj.ku*(((xu-xu0)-(xb-xb0)) + BCMobj.NonLinMode*BCMobj.etau*((xu-xu0)-(xb-xb0))^3); % [N] Elastic force in the upper mass
  Fkl = - BCMobj.kl*(((xl-xl0)-(xb-xb0)) + BCMobj.NonLinMode*BCMobj.etal*((xl-xl0)-(xb-xb0))^3); % [N] Elastic force in the lower mass
  Fkb = - BCMobj.kb*((xb-xb0) + BCMobj.NonLinMode*BCMobj.etab*(xb-xb0)^3); % [N] Elastic force in the body mass
  Fkc = - BCMobj.kc*((xl-xl0) - (xu-xu0)); % [N] Coupling force between the upper and lower mases

end