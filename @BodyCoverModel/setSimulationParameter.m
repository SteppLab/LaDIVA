%%
% setSimulationParameter: Function for setting the simulation parameters.
%
% Structure 1: setSimulationParameter(BCMobj,fs)
% where
% BCMObj: is an object from BodyCoverModel (handle) class,
% fs: is the sampling frequency in Hertz (real numbre, >=1).
%
% Structure 2: setSimulationParameter(BCMobj,Ts)
% where
% BCMObj: is an object from BodyCoverModel (handle) class,
% Ts: is the sampling period in seconds (real numbre, <1).
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function setSimulationParameter(BCMobj,Param)
  if isnumeric(Param)&&(Param>0)
    if (Param>=1) % Param is sampling frequency in Hertz
      BCMobj.fs = Param;
      BCMobj.Ts = 1/Param;
    else % Param is sampling period in seconds
      BCMobj.Ts = Param;
      BCMobj.fs = 1/Param;
    end
    BCMobj.SimParamOK = true;
  else
    error('Incorrect input arguments! See function structures.')
  end
end