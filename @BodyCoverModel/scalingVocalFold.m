%%
% BCMSimulate: Function for scaling the volume of the vocal fold according 
% to the rule: VFvol_out = scale * VFvol_rest;
%
% Structure: BCMSimulate(BCMobj, scale)
% where
%
% BCMObj: is an object from BodyCoverModel (handle) class,
% scale: is the scaling factor.
%
% Coded by Gabriel Alzamendi, July 2019.
% Based on previous code by Matías Zañartu and Gabriel Galindo.
function scalingVocalFold(BCMobj,varargin)
  if nargin==1
      scale = 1.0;
  elseif (nargin==2)&&isnumeric(varargin{1})
      if varargin{1}>0
        scale = varargin{1};
      else
        error('Parameter ''scale'' should be positive!')
      end 
  else
    error('Incorrect input arguments! See function structure.')
  end
  
  % Warning message for unusually large vocal fold volume
  if scale>2
    warning('Resulting Vocal Fold volume is unusually large (larger than twice the volume at rest)!')
  end
  
  % Scaling
  BCMobj.scale = scale;
  scale_3root = nthroot(scale,3); % cube root of scale factor
  BCMobj.Lg0 = scale_3root*BCMobj.Lg_init; % [m] Vocal fold length at rest
  BCMobj.Tg0 = scale_3root*BCMobj.Tg_init;% [m] Vocal fold thickness at rest
  BCMobj.DMuc0 = scale_3root*BCMobj.DMuc_init;% [m] Mucosa depth at rest
  BCMobj.DLig0 = scale_3root*BCMobj.DLig_init;% [m] Ligament depth at rest
  BCMobj.DMus0 = scale_3root*BCMobj.DMus_init;% [m] TA muscle depth at rest

end