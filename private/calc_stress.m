function p = calc_stress(p)
% Stress calculations for the three-mass model based on Titze and Story
% 2002
% Author: Brad Story
%edits MCler 10272015 for efficiency

ep1 = [-.5;-.5;-.5];
%strs = [0,0,0];
ep2 = [-.35; 0.0; -.05]; %zeros(3,1);
sig0 = [5000; 4000; 10000];
sig2 = [300000; 13930; 15000];
Ce = [4.4; 17; 6.5];

strs = (-sig0./ep1).*(ones(3,1) * p.eps - ep1);

for i=1:3
  if(p.eps > ep2(i) )
    strs(i) = strs (i) + sig2(i)*( exp(Ce(i)*(p.eps-ep2(i)))-1-Ce(i)*(p.eps-ep2(i)));
  end;  
end;

p.sigmuc = strs(1);
p.sigl = strs(2);
p.sigp = strs(3);
