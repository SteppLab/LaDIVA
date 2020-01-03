function [y] = rk3m(neq,y,time,dt,p)
%
% Fourth order Runge-Kutta algorithm for Three-mass model */
% Brad H. Story    6-1-98            */
% Mcler 10-27-15 for efficiency
% yo = [];
% dyo = [];
% k1 = [];
% k2 =[];
% k3 = [];
% k4 = [];

% % 	for i=1:neq
% % 		yo(i) = y(i);
% % 	end;
% % 	
yo = y; 
	tt = time;
	dyo = eom_3m(yo,dyo,tt,p);
	
% 	for i=1:neq
% 		k1(i) = dt * dyo(i);
% 		yo(i) = y(i) + k1(i)/2.0;
% 	end;

k1 = dt*dyo;
yo = y + k1/2;

	
	tt = time + 0.5 * dt;
	dyo =eom_3m(yo,dyo,tt,p);

% 	for i=1:neq
%         k2(i) = dt * dyo(i);
%         yo(i) = y(i) + k2(i)/2.0;
%     end;

k2 = dt*dyo;
yo = y + k2/2;


	tt = time + 0.5 * dt;
  dyo =eom_3m(yo,dyo,tt,p);

% 	for i=1:neq
%     k3(i) = dt * dyo(i);
%     yo(i) = y(i) + k3(i)/2.0;
% 	end;

k3 = dt * dyo;
yo = y+k3/2;

	tt = time + dt;
  dyo =eom_3m(yo,dyo,tt,p);

% % 	for i=1:neq
% %     k4(i) = dt * dyo(i);
% %   end;

k4 = dt * dyo;

% 	for i=1:neq
%     y(i) = y(i) + k1(i)/6.0 + k2(i)/3.0 + k3(i)/3.0 + k4(i)/6.0;
% 	end;

y = y + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0;

