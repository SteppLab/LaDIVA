function test_LeTalker_run_continuously_with_noise()
% %     load('..\happy_v3.mat');
% %     %with noise on
% %     while time<(ndata+1)*dt;
% %         t0=floor(time/dt);
% %         t1=(time-t0*dt)/dt;
% %     
% %         if time == 0; restartSimulation = 1; else restartSimulation=0; end
% %         [all_pulses, new_pulse, ~] = LeTalker_run_continuously(Art(11,min(ndata,1+t0)), .1, Art(7,min(ndata,1+t0)), Art(13,min(ndata,1+t0)),restartSimulation); %(.2,.2,7820, .01, restartSimulation);%(Art(11), .1, Art(7), Art(13)); %pass function CT activation, TA activation (doesn't change yet),  lung pressure -1 to 1, and voicing -1 to 1
% % 
% %         time=time+numel(new_pulse)/11025;
% %     end
    figure(99);
    subplot(2,1,1)
    load('letalker_with_noise.mat');
    plot(all_pulses);
    with = all_pulses;
    subplot(2,1,2)
    load('letalker_NO_noise.mat')
    plot(all_pulses);
    no=all_pulses;
    
    figure(44);
    plot(with, no,'.');
    sum(with==no)/length(with) % if this is 1, these are equivalent
end