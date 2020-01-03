%Determines good targets for DIVA_GUI
%mcler 12/21/2015
% ASSUMES THAT DIVA GOES FROM 0 to 1; if from -1 to 1, change y to -1 on line 41

% % Table 1. LeTalker Inputs and Outputs
% % DIVA Art #	Input (to LeTalker)	Range  	Output (calculated from LeTalker glottal pressure output) 	Range
% % (11)	|   CT tension                  0 - 1           |    F0                          0; 68 - 613
% % (12)	|   Lung pressure               2000 - 20000    |    Intensity                   1.7 - 2384
% % (13)	|   VF prephonatory distance	0 - 0.03        |    Harmonic-to-noise ratio     -.0087 - 52.4

%load('PermutationsOut_19-Nov-2015_02-18-58 9hr_overnight.mat')


%% find outputs with F0 from 100-125
ranges = [100 125; 125 175; 175 250];
for r = 1:length(ranges)
    
    idx = find( (output(:,1)>ranges(r,1) & (output(:,1)<ranges(r,2))));
    
    for v=1:2 %variables
        if v==1
            %intensity
            str = 'int';
            a = 1.7;
            b = 2384;
            outRange = [min(output(idx,2)) max(output(idx,2))];
        else
            %HNR
            str ='HNR';
            a = -.0087;
            b = 52.4;
            outRange = [min(output(idx,3)) max(output(idx,3))];
        end
        
        str2 = ['min'; 'max'];
        for m=1:2 %min max
            
            c = outRange(m);
            %a = 1.7;
            %b = 2384;
            y = 0;
            z = 1;
            x = (c - a) * (z - y) / (b - a) + y; 
            disp([str  ' '  str2(m,:)  ': '  num2str(x)]); 
        end
    end
        
end

%% Get good inputs for output F0
ranges = [100 125; 125 250; 225 250];
ranges = [104; 219; 104];
load('fake_input.mat');
t = size(fake,2);
block = ceil(t/3);
col = 'rgb';
figure(12); hold on;
for r = 1:length(ranges)
    
    idx = find( (output(:,1)==ranges(r) ));%& output(:,3)>10 & output(:,2)>800 & output(:,2)<1000 ));
    subplot(2,1,1); hold on;plot(input(idx,2),[col(r) 'o'])
    subplot(2,1,2); hold on;plot(input(idx,3),[col(r) 'o'])

    input_ = input(idx,:);
    row = input_(ceil(length(input_)/2),:);
    rescaled = row;
    
    %The general case (when you have a value c between a and b and you want a value x between y and z), x is calculated as follows:
    %CT 
    a = 0;
    b = 1;
    c = row(1);
    y = 0;
    z = 1;
    x = (c - a) * (z - y) / (b - a) + y; 
    rescaled(1) = x;
    
    %pressure
    c = row(2);
    a= 2000;
    b= 20000;
    y = 0;
    z = 1;
    x = (c - a) * (z - y) / (b - a) + y; 
    rescaled(2) = x;
    
    %Voicing
    c = row(3);
    a= .03;
    b= 0;
    y = 0;
    z = 1;
    x = (c - a) * (z - y) / (b - a) + y; 
    rescaled(3) = x;    
    
    disp([num2str(r) ' ct: ' num2str(row(1)) ' pressure: ' num2str(row(2)) ' voicing: ' num2str(row(3))])
    disp([num2str(r) ' ct: ' num2str(rescaled(1)) ' pressure: ' num2str(rescaled(2)) ' voicing: ' num2str(rescaled(3))])

    fake(11:13,block*r-block+1:block*r) = repmat(rescaled',[1 block]); 
end


%fake(11:13,3) = 0.0001;
fake(11,:) = linspace(0,.8,length(fake(11,:))); %8000; 
fake(12,:) = .2; 
fake(13,:) = .98; %repmat(rescaled',[1 block]); 
save('fake_input_redux.mat','fake');

%CT IN = 0, .4, .8
%LeTalker (rescaled) = .5, .7, .9

%pressure IN = .2
%LeTalker (rescaled) = 12800

%Voicing IN = .98
% LeTalker (rescaled) = .0003;

idx = find( (input(:,1)==.9 )& (input(:,2)==12500 ) & (input(:,3)==0.0005 ))
output(idx,:)


%F0 %Int %HNR
%143.0000  897.5479   11.1903
%197.0000  547.3176   13.2879
%373.0000  199.5766   22.8280

%now rescale these to 0 to 1
%143.0000  897.5479->0.3760   11.1903-> 0.2137
%197.0000  547.3176->0.2290   13.2879-> 0.2537
%373.0000  199.5766->0.0831   22.8280-> 0.4357

% 130  160  340
% .3      .2      0
% .15 .2  .35
% 
% 150  220  400
% .4       .3        .12
% .25 .3 .45

%
c =  22.8280
a= -.0087;
b=  52.4;
y = 0;
z = 1;
x = (c - a) * (z - y) / (b - a) + y 


rescaled(3) = x;    
    
            str = 'int';
            a = 1.7;
            b = 2384;
            outRange = [min(output(idx,2)) max(output(idx,2))];
        else
            %HNR
            str ='HNR';
            a = -.0087;
            b = 52.4;
            outRange = [min(output(idx,3)) max(output(idx,3))];
        end
        
        str2 = ['min'; 'max'];
        for m=1:2 %min max
            
            c = outRange(m);
            %a = 1.7;
            %b = 2384;
            y = 0;
            z = 1;
            x = (c - a) * (z - y) / (b - a) + y; 
            disp([str  ' '  str2(m,:)  ': '  num2str(x)]); 

subplot(2,1,1); hold on;plot(input(idx,2),[col(r) 'o'])
subplot(2,1,2); hold on;plot(input(idx,3),[col(r) 'o'])

input_ = input(idx,:);
row = input_(ceil(length(input_)/2),:);
rescaled = row;
    




%for i=1:4

    %load('U:\Meredith Cler\Analysis-Code\DIVA + StoryTitze\Graphs for paper\PermutationsOut_19-Nov-2015_02-18-58 9hr_overnight.mat')
    %The general case (when you have a value c between a and b and you want a value x between y and z), x is calculated as follows:
    %x := (c - a) * (z - y) / (b - a) + y

%end


%pasting output here
% int min: 0.024906
% int max: 0.99988
% HNR min: 0.10522
% HNR max: 0.22547

% int min: 0.0047323
% int max: 0.86973
% HNR min: 0.019694
% HNR max: 0.5133

% int min: 0.0030129
% int max: 0.25461
% HNR min: 0.027261
% HNR max: 0.6652

% 100  125  300
% 0.02  .004  .003
% 0.1 .01  .02
% 
% 
% 125  275  350
% 1        0.86        0.25
% 0.22        0.51        0.66




