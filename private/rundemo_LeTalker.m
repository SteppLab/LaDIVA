%Run Code as Demonstration 
%B. Story
%10.25.11

clear all
Fs = 44100;
InitializeLeTalker
[p, r] = LeTalker(p,c,4000,Fs,[],[],[]);
PlotLeTalkerWaveforms(r,p,Fs,1,2);