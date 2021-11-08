if 1
    x=convn(randn(13,1000),conn_hanning(51)'/sum(conn_hanning(51)),'valid')*1.25;
    x=x+fliplr(x);
    x(11,:)=0;
    x(12,:)=.5;
    x(13,:)=x(13,:)+.5;
end

MODE=1;
fs=11025;
s=diva_synth(x,'sound');
[Aud,Som,Outline]=diva_synth(x,'AudSom');

if MODE==1, 
    soundsc(s,fs);
    K=22;
else
    fps=20;
    K=1/20/.005;
    defs_videowriteframerate=fps;
    if isempty(which('VideoWriter')), uiwait(errordlg('Sorry. VideoWriter functionality only supported on newer Matlab versions')); return; end
    videoformats={'*.avi','Motion JPEG AVI (*.avi)';'*.mj2','Motion JPEG 2000 (*.mj2)';'*.mp4;*.m4v','MPEG-4 (*.mp4;*.m4v)';'*.avi','Uncompressed AVI (*.avi)'; '*.avi','Indexed AVI (*.avi)'; '*.avi','Grayscale AVI (*.avi)'};
    [filename, pathname,filterindex]=uiputfile(videoformats,'Save video as','video01.avi');
    if isequal(filename,0), return; end
    objvideo = VideoWriter(fullfile(pathname,filename),regexprep(videoformats{filterindex,2},'\s*\(.*$',''));
    set(objvideo,'FrameRate',defs_videowriteframerate);
    open(objvideo);
end    
for time=1:K:size(x,2)
    subplot(421);
    h=bar(x(:,time));
    set(h,'facecolor','r');
    set(gca,'ylim',[-1 1],'xlim',[0 size(x,1)+1]);
    axis off;
    
    subplot(423);
    h=plot(s(max(1,min(length(s),round((time-.5)*.005*fs)+(-500:500)))),'b-');
    axis off;
    set(gca,'ylim',[-.5 .5],'xlim',[1,1001]);
    
    subplot(425);
    h=bar(Aud(:,time));
    set(h,'facecolor','b');
    set(gca,'ylim',[0 3000],'xlim',[0 size(Aud,1)+1]);
    axis off;
    
    subplot(427);
    h=bar(Som(:,time));
    set(h,'facecolor','b');
    set(gca,'ylim',[-1 1],'xlim',[0 size(Som,1)+1]);
    axis off;
    
    subplot(122);
    plot(Outline(:,time),'k','linewidth',2);
    axis equal off;
    
    if MODE==1
        drawnow;
    else
        currFrame=getframe(gcf);
        writeVideo(objvideo,currFrame);
    end
end
if MODE==1
else
    close(objvideo);
    audiowrite(regexprep(fullfile(pathname,filename),'\..*$','.wav'),s,fs);
end