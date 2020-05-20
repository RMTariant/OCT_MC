function savepic(jfig,scale,sname)
%function savepic(jfig,scale,sname)
% Saves figure(jfig) as sname as either a .png or .jpg file.
% scale = size of final figure. Not WYSIWYG ! 
%   Typically use [4 3] for standard single figure. 
%   Use [4 7] for two vertical subplots: subplot(2,1,1) & subplot(2,1,2).
%   Font sizes should be reduced when using savepic.m.
%   Text is usually shifted right. 
%   IN SUMMARY, you have to check how the final figures look.
%   One needs to learn how to use to it, but it is worth it.
%
% Examples:
%   savepic(1,[4 3],'singleplot.png')
%   savepic(1,[4 7],'verticalsubplots.jpg')
%
% Currently set to generate a 600 dpi figure (-r600). Adjustable.
%
% s.l.jacques, June 1,2017

if strcmp(sname(end-3:end),'png')
    CH = 1;
else
    CH = 2;
end

switch CH
    case 1  % png 
        set(figure(jfig),'PaperPosition',[0 0 scale],'PaperSize',scale)
        print('-dpng','-r600',sname)
        fprintf('%s saved\n',sname)

    case 2 % jpeg
        set(figure(jfig),'PaperPosition',[0 0 scale],'PaperSize',scale)
        print('-djpeg','-r600',sname)
        fprintf('%s saved\n',sname)
end
