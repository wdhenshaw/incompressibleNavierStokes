% ------------------------------------------------------
%   Open and position a figure 
%
% Input: 
%    fWidth, fHeight : figure width and height, e.g. 400, 350 
%
% ------------------------------------------------------
function [fh] = openFig(num,fWidth,fHeight)

windowTop=1150; windowLeft=10; 

% ******** FIX ME  *********
numFigY=2; % number of figures in vertical direction
ny = mod(num-1,numFigY);
nx = floor( (num-1)/numFigY ); 

topExtraHeight=round(75*(fHeight/450));   % accounts for height of toolbar etc. on the top of the figure
fyShift=fHeight*1.2 + topExtraHeight; 
fxShift=fWidth*1.05;

fh=figure(num);
% fh.Position
% fh.OuterPosition   % [left bottom width height]

% windowWidth=fh.OuterPosition(3);
% windowHeight=fh.OuterPosition(4);
fxShift=fWidth + fh.OuterPosition(3)-fh.Position(3) + 10; 
fyShift=fHeight+ fh.OuterPosition(4)-fh.Position(4) + 50;

set(fh,'Position',[windowLeft+nx*fxShift windowTop-ny*fyShift fWidth fHeight])

fxShift=fWidth + fh.OuterPosition(3)-fh.Position(3) + 10; 
fyShift=fHeight+ fh.OuterPosition(4)-fh.Position(4) + 50;

set(fh,'Position',[windowLeft+nx*fxShift windowTop-ny*fyShift fWidth fHeight])


end
