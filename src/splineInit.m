%
% Initialize the spline 
% 
% INPUT: 
%  x(:),y(:) : data points for the spline
%  conds     : end conditions
%     conds = 'clamped';
%     conds = 'periodic';
%     conds = 'complete';
%     conds = 'not-a-knot';
%     conds = 'second';    % can input 2nd derivatives as extra parameters, zero by default **FINISH ME**
% OUTPUT
%   cs : structure holding the spline data
%
% USAGE:
%     cs = splineInit( xs,ys,conds );
% 
%     m=201;
%     x=linspace(a-.1,b+.1,m)';
% 
%     [y,yx,yxx] = spline(cs , x);
%
function cs = splineInit( xs,ys,conds )

    cs.pp   = csape( xs,ys,conds );
    cs.ppx  = fnder(cs.pp, 1);   % 
    cs.ppxx = fnder(cs.pp, 2); 

  return
end