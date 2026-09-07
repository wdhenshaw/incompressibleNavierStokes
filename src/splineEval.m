%
% Evaluate a spline constructed by splineInit.
% Optionally return the first and second derivatives
%
% INPUT:
%   cs : structure returned by splineInit
%   x(i) : evale spline at these points 
% OUTPUT
%   y(:) : spline evaluated at points x(i)
%   yx(i) : optinally compute the first derivative
%   yxx(i) : optinally compute the second derivative
%
function [y,yx,yxx] = splineEval( cs,x )

  y = fnval(cs.pp,x);

  if( nargout>1 )
    yx = fnval(cs.ppx,x);
  end 
  if( nargout>2 )
    yxx = fnval(cs.ppxx,x);
  end

  return
end