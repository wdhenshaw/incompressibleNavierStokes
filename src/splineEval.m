%
% Evaluate the spline y=S(x) and its derivative yx=S'(x)  
%
% xs,ys,coeff (input) : spline data from splineCoeff
% x (input) : evaluate the spline at these points
% y (output) : array of evaluated spline values 
% yx (output) : array of spline first derivatives
% 
function [y,yx ] = splineEval( xs,ys,coeff, x )

n=length(xs);
m=length(x);   % evaluate this many points

for i=1:m % loop over points 

  % find the interval that holds xs(i) 
  % if x(i) <= xs(1) use interval 1
  % if x(i) >= xs(n) use interval n-1 
  % -- this search is could be made faster --
  is=1; % holds spline interval
  for j=1:n-1
   is=j; 
   if( x(i) <= xs(j+1) )
     break;
   end;
  end; 

  % fprintf('i=%d: x(i)=%9.2e is evaluated from interval %d, knots=[%9.3e,%9.3e]\n',i,x(i),is,xs(is),xs(is+1))

  dx=x(i)-xs(is); 
  y(i) = ys(is) + dx*( coeff(is,1) + dx*( coeff(is,2) + dx*coeff(is,3) ) );

  yx(i) = coeff(is,1) + dx*( 2.*coeff(is,2) + 3.*dx*coeff(is,3) );
end;


