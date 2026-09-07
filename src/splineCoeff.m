%
% Calculate spline coefficients (see textbook)
%
% xs,ys (input) : x and y coordinates of the spline data
% bcLeft, gLeft (input): boundary condition option and value on the left side, x=a. 
%    bcLeft = 0 : natural boundary condition (gLeft not used)
%    bcLeft = 1 : clamped boundary condition, set S'(a)=gLeft
%    bcLeft = 2 : curvature boundary condition, set S''(a)=gLeft
% bcRight, gRight (input): boundary condition option and value on the right side, x=b.
%    bcRight = 0 : natural boundary condition (gRight not used)
%    bcRight = 1 : clamped boundary condition, set S'(b)=gRight
%    bcRight = 2 : curvature boundary condition, set S''(b)=gRight
% coeff (output) : spline coefficients
%   
function coeff = splineCoeff( xs,ys, bcLeft,gLeft, bcRight,gRight )

  n=length(xs);
  A=sparse(n,n);
  r=zeros(n,1);

  for i=1:n-1
    dx(i)=xs(i+1)-xs(i); dy(i)=ys(i+1)-ys(i);
  end; 

  % Load the matrix and rhs
  for i=2:n-1
   A(i,i-1:i+1)= [ dx(i-1), 2*(dx(i-1)+dx(i)), dx(i)];

   r(i)=3.*( dy(i)/dx(i)-dy(i-1)/dx(i-1)); 
  end;

  % Set boundary conditions
  if( bcLeft==0 )
    A(1,1)=1;     % natural 
  elseif( bcLeft==1 )  
    A(1,1:2)=[2*dx(1), dx(1)]; r(1)=3.*(dy(1)/dx(1)-gLeft); % clamped 
  elseif( bcLeft==2 )
    A(1,1)=2; r(1)=gLeft;  % curvature
  else
    error('splineCoeff:ERROR: unknown bcLeft');
  end; 

  if( bcRight==0 )
    A(n,n)=1;     % natural 
  elseif( bcRight==1 )  
    A(n,n-1:n)=[dx(n-1), 2*dx(n-1)]; r(n)=3.*(gRight-dy(n-1)/dx(n-1)); % clamped 
  elseif( bcRight==2 )
    A(n,n)=2; r(n)=gRight;
  else
    error('splineCoeff:ERROR: unknown bcRight');
  end; 

  coeff=zeros(n,3);
  coeff(:,2)=A\r; % for for "c"

  for i=1:n-1
    coeff(i,3)=(coeff(i+1,2)-coeff(i,2))/(3.*dx(i));
    coeff(i,1)=dy(i)/dx(i)-dx(i)*(2.*coeff(i,2)+coeff(i+1,2))/3.;
  end;
  coeff=coeff(1:n-1,1:3); 

return
end
