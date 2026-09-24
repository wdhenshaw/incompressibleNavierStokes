%
% Evaluate the derivative of polynomial TZ function 
%
% ntd,nxd,nyd,nzd (input): Specify the derivative to compute by indicating the order
%   of each partial derivative. 
%    ntd : number of time derivatives (order of the time derivative).
%    nxd : number of x derivatives (order of the x derivative).
%    nyd : number of y derivatives (order of the y derivative).
% 
% x,y (input) : arrays of x and y values
% t (input) : time
%
% NOTE: call initPoly( par ) FIRST   if you have not supplied the coefficients par.xpoly and par.tpoly
%
function result = polyTZ( ntd,nxd,nyd, x,y,t,ic,par )

   % if( nargin<8 ) d=1; end 

   degreex    = par.degreex;
   degreet    = par.degreet;
   nc         = par.nc; 
   
   if( ~isfield(par,'xpoly') )
     error('You should call initPolyTZ first!');
   end
   

  nx = size(x,1); ny=size(x,2); 

  result=zeros(nx,ny);
  
  if( nxd>degreex || nyd>degreex || ntd>degreet )
    return;
  end 

  tFactorial = factorial(ntd); 
  xFactorial = factorial(nxd);
  yFactorial = factorial(nyd); 


  % vector version 
  I1=1:nx; I2=1:ny; 
  yPow=yFactorial;
  for( iy=nyd:degreex )    % base 0 
    xPow=xFactorial;
    for( ix=nxd:degreex )  % base 0 

      result(I1,I2) = result(I1,I2) + par.xpoly(ix+1,iy+1,ic).*xPow.*yPow;
      
      if( ix<degreex )
        xPow = xPow.*x(I1,I2)*((ix+1)/(ix-nxd+1));
      end 
    end 
    if( iy<degreex )
      yPow = yPow.*y(I1,I2)*((iy+1)/(iy-nyd+1));
    end 
  end 
  
  tPow = tFactorial;
  pt = 0; 
  for( it=ntd:degreet )

    pt = pt + par.tpoly(it+1,ic)*tPow;

    if( it<degreet )
      tPow = tPow * (it+1)*t/(it-ntd+1);
    end 
  end

  result = result.*pt; 


return
end
