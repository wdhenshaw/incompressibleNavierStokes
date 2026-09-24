%
% Evaluate the derivative of trigonometric TZ function 
%      ---- ue = sin(kx*(x+phix)).*sin(ky*(y+phiy))*sin(kt*(t+phit)) -----
%
% ntd,nxd,nyd,nzd (input): Specify the derivative to compute by indicating the order
%   of each partial derivative. 
%    ntd : number of time derivatives (order of the time derivative).
%    nxd : number of x derivatives (order of the x derivative).
%    nyd : number of y derivatives (order of the y derivative).
% 
% x,y (input) : arrays of x and y values
% t (input) : time
% ic (input) : component number
%
% Parameter input:
%   par.nc                : number of components
%   par.kx(nc) : 
%   par.ky(nc) : 
%   par.kt(nc) : 
%   par.phix(nc) : phase in x
%   par.phiy(nc) : phase in y
%   par.phit(nc) : phase in time
%
function result = trigTZ( ntd,nxd,nyd, x,y,t,ic,par )

  % if( nargin<8 ) d=1; end 

  % kx         = par.ky;
  % ky         = par.kx;
  % kt         = par.kt;
  nc         = par.nc; 

  %  if( ~isfield(par,'phit') )
  %    par.phit = zeros(nc,numDomains); 
  %  end


  % printArray(kx,'trigTZ: kx','%f5.2 ');
  
  nx = size(x,1); ny=size(x,2); 
  result=zeros(nx,ny);

  I1=1:nx; I2=1:ny; 

  xsign = 1-2*mod(floor(nxd/2),2);  % +1, +1, -1, -1, ...
  ysign = 1-2*mod(floor(nyd/2),2);
  tsign = 1-2*mod(floor(ntd/2),2);
  
  % fprintf('trigTZ: ntd=%d, nxd=%d, nyd=%d, xsign=%g, ysign=%g, tsign=%g\n',ntd,nxd,nyd,xsign,ysign,tsign);
  
  % ---- ue = sin(kx*x).*sin(ky*y)*sin(kt*t) -----
  kx = par.kx(ic); %  2*pi;   % kxv(ic);
  ky = par.ky(ic); % 2.5*pi; % kyv(ic);
  kt = par.kt(ic); % 3.3*pi; % ktv(ic);
  phit = par.phit(ic); % phase 
  phix = par.phix(ic);
  phiy = par.phiy(ic); 
  if( mod(ntd,2)==0 )
    timeFactor = tsign*(kt^ntd)*sin(kt*(t+phit));
  else
    timeFactor = tsign*(kt^ntd)*cos(kt*(t+phit));
  end 
  if( mod(nxd,2)==0 )
    xFactor = ( timeFactor*xsign*(kx^nxd) )*sin(kx*(x+phix));
  else
    xFactor = ( timeFactor*xsign*(kx^nxd) )*cos(kx*(x+phix));
  end 
  if( mod(nyd,2)==0 )
    yFactor = (ysign *(ky^nyd))*sin(ky*(y+phiy));
  else
    yFactor = (ysign *(ky^nyd))*cos(ky*(y+phiy));
  end 
  result(I1,I2) = xFactor.*yFactor;


return
end
