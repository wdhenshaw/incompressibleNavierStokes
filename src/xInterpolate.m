%
% Interpolate a point (xi,yi) using linear interpolation
%
function uiv = xInterpolate( xi,yi,u,v,par )

   offset = 1 + par.numGhost;

   ix = floor((xi-par.xa)/par.dx) + offset; % closest grid point less than xi
   iy = floor((yi-par.ya)/par.dy) + offset; % closest grid point less than yi

   ix = max(1,min(ix,par.Ngx-1)); 
   iy = max(1,min(iy,par.Ngy-1)); 

   alphax = (xi-par.xa)/par.dx - (ix-offset);  % The point xi is this fraction through the cell [ix,ix+1]
   alphay = (yi-par.ya)/par.dy - (iy-offset);

   if( alphax<0 || alphax>1 ) fprintf('xInterpolate:ERROR alphax=%g is outside [0,1]\n',alphax); pause; end
   if( alphay<0 || alphay>1 ) 
    fprintf('xInterpolate:ERROR alphay=%g is outside [0,1], yi=%9.3e\n',alphay,yi); 
    pause; 
   end

   if( par.idebug>1 )
     fprintf('xInterpolate: [xi,yi]=[%9.3e,%9.3e] [ix,iy]=[%3d,%3d], [x,y]=[%9.3e,%9.3e] [alphax,alphay]=[%5.3f,%5.3f]\n',...
              xi,yi,ix,iy,par.x(ix,iy,1),par.x(ix,iy,2),alphax,alphay);
   end

   uiv(1) = (1-alphay)*( (1-alphax)*u(ix,iy  ) + alphax*u(ix+1,iy  ) ) ...
              +alphay *( (1-alphax)*u(ix,iy+1) + alphax*u(ix+1,iy+1) );

   uiv(2) = (1-alphay)*( (1-alphax)*v(ix,iy  ) + alphax*v(ix+1,iy  ) ) ...
              +alphay *( (1-alphax)*v(ix,iy+1) + alphax*v(ix+1,iy+1) );
  return
end