%
%  Return the partial derivatives of u and v 
%   ux(1:2,1:2)      : ux(i,j) = partial u_i / partial x_j
%   uxx(1:2,1:2,1:2) : uxx(i,j,k) = partial^2 u_i / p_j p_k 
% 
function [ux,uxx] = getDerivatives( i1,i2,un,vn,par )

  dr = par.dr(1);
  ds = par.dr(2);  

  % declare operators 
  % --- Difference Operators ---
  % derivatives of entries in the Jacobian matrix
  DJzr = @(u,I1,I2,m1,m2) ( u(I1+1,I2,m1,m2) -u(I1-1,I2,m1,m2) )/(2*dr);
  DJzs = @(u,I1,I2,m1,m2) ( u(I1,I2+1,m1,m2) -u(I1,I2-1,m1,m2) )/(2*ds);

  Dr2 = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )/(2.*dr);   % u.r to second order 
  Ds2 = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )/(2.*ds);   % u.s

  Drr2 = @(u,I1,I2) ( u(I1+1,I2) -2*u(I1,I2) +u(I1-1,I2) )/(dr^2);                             % u.rr to second order
  Dss2 = @(u,I1,I2) ( u(I1,I2+1) -2*u(I1,I2) +u(I1,I2-1) )/(ds^2);                             % u.ss
  Drs2 = @(u,I1,I2) ( u(I1+1,I2+1) - u(I1-1,I2+1) - u(I1+1,I2-1) + u(I1-1,I2-1) )/(4*dr*ds);   % u.rs
  
  rx = par.rx(i1,i2,1,1);
  ry = par.rx(i1,i2,1,2);
  sx = par.rx(i1,i2,2,1);
  sy = par.rx(i1,i2,2,2);

  rxr = DJzr(par.rx,i1,i2,1,1);
  ryr = DJzr(par.rx,i1,i2,1,2);
  sxr = DJzr(par.rx,i1,i2,2,1);
  syr = DJzr(par.rx,i1,i2,2,2);

  rxs = DJzs(par.rx,i1,i2,1,1);
  rys = DJzs(par.rx,i1,i2,1,2);
  sxs = DJzs(par.rx,i1,i2,2,1);
  sys = DJzs(par.rx,i1,i2,2,2);        

  rxx = rx*rxr + sx*rxs;
  ryy = ry*ryr + sy*rys;
  sxx = rx*sxr + sx*sxs;
  syy = ry*syr + sy*sys;

  rxy = ry*rxr + sy*rxs;
  sxy = ry*sxr + sy*sxs;

  ur = Dr2(un,i1,i2); us = Ds2(un,i1,i2);
  vr = Dr2(vn,i1,i2); vs = Ds2(vn,i1,i2);

  urr = Drr2(un,i1,i2); urs = Drs2(un,i1,i2); uss = Dss2(un,i1,i2);
  vrr = Drr2(vn,i1,i2); vrs = Drs2(vn,i1,i2); vss = Dss2(vn,i1,i2);

  ux(1,1) = rx*ur + sx*us;  % u.x 
  ux(2,1) = rx*vr + sx*vs;  % v.x 

  ux(1,2) = ry*ur + sy*us;  % u.y 
  ux(2,2) = ry*vr + sy*vs;  % v.y 

  uxx(1,1,1) = (rx^2)*urr + 2*(rx*sx)*urs + (sx^2)*uss + (rxx)*ur + (sxx)*us;         % u.xx 
  uxx(2,1,1) = (rx^2)*vrr + 2*(rx*sx)*vrs + (sx^2)*vss + (rxx)*vr + (sxx)*vs;         % v.xx 

  uxx(1,2,2) = (ry^2)*urr + 2*(ry*sy)*urs + (sy^2)*uss + (ryy)*ur + (syy)*us;         % uyy 
  uxx(2,2,2) = (ry^2)*vrr + 2*(ry*sy)*vrs + (sy^2)*vss + (ryy)*vr + (syy)*vs;         % vyy 
 
  uxx(1,1,2) = (rx*ry)*urr + (rx*sy+ry*sx)*urs + (sx*sy)*uss + (rxy)*ur + (sxy)*us;   % u.xy
  uxx(2,1,2) = (rx*ry)*vrr + (rx*sy+ry*sx)*vrs + (sx*sy)*vss + (rxy)*vr + (sxy)*vs;   % v.xy

  uxx(:,2,1) = uxx(:,1,2); % u.yx = u.xy 

  return
end