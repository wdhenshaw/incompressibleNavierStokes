% declare operators 
% --- Difference Operators ---
Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y

DpxDmx = @(u,I1,I2) ( u(I1+1,I2) -2.*u(I1,I2) +u(I1-1,I2) )*(1./dx^2);  % u.xx 
DpyDmy = @(u,I1,I2) ( u(I1,I2+1) -2.*u(I1,I2) +u(I1,I2-1) )*(1./dy^2);  % u.yy

DzxDzy = @(u,I1,I2) ( u(I1+1,I2+1) - u(I1-1,I2+1) - u(I1+1,I2-1) +u(I1-1,I2-1) )*(1./(4.*dx*dy)); % u.xy 

% Define extrapolations: (is1=+1/-1 and is2=+1/-1 defines the direction ("shift") of extrapolation)
extrap3 = @(u,I1,I2,is1,is2) (3.*u(I1+is1,I2+is2) - 3.*u(I1+2*is1,I2+2*is2) + u(I1+3*is1,I2+3*is2));   % 3rd-order extrapolation

% convert (ix,iy) to equation number in the matrix:  
eqn = @(ix,iy)  1 + ix-1 + Ngx*( iy-1 );
