%
%  Test the trig TZ function 
%
  clear global; clf;
  
  rainbowMap = getRainbow();
  
  ax=0; bx=1;
  ay=0; by=1;
  
  Nx=20; Ny=20; 
  Ngx=Nx+1; Ngy=Nx+1; 
  dx=(bx-ax)/Nx; dy=(by-ay)/Ny;  % grid spacing
  iax=1; iay=1;
  
  % -- form the 2D grid points ---
  x = zeros(Ngx,Ngy); y = zeros(Ngx,Ngy); % grid 
  for( iy=1:Ngy )
    for( ix=1:Ngx )
      x(ix,iy)=ax + (ix-iax)*dx; 
      y(ix,iy)=ay + (iy-iay)*dy; 
    end;
  end;

  nc = 2; 
  par.nc=nc;
  par.degreex =2;
  par.degreet =2;
  par.numDomains = 1; % not used in this version 

  numDomains=par.numDomains; 
  par.kx=zeros(nc,numDomains);  
  par.ky=zeros(nc,numDomains);  
  par.kt=zeros(nc,numDomains);
  par.phit=zeros(nc,numDomains);  % phase in t
  par.phix=zeros(nc,numDomains);  % phase in x
  par.phiy=zeros(nc,numDomains);  % phase in y
  for( ic=1:nc )
    for( d=1:numDomains )
      par.kx(ic,d) = (2+ic+d/2)*pi;
      par.ky(ic,d) = (3+1.5*ic+d/2)*pi;
      par.kt(ic,d) = (3-ic+d)*pi;
    end
  end 

  I1=1:Ngx; I2=1:Ngy; 
  
  t=.1;
  ntd=0; nxd=0; nyd=0; 
  u = zeros(Ngx,Ngy,nc);
  for( ic=1:nc )
    subplot(nc,1,ic);  
    u(:,:,ic) = trigTZ( ntd,nxd,nyd, x,y,t,ic,par );
    surf( x,y,u(I1,I2,ic) );
    title(sprintf('trigTZ, component %d',ic)); xlabel('x'); ylabel('y'); 
    % view(0,90);  % top view
    colormap(rainbowMap);  colorbar; hold on;
    shading interp;
  end
  hold off; 

  if( nc==2 )
    % resize the window 
    xwidth = 400;
    ywidth = 600;
    
    pos = get(gcf,'position');  
    pos(3) = xwidth;
    pos(4) = ywidth;
    set(gcf,'position',pos);
  
  end


 delta = sqrt(eps); 

 ic=1; 
 uex = @(x,y,t) trigTZ( 0,1,0, x,y,t,ic,par );
 uey = @(x,y,t) trigTZ( 0,0,1, x,y,t,ic,par );

 % -- check a derivative ---
 nxd=0;
 % ux = trigTZ( ntd,nxd,nyd, x,y,t,par );
 ux = uex(x,y,t);

 nxd=0;
 uxp = trigTZ( ntd,nxd,nyd, x+delta,y,t,ic,par );
 uxm = trigTZ( ntd,nxd,nyd, x-delta,y,t,ic,par );
 uxd = (uxp-uxm)/(2*delta);

 maxErr = max(max(max(abs(ux-uxd))));
 fprintf('Max error in ux  = %8.1e\n',maxErr); 

 uy = uey(x,y,t);
 uyp = trigTZ( ntd,nxd,nyd, x,y+delta,t,ic,par );
 uym = trigTZ( ntd,nxd,nyd, x,y-delta,t,ic,par );
 uyd = (uyp-uym)/(2*delta);

 maxErr = max(max(max(abs(uy-uyd))));
 fprintf('Max error in uy  = %8.1e\n',maxErr); 

 % -- check a derivative ---
 nyd=2;
 uyy = trigTZ( ntd,nxd,nyd, x,y,t,ic,par );
 nyd=0;

 delta2 = eps^(1/3); 
 uyp = trigTZ( ntd,nxd,nyd, x,y+delta2,t,ic,par );
 u   = trigTZ( ntd,nxd,nyd, x,y       ,t,ic,par );
 uym = trigTZ( ntd,nxd,nyd, x,y-delta2,t,ic,par );
 uyyd= (uyp-2*u+uym)/(delta2^2);

 maxErr = max(max(max(abs(uyy-uyyd))));
 fprintf('Max error in uyy = %8.1e\n',maxErr); 



 % -- check a derivative ---
 ntd=1;
 ut = trigTZ( ntd,nxd,nyd, x,y,t,ic,par );

 ntd=0;
 utp = trigTZ( ntd,nxd,nyd, x,y,t+delta,ic,par );
 utm = trigTZ( ntd,nxd,nyd, x,y,t-delta,ic,par );
 utd = (utp-utm)/(2*delta);

 maxErr = max(max(max(abs(ut-utd))));
 fprintf('Max error in ut  = %8.1e\n',maxErr); 


