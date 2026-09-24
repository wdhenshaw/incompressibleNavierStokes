%
%  Test the polynomial TZ function 
%
% Examples:
%    testPoly
%    testPoly -degreex=1 -degreet=2 -userCoeff=1
%

function testPoly(varargin)

  clearvars -except varargin;
  clear global;

  % --- Clear all open figures ----
  clearOpenFigures(1:2);
  rainbowMap = getRainbow();
  
  ax=0; bx=1;
  ay=0; by=1;

  par.nc              = 2;
  par.degreex         = 2;
  par.degreet         = 2;
  par.userCoeff       = 0; % 1 = use user supplied coeff
  % par.numDomains = 1; 

  par.echo=0;
  for i = 1 : nargin
    line = varargin{i};
    par = assignCommandLineOption( line, par, par.echo );    
  end 

  
  Nx=10; Ny=10; 
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

  nc = par.nc; 


  % Initialize the polynomial TZ function
  if( par.userCoeff==0 )
    % generate default coefficients 
    par = initPolyTZ( par );
  else
    % user defined coefficients
    par.xpoly = zeros(par.degreex+1,par.degreex+1,nc);
    par.xpoly(:,:,1)=1;
    par.xpoly(:,:,2)=1/2;

    par.tpoly = zeros(par.degreet+1,nc);
    par.tpoly(:,1)=0.5;
    par.tpoly(:,2)=0.25;

  end

  I1=1:Ngx; I2=1:Ngy; 
  
  t=.1;
  ntd=0; nxd=0; nyd=0; 
  u = zeros(Ngx,Ngy,nc);
  for( ic=1:nc )
    subplot(nc,1,ic);  
    u(:,:,ic ) = polyTZ( ntd,nxd,nyd, x,y,t,ic,par );
    surf( x,y,u(I1,I2,ic) );
    title(sprintf('polyTZ, component %d',ic)); xlabel('x'); ylabel('y'); 
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


   % define an anonymous function to eval a derivative
   ic=1; 
   uex = @(x,y,t) polyTZ( 0,1,0, x,y,t,ic,par );

   % -- check a derivative ---
   nxd=1;
   % ux = polyTZ( ntd,nxd,nyd, x,y,t,par );
   ux = uex(x,y,t);

   nxd=0;
   delta = eps^(1/3); 
   uxp = polyTZ( ntd,nxd,nyd, x+delta,y,t,ic,par );
   uxm = polyTZ( ntd,nxd,nyd, x-delta,y,t,ic,par );
   uxd = (uxp-uxm)/(2*delta);

   maxErr = max(max(max(abs(ux-uxd))));
   fprintf('Max error in ux  = %8.1e\n',maxErr); 

   % -- check a derivative ---
   nyd=2;
   uyy = polyTZ( ntd,nxd,nyd, x,y,t,ic,par );

   nyd=0;
   delta = eps^(1/4);
   uyp = polyTZ( ntd,nxd,nyd, x,y+delta,t,ic,par );
   uz  = polyTZ( ntd,nxd,nyd, x,y      ,t,ic,par );
   uym = polyTZ( ntd,nxd,nyd, x,y-delta,t,ic,par );
   uyyd= (uyp-2*uz+uym)/(delta^2);

   maxErr = max(max(max(abs(uyy-uyyd))));
   fprintf('Max error in uyy = %8.1e\n',maxErr); 



   % -- check a derivative ---
   ntd=1;
   ut = polyTZ( ntd,nxd,nyd, x,y,t,ic,par );

   ntd=0;
   delta = eps^(1/3);
   utp = polyTZ( ntd,nxd,nyd, x,y,t+delta,ic,par );
   utm = polyTZ( ntd,nxd,nyd, x,y,t-delta,ic,par );
   utd = (utp-utm)/(2*delta);

   maxErr = max(max(max(abs(ut-utd))));
   fprintf('Max error in ut  = %8.1e\n',maxErr); 


  return
end