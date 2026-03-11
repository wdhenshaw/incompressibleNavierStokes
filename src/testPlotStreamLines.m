%
% Test function for plotting streamlines
%
function testPlotStreamLines( varargin )

 clearvars -except varargin;
   % --- Clear all open figures ----
  clearOpenFigures(1:7);

  addpath(genpath(pwd)); % allow matlab to find files in subfolders

  par.Nx = 50;
  par.Ny = -1; % by default set to par.Nx
  par.xa=0.; par.xb=1.; par.ya=0.; par.yb=1.;  % space interval interval
  par.numGhost=1;
  par.kx=1; 
  par.ky=1;
  par.echo=0;
  par.idebug=1;

  % load the rainbow colour table
  rainbow;


  % --- read command line args ---
  for i = 1 : nargin
    line = varargin{i};

    % New way for any variable in par
    par = assignCommandLineOption( line, par, par.echo );

  end

  if( par.Ny==-1 ) par.Ny=par.Nx; end

  par.kx = par.kx*2*pi;
  par.ky = par.ky*2*pi;

  ax = par.xa; bx = par.xb;
  ay = par.ya; by = par.yb;

  Nx = par.Nx;
  Ny = par.Ny;

  numGhost = par.numGhost;



  % --- Setup the grid ---
  dx=(bx-ax)/Nx; dy=(by-ay)/Ny;  % grid spacing 
  
  iax=1+numGhost;    iay=1+numGhost;       % index of boundary point at x=xa, y=ya
  ibx=iax+Nx;        iby=iay+Ny;           % index of boundary point at x=xb, y=yb
  Ngx=ibx+numGhost;  Ngy=iby+numGhost;     % number of grid points in x and y
  Ng = Ngx*Ngy; 
  
  % grid index range: gid(side,axis)
  gid(1,1)=iax; gid(2,1)=ibx; gid(1,2)=iay; gid(2,2)=iby; 


  par.dx   = dx;
  par.dy   = dy;

  par.gid = gid;
  % par.iax = iax; par.ibx= ibx;
  % par.iay = iay; par.iby= iby;
  par.Ngx = Ngx;
  par.Ngy = Ngy;

  % -- form the 2D grid points ---
  par.x = zeros(Ngx,Ngy,2); 
  for( iy=1:Ngy )
    for( ix=1:Ngx )
      par.x(ix,iy,1)=ax + (ix-iax)*dx; 
      par.x(ix,iy,2)=ay + (iy-iay)*dy; 
    end
  end

  u = zeros(Ngx,Ngy);
  v = zeros(Ngx,Ngy);

  % ue.x + ve.y = 0 if kx=ky
  ue = @(x,y) sin(par.kx*x).*sin(par.ky*y);
  ve = @(x,y) cos(par.kx*x).*cos(par.ky*y);

  u = ue(par.x(:,:,1),par.x(:,:,2));
  v = ve(par.x(:,:,1),par.x(:,:,2));


 % Test interpolation routine:
  if( 1==0 )
    Ni = 4;
    xiv = linspace(par.xa,par.xb,Ni)';
    yiv = linspace(par.ya,par.yb,Ni)';
    maxErr=0;
    for( ky=1:Ni )
    for( kx=1:Ni )
      xi = xiv(kx); 
      yi = yiv(ky);
      uiv = xInterpolate( xi,yi,u,v,par );
      ui = uiv(1); vi=uiv(2);

      uierr = abs(ue(xi,yi) - ui);
      vierr = abs(ve(xi,yi) - vi);
      if( par.idebug>0 )
        fprintf('[Nx,Ny]=[%3d,%3d] [ui,vi]=[%10.3e,%10.3e] err=[%9.2e,%9.2e]\n',Nx,Ny,ui,vi,uierr,vierr);
      end
      maxErr=max(maxErr,max(uierr,vierr));
    end
    end
    fprintf('[Nx,Ny]=[%3d,%3d] maxErr=%9.3e\n',Nx,Ny,maxErr);
  end 

  if( 1==0 )
    % draw a coloured line
    figure(3)
    x = linspace(1,10,50);
    y = sin(x);
    y(end) = NaN;
    c = y;
    patch(x,y,c,'EdgeColor','interp','LineWidth',3,'LineJoin','round');
    colormap(par.rainbowMap); colorbar; 
    title('test draw a coloured line');

  end

  if( 1==1 )
    figure(1)
    tl =tiledlayout('flow','TileSpacing','Compact');

    [I1,I2]=getIndex(par.gid);

    nexttile;
    surf( par.x(I1,I2,1),par.x(I1,I2,2),u(I1,I2) ); hold on
    contour3( par.x(I1,I2,1),par.x(I1,I2,2),u(I1,I2),'k-' ); 
    colormap(par.rainbowMap); colorbar; shading interp; 
    xlabel('x'); ylabel('y'); view(0,90);  % top view 
    title('u')
    setAspectRatio();

    nexttile;
    surf( par.x(I1,I2,1),par.x(I1,I2,2),v(I1,I2) ); hold on
    contour3( par.x(I1,I2,1),par.x(I1,I2,2),v(I1,I2),'k-' ); 
    colormap(par.rainbowMap); colorbar; shading interp; 
    xlabel('x'); ylabel('y'); view(0,90);  % top view 
    title('v')
    setAspectRatio();

    hold off
  end

  figure(2)
  par = plotStreamLines( u,v,par );

  % % streamline mask:
  % par.nxg=ceil(par.Nx/2);
  % par.nyg=ceil(par.Ny/2);
  % par.maskForStreamLines = zeros(par.nxg,par.nyg);

  % par.uMin = min( min(min(u)), min(min(v)) );
  % par.uMax = max( max(max(u)), max(max(v)) );

  % par.streamLineStoppingTolerance=1e-3; 



  % figure(2);

  % xa = par.xa;
  % xb = par.xb;
  % ya = par.ya;
  % yb = par.yb;

  % xba=xb-xa;
  % yba=yb-ya;  
  % nxg = par.nxg;
  % nyg = par.nyg;  

  % % plot a fake line off screen
  % h = plot( [xa-1,xa-1], [ya,yb], 'k-' );

  % size = max(xb-xa,yb-ya); % over-all domain size : determines size of arrows 

  % hold on;
  % for( i=1:par.nxg )
  %   for( j=1:par.nyg )
  %     if( par.maskForStreamLines(i,j)==0 )  % no streamline has passed through this point
  %       xtp=xa+xba*(i+.5)/nxg;  % starting point for streamline
  %       ytp=ya+yba*(j+.5)/nyg;
  %       % first integrate backwards in time from this spot
  %       par.cflStreamLine = -.5; % plot in backward direction
  %       [xp,yp,cp,par] = drawAStreamLine( xtp,ytp, u,v,par );  

  %       % yp(end)=NaN;
  %       % fprintf('xp=%d, yp=%d, cp=%d\n',length(xp),length(yp),length(cp));
  %       % pause
  %       if( length(xp)>3 )
  %         patch([xp,xp(end)],[yp,NaN],[cp,cp(end)],'EdgeColor','interp','LineWidth',2,'LineJoin','round');
   
  %         if( length(xp)>10 )
  %           plotArrows(h, xp, yp, 'number', 1, 'color', 'k', 'LineWidth', 1,'scale', 0.5, 'ratio', 'equal', 'size',size );
  %         end
  %       end
  %       % plot( xp,yp,'k-'); 
  %       hold on

  %       % ---  Now plot the streamline in the forward direction ---
  %       par.cflStreamLine = .5; % plot in forward direction
  %       [xp,yp,cp,par] = drawAStreamLine( xtp,ytp, u,v,par ); 
  %       % fprintf('xp=%d, yp=%d, cp=%d\n',length(xp),length(yp),length(cp));
  %       if( length(xp)>3 )
  %         patch([xp,xp(end)],[yp,NaN],[cp,cp(end)],'EdgeColor','interp','LineWidth',2,'LineJoin','round');

  %         if( length(xp)>10 )
  %           % h = plot( [xa-1,xb-1], [ya,yb], 'k-' ); % plot a fake line off screen
  %           plotArrows(h,flip(xp), flip(yp), 'number', 1, 'color', 'k', 'LineWidth', 1,'scale', 0.4, 'ratio', 'equal', 'size',size ); 
  %         end         
  %       end       
  %       % plot( xp,yp,'k-');
  %     end
  %   end
  % end

  % hold off;
  % colormap(par.rainbowMap);  colorbar; % shading interp; 
  % xlabel('x'); ylabel('y');
  % xlim([xa,xb]);
  % ylim([ya,yb]);
  % title('streamlines');


  if( 1==0 )

  figure(4)
   t = [0:0.01:20];
   x = t.*cos(t);
   y = t.*sin(t);
   arrowPlot(x, y, 'number', 10, 'LineWidth', 2 );

   figure(5)
   t = [0:0.01:20];
   x = t.*cos(t);
   y = t.*sin(t);
   arrowPlot(x, y, 'number', 5, 'color', 'r', 'LineWidth', 1, 'scale', 0.8, 'ratio', 'equal');

end


  % xtp = .35; ytp=.3; 
  % par.cflStreamLine = .5; % plot in forward direction
  % [xp,yp,cp,par] = drawAStreamLine( xtp,ytp, u,v,par );  
  % plot( xp,yp,'k-'); hold on

  % par.cflStreamLine = -.5; % plot in backward direction
  % [xp,yp,cp,par] = drawAStreamLine( xtp,ytp, u,v,par );  
  % plot( xp,yp,'k-'); hold on

  % xtp = .65; ytp=.7; 
  % par.cflStreamLine = .5; % plot in forward direction
  % [xp,yp,cp,par] = drawAStreamLine( xtp,ytp, u,v,par );  
  % plot( xp,yp,'k-'); hold on

  % par.cflStreamLine = -.5; % plot in backward direction
  % [xp,yp,cp,par] = drawAStreamLine( xtp,ytp, u,v,par );  
  % plot( xp,yp,'k-'); hold on

  % hold off;
  % title('streamline');

  % grid on;




 return
end