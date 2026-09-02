%
%  Construct the grid
% 
function [par] = setupGrid( par )


  % --- Setup the grid ---
  if( par.Nx<0 )
    if( strcmp(par.map,'Annulus') )
      % choose number of grid points on the Annulus to have nearly equal grid spacings
      Nr = par.N0; % Nr 
      deltaR = (par.outerRadius-par.innerRadius)/Nr;
      radiusAverage = .25*par.outerRadius + .75*par.innerRadius;
      Ntheta = ceil( 2*pi*radiusAverage/deltaR );
      par.Nx = Ntheta;
      par.Ny = Nr;
    else
      par.Nx = par.N0;
    end

  end
  if( par.Ny<0 )
    % if( strcmp(par.map,'Cartesian') || strcmp(par.map,'Rectangle)' )
    par.Ny=ceil( par.Nx * (par.yb-par.ya)/(par.xb-par.xa) );   % number of space intervals

  end

  par.dx=(par.xb-par.xa)/par.Nx;  % grid spacing 
  par.dy=(par.yb-par.ya)/par.Ny;  % grid spacing 

  par.dr(1) = par.dx;  % this will be redefined for a curvlinear grid
  par.dr(2) = par.dy;
  
  numGhost = par.orderInSpace/2;         % number of ghost points
  % numGhost = 0;%
  
  iax=1+numGhost;    iay=1+numGhost;       % index of boundary point at x=xa, y=ya
  ibx=iax+par.Nx;    iby=iay+par.Ny;       % index of boundary point at x=xb, y=yb
  Ngx=ibx+numGhost;  Ngy=iby+numGhost;     % number of grid points in x and y
  Ng = Ngx*Ngy; 
  
  % grid index range: par.gid(side,axis)
  par.gid(1,1)=iax; par.gid(2,1)=ibx; par.gid(1,2)=iay; par.gid(2,2)=iby; 

  par.dim(1,1)=1; par.dim(2,1)=Ngx;
  par.dim(1,2)=1; par.dim(2,2)=Ngy;

  if( par.idebug>1 )
    fprintf('dx=%g, dy=%g, par.gid=[%d,%d]x[%d,%d], dim=[%d,%d]x[%d,%d]\n',...
             par.dx,par.dy,par.gid(1,1),par.gid(2,1),par.gid(1,2),par.gid(2,2), par.dim(1,1),par.dim(2,1),par.dim(1,2),par.dim(2,2));
  end

  par.Ngx = Ngx;
  par.Ngy = Ngy;

  par.numGhost = numGhost; 

  % -- form the 2D grid points ---

  if( strcmp(par.map,'Cartesian') && par.gridMotion==par.noMotion)
    par.isCartesian =1;

    par.x = zeros(par.Ngx,par.Ngy,2); 
    for( iy=1:par.Ngy )
      for( ix=1:par.Ngx )
        par.x(ix,iy,1)=par.xa + (ix-iax)*par.dx; 
        par.x(ix,iy,2)=par.ya + (iy-iay)*par.dy; 
      end
    end
  else
    % evaluate a curvilinear grid 
    par.isCartesian = 0;
    par = evalMap( par );
  end

  % % -- new way for moving grids : store grid and metrics in the grid function
  % for igf=1:par.numberOfGridFunctions
  %   gf{igf}.x = par.x(:,:,:);
  %   if( ~strcmp(par.map,'Cartesian') )
  %     gf{igf}.rx = par.rx(:,:,:,:);
  %   end
  % end
  


  return
end