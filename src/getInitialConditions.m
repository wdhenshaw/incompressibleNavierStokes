%
%  Assign initial conditions
%
function [gf,par] = getInitialConditions( t,gf,cur, par )

  if( strcmp(par.ic,'default') )
    
    gf{cur}.u(:,:) = par.u0(par.x(:,:,1),par.x(:,:,2)); 
    gf{cur}.v(:,:) = par.v0(par.x(:,:,1),par.x(:,:,2)); 

  elseif( strcmp(par.ic,'zero') )

    gf{cur}.u(:,:)=0.;
    gf{cur}.v(:,:)=0.;

  elseif( strcmp(par.ic,'constant') )
  	
   gf{cur}.u(:,:) = par.uic;
   gf{cur}.v(:,:) = par.vic;

  elseif( strcmp(par.ic,'shear') )
    
   % shear flow 
   par.plotErrors=0;

   % un = zeros(par.Ngx,par.Ngy);
   % vn = zeros(par.Ngx,par.Ngy);

   beta = par.shearBeta;
   ym = 0.5*(par.ya + par.yb); % center in y 
   gf{cur}.u(:,:) = tanh( beta*( par.x(:,:,2)-ym ));

   % perturbation for v 
   delta=par.shearDeltav;
   gf{cur}.v(:,:) = delta*cos( 2*pi*(par.x(:,:,2)-ym) ).*sin( par.kx*par.x(:,:,1));

  % elseif( strcmp(par.ic,'freeSurface') )

  %   gf{cur}.u(:,:)=0.;
  %   gf{cur}.v(:,:)=0.;

    

   elseif ( strcmp(par.ic, 'perturbedPoiseuille') )
		
	ampu = (par.pressureInflowValue - par.pOutflow) / ( 2 * par.mu );

	epsilon = par.perturbation;

	uPert =   cos( par.kx * par.x(:,:,1) ) .* sin( par.kx * par.x(:,:,2)) ;
	vPert = - sin( par.kx * par.x(:,:,1) ) .* cos( par.kx * par.x(:,:,2)) ;

	% uPert = ( (cos( par.kx*par.x(:,:,1))).^3 + cos( par.kx*( par.x(:,:,1) + par.x(:, :, 2))) ) .* sin( par.ky*(par.x(:,:,2)) ) ;
	% vPert = sin( par.kx*(par.x(:,:,1)) ) .* sin( par.kx*par.x(:,:,2)) ;

	  gf{cur}.u(:,:)= ampu * ( 1 - par.x(:,:,2) ) .* par.x(:,:,2) + epsilon * uPert;
      gf{cur}.v(:,:)= epsilon * vPert;
  
  else
  	fprintf('getInitialConditions: unknown initial condition ic=[%s]\n',par.ic);
  	pause;
  end

  % if( strcmp(par.map,'freeSurface') )  

  %   % -- initialize free surface ---
  %   i1a = par.gid(1,1);
  %   i2a = par.gid(1,2);
  %   [I1,I2]=getIndex(par.dim);
  %   rv = zeros(par.Ngx,1);
  %   for i1=I1
  %     rv(i1) = (i1-i1a)*par.dr(1);
  %   end 

  %   omegat=pi;
  %   ampt=         0.1*cos(omegat*t);
  %   amptt=-omegat*0.1*sin(omegat*t);

  %   yv  =  ampt*sin(par.kx*rv); %% HARD CODE FOR NOW 


  % end 



return
end
