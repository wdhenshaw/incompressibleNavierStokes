function [un,vn,par] = getInitialConditions( t,un,vn,par )

  if( strcmp(par.ic,'default') )
    
    un(:,:) = par.u0(par.x(:,:,1),par.x(:,:,2)); 
    vn(:,:) = par.v0(par.x(:,:,1),par.x(:,:,2)); 

  elseif( strcmp(par.ic,'zero') )

    un(:,:)=0.;
    vn(:,:)=0.;

  elseif( strcmp(par.ic,'constant') )
  	
   un(:,:) = par.uic;
   vn(:,:) = par.vic;

  elseif( strcmp(par.ic,'shear') )
    
   % shear flow 
   par.plotErrors=0;

   % un = zeros(par.Ngx,par.Ngy);
   % vn = zeros(par.Ngx,par.Ngy);

   beta = par.shearBeta;
   ym = 0.5*(par.ya + par.yb); % center in y 
   un(:,:) = tanh( beta*( par.x(:,:,2)-ym ));

   % perturbation for v 
   delta=par.shearDeltav;
   vn(:,:) = delta*cos( 2*pi*(par.x(:,:,2)-ym) ).*sin( par.kx*par.x(:,:,1));

  else
  	fprintf('getInitialCOnditions: unknown initial condition ic=[%s]\n',par.ic);
  	pause;
  end
return
end