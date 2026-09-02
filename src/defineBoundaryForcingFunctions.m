  function par = defineBoundaryForcingFunctions( par )

  % For testing, choose BC's and IC to match the true solution

  if( strcmp(par.ic,'shear') )
    % shear flow 
    par.guax = @(x,y,t)  0.;  % BC RHS at x=ax
    par.gubx = @(x,y,t)  0;   % BC RHS at x=bx
    par.guay = @(x,y,t) -1.;  % BC RHS at y=ay
    par.guby = @(x,y,t) +1.;  % BC RHS at y=by

    par.gvax = @(x,y,t) 0.;  % BC RHS at x=ax
    par.gvbx = @(x,y,t) 0.;  % BC RHS at x=bx
    par.gvay = @(x,y,t) 0.;  % BC RHS at y=ay
    par.gvby = @(x,y,t) 0.;  % BC RHS at y=by

     par.pgax = @(x,y,t) 0.; 
     par.pgbx = @(x,y,t) 0.; 
     par.pgay = @(x,y,t) 0.; 
     par.pgby = @(x,y,t) 0.; 

  elseif( ~strcmp(par.ms,'none')            || ...
          ~strcmp(par.knownSolution,'none')        )
    par.guax = @(x,y,t) par.ue(x,y,t);  % BC RHS at x=ax
    par.gubx = @(x,y,t) par.ue(x,y,t);  % BC RHS at x=bx
    par.guay = @(x,y,t) par.ue(x,y,t);  % BC RHS at y=ay
    par.guby = @(x,y,t) par.ue(x,y,t);  % BC RHS at y=by
    par.u0   = @(x,y) par.ue(x,y,0);    % IC function

    par.gvax = @(x,y,t) par.ve(x,y,t);  % BC RHS at x=ax
    par.gvbx = @(x,y,t) par.ve(x,y,t);  % BC RHS at x=bx
    par.gvay = @(x,y,t) par.ve(x,y,t);  % BC RHS at y=ay
    par.gvby = @(x,y,t) par.ve(x,y,t);  % BC RHS at y=by
    par.v0   = @(x,y) par.ve(x,y,0);    % IC function

    % -- Assign function handles for pressure BC's ---
    if( par.bc(1,1)==par.dirichlet || par.bc(1,1)==par.periodic || par.bc(1,1)==par.pressureInflow || par.bc(1,1)==par.traction )
      par.pgax = @(x,y,t) par.pe(x,y,t);   % p: Dirichlet BC RHS at x=ax
    elseif( par.bc(1,1)==par.noSlipWall || par.bc(1,1)==par.inflow || par.bc(1,1)==par.slipWall )
      par.pgax = @(x,y,t) -par.pex(x,y,t);  % p: Neumann BC RHS at x=ax
    else
      fprintf('Define BC forcing: finish me...\n'); pause; 
    end 

    if( par.bc(2,1)==par.dirichlet || par.bc(2,1)==par.periodic || par.bc(2,1)==par.pressureInflow || par.bc(2,1)==par.traction  )
      par.pgbx = @(x,y,t) par.pe(x,y,t);   % p: Dirichlet BC RHS at x=bx
    elseif( par.bc(2,1)==par.noSlipWall || par.bc(2,1)==par.inflow || par.bc(2,1)==par.slipWall )
      par.pgbx = @(x,y,t) par.pex(x,y,t);   % p: Neumann BC RHS at x=bx
    elseif( par.bc(2,1)==par.outflow )
      par.pgbx = @(x,y,t) par.outflowPressureCoeffp*par.pe(x,y,t) + par.outflowPressureCoeffpn*par.pex(x,y,t);   % p: mixed BC 
    else
      fprintf('Define BC forcing: finish me...\n'); pause; 
    end 

    if( par.bc(1,2)==par.dirichlet || par.bc(1,2)==par.periodic || par.bc(1,2)==par.pressureInflow || par.bc(1,2)==par.traction  )
      par.pgay = @(x,y,t) par.pe(x,y,t);   % p: Dirichlet BC RHS at y=ay
    elseif( par.bc(1,2)==par.noSlipWall || par.bc(1,2)==par.inflow || par.bc(1,2)==par.slipWall)
      par.pgay = @(x,y,t) -par.pey(x,y,t);  % p: Neumann BC RHS at y=ay
    else
      fprintf('Define BC forcing: finish me...\n'); pause; 
    end 

    if( par.bc(2,2)==par.dirichlet || par.bc(2,2)==par.periodic || par.bc(2,2)==par.pressureInflow || par.bc(2,2)==par.traction )
      par.pgby = @(x,y,t) par.pe(x,y,t);   % p: Dirichlet BC RHS at y=by
    elseif( par.bc(2,2)==par.noSlipWall || par.bc(2,2)==par.inflow || par.bc(2,2)==par.slipWall)
      par.pgby = @(x,y,t) par.pey(x,y,t);   % p: Neumann BC RHS at y=by
    else
      fprintf('Define BC forcing: finish me...\n'); pause; 
    end 

  else 
    % default case
    par.guax = @(x,y,t) par.uInflow;  % BC RHS at x=ax
    par.gubx = @(x,y,t) 0.;  % BC RHS at x=bx
    par.guay = @(x,y,t) 0.;  % BC RHS at y=ay
    par.guby = @(x,y,t) 0.;  % BC RHS at y=by
    % par.u0   = @(x,y) 0.;    % IC function

    par.gvax = @(x,y,t) 0.;  % BC RHS at x=ax
    par.gvbx = @(x,y,t) 0.;  % BC RHS at x=bx
    par.gvay = @(x,y,t) 0.;  % BC RHS at y=ay
    par.gvby = @(x,y,t) 0.;  % BC RHS at y=by
    % par.v0   = @(x,y) 0.;    % IC function   

    par.pgax = @(x,y,t) 0.;   % p: Dirichlet BC RHS at x=ax
    par.pgbx = @(x,y,t) 0.;   % p: Dirichlet BC RHS at x=bx
    par.pgay = @(x,y,t) 0.;   % p: Dirichlet BC RHS at y=ay
    par.pgby = @(x,y,t) 0.;   % p: Dirichlet BC RHS at y=by

  end
  % pgbx = @(x,y,t) pe2(x,y,t);  % BC RHS at x=bx
  % pgay = @(x,y,t) pe2(x,y,t);  % BC RHS at y=ay
  % pgby = @(x,y,t) pe2(x,y,t);  % BC RHS at y=by

  % create a cell array of function handles for the BC functions
  par.gu = { par.guax, par.gubx, par.guay, par.guby }; 
  par.gv = { par.gvax, par.gvbx, par.gvay, par.gvby }; 
  par.gp = { par.pgax, par.pgbx, par.pgay, par.pgby }; 
