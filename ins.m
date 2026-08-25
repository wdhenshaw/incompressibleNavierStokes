%
%  --- Solve the Incompressible Navier-Stokes Equations ---
%       Velocity-Pressure (Pressure-Poisson) Formulation
%
%    u_t  + u*u.x + v*u.y + p.x = nu*( u.xx + u.yy ) + fu(x,y,t) ,    x in [xa,xb], y in [ya,yb]
%    v_t  + u*v.x + v*v.y + p.y = nu*( v.xx + v.yy ) + fv(x,y,t) ,    0 < t <= tFinal
%    p.xx + p.yy = -( u.x*u.x + 2*v.x*u.y + v.y*v.y ) + fp(x,y,t)
% 
%
% For usage run ins with no arguments (or see fprintf statements below)
%
function outPar = ins(varargin)

  % clear; clf; 

  clearvars -except varargin;
   % --- Clear all open figures ----
  clearOpenFigures(1:7);

  addpath(genpath(pwd)); % allow matlab to find files in subfolders

  % % ------ Restrict the number of threads so that cputime gives accurate answers ----------
  % maxThreads=maxNumCompThreads;
  % numThreads=1;
  % maxNumCompThreads(par.numThreads);   

  % clear; clf;  
  fontSize=16; lineWidth=2; markerSize=8; 
  set(0,'DefaultLineMarkerSize',markerSize);
  set(0,'DefaultLineLineWidth',lineWidth);
  set(0,'DefaultAxesFontSize',fontSize);

  % load the rainbow colour table
  rainbow;


  par.periodic=-1; par.dirichlet=1; par.noSlipWall=2; par.slipWall=3; par.inflow=4; par.outflow=5; par.pressureInflow=6;

  par.uc=1; par.vc=2; par.pc=3;    % component numbers

  par.numThreads=1;                % max number of threads Matlab is allowed to use 
  par.nd=2;                        % number of space dimensions
  par.tf=.5;                       % final time 
  par.tp=.1;                       % times to plot
  par.movieMode=0;                 % 1=run movie
  par.dtMax = 1e10;                % max dt (usually for implicit time stepping)
  par.cfl=.9; 
  par.ts = 'ab2';                  % time-stepping 
  par.orderInSpace=2; 
  par.checkTimeStep = 10000;       % check the time-step every this many steps
  par.savePlots = 0;               % 1 = save plots
  par.figDir    = 'fig';           % figure directory
  par.plotName  = 'ins';           % for plot name 

  par.map = 'Cartesian';           % 'Cartesian', 'Rectangle', 'Annulus', 'TFI', ...


  par.checkFileName = 'ins.check'; % name of the check file

  par.useOptFill=1;                % use optimized fill method for matrices
  par.useNew    =1;                % use new re-organized functions

  par.xa=0.; par.xb=1.;            % spatial dimensions
  par.ya=0.; par.yb=1.;            % spatial dimensions

  par.orderInSpace=2; 
  par.idebug=0;                    % set to 1 for debugging 
  par.plotOption=1;                % set to 1 for plotting
  par.computeErrors=0;           
  par.plotErrors=1;  
  par.plotSolutionOnGhost=0;       % 1 = plot solution and errors on ghost points   
  par.plotEveryStep=0;             % 1 = plot very step for debugging
  par.plotGrid=0;                  % 1 = plot grid
  par.nu=.1;                       % coefficient of diffusion
  par.rho=1;                       % FIX ME for rho .ne. 1
  par.kx=1.;                       % x-wave number in the IC and exact solution (scaled by 2 pi below)
  par.ky=1.;                       % y-wave number in the IC and exact solution (scaled by 2 pi below)
  par.kt=.5;                       % for TZ manufactured solution (scaled by 2 pi below)
  par.degreex = 2;
  par.degreet = 2;
  par.tzScale=0;                   % Trig TZ: 0=scale yb 1/sqrt(kx^2+ky^2) , 1=scale TZ yb 1/(kx^2+ky^2) 
  par.xa=0.; par.xb=1.;            % space interval interval
  par.ya=0.; par.yb=1.;            % space interval interval
  par.knownSolution='none';        % known solution, if any 
  par.ms = 'none';                 % manufactured solution, trig or poly
  par.degreex=2; degreet=2;        % degree of poly MS

  par.N0=10;                       % grid points in x and y if Nx and Ny are not set
  par.Nx = -1; 
  par.Ny = -1;

  par.shade='faceted'; 
  par.cdv=1.;                      % coefficient of divergence damping 
  par.ms = 'none';  
  par.bcs='nnnn';

  par.ic = 'default';              % initial condition
  par.shearBeta = 40;              % parameter in shear flow IC u = tanh(beta*(y-ym))
  par.shearDeltav=1e-2;            % amplitude of perturbation in v for shear flow IV

  par.uic = 1;                     % constant initial condition values
  par.vic = 0; 


  par.ad   = 0;                  % set to 1 to turn on artificial dissipation
  par.ad21 = .1;                 % coeff of linear AD
  par.ad22 = .1;                 % coeff of non-linear AD

  % Outflow BC for p is a0*p + a1*p.n = a0*pOutflow 
  par.outflowPressureCoeffp =1;  % a0 
  par.outflowPressureCoeffpn=1;  % a1
  par.pOutflow              =0;   

  par.pressureInflowValue   =1; % value for pressInflow BC 

  par.uInflow=1;  

  % Annulus map:
  par.x0 = 0;  % centre
  par.y0 = 0;
  par.startAngle  =0.; % angle variable on [0,1]
  par.endAngle    =1.; 
  par.innerRadius =0.5;
  par.outerRadius =1.;


  par.echo = 0;

  % --- read command line args ---
  for i = 1 : nargin
    line = varargin{i};

    % Read command line arguments for any entry in the "par" class
    par = assignCommandLineOption( line, par, par.echo );

  end

  if( nargin==0 )
    fprintf('Usage\n');
    fprintf(' ins -nu=<f> -tf=<f> -tp=<f> -ts=[ab2|pc2|im2] -cfl=<f> -bcs=<s> -ms=[none|poly|trig] -ic=[default|constant] ...\n');
    fprintf('     -idebug=<i> -known=[none|Poiseuille|TaylorGreen] -plotOption=<i> -kx=<f> -ky=<f> -kt=<f} -degreex=<i> -degreet=<i> ...\n');
    fprintf('     -ad=<i> -ad21=<f> -ad22=<f> -movieMode=[0|1] -savePlots=[0|1] -checkFileName=<s> -map=[Cartesian|Rectangle|Annulus|TFI] \n')
    fprintf('where:')
    fprintf(' bcs : list of four letters, d=Dirichlet, n=no-slip wall, s=slip wall, i=inflow, I=pressure inflow, o=outflow\n');
    fprintf('     : examples -bcs=ions : left=i, right=o, bottom=n, top=s\n');
    fprintf(' degreex,degreet : degrees of polynomial manufactured solution\n');

    outpar=0;
    return
  end

  % ------ Restrict the number of threads so that cputime gives accurate answers ----------
  maxThreads=maxNumCompThreads;
  maxNumCompThreads(par.numThreads);  


  par.savePlotThisStep=0; % set to 1 when we want to save plots

  par.tFinal = par.tf;
  N0 = par.N0;
  nu = par.nu;
  % xa = par.xa; xb = par.xb;
  % ya = par.ya; yb = par.yb;

  % set to 1 if implicit solvers for the velocity component are different:
  par.multipleImplicitSolversNeeded=0; 
  numbc = strlength(par.bcs); 
  par.bcLabel=''; 
  for dir=1:par.nd
    for side=1:2
      m= side + 2*(dir-1); 
      if( m <= numbc )
        if( par.bcs(m)=='d' )
          par.bc(side,dir)=par.dirichlet;
          par.bcLabel=strcat(par.bcLabel,'d');
        elseif( par.bcs(m)=='n' )
          par.bc(side,dir)=par.noSlipWall;
          par.bcLabel=strcat(par.bcLabel,'n');
        elseif( par.bcs(m)=='s' )
          par.bc(side,dir)=par.slipWall;
          par.bcLabel=strcat(par.bcLabel,'s');
          par.multipleImplicitSolversNeeded=1; 
        elseif( par.bcs(m)=='I' )
          par.bc(side,dir)=par.pressureInflow;
          par.bcLabel=strcat(par.bcLabel,'I');   
          par.multipleImplicitSolversNeeded=1;                     
        elseif( par.bcs(m)=='i' )
          par.bc(side,dir)=par.inflow;
          par.bcLabel=strcat(par.bcLabel,'i');
        elseif( par.bcs(m)=='o' )
          par.bc(side,dir)=par.outflow;
          par.bcLabel=strcat(par.bcLabel,'o');                     
        elseif( par.bcs(m)=='p' )
          par.bc(side,dir)=par.periodic;
          par.bcLabel=strcat(par.bcLabel,'p');
        else
          fprintf('ERROR: Unknown bcn=[%s]\n',par.bcs);
          error('error');
        end
      end
    end
  end  

  if( ~strcmp(par.ms,'none') && ~strcmp(par.map,'Cartesian') )
    % For testing : Turn off divergence damping for Manufactured solutions and non-Cartesian grids 
    par.cdv=0;
  end


  par.kx = par.kx*2*pi;
  par.ky = par.ky*2*pi;
  par.kt = par.kt*2*pi;

  % timings:
  par.cpuSetup           = 0;
  par.cpuGetUt           = 0;
  par.cpuPressure        = 0; % pressure solves (not factor)
  par.cpuBC              = 0;
  par.cpuImplicit        = 0; % implicit time-stepping
  par.cpuFactorImpMatrix = 0; % time to factor the implicit matrices
  par.cpuPlot            = 0;

  % --- define the manufactured solution and forcing functions ---
  par = defineManufacturedSolution( par );

  % --- define any known solution ---
  par = defineKnownSolution( par );

  % Define guax, gubx, ...
  par = defineBoundaryForcingFunctions( par );

  outPar.maxErr = zeros(4,1);

  %  --- Setup the grid ---
  par = setupGrid( par );
  
  if( par.idebug> 0 )
    fprintf('-----------------Incompressible Navier-Stokes ---------------------------------------------\n');
    fprintf(' ts=%s, tFinal=%g, nu=%g, cfl=%g, cdv=%g, knownSolution=%s, N0=%d, idebug=%d numThreads=%d, maxThreads=%d\n',...
             par.ts,par.tFinal,par.nu,par.cfl,par.cdv,par.knownSolution,N0,par.idebug,par.numThreads,maxThreads);
    fprintf(' map=%s\n',par.map);
    fprintf(' manufactured solution ms=%s, degreex=%d, degreet=%d, [kx,ky,kt]=[%g,%g,%g]*2*pi, tzScale=%d\n',...
             par.ms,par.degreex,par.degreet,par.kx/(2*pi),par.ky/(2*pi),par.kt/(2*pi),par.tzScale);
    fprintf(' par.bcLabel=%s, par.bc=[%d,%d,%d,%d] par.gid=[%d,%d,%d,%d]\n',par.bcLabel,par.bc(1,1),par.bc(2,1),par.bc(1,2),par.bc(2,2),...
             par.gid(1,1),par.gid(2,1),par.gid(1,2),par.gid(2,2) );
    fprintf(' useOptFill = %d (use optimized fill method for matrices)\n',par.useOptFill)
    fprintf(' ad=%d, ad21=%g, ad22=%g (ad=1 : add artificial dissipation)\n',par.ad,par.ad21,par.ad22);
    % fprintf(' useNew=%d (1=use new re-organized functions)\n',par.useNew);
    fprintf(' plotEveryStep=%d, plotOption=%d\n',par.plotEveryStep,par.plotOption);
    fprintf(' outflow: outflowPressureCoeffp=%g, outflowPressureCoeffpn=%g, pOutflow=%g\n',...
              par.outflowPressureCoeffp,par.outflowPressureCoeffpn,par.pOutflow);
    fprintf(' ic=%s\n',par.ic);
    fprintf('-------------------------------------------------------------------------------------------\n');

  end 


  % allocate space for the solution 
  Ngx = par.Ngx; Ngy=par.Ngy;
  un   = zeros(Ngx,Ngy);   % holds U_i^n
  vn   = zeros(Ngx,Ngy);   % holds V_i^n
  pn   = zeros(Ngx,Ngy);   % holds P_i^n

  unp1 = zeros(Ngx,Ngy);   % holds U_i^n+1
  vnp1 = zeros(Ngx,Ngy);   % holds V_i^n+1
  pnp1 = zeros(Ngx,Ngy);   % holds P_i^n+1
  
  ut   = zeros(Ngx,Ngy);   % holds u.t
  vt   = zeros(Ngx,Ngy);   % holds v.t

  par.cpuTotal = cputime; % start total timing here

  % fprintf('ims: START un=[%d,%d] vn=[%d,%d]\n',...
  %        size(un,1),size(un,2), size(vn,1),size(vn,2) );


  par.cpuSetup = par.cpuTotal;
  % --- Initial conditions ---
  t=0.; 
  [un,vn,par] = getInitialConditions( t,un,vn,par );


   % fprintf('AFTER getIC un=[%d,%d] vn=[%d,%d]\n',...
   %       size(un,1),size(un,2), size(vn,1),size(vn,2) );

  par.dtOld = -1;
  [dt,par] = getTimeStep( 0,un,vn,par );
  Nt = max(2,ceil(par.tFinal/dt));      % estimated number of time-steps 
  dt = par.tFinal/Nt;                   % adjust dt to reach tFinal exactly

  % coefficients in the Adams-Bashhforth scheme: 
  dtOld = dt;
  par.ab1=  dt*(1.+dt/(2.*dtOld));  % becomes 1.5*dt  if dt==dtOld
  par.ab2= -dt*    dt/(2.*dtOld);   %         -.5*dt      


  % fprintf('>>> nu*dt/dx^2 = %10.2e',nu*dt/dx^2);  
  % dt = cfl*.5*sqrt( dx^2 + dy^2);  % time step (adjusted below)  

  factorMatrix=1; 
  [pn,par] = pressureEquation( t,un,vn,dt,factorMatrix,par ); 
  factorMatrix=0; 

  nuScaleFactor=1.; % scale factor of nu*Delta( ) in getUt


  if( par.plotOption>2 )
    % ---- plot initial conditions ----
    par.step=1;
    par = plotSolution( t,un,vn,pn, par);
    fprintf('Plot initial condition and pause...\n');
    pause 
  end
  
  if( strcmp(par.ts,'ab2') || strcmp(par.ts,'pc2') )
    % -- evaluate du/dt at t=-dt for some schemes ---
    t=-dt;
    % do this for now: 
    ut = par.uet(par.x(:,:,1),par.x(:,:,2),t);
    vt = par.vet(par.x(:,:,1),par.x(:,:,2),t);
    % unp1 = ue(x,y,t);
    % vnp1 = ve(x,y,t);
    % pnp1 = pressureEquation( t,un,vn,dt,factorMatrix ); 
    % [ut,vt] = getUt( t,unp1,vnp1,pnp1,I1,I2,nuScaleFactor );
  end
  if( strcmp(par.ts,'im2') )
    % implicit scheme: ut holds explicit part of the operator
    t=-dt; 
    unp1 = par.ue(par.x(:,:,1),par.x(:,:,2),t);
    vnp1 = par.ve(par.x(:,:,1),par.x(:,:,2),t);
    pnp1 = par.pe(par.x(:,:,1),par.x(:,:,2),t); 
    nuScaleFactor=0.; % leave off viscous terms
    [ut,vt,par] = getUt( t,unp1,vnp1,pnp1,nuScaleFactor,par );
    nuScaleFactor=1;
    
  end

  t=0.;     

  % --- Start time-stepping loop ---

  if( strcmp(par.ts,'im2') )
    % ----- Form the implicit matrix ------
    par = formImplicitTimeSteppingMatrix( dt, par );

  end 

  par.cpuSetup = cputime - par.cpuSetup;
  par.cpuGetUt    = 0;
  par.cpuPressure = 0; % pressure solves (not factor)
  par.cpuBC       = 0;
  par.cpuImplicit = 0; % implicit time-stepping

  % --- Difference Operators ---
  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y

  DpxDmx = @(u,I1,I2) ( u(I1+1,I2) -2.*u(I1,I2) +u(I1-1,I2) )*(1./dx^2);  % u.xx 
  DpyDmy = @(u,I1,I2) ( u(I1,I2+1) -2.*u(I1,I2) +u(I1,I2-1) )*(1./dy^2);  % u.yy

  DzxDzy = @(u,I1,I2) ( u(I1+1,I2+1) - u(I1-1,I2+1) - u(I1+1,I2-1) +u(I1-1,I2-1) )*(1./(4.*dx*dy)); % u.xy   
   
  [I1,I2] = getIndex( par.gid );

  nextTimeToPlot = par.tp;

  maxNumberOfSteps = 1000000;
  tnp1=0; 
  for( n=1:maxNumberOfSteps )
  
    t    = tnp1;      % current time 
    tnp1 = t+dt;      % new time

    par.step = n;     % for titles
    
    if( strcmp(par.ts,'ab2') )
      % --- Adams-Bashforth order 2 ---

      par.nuScaleFactor=nuScaleFactor;
      [unp1,vnp1,pnp1,par,ut,vt] = advanceAdams( t,dt, un,vn,pn,unp1,vnp1,pnp1, ut,vt,par );


      % copy over for next step 
      un=unp1;
      vn=vnp1; 
      pn=pnp1; 

    elseif( strcmp(par.ts,'pc2') )
      % ---PC2 = AB2 + AM2 predictor corrector order 2 ---


      par.nuScaleFactor=nuScaleFactor;
      [unp1,vnp1,pnp1,par,ut,vt] = advancePC( t,dt, un,vn,pn,unp1,vnp1,pnp1, ut,vt,par );

      % copy over for next step 
      un=unp1;
      vn=vnp1; 
      pn=pnp1; 

    elseif( strcmp(par.ts,'im2') )

      % -- IMEX Implicit-Explicit scheme ---
      par.nuScaleFactor=nuScaleFactor;
      [unp1,vnp1,pnp1,par,ut,vt] = advanceIM( t,dt, un,vn,pn,unp1,vnp1,pnp1, ut,vt,par );

      % copy over for next step 
      un=unp1;
      vn=vnp1; 
      pn=pnp1;

    else
      fprintf('ERROR: unknown ts=%s\n',par.ts);
      pause;
      pause; 
    end

    
    if(  par.plotEveryStep || tnp1 >= nextTimeToPlot-.5*dt ) 
      nextTimeToPlot = nextTimeToPlot + par.tp;

      errorsComputed=0; 
      if( par.plotEveryStep || par.movieMode || (par.plotOption>0 && mod(floor(par.plotOption/2),2)==1) )
        par.savePlotThisStep=1;

        par = plotSolution( tnp1,un,vn,pn, par);

        if( par.movieMode==0 )   pause; else drawnow; end 
        errorsComputed=1;
      end
      cpuCurrent = cputime-par.cpuTotal;
      if( par.idebug ) 
        [maxDivU,maxGradU,par] = getMaxDivergence( un,vn,par );
        divOverGrad = maxDivU/max(maxGradU,1e-10);

        if( par.computeErrors  )
          if( errorsComputed )
            % These next were computed in plotSolution:
            pErrMax = par.maxErr(1);
            uErrMax = par.maxErr(2);
            vErrMax = par.maxErr(3);
            divMax  = par.maxErr(4);
          else
            [maxErr,perr,uerr,verr,div] = getErrors( tnp1,un,vn,pn,par );
            pErrMax = maxErr(1);
            uErrMax = maxErr(2);
            vErrMax = maxErr(3);
            divMax  = maxErr(4);
          end
          fprintf('%s: t=%9.3e step=%6d dt=%9.3e err-[u,v,p]=[%8.2e,%8.2e,%8.2e] div/grad=%9.2e, cpu=%9.2e(s)\n',...
               par.ts,tnp1,par.step,dt,uErrMax,vErrMax,pErrMax, divOverGrad,cpuCurrent); 
        else
          fprintf('%s: t=%9.3e step=%6d dt=%9.3e div/grad=%9.2e, cpu=%9.2e(s)\n',par.ts,tnp1,par.step,dt,divOverGrad,cpuCurrent); 
        end
      end

    end 

    if( tnp1 > par.tFinal-.5*dt )
      break;
    end

    if( mod(n,par.checkTimeStep)==0 )
      [dtNew,par] = getTimeStep( n,un,vn,par );

      numSteps = max(1,ceil((par.tFinal-tnp1)/dtNew));     % estimated number of time-steps 
      dtNew = (par.tFinal-tnp1)/numSteps;                     % adjust dt to reach tFinal exactly

      dtDiff = abs(dtNew-dt)/dt; % relative change 
      if( dtDiff > 0.1 )

        fprintf('Change the time-step: step=%6d, dt=%10.4e, dtNew=%10.4e, numStepsRemaining=%d (rel-dtdiff=%9.2e)\n',n,dt,dtNew,numSteps,dtDiff);
        dtOld = dt;
        dt = dtNew;
        par.ab1=  dt*(1.+dt/(2.*dtOld));  % becomes 1.5*dt  if dt==dtOld
        par.ab2= -dt*    dt/(2.*dtOld);   %         -.5*dt        

        if( strcmp(par.ts,'im2') )
          fprintf('time-step has changed, refactor the implicit matrix...\n');
          par = formImplicitTimeSteppingMatrix( dt, par );
        end
      else
        if( par.idebug > 1 )
          fprintf('Do NOT change the time-step: rel-dt-diff=%9.2e\n',dtDiff);
        end
      end 
    else
      par.ab1 = 1.5*dt;
      par.ab2 = -.5*dt;
    end   
  
  end
  % --- End time-stepping loop ---
  par.cpuTotal = cputime-par.cpuTotal;

  Nt = n; % actual number of time-steps

  if( abs(tnp1-par.tFinal)/par.tFinal > 1e-6*dt )
    fprintf('ERROR: tnp1=%16.10e is not equal to tFinal=%16.10e\n',tnp1,tFinal);
  end

  % plot results 
  if par.plotOption>1
    par.savePlotThisStep=1;
    par = plotSolution( tnp1,un,vn,pn, par);
  end

  if( par.computeErrors )
    maxErr = getErrors( tnp1,un,vn,pn,par );
    pErrMax = maxErr(1);
    uErrMax = maxErr(2);
    vErrMax = maxErr(3);
    divMax  = maxErr(4);
    if( par.plotOption>=-0 )
      fprintf('%s: t=%8.2e: Nx=%3d Ny=%3d Nt=%5d dt=%8.2e max-Err(p,u,v)=(%8.2e,%8.2e,%8.2e) cpu=%8.2e(s)\n',par.ts,tnp1,par.Nx,par.Ny,Nt,dt,pErrMax,uErrMax,vErrMax,par.cpuTotal);
    end
  else
    maxErr = zeros(4,1);
    if( par.plotOption>=-0 )
      fprintf('%s: t=%8.2e Nx=%3d Ny=%3d Nt=%5d dt=%8.2e cpu=%8.2e(s)\n',par.ts,tnp1,par.Nx,par.Ny,Nt,dt,par.cpuTotal);
    end
  end


  if( par.plotOption~=-1 )
    fprintf('  ------------- TIMINGS %s Nx=%3d Ny=%3d Nt=%5d --------------------\n',par.ts,par.Nx,par.Ny,Nt)
    fprintf('                  cpu (s)     %%   \n')
    fprintf('Total             %8.2e  %5.1f\n',par.cpuTotal,              par.cpuTotal/par.cpuTotal*100);
    fprintf('  setup           %8.2e  %5.1f\n',par.cpuSetup,              par.cpuSetup/par.cpuTotal*100);
    fprintf('  getUt           %8.2e  %5.1f\n',par.cpuGetUt,              par.cpuGetUt/par.cpuTotal*100);
    fprintf('  implicit solve  %8.2e  %5.1f\n',par.cpuImplicit,           par.cpuImplicit/par.cpuTotal*100);
    fprintf('  implicit factor %8.2e  %5.1f\n',par.cpuFactorImpMatrix,    par.cpuFactorImpMatrix/par.cpuTotal*100);
    fprintf('  pressure solve  %8.2e  %5.1f\n',par.cpuPressure,           par.cpuPressure/par.cpuTotal*100);
    fprintf('  bc              %8.2e  %5.1f\n',par.cpuBC,                 par.cpuBC/par.cpuTotal*100);
    fprintf('  plotting        %8.2e  %5.1f\n',par.cpuPlot,               par.cpuPlot/par.cpuTotal*100);
  end   

  % fprintf('%s: t=%8.2e CPU=%7.1e ',ts,tnp1,cpu);
  % if( m>1 )
  %   fprintf(' ratios(u,v,p,div)=(%4.2f,%4.2f,%4.2f,%4.2f) rates=(%4.2f,%4.2f,%4.2f,%4.2f)\n',...
  %                    erru(m-1)/erru(m),errv(m-1)/errv(m),errp(m-1)/errp(m),errDiv(m-1)/errDiv(m), ...
  %                    log2(erru(m-1)/erru(m)),log2(errv(m-1)/errv(m)),log2(errp(m-1)/errp(m)),log2(errDiv(m-1)/errDiv(m))); 
  % end
  
  par.maxErr = maxErr;
  par.Nt = Nt;
  par.dt = dt;
  par.uNorm(1) = max(abs(un), [], "all");
  par.uNorm(2) = max(abs(vn), [], "all");
  par.uNorm(3) = max(abs(pn), [], "all");
  [maxDivU,maxGradU,par] = getMaxDivergence( un,vn,par );
  par.maxDivU = maxDivU;
  par.maxGradU = maxGradU;
  writeCheckFile( par )

  outPar.maxErr = maxErr;
  outPar.cpu = par.cpuTotal;  
end


% --- Utility functions ---

% -----------------------------------------------------------------------
% convert a bc name (e.g. 'noSlipWall') to the corresponding integer flag
% -----------------------------------------------------------------------
function [ bcNumber ] = bcNameToNumber( bcName )

 globalDeclarations;

 if( strcmp(bcName,'periodic') )
   bcNumber=periodic;
 elseif( strcmp(bcName,'dirichlet') )
   bcNumber=dirichlet;
 elseif( strcmp(bcName,'noSlipWall') )
   bcNumber=noSlipWall;
 elseif( strcmp(bcName,'slipWall') )
   bcNumber=slipWall;
 elseif( strcmp(bcName,'inflow') )
   bcNumber=inflow;
 elseif( strcmp(bcName,'outflow') )
   bcNumber=outflow;
 else
  fprintf('bcNameToNumber:ERROR: unknown bcName=[%s]\n',bcName); pause; 
 end
end



% -----------------------------------------------------------------------
% Function to compute modulus with base 1 
% -----------------------------------------------------------------------
function [ mm ] = myMod( m,numberOfTimeLevels )
  mm = mod( m-1,numberOfTimeLevels)+1; 
end

% -----------------------------------------------------------------------
% Function getReal: read a command line argument for a real variable
% -----------------------------------------------------------------------
function [ val ] = getReal( line,name,val)
 % fprintf('getReal: val=%g line=[%s] name=[%s]\n',val,line,name);
 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%e',name)); 
   % fprintf('getReal: scan for val=%g\n',val);
 end
end

% -----------------------------------------------------------------------
% Function getInt: read a command line argument for an integer variable
% -----------------------------------------------------------------------
function [ val ] = getInt( line,name,val)
 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%d',name)); 
 end
end

% -----------------------------------------------------------------------
% Function getString: read a command line argument for a string variable
% -----------------------------------------------------------------------
function [ val ] = getString( line,name,val)
 if( strncmp(line,strcat(name,'='),length(name)+1) )
   val = sscanf(line,sprintf('%s=%%s',name)); 
 end
end
