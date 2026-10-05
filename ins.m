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
  addpath('src/'); addpath('matlabUtilities/'); % Add src files and matlabUtils

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

  % load the rainbow colour table & load to par
  % NOTE: Not assigned by user, just assigning variables
  rainbow;
  par.rainbowMap = rainbowMap;

  % ---- Boundary Condition Types ---
  par.periodic       =-1; 
  par.dirichlet      = 1; 
  par.noSlipWall     = 2; 
  par.slipWall       = 3; 
  par.inflow         = 4; 
  par.outflow        = 5; 
  par.pressureInflow = 6;
  par.traction       = 7;
 
  % ---- grid motion types ----
  par.noMotion          =0; 
  par.translate         =1; 
  par.rotate            =2;
  par.deform            =3;
  par.freeSurfaceMotion =4; 


  par.uc=1; par.vc=2; par.pc=3;    % component numbers
  par.nd=2;                        % number of space dimensions

  % WARN: To be implemented
  par.orderInSpace=2;
  
  % TODO: Fix for rho n not equal 1
  par.rho     = 1;
  
  % STEP: Define defaults for user-defined parameters

% Solution Description parameters

  % The three separate ways + par.bcs for defining an IBVP for the code to solve
  par.knownSolution='none';        % known solution, if any 
  par.ms = 'none';                 % manufactured solution, trig or poly
  par.ic = 'default';              % initial condition [default|constant|zero|shear]
  par.bcs='dddd';

  par.degreex=2; par.degreet=2;    % degree of poly MS % WARN: Maybe case specific?
  par.uic = 1;                     % constant initial condition values % WARN: Maybe case specific?
  par.vic = 0; 

  %BC Term parameters
  % Outflow BC for p is a0*p + a1*p.n = a0*pOutflow 
  par.outflowPressureCoeffp =1;  % a0 
  par.outflowPressureCoeffpn=1;  % a1
  par.pOutflow              =0;   
  par.pressureInflowValue   =1; % value for pressInflow BC 
  par.uInflow=1;  

  % Physical parameters
  par.nu      =.1;                 % coefficient of diffusion
  par.gravity = 0;                 % acceleration due to gravity is [0,par.gravity] NOTE: Negative!
  par.gamma   = 0;                 % coefficient of surface tension
  par.kx=1.;                       % x-wave number in the IC and exact solution (scaled by 2 pi below)
  par.ky=1.;                       % y-wave number in the IC and exact solution (scaled by 2 pi below)
  par.kt=.5;                       % for TZ manufactured solution (scaled by 2 pi below)
  par.tzScale=0;                   % Trig TZ: 0=scale yb 1/sqrt(kx^2+ky^2) , 1=scale TZ yb 1/(kx^2+ky^2) 
  par.perturbation=1e-3;           % Parameter for any eprturbation in the problem % WARN: Might be problem specific?

  % Time & Timestepping related parameters
  par.tf=.5;                       % final time 
  par.dtMax = 1e10;                % max dt (usually for implicit time stepping)
  par.cfl=.9; 
  par.ts = 'ab2';                  % time-stepping scheme
  par.checkTimeStep = 10000;       % check the time-step every this many steps

  % Grid / Geometry Parameters
  par.xa=0.; par.xb=1.;            % spatial dimensions
  par.ya=0.; par.yb=1.;            % spatial dimensions
  par.N0=10;                       % grid points in x and y if Nx and Ny are not set
  par.Nx = -1; 
  par.Ny = -1;
  par.map = 'Cartesian';           % 'Cartesian', 'Rectangle', 'Annulus', 'TFI', 'rotatedSquare', 'freeSurface', ...
  par.icfs  = 'sine';              % free surface initial condition ['sine','gaussian','cos']
  par.ampfs = 0.05;                % initial free surface amplitude

  % Motion Parameters
  par.motion='none';   % [none|translate|rotate|deform]
  par.transVect = [1,1]; % direction of the tranlate motion % WARN: Maybe case specific?
  par.numberOfTimeLevels=3; % ---- For moving grids ---
  par.numberOfGridFunctions=par.numberOfTimeLevels;  
  par.predictGrid=0; % par.gridMotionOption is passed to getGrid: 
  par.correctGrid=1;
  par.gridMotionOption=par.predictGrid;

  % Optimization & Computation parameters
  par.numThreads=1;                % max number of threads Matlab is allowed to use 
  par.useOptFill=1;                % use optimized fill method for matrices
  par.useNew    =1;                % use new re-organized functions
  par.combinedImplicitSolverNeeded =0; % set to 1 if implicit solve couples u and v

  % Artifical Terms Parameters
  par.cdv=1.;                      % coefficient of divergence damping 
  par.ad   = 0;                    % set to 1 to turn on artificial dissipation
  par.ad21 = .1;                   % coeff of linear AD
  par.ad22 = .1;                   % coeff of non-linear AD

  % Printing parameters
  par.idebug=0;                    % set to 1 for debugging TODO: Add more about other debugging types!!!
  par.computeErrors=0;           
  par.verbose=0; % Change to increase amount outputed                % TODO: Implement
  par.echo = 1;   % Parameter involved in assignCommandLineArguments
  par.checkFileName = 'ins.check'; % name of the check file

  % Plotting parameters
  par.plotOption=1;                % set to 1 for plotting
  par.movieMode=0;                 % 1=run movie 
  par.savePlots = 0;               % 1 = save plots
  par.tp=.1;                       % times to plot
  par.plotErrors=1;  
  par.plotSolutionOnGhost=0;       % 1 = plot solution and errors on ghost points   
  par.plotEveryStep=0;             % 1 = plot very step for debugging
  par.plotGrid=0;                  % 1 = plot grid
  par.plotVorticity=0;             % 1 = plot vorticity and streamlines
  par.plotAspectRatio=-1;          % plot aspect ratio
  par.shade='faceted'; 
  par.figDir    = 'fig';           % figure directory
  par.plotName  = 'ins';           % for plot name 

  % Case-specific parameters:
  % IMP: Need to migrate these to case specific things
  par.shearBeta = 40;              % parameter in shear flow IC u = tanh(beta*(y-ym))
  par.shearDeltav=1e-2;            % amplitude of perturbation in v for shear flow IV

  % Annulus map:
  par.x0 = 0;  % centre
  par.y0 = 0;
  par.startAngle  =0.; % angle variable on [0,1]
  par.endAngle    =1.; 
  par.innerRadius =0.5;
  par.outerRadius =1.;

  % Free surface gaussian
  par.betag = 10.;                 % gaussian parameter affecting how narrow the gaussian is 
  par.x0g   = 0.5;                 % centre for the gaussian

  % NOTE: Default par arguments end here!!!!

  % STEP: Prepare user defined arguments for the code

  % --- read command line args ---
  for i = 1 : nargin
    line = varargin{i};

    % Read command line arguments for any entry in the "par" class
    par = assignCommandLineOption( line, par, par.echo );

  end

  if( nargin==0 )
    fprintf('Usage\n');
    fprintf(' ins -nu=<f> -gamma=<f> -tf=<f> -tp=<f> -ts=[ab2|pc2|im2] -cfl=<f> -bcs=<s> -ms=[none|poly|trig]  -gravity-<f> ...\n');
    fprintf('     -idebug=<i> -known=[none|Poiseuille|TaylorGreen] -plotOption=<i> -kx=<f> -ky=<f> -kt=<f} -degreex=<i> -degreet=<i> ...\n');
    fprintf('     -ad=<i> -ad21=<f> -ad22=<f> -movieMode=[0|1] -savePlots=[0|1] -checkFileName=<s> \n');
    fprintf('     -map=[Cartesian|Rectangle|rotatedSquare|Annulus|TFI|freeSurface] \n')
    fprintf('     -motion=[none|translate|rotate|deform]\n');
    fprintf('     -ic=[default|constant|zero|shear]\n');
    fprintf('     -icfs=[sine|gaussian]\n');
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

  par.mu = par.nu*par.rho; 
  par.tFinal = par.tf;

  par.gravityVector = [0; par.gravity];

  N0 = par.N0;
  nu = par.nu;
  % xa = par.xa; xb = par.xb;
  % ya = par.ya; yb = par.yb;

  if( strcmp(par.motion,'none') )
    par.gridMotion=par.noMotion;
  elseif( strcmp(par.motion,'translate') )
    par.gridMotion=par.translate;
  elseif( strcmp(par.motion,'rotate') )
    par.gridMotion=par.rotate;
  elseif( strcmp(par.motion,'deform') )
    par.gridMotion=par.deform;    
  elseif( strcmp(par.motion,'freeSurfaceMotion') )
    par.gridMotion=par.freeSurfaceMotion;     
  else
    error('Unknown motion');
  end

  if( strcmp(par.map,'Cartesian') && par.gridMotion==par.noMotion)
    par.isCartesian =1;
  else
    par.isCartesian =0;
  end  

  % par.isCartesian =1; % -----************

  % multipleImplicitSolversNeeded : set to 1 if implicit solvers for the velocity component are different:
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
          if( par.isCartesian )
            par.multipleImplicitSolversNeeded=1; 
          else
            % slipWall BC couples u and v on the boundary and ghost 
            par.combinedImplicitSolverNeeded=1; 
          end
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
        elseif( par.bcs(m)=='t' )
          par.bc(side,dir)=par.traction;
          par.bcLabel=strcat(par.bcLabel,'t');  
          % traction BC couples u and v on the boundary and ghost 
          par.combinedImplicitSolverNeeded=1;     
        else
          fprintf('ERROR: Unknown bcn=[%s]\n',par.bcs);
          error('error');
        end
      end
    end
  end  
  if( par.combinedImplicitSolverNeeded ) par.multipleImplicitSolversNeeded=0; end 

  if( strcmp(par.map,'tfi') )
    par.map='TFI';
  end
  if( strcmp(par.map,'rectangle') )
    par.map='Rectangle';
  end  
  
  if( (~strcmp(par.ms,'none') && ~strcmp(par.map,'Cartesian')) ) % || strcmp(par.ms,'poly')  )
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

  % NOTE: If solution is following Manufactured / Known solutions, that is defined here!

  % --- define the manufactured solution and forcing functions ---
  par = defineManufacturedSolution( par );

  % --- define any known solution ---
  par = defineKnownSolution( par );

  % Define guax, gubx, ...
  % NOTE: Depending on solution type, define forcing functions (ms, known, ic)
  par = defineBoundaryForcingFunctions( par );

  outPar.maxErr = zeros(4,1);

  %  --- Setup the grid ---
  [par] = setupGrid( par );
  
  if( par.idebug> 0 )
    fprintf('-----------------Incompressible Navier-Stokes ---------------------------------------------\n');
    fprintf(' ts=%s, tFinal=%g, nu=%g, gamma=%g, gravity=[%g,%g] cfl=%g, cdv=%g, knownSolution=%s, N0=%d, idebug=%d numThreads=%d, maxThreads=%d\n',...
             par.ts,par.tFinal,par.nu,par.gamma,par.gravityVector(1),par.gravityVector(2),par.cfl,par.cdv,par.knownSolution,N0,par.idebug,par.numThreads,maxThreads);
    fprintf(' map=%s, isCartesian=%d, motion=%s (gridMotion=%d)\n',par.map,par.isCartesian,par.motion,par.gridMotion);
    fprintf(' manufactured solution ms=%s, degreex=%d, degreet=%d, [kx,ky,kt]=[%g,%g,%g]*2*pi, tzScale=%d\n',...
             par.ms,par.degreex,par.degreet,par.kx/(2*pi),par.ky/(2*pi),par.kt/(2*pi),par.tzScale);
    fprintf(' par.bcLabel=%s, par.bc=[%d,%d,%d,%d] par.gid=[%d,%d,%d,%d]\n',par.bcLabel,par.bc(1,1),par.bc(2,1),par.bc(1,2),par.bc(2,2),...
             par.gid(1,1),par.gid(2,1),par.gid(1,2),par.gid(2,2) );
    fprintf(' useOptFill = %d (use optimized fill method for matrices)\n',par.useOptFill)
    fprintf(' ad=%d, ad21=%g, ad22=%g (ad=1 : add artificial dissipation)\n',par.ad,par.ad21,par.ad22);
    % fprintf(' useGridFunctions=%d (1=use new grid functions)\n',par.useGridFunctions);
    fprintf(' plotEveryStep=%d, plotOption=%d\n',par.plotEveryStep,par.plotOption);
    fprintf(' outflow: outflowPressureCoeffp=%g, outflowPressureCoeffpn=%g, pOutflow=%g\n',...
              par.outflowPressureCoeffp,par.outflowPressureCoeffpn,par.pOutflow);
    fprintf(' implicitTimeStep: multipleImplicitSolversNeeded=%d combinedImplicitSolverNeeded=%d\n',par.multipleImplicitSolversNeeded,par.combinedImplicitSolverNeeded);
    fprintf(' ic=%s (initial condition), icfs=%s (free surface initial condition)\n',par.ic,par.icfs);
    fprintf(' freeSurface: ampfs=%g\n',par.ampfs);
    fprintf('-------------------------------------------------------------------------------------------\n');

  end 

  % STEP: Start solve

  % --- allocate space for the solution ---
  Ngx = par.Ngx; Ngy=par.Ngy;
  for igf=1:par.numberOfGridFunctions
    gf{igf}.u   = zeros(Ngx,Ngy);
    gf{igf}.v   = zeros(Ngx,Ngy);
    gf{igf}.p   = zeros(Ngx,Ngy);
    gf{igf}.eta = zeros(Ngx,par.nd); % free surface coordinates 
  end
  
  ut   = zeros(Ngx,Ngy);   % holds u.t
  vt   = zeros(Ngx,Ngy);   % holds v.t

  par.cpuTotal = cputime; % start total timing here

  % fprintf('ims: START un=[%d,%d] vn=[%d,%d]\n',...
  %        size(un,1),size(un,2), size(vn,1),size(vn,2) );


  par.cpuSetup = par.cpuTotal;
  % --- Initial conditions ---
  t=0.; 
  cur =1;  % cuurent solution index into gf{}
  next=mod(cur  +par.numberOfGridFunctions,par.numberOfGridFunctions)+1; % = 2 
  prev=mod(cur-2+par.numberOfGridFunctions,par.numberOfGridFunctions)+1; % = 3
  % -- get the grid and grid velocity ---
  [gf,par] = getGrid( t, gf,cur, par ); 

  % --- assign initial conditions ---  
  [gf,par] = getInitialConditions( t,gf,cur,par );

   % fprintf('AFTER getIC un=[%d,%d] vn=[%d,%d]\n',...
   %       size(un,1),size(un,2), size(vn,1),size(vn,2) );

  par.dtOld = -1;
  [dt,par] = getTimeStep( 0,gf{cur}.u,gf{cur}.v, gf,cur, par );
  Nt = max(2,ceil(par.tFinal/dt));      % estimated number of time-steps 
  dt = par.tFinal/Nt;                   % adjust dt to reach tFinal exactly
  par.dt=dt;

  % coefficients in the Adams-Bashhforth scheme: 
  dtOld = dt;
  par.ab1=  dt*(1.+dt/(2.*dtOld));  % becomes 1.5*dt  if dt==dtOld
  par.ab2= -dt*    dt/(2.*dtOld);   %         -.5*dt      


  % fprintf('>>> nu*dt/dx^2 = %10.2e',nu*dt/dx^2);  
  % dt = cfl*.5*sqrt( dx^2 + dy^2);  % time step (adjusted below)  

  par.factorPressureMatrix=1; 
  [gf{cur}.p,par] = pressureEquation( t,gf{cur}.u,gf{cur}.v,dt, gf,cur, par ); 

  nuScaleFactor=1.; % scale factor of nu*Delta( ) in getUt


  if( par.plotOption>2 )
    % ---- plot initial conditions ----
    par.step=1;
    par = plotSolution( t,gf{cur}.u,gf{cur}.v,gf{cur}.p, gf,cur, par);
    fprintf('Plot initial condition and pause...\n');
    pause 
  end
  
  if( strcmp(par.ts,'ab2') || strcmp(par.ts,'pc2') )
    % -- evaluate du/dt at t=-dt for some schemes ---

    % CHECK ME FOR MOVING GRIDS

    t=-dt;

    [gf,par] = getGrid( t, gf,prev, par ); 

    if( ~strcmp(par.ms,'none') || ~strcmp(par.knownSolution,'none') )
      gf{prev}.u = par.ue(gf{prev}.x(:,:,1),gf{prev}.x(:,:,2),t);
      gf{prev}.v = par.ve(gf{prev}.x(:,:,1),gf{prev}.x(:,:,2),t);
      gf{prev}.p = par.pe(gf{prev}.x(:,:,1),gf{prev}.x(:,:,2),t); 


      if( par.gridMotion==par.noMotion )
        ut = par.uet(gf{prev}.x(:,:,1),gf{prev}.x(:,:,2),t);
        vt = par.vet(gf{prev}.x(:,:,1),gf{prev}.x(:,:,2),t);
      else
        nuScaleFactor=1;
        [ut,vt,par] = getUt( t,gf{prev}.u,gf{prev}.v,gf{prev}.p, gf,prev, nuScaleFactor,par );
      end 
    else
      % do this for now:
      gf{prev}.u = gf{cur}.u;
      gf{prev}.v = gf{cur}.v;
      gf{prev}.p = gf{cur}.p;
      ut(:,:)=0;
      vt(:,:)=0;
    end
 
  end
  if( strcmp(par.ts,'im2') )
    % implicit scheme: ut holds explicit part of the operator
    t=-dt; 

    [gf,par] = getGrid( t, gf,prev, par );    

    if( ~strcmp(par.ms,'none') || ~strcmp(par.knownSolution,'none') )
      gf{prev}.u = par.ue(gf{prev}.x(:,:,1),gf{prev}.x(:,:,2),t);
      gf{prev}.v = par.ve(gf{prev}.x(:,:,1),gf{prev}.x(:,:,2),t);
      gf{prev}.p = par.pe(gf{prev}.x(:,:,1),gf{prev}.x(:,:,2),t); 

      % gf{prev}.u = par.ue(par.x(:,:,1),par.x(:,:,2),t);
      % gf{prev}.v = par.ve(par.x(:,:,1),par.x(:,:,2),t);
      % gf{prev}.p = par.pe(par.x(:,:,1),par.x(:,:,2),t); 


      nuScaleFactor=0.; % leave off viscous terms
      [ut,vt,par] = getUt( t,gf{prev}.u,gf{prev}.v,gf{prev}.p, gf,prev, nuScaleFactor,par );
      nuScaleFactor=1; 
    else
     % do this for now:
      gf{prev}.u = gf{cur}.u;
      gf{prev}.v = gf{cur}.v;
      gf{prev}.p = gf{cur}.p;     
      ut(:,:)=0;
      vt(:,:)=0;
    end                
    
  end

  t=0.;     

  % STEP: --- Start time-stepping loop ---
  par.factorImplicitMatrix=1;

  % if( strcmp(par.ts,'im2') )
  %   % ----- Form the implicit matrix ------
  %   par = formImplicitTimeSteppingMatrix( t,dt, gf,cur, par );

  % end 

  par.cpuSetup = cputime - par.cpuSetup;
  par.cpuGetUt    = 0;
  par.cpuPressure = 0; % pressure solves (not factor)
  par.cpuBC       = 0;
  par.cpuImplicit = 0; % implicit time-stepping

   
  [I1,I2] = getIndex( par.gid );

  nextTimeToPlot = par.tp;

  maxNumberOfSteps = 1000000;
  tnp1=0; 
  % -------------------------------------------------
  % ---------------- START TIME STEPS ---------------
  % -------------------------------------------------
  for( n=1:maxNumberOfSteps )
  
    t    = tnp1;      % current time 
    tnp1 = t+dt;      % new time

    par.step = n;     % for titles

    if( strcmp(par.ts,'ab2') )
      % --- Adams-Bashforth order 2 ---
      par.nuScaleFactor=nuScaleFactor;
      [gf,par,ut,vt] = advanceAdams( t,dt, gf,cur, ut,vt,par );

    elseif( strcmp(par.ts,'pc2') )
      % ---PC2 = AB2 + AM2 predictor corrector order 2 ---
      par.nuScaleFactor=nuScaleFactor;
      [gf,par,ut,vt] = advancePC( t,dt, gf,cur, ut,vt,par );

    elseif( strcmp(par.ts,'im2') )

      % -- IMEX Implicit-Explicit scheme ---
      par.nuScaleFactor=nuScaleFactor;
      [gf,par,ut,vt] = advanceIM( t,dt, gf,cur, ut,vt,par );

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

        par = plotSolution( tnp1,gf{next}.u,gf{next}.v,gf{next}.p, gf,next, par);

        if( par.movieMode==0 )   pause; else drawnow; end 
        errorsComputed=1;
      end
      cpuCurrent = cputime-par.cpuTotal;
      if( par.idebug ) 
        [maxDivU,maxGradU,par] = getMaxDivergence( gf{next}.u,gf{next}.v,par );
        divOverGrad = maxDivU/max(maxGradU,1e-10);

        if( par.computeErrors  )
          if( errorsComputed )
            % These next were computed in plotSolution:
            pErrMax = par.maxErr(1);
            uErrMax = par.maxErr(2);
            vErrMax = par.maxErr(3);
            divMax  = par.maxErr(4);
          else
            [maxErr,perr,uerr,verr,div] = getErrors( tnp1,gf{next}.u,gf{next}.v,gf{next}.p, gf,next, par );
            pErrMax = maxErr(1);
            uErrMax = maxErr(2);
            vErrMax = maxErr(3);
            divMax  = maxErr(4);
          end
          fprintf('%s: t=%9.3e step=%6d Nx=%3d dt=%9.3e err-[p,u,v]=[%8.2e,%8.2e,%8.2e] div/grad=%9.2e (grad=%8.2e), cpu=%9.2e(s)\n',...
               par.ts,tnp1,par.step,par.Nx,dt,pErrMax,uErrMax,vErrMax, divOverGrad,maxGradU,cpuCurrent); 
        else
          fprintf('%s: t=%9.3e step=%6d Nx=%3d dt=%9.3e div/grad=%9.2e (grad=%8.2e) cpu=%9.2e(s)\n',par.ts,tnp1,par.step,par.Nx,dt,divOverGrad,maxGradU,cpuCurrent); 
        end
      end

    end 

    % update grid function counters 
    prev = mod(prev,par.numberOfGridFunctions)+1;
    cur  = mod(cur ,par.numberOfGridFunctions)+1;
    next = mod(next,par.numberOfGridFunctions)+1;

    if( tnp1 > par.tFinal-.5*dt )
      break;
    end

    if( mod(n,par.checkTimeStep)==0 )
      [dtNew,par] = getTimeStep( n,gf{cur}.u,gf{cur}.v, gf,cur, par );

      numSteps = max(1,ceil((par.tFinal-tnp1)/dtNew));     % estimated number of time-steps 
      dtNew = (par.tFinal-tnp1)/numSteps;                     % adjust dt to reach tFinal exactly

      dtDiff = abs(dtNew-dt)/dt; % relative change 
      if( dtDiff > 0.1 )

        fprintf('Change the time-step: step=%6d, dt=%10.4e, dtNew=%10.4e, numStepsRemaining=%d (rel-dtdiff=%9.2e)\n',n,dt,dtNew,numSteps,dtDiff);
        dtOld = dt;
        dt = dtNew;
        par.dt=dt;
        % Adjust AB2 coefficients for a variable time step:
        par.ab1=  dt*(1.+dt/(2.*dtOld));  % becomes 1.5*dt  if dt==dtOld
        par.ab2= -dt*    dt/(2.*dtOld);   %         -.5*dt        

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
    par = plotSolution( tnp1,gf{cur}.u,gf{cur}.v,gf{cur}.p, gf,cur, par);
  end

  if( par.computeErrors )
    maxErr = getErrors( tnp1,gf{cur}.u,gf{cur}.v,gf{cur}.p, gf,cur, par );
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

  par.maxErr = maxErr;
  par.Nt = Nt;
  par.dt = dt;
  par.uNorm(1) = max(abs(gf{cur}.u), [], "all");
  par.uNorm(2) = max(abs(gf{cur}.v), [], "all");
  par.uNorm(3) = max(abs(gf{cur}.p), [], "all");
  [maxDivU,maxGradU,par] = getMaxDivergence( gf{cur}.u,gf{cur}.v, par );
  par.maxDivU = maxDivU;
  par.maxGradU = maxGradU;
  writeCheckFile( par )

  % return output parameters for runConvergence
  outPar.maxErr = maxErr;
  outPar.cpu    = par.cpuTotal; 
  outPar.par    = par; 
end



