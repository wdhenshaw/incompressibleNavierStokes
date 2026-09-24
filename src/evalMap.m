%
%  Evaluate the mapping grid points and metric derivatives
% 
%  par.x(i1,i2,1:2) (output) : grid points 
%  par.rx(i1,i2,1:2,1:2) (output) : metric derivatives dr_i/dx_j
%
function par = evalMap( par )


  r = zeros(par.Ngx,par.Ngy,2);

  par.dr(1) = 1./par.Nx;
  par.dr(2) = 1./par.Ny;

  [I1,I2]=getIndex(par.dim);

  i1a = par.gid(1,1);
  i2a = par.gid(1,2);

  for i2=I2
  for i1=I1
    r(i1,i2,1) = (i1-i1a)*par.dr(1);
    r(i1,i2,2) = (i2-i2a)*par.dr(2);
  end
  end

  par.x  = zeros(par.Ngx,par.Ngy,par.nd);
  par.rx = zeros(par.Ngx,par.Ngy,par.nd,par.nd);

  if( strcmp(par.map,'Rectangle' ) || strcmp(par.map,'Cartesian') )
    
    % for testing make a rectangle that pretends to vbe a more general curvilinear grid

    for i2=I2
    for i1=I1
      par.x(i1,i2,1)=par.xa + (i1-i1a)*par.dx; 
      par.x(i1,i2,2)=par.ya + (i2-i2a)*par.dy; 

      par.rx(i1,i2,1,1) = par.dr(1)/par.dx;
      par.rx(i1,i2,1,2) = 0.;
      par.rx(i1,i2,2,1) = 0.;
      par.rx(i1,i2,2,2) = par.dr(2)/par.dy;
    end
    end

  elseif( strcmp(par.map,'rotatedSquare' )  )
    % counter-clockwise rotation about the origin
    %     xv = R xv0
    %     R = [ cos(angle) -sin(angle) ]
    %         [ sin(angle)  cos(angle) ]    
    
    angle=pi/4;
    ca = cos(angle);
    sa = sin(angle);

    for i2=I2
    for i1=I1
      x0 = par.xa + r(i1,i2,1)*(par.xb-par.xa);
      y0 = par.ya + r(i1,i2,2)*(par.yb-par.ya);
      xr0 = (par.xb-par.xa);
      xs0 = 0;
      yr0 = 0;
      ys0 = (par.yb-par.ya);

      par.x(i1,i2,1) = ca*x0 - sa*y0;
      par.x(i1,i2,2) = sa*x0 + ca*y0;

      xr(1,1) = ca*xr0 - sa*yr0;
      xr(1,2) = ca*xs0 - sa*ys0;

      xr(2,1) = sa*xr0 + ca*yr0; 
      xr(2,2) = sa*xs0 + ca*ys0; 

      rx = inv(xr);

      par.rx(i1,i2,1,1) = rx(1,1);
      par.rx(i1,i2,1,2) = rx(1,2);
      par.rx(i1,i2,2,1) = rx(2,1);
      par.rx(i1,i2,2,2) = rx(2,2);
    end
    end

  elseif( strcmp(par.map,'Annulus' ) )


    scale= 2*pi*(par.endAngle-par.startAngle);
    rad  = par.outerRadius-par.innerRadius;

    for i2=I2
    for i1=I1
      r1 = r(i1,i2,1); % angle variable
      r2 = r(i1,i2,2); % radial variable

       angle = scale*r1 + 2*pi*par.startAngle;
       radius= rad*r2 + par.innerRadius;
       sina=sin(angle);
       cosa=cos(angle);
     
       par.x(i1,i2,1)=radius*cosa+par.x0;
       par.x(i1,i2,2)=radius*sina+par.y0;

       xr(1,1)=-radius*scale*sina;
       xr(2,1)= radius*scale*cosa;
       xr(1,2)= rad*cosa;
       xr(2,2)= rad*sina; 
   
       rx = inv(xr);

       for m2=1:par.nd
       for m1=1:par.nd
         par.rx(i1,i2,m1,m2) = rx(m1,m2);
       end
       end

    end
    end

  elseif( strcmp(par.map,'TFI' ) )

    % Transfinite interpolation mapping 

    % --- Hard code bottom and top curves for now ---

    % --- bottom curve ---
    ampb=0.05;
    curveb =  @(r) ampb*sin(2*pi*r);
    curvebr = @(r) (2*pi*ampb)*cos(2*pi*r);

    % --- top curve ---
    ampt=0.1;
    cty=0.5; 
    curvet =  @(r) cty + ampt*sin(2*pi*r);
    curvetr = @(r) (2*pi*ampt)*cos(2*pi*r);



    for i2=I2
    for i1=I1
      r1 = r(i1,i2,1); 
      r2 = r(i1,i2,2); 


       par.x(i1,i2,1)=r1; 
       par.x(i1,i2,2)=(1-r2)*curveb(r1) + r2*curvet(r1);

       xr(1,1)= 1;
       xr(2,1)= (1-r2)*curvebr(r1) + r2*curvetr(r1);
       xr(1,2)= 0; 
       xr(2,2)= -curveb(r1) + curvet(r1);
   
       rx = inv(xr);

       for m2=1:par.nd
       for m1=1:par.nd
         par.rx(i1,i2,m1,m2) = rx(m1,m2);
       end
       end

    end
    end

  elseif( strcmp(par.map,'freeSurface') )

    % freeSurface : TFI with top curve defined by data points 

    % --- bottom curve ---
    cby=-0.5;
    ampb=0.0;
    curveb =  @(r) cby + ampb*sin(2*pi*r);
    curvebr = @(r) (2*pi*ampb)*cos(2*pi*r);

    % --- top curve ---
    rv = r(:,1,1); 
    ampt= par.ampfs; % 0.05;
    xv = par.xa + (par.xb-par.xa)*rv;
    if( strcmp(par.icfs,'sine') )
      yv = ampt*sin(par.kx*rv);
    elseif( strcmp(par.icfs,'gaussian') )
      yv = ampt*exp( - (par.betag*(xv-par.x0g)).^2 );
	elseif (strcmp(par.icfs,'cos'))
	  yv = ampt*cos(par.kx*xv);
    else
      fprintf('\n getGrid:ERROR: unknown free surface initial condition, icsf=[%s]\n',par.icsf);
    end    
    % yv = ampt*sin(par.kx*rv);  %% HARD CODE FOR NOW -- MUST MATCH VALUES IN getGrid

    conds='second'; % second derivative zero ('natural')
    if( par.bc(1,1)==par.periodic )
      conds='periodic';
    end 
    par.cs = splineInit( rv,yv,conds );
    [yv0,yvr] = splineEval( par.cs, rv ); 


    for i2=I2
    for i1=I1
      r1 = r(i1,i2,1); 
      r2 = r(i1,i2,2); 


       par.x(i1,i2,1)=r1; 
       par.x(i1,i2,2)=(1-r2)*curveb(r1) + r2*yv(i1);

       xr(1,1)= 1;
       xr(2,1)= (1-r2)*curvebr(r1) + r2*yvr(i1);
       xr(1,2)= 0; 
       xr(2,2)= -curveb(r1) + yv(i1);
   
       rx = inv(xr);

       for m2=1:par.nd
       for m1=1:par.nd
         par.rx(i1,i2,m1,m2) = rx(m1,m2);
       end
       end

    end
    end



  else
    fprintf('evalMap: ERROR: unknown map=[%s]\n',par.map);
    error('evalMap:error');
  end


  if( par.plotGrid )

    % --- plot the grid points --

    figure(5)
    for i2=I2
      if( i2<par.gid(1,2) || i2>par.gid(2,2) )
        plot( par.x(I1,i2,1),par.x(I1,i2,2),'b--','LineWidth',1); % ghost line
      else
        plot( par.x(I1,i2,1),par.x(I1,i2,2),'b-'); 
      end
      hold on;
    end
    for i1=I1
      if( i1<par.gid(1,1) || i1>par.gid(2,1) )
        plot( par.x(i1,I2,1),par.x(i1,I2,2),'k--','LineWidth',1);
      else
        plot( par.x(i1,I2,1),par.x(i1,I2,2),'k-');
      end
    end

    title(sprintf('Grid: map=%s',par.map)); xlabel('x'); ylabel('y');

    % set plot bounds
    xMin = min(par.x(:,:,1), [], "all");
    xMax = max(par.x(:,:,1), [], "all");

    yMin = min(par.x(:,:,2), [], "all");
    yMax = max(par.x(:,:,2), [], "all");

    delta=.1; 
    xPlotMin = xMin- delta*(xMax-xMin);
    xPlotMax = xMax+ delta*(xMax-xMin);

    yPlotMin = yMin- delta*(yMax-yMin);
    yPlotMax = yMax+ delta*(yMax-yMin);

    xlim([xPlotMin,xPlotMax]); ylim([yPlotMin,yPlotMax]);
    grid on;
    axis equal;

    hold off;
    % pause

  end


  return
end
