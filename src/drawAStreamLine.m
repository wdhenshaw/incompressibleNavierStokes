
% =======================================================================
%         drawAStreamLine
%  
%  Integrate the streamline whose initial position is (xtp,ytp)
%  Mark the cells in the array par.maskForStreamLines(i,j) i=1,nxg j=1,nyg
%  which the streamline passes through and stop plotting the
%  streamline if the streamline enters a cell with par.maskForStreamLines(i,j) >= 2
%  Plot arrows on the streamline when the streamline enters a
%  cell (i,j) satisfying mod(i,lax)=0 and mod(j,lay)=0 where
%  (lax,lay) are assign below
%   Colour the contours according to the value of u**2+v**2.
% 
% This file started from Overuture/ogshow/streamLinesOpt.C
% =======================================================================
function [xp,yp,cp,par] = drawAStreamLine( xtp,ytp, u,v,par )



% void PlotIt:: 
% drawAStreamLine(GenericGraphicsInterface &gi, const GridCollection & gc, 
%                 const realGridCollectionFunction & uv, 
%                 int *componentsToInterpolate,
%                 IntegerArray & par.maskForStreamLines_,
%                 real arrowSize,
%                 GraphicsParameters & psp,
%                 real *uip, int *indexGuessp,
%                 real & xa, real &ya, real & xb, real &yb, real &xba, real &yba, 
%                 real &uMin, real &uMax, real &cfl, int &nrsmx,
%                 int & nxg, int & nyg, real &xtp, real &ytp, int & intopt )

  lax=5; lay=5;  % for plotting arrows

  xa = par.xa;
  xb = par.xb;
  ya = par.ya;
  yb = par.yb;

  xba=xb-xa;
  yba=yb-ya;

  uMin  = par.uMin;
  uMax  = par.uMax;
  nrsmx = max(par.Nx-1,par.Ny-1);

  epsU = par.streamLineStoppingTolerance*uMax;

  % for streamline mask:
  nxg = par.nxg;
  nyg = par.nyg;

  cfl = par.cflStreamLine;

  dxmx=max(xba,yba);

  uvfact=1./(uMax-uMin);     % normalization factor for colour table

  nt=nrsmx*2;   % maximum number of time steps

  t=0.;
  cfla=abs(cfl)*dxmx/nrsmx;    % *** could use dxmx and nrsmx from a particular grid?? ****
  dtmx=cfla*2.;
  dtMax = 100.*cfla/uvfact;
  
  xgf=nxg/xba; ygf=nyg/yba;

  x=xtp; y=ytp; 
  % xi=x; yi=y;

  % xvp = zeros(2,1);
  if( par.idebug>1 ) fprintf('Entering drawAStreamLine (xtp,ytp)=(%g,%g) uMin=%9.2e uMax=%9.2e epsU=%9.3e idebug=%d...\n',xtp,ytp,uMin,uMax,epsU,par.idebug); end


  if( par.idebug>1 ) fprintf("***Start a line x=%e, y=%e, nt=%d \n",x,y,nt); end


  ixg0=-1;              
  iyg0=-1; 
  ixg1=floor((x-xa)*xgf)+1; % starting cell
  iyg1=floor((y-ya)*ygf)+1;

  dpath=0.;
  
  uiv =xInterpolate( x,y,u,v,par );

  it=1; 
  % uValue=(sqrt( uiv(1)^2 + uiv(2)^2 )-uMin)*uvfact;
  uValue=sqrt( uiv(1)^2 + uiv(2)^2 );
  xp(it) = x; 
  yp(it) = y; 
  cp(it) = uValue; % uValue is in [0,1]  

  % ...........Take time steps it=2,3,...,nt  (t=dt,2dt,...,(nt-1)*dt)
  for( it=2:nt )
    if( par.idebug>1 ) fprintf('it=%d...\n',it); end

    xi=x;   yi=y;
    uiabs=abs(uiv(1))+abs(uiv(2));
    if( uiabs < epsU )  % flow is too slow here to move anywhere
      % if( it==2 ) printf("******** Exiting, flow too slow, it=%i *************\n",it);
      if( par.idebug>1 ) fprintf('******** Exiting, flow too slow, it=%i *************\n',it); end
      break
    end 

    dt= cfla/uiabs;
    if( cfl<0 ) dt = -dt; end 

    t = t+dt;
    xs=x+dt*uiv(1);
    ys=y+dt*uiv(2);

    if( xs<xa || xs>xb || ys<ya || ys>yb )
      if( par.idebug>1 ) fprintf("**Line leaves region xs=%e, ys=%e, dt=%e \n",xs,ys,dt); end
      break;   %      ....outside plotting bounds
    end

    uOld=uiv(1); vOld=uiv(2);
      
    % interpolate velocity (UI(0),UI(1)) at (xs,ys)
    uiv =xInterpolate( xs,ys,u,v,par );

    % trapezoidal rule
    uiabs=max(abs(uiv(1))+abs(uiv(2)),epsU);
    dt = min( abs(dt), cfla/uiabs);  
    if( cfl<0 ) dt = -dt; end 


    x=x + .5*dt*(uOld+uiv(1));
    y=y + .5*dt*(vOld+uiv(2));
      
    % if( x<xa || x>xb || y<ya || y>yb )
    %   if( par.idebug ) fprintf("**Line leaves region xs=%e, ys=%e, dt=%e \n",xs,ys,dt); end
    %   break;   %      ....outside plotting bounds
    % end

    if( par.idebug>1 )
      fprintf(" it=%3d: x=%10.2e, y=%10.2e, ui(0,0)=%10.2e ui(0,1)=%10.2e, dt=%10.2e, cfl=%4.1f \n",it,x,y,uiv(1),uiv(2),dt,cfl);
    end

    % uValue=(sqrt( uiv(1)^2 + uiv(2)^2 )-uMin)*uvfact;
    uValue=sqrt( uiv(1)^2 + uiv(2)^2 );

    par.maxSpeed = max(par.maxSpeed,uValue);

    % save point and normalized colour value
    xp(it) = x; 
    yp(it) = y; 
    cp(it) = uValue; % uValue is in [0,1]

    %  glVertex2(XSCALE(x),YSCALE(y));

    %  mark cell as being passed through
    ixg=floor((x-xa)*xgf) + 1;
    iyg=floor((y-ya)*ygf) + 1;
      
    if( ixg>=1 && ixg<=nxg && iyg>=1 && iyg<=nyg )

      if( ixg~=ixg0 || iyg~=iyg0 )

        % only 2 streamlines allowed per cell
        if( par.maskForStreamLines(ixg,iyg)>=2 )
        	if( par.idebug>1 ) fprintf('stop streamline, only 2 streamlines per cell\n'); end
          break;
        end
        if( par.idebug>1 ) fprintf('mark streamline mask (ixg,iyg)=(%3d,%3d)\n',ixg,iyg); end
        par.maskForStreamLines(ixg,iyg)=par.maskForStreamLines(ixg,iyg)+1;
        ixg0=ixg;
        iyg0=iyg;
        if( par.maskForStreamLines(ixg,iyg)==1 &&  mod(ixg,lax)==0 && mod(iyg,lay)==0 )
         
          % draw arrow ...  
          par.numArrows=par.numArrows+1;
          par.arrows(par.numArrows,1) = x;
          par.arrows(par.numArrows,2) = y;
          
          if( cfl>0 ) asign=1.; else asign=-1; end 
          if( it>1 )
	          dxp = asign*(xp(it)-xp(it-1));
	          dyp = asign*(yp(it)-yp(it-1));
	        else
	        	dxp =1; dyp=1;
	        end
          theta = atan2(dyp,dxp);
          par.arrows(par.numArrows,3) = theta;

        end % end draw arrow 

        %   try and check for closed loops, since we allow for 2 lines per cell we want
        %  to prevent drawing a closed loop twice
        dpath = dpath + abs(x-xi)+abs(y-yi);
      end

      if( it>25 && ixg==ixg1 && iyg==iyg1 && dpath > 0.05*dxmx )
      	if( par.idebug>1 ) fprintf('stop streamline, loop?\n'); end
        break;
      end

    end 

  end % end for it

  xtp=xi; % return last point
  ytp=yi;

  return
end



% % #define XV(i1) xvp[i1]
% %   XV(0)=x; XV(1)=y;

% % #define indexGuess(i) indexGuessp[i]
% % #define UI(i1) uip[i1]
  
% %     int * par.maskForStreamLinesp = par.maskForStreamLines_.Array_Descriptor.Array_View_Pointer1;
% %     const int par.maskForStreamLinesDim0=par.maskForStreamLines_.getRawDataSize(0);
% % #define par.maskForStreamLines(i0,i1) par.maskForStreamLinesp[i0+par.maskForStreamLinesDim0*(i1)]
  

%   int notOk = xInterpolateOpt(2,componentsToInterpolate,xvp,indexGuessp,uip,uv,gc,intopt );

%   % starting guess for next point is the start position of this point:
%   indexGuess0=indexGuess(0),indexGuess1=indexGuess(1),indexGuess3=indexGuess(3); % save


%   % GL_GraphicsInterface & gigl = (GL_GraphicsInterface &)gi;
%   % RealArray & gb = gigl.globalBound[0];
%   % real *rotationCenter = gigl.rotationCenter[0];

%   bool globalBoundSet = gb(0,0)< 1.e100;
%   const real aspectRatio = gi.getKeepAspectRatio() || !globalBoundSet ? 1. : (gb(1,1)-gb(0,1))/max(REAL_MIN,gb(1,0)-gb(0,0));
%    //  printF("drawAStreamLine: aspectRatio=%6.3f globalBoundSet=%i\n",aspectRatio,(int)globalBoundSet);

%   const bool glLines=false ; // line strip is a bit faster. true;
%   const real epsU = psp.streamLineStoppingTolerance*uMax;   // stop streamlines when |u|+|v| < epsU )
%   const real size=arrowSize*min(1.,40./max(nxg,nyg));   //  size for arrows
  
%   if( notOk==0 ) // This means we can interpolate
%   {
%     if( debug) printf("***Start a line x=%e, y=%e \n",x,y);
    
% //     glLineWidth(psp.size(GraphicsParameters::streamLineWidth)*psp.size(GraphicsParameters::lineWidth)*
% //              gi.getLineWidthScaleFactor());
%     if( glLines )
%       glBegin(GL_LINES); 
%     else
%       glBegin(GL_LINE_STRIP);
    
%     gi.setColourFromTable( (sqrt( SQR(UI(0))+SQR(UI(1)) )-uMin)*uvfact,psp);
%     glVertex2(XSCALE(x),YSCALE(y));

%     int ixg0=-1, iyg0=-1, ixg1=int((x-xa)*xgf), iyg1=int((y-ya)*ygf);
%     real dpath=0.;
    
%     int index=-1, index2;  // for colour table
%     // ...........Take time steps it=2,3,...,nt  (t=dt,2dt,...,(nt-1)*dt)
%     for( int it=2; it<=nt; it++ )
%     {
%       xi=x;   yi=y;
%       real uiabs=fabs(UI(0))+fabs(UI(1));
%       if( uiabs < epsU )  // flow is too slow here to move anywhere
%       {
%         // if( it==2 ) printf("******** Exiting, flow too slow, it=%i *************\n",it);
%         break;
%       }
      
%       // real dt= min(dtmx,cfla)/uiabs; // *wdh* 080131
%       real dt= cfla/uiabs;
%       dt= cfl > 0. ? dt : -dt;
%       t+=dt;
%       real xs=x+dt*UI(0);
%       real ys=y+dt*UI(1);
%       if( xs<xa || xs>xb || ys<ya || ys>yb )
%       {
%         if( debug ) printf("**Line leaves region xs=%e, ys=%e, dt=%e \n",xs,ys,dt);
%         break;  //     ....outside plotting bounds
%       }
%       real uOld=UI(0), vOld=UI(1);
      
%       //    interpolate velocity (UI(0),UI(1)) at (xs,ys)

%       XV(0)=xs; XV(1)=ys;
%       // notOk = xInterpolate(numberOfPointsToInterpolate,componentsToInterpolate,xv,indexGuess,
%       //                   ui,uv,gc,intopt); 
%       notOk = xInterpolateOpt(2,componentsToInterpolate,xvp,indexGuessp,uip,uv,gc,intopt );
      
%       if( notOk!=0 ) // this means we cannot interpolate
%       {
%         if( debug ) printf("**Line ends (notOk) xs=%e, ys=%e, dt=%e \n",xs,ys,dt);
%         break;   //  ....unable to interpolate
%       }
      
%       // trapezoidal rule
%       uiabs=max(fabs(UI(0))+fabs(UI(1)),epsU);
%       // *wdh* 080131 : do not recompute dt here -- but what if uiabs increases a lot ?
%       // dt=min(dtmx,cfla)/uiabs; // *wdh* 080131 

%       // dt = cfla/uiabs;  
%       dt = min( fabs(dt), cfla/uiabs);  

%       dt= cfl > 0. ? dt : -dt;
% //       x=.5*(xs + x+dt*UI(0));
% //       y=.5*(ys + y+dt*UI(1));
%       x=x + .5*dt*(uOld+UI(0));
%       y=y + .5*dt*(vOld+UI(1));
      
%       if( debug )
%         printf(" x=%e, y=%e, ui(0,0)=%e ui(0,1)=%e, dt=%e, cfl=%4.1f \n",x,y,UI(0),UI(1),dt,cfl);

%       real uValue=(sqrt( SQR(UI(0))+SQR(UI(1)) )-uMin)*uvfact;
%       if( psp.colourTable==GraphicsParameters::rainbow )
%       {
%         index2 = min(max(int(uValue*255+.5),0),255);
%         if( index2!=index )
%         {
%           index=index2;
%           glColor3f(colourTable[index][0]/255.,colourTable[index][1]/255.,colourTable[index][2]/255.);
%         }
%       }
%       else 
%         gi.setColourFromTable( uValue,psp);

%       glVertex2(XSCALE(x),YSCALE(y));
%       if( glLines ) glVertex2(XSCALE(x),YSCALE(y));  // what is this ?

%       //  mark cell as being passed through
%       int ixg=int((x-xa)*xgf);
%       int iyg=int((y-ya)*ygf);
      
%       if( ixg>=0 && ixg<nxg && iyg>=0 && iyg<nyg )
%       {
%         if( ixg!=ixg0 || iyg!=iyg0 )
%         {
%           // only 2 streamlines allowed per cell
%           if( par.maskForStreamLines(ixg,iyg)>=2 )
%             break;
%           par.maskForStreamLines(ixg,iyg)++;
%           ixg0=ixg;
%           iyg0=iyg;
%           if( par.maskForStreamLines(ixg,iyg)==1 &&  (ixg % lax)==0 && (iyg % lay)==0 )
%           {
%             // draw arrow ...   ** it did not make a big difference in cpu to turn these arrows off 
%             glEnd();   // end line so we can plot the arrow

%             // real angle=atan2((double)UI(1),(double)UI(0))*180./Pi+90.;
%             real xScale=XSCALE(x), yScale=YSCALE(y);
%             real angle=atan2(double(UI(1)*psp.yScaleFactor),double(UI(0)*psp.xScaleFactor*aspectRatio))*180./Pi+90.;
%             gi.setColour(GenericGraphicsInterface::textColour);  // arrow colour
%             gi.xLabel("V",xScale,yScale,size,0,angle,psp);

%             glLineWidth(psp.size(GraphicsParameters::streamLineWidth)*
%                         psp.size(GraphicsParameters::lineWidth)*
%                         gi.getLineWidthScaleFactor());
%             if( glLines )
%               glBegin(GL_LINES);   // restart the line
%             else
%               glBegin(GL_LINE_STRIP);   // restart the line

%             real uValue=(sqrt( SQR(UI(0))+SQR(UI(1)) )-uMin)*uvfact;
%             if( psp.colourTable==GraphicsParameters::rainbow )
%             {
%               index = min(max(int(uValue*255+.5),0),255);
%               glColor3f(colourTable[index][0]/255.,colourTable[index][1]/255.,colourTable[index][2]/255.);
%             }
%             else 
%               gi.setColourFromTable( uValue,psp);

%             glVertex2(xScale,yScale);
%           }
%         }
%       }
%       //  try and check for closed loops, since we allow for 2 lines per cell we want
%       // to prevent drawing a closed loop twice
%       dpath+=fabs(x-xi)+fabs(y-yi);
%       if( it>25 && ixg==ixg1 && iyg==iyg1 && dpath > 0.05*dxmx )
%         break;
%     }
%     glEnd();
%   }
  
%   xtp=xi; 
%   ytp=yi;

%   // starting guess for next point is the start position of this point:
%   indexGuess(0)=indexGuess0; indexGuess(1)=indexGuess1; indexGuess(3)=indexGuess3;   

% }


