%
%  Test routine for spline functions
%
clear all; % clear all variables 
set(gca,'FontSize',16); % set font-size for labels

par.savePlots=0; 

clearOpenFigures(1:7);

if( 1==1 )
  % --- Check the curvature routine ---

  % spline for a circle

  ns=21; 
  a=0; b=2*pi;
  thetas=linspace(a,b,ns)'; 

  rada=1;
  radb=2;

  xs = rada*cos(thetas);
  ys = radb*sin(thetas);

  % create the cubic splines 
  xSpline = csape(thetas,xs,'periodic');
  ySpline = csape(thetas,ys,'periodic');


  % evaluate the spline 
  m=101;
  thetav = linspace(a,b,m)';

  xv = fnval(xSpline,thetav);
  yv = fnval(ySpline,thetav);

  figure(1);
  plot( xv,yv,'b-', xs,ys,'o');
  grid on; xlabel('x'); ylabel('y');
  title(sprintf('spline: ellipse (ra,rb)=(%g,%g)',rada,radb));
  axis equal;

  xSpliner  = fnder(xSpline, 1);   %  spline for first deriv 
  xSplinerr = fnder(xSpline, 2);

  ySpliner  = fnder(ySpline, 1);   
  ySplinerr = fnder(ySpline, 2);   


  I = 1:m;
  yvr  = zeros(m,2);
  yvrr = zeros(m,2);
  yvr(I,1) = fnval(xSpliner,thetav);
  yvr(I,2) = fnval(ySpliner,thetav);

  yvrr(I,1) = fnval(xSplinerr,thetav);
  yvrr(I,2) = fnval(ySplinerr,thetav);

  % --- evaluate the curvature 
  %      kappa(i) =  ( nv(i)^T yv.rr(i) )/( yv_r(i) ^T  yv_r(i) )
  % 
  % normal : nv(i) = [-yv(i,2),yv(i,1)]/( sqrt(yv(i,1)^2 + yv(i,2)^2 )
  nv=zeros(m,2);
  nv(I,1) = -yvr(I,2);
  nv(I,2) =  yvr(I,1);

  kappa=zeros(m,1);
  kappa = ( nv(I,1).*yvrr(I,1) + nv(I,2).*yvrr(I,2) )./( (yvr(I,1).*yvr(I,1) + yvr(I,2).*yvr(I,2)).^(3/2)  );

  % true curvature of an ellipse: 
  kappaTrue = zeros(m,1);
  kappaTrue = (rada*radb)./( ( rada^2*sin(thetav).^2 + radb^2*cos(thetav).^2 ).^(3/2) );
  figure(2);
  plot( thetav,kappa,'r-', 'LineWidth', 2); hold on;
  plot( thetav,kappaTrue,'k-',  'LineWidth', 1); hold off;
  xlabel('\theta'); 
  legend('\kappa (computed)', '\kappa (true)');
  grid on;
  title(sprintf('curvature: ellipse (ra,rb)=(%g,%g)',rada,radb));
  axis equal;

  maxErr = max(abs(kappa-kappaTrue));
  fprintf('Curvature: max-err=%9.2e (rada=%g, radb=%g, ns=%d, m=%d)\n',maxErr,rada,radb,ns,m);

  figure(3);
  plot( thetav,kappa-kappaTrue,'r-', 'LineWidth', 2);
  xlabel('\theta'); 
  legend('\kappa - \kappa (true)');
  grid on;
  title(sprintf('curvature error: ellipse (ra,rb)=(%g,%g)',rada,radb));

  pause

end


if( 0==1 )
  % test: (exact for clamped and curvature BC's)
  f = @(x) 1+x+x.^2+x.^3; 
  fx = @(x) 1. + 2.*x+ 3.*x.^2; 
  fxx = @(x) 2. + 6.*x; 
else
  f =  @(x) exp(sin(pi*x).^2); 
  fx =  @(x) pi*sin(2*pi*x).*f(x); 
  fxx =  @(x) 2.*pi^2*cos(2*pi*x).*f(x) + (pi*sin(2*pi*x)).^2.*f(x);
end;


n=9; 
a=-1; b=1.; 
xs=linspace(a,b,n)'; % spline data 
ys=f(xs);

if( 1==1 )
  % TEST the matlab routine csape 

  conds{1} = 'clamped';
  conds{2} = 'periodic';
  conds{3} = 'complete';
  conds{4} = 'not-a-knot';
  conds{5} = 'second';    % can input 2nd derivatives as extra parameters, zero by default 

  for icase=1:5
    label=conds{icase};

    pp = csape(xs,ys,conds{icase});
    % evaluate the spline 
    m=201;
    x=linspace(a-.1,b+.1,m)';

    y = fnval(pp,x);

    ff = f(x); 
    plot( x,ff,'r-', x,y,'b-', xs,ys,'o', 'LineWidth',2 );
    title(sprintf('Spline (%s)',label)); 
    xlabel('x'); ylabel('y');  grid on;
    legend('f','spline','knots');
    xlim([a-.1,b+.1]);
    pause

    %   Compute the first derivative
    ppx = fnder(pp, 1); 
    yx = fnval(ppx,x);

    ffx=fx(x); 

    plot( x,ffx,'r-', x,yx,'b-', 'LineWidth',2 );
    title(sprintf('Spline first derivative (%s)',label)); 
    xlabel('x'); ylabel('y'); grid on;
    legend('fx','Sx');
    xlim([a-.1,b+.1]);
    pause

  end 
end




% ----------- OLD WAY -------


% for icase=1:3

%   if( icase==1 )
%     label='clamped'; 
%     bcLeft =1; gLeft =fx(a);  % clamped
%     bcRight=1; gRight=fx(b);  % clamped 
%   elseif( icase==2 )
%     label='curvature'; 
%     bcLeft =2; gLeft =fxx(a);  % 2nd-derivatve
%     bcRight=2; gRight=fxx(b);  % 2nd-derivative 
%   else
%     label='natural'; 
%     bcLeft =0; gLeft =0; % natural 
%     bcRight=0; gRight=0;  % natural
%   end; 

%   coeff = splineCoeff( xs,ys, bcLeft,gLeft, bcRight,gRight );
  
%   % evaluate the spline 
%   m=201;
%   x=linspace(a-.1,b+.1,m);
%   [y,yx] = splineEval( xs,ys,coeff, x );
  
  
%   ff = f(x); 
  
%   plot( x,ff,'r-', x,y,'b-', xs,ys,'o', 'LineWidth',2 );
%   title(sprintf('Spline (%s)',label)); 
%   xlabel('x'); ylabel('y'); 
%   legend('f','spline','knots');
%   xlim([a-.1,b+.1]);
%   if( par.savePlots )
%     print('-depsc2',sprintf('splineCase%s.eps',label));
%   end 

%   pause;

%   ffx=fx(x); 

%   plot( x,ffx,'r-', x,yx,'b-', 'LineWidth',2 );
%   title(sprintf('Spline first derivative (%s)',label)); 
%   xlabel('x'); ylabel('y'); 
%   legend('fx','Sx');
%   xlim([a-.1,b+.1]);

%   if( par.savePlots )
%     print('-depsc2',sprintf('splineDerivativeCase%s.eps',label));
%   end
%   pause;

%   err = ff-y;
%   errx = ffx-yx;
%   plot( x,err,'r-', x,errx/10,'b-', 'LineWidth',2 );
%   title(sprintf('Error in the spline and derivative (%s)',label));
%   xlabel('x'); ylabel('err'); 
%   legend('error in f','error in fx/10');
%   xlim([a-.1,b+.1]);
%   if( par.savePlots )
%     print('-depsc2',sprintf('splineErrorCase%s.eps',label));
%   end
%   pause
% end; % end icase
