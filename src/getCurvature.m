%
% Compute the curvature of the free-surface
%
%    If the free-surface curve is defined by yv(r) then the curvature is
%
%      kappa(i) =  ( nv(i)^T yv.rr(i) )/( yv_r(i) ^T  yv_r(i) )
%            
%
function kappa = getCurvature( t, gf,cur, par )


  if( par.idebug>0 )
  	fprintf('*** getCurvature called for cur=%d ***\n',cur);
  end
  % free surface positions: 
  % gf{cur}.eta(I1,1) = xv(I1);
  % gf{cur}.eta(I1,2) = yv(I1);  

  %  gf{cur}.xcs : x-spline
  %  gf{cur}.ycs : y-spline

  % xSpliner  = fnder(gf{cur}.xcs, 1);   %  spline for first deriv 
  % xSplinerr = fnder(gf{cur}.xcs, 2);

  % ySpliner  = fnder(gf{cur}.ycs, 1);   
  % ySplinerr = fnder(gf{cur}.ycs, 2); 

  I = 1:par.Ngx;
  yvr  = zeros(par.Ngx,2);
  yvrr = zeros(par.Ngx,2);
  yvr(I,1) = fnval(gf{cur}.xcs.ppx,gf{cur}.rvs); % x.r = first derivative
  yvr(I,2) = fnval(gf{cur}.ycs.ppx,gf{cur}.rvs); % y.r 

  yvrr(I,1) = fnval(gf{cur}.xcs.ppxx,gf{cur}.rvs);  % x.rr second deriv
  yvrr(I,2) = fnval(gf{cur}.ycs.ppxx,gf{cur}.rvs);  % y.rr  

  nv=zeros(par.Ngx,2);
  nv(I,1) = -yvr(I,2);
  nv(I,2) =  yvr(I,1);

  kappa=zeros(par.Ngx,1);
  kappa = ( nv(I,1).*yvrr(I,1) + nv(I,2).*yvrr(I,2) )./( (yvr(I,1).*yvr(I,1) + yvr(I,2).*yvr(I,2)).^(3/2)  );

  if( 1==1 && t<= 5*par.dt )
  	figure(7);
  	plot(gf{cur}.rvs,kappa,'-x');
  	xlabel('r');
  	title(sprintf('Curvature, t=%9.2e, \\gamma=%9.2e, cur=%d\n',t, par.gamma,cur));
  	grid on;
  	pause
  end


	return
end