%
%  Compute the true solution and the errors
%

function [maxErr,perr,uerr,verr,div] = getErrors( tn,un,vn,pn,par )

  dx  = par.dx;
  dy  = par.dy;  
  dr = par.dr(1);
  ds = par.dr(2);  

	[I1,I2] = getIndex( par.gid );

  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y 

  Dr2 = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )/(2.*dr);   % u.r to second order 
  Ds2 = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )/(2.*ds);   % u.s   	

	uTrue = par.ue(par.x(:,:,1),par.x(:,:,2),tn);  % eval exact solution
	vTrue = par.ve(par.x(:,:,1),par.x(:,:,2),tn);  % eval exact solution
	pTrue = par.pe(par.x(:,:,1),par.x(:,:,2),tn);  % eval exact solution
	uerr = un-uTrue;
	verr = vn-vTrue;
	perr = pn-pTrue;
	uErrMax = max(max(abs(un(I1,I2)-uTrue(I1,I2))));  % max-norm error
	vErrMax = max(max(abs(vn(I1,I2)-vTrue(I1,I2))));  % max-norm error
	pErrMax = max(max(abs(pn(I1,I2)-pTrue(I1,I2))));  % max-norm error

  div=zeros(par.Ngx,par.Ngy);

  if( par.isCartesian )

	  div(I1,I2) = Dzx(un,I1,I2)+Dzy(vn,I1,I2);

	else

    % --- curvilinear ---

    for( i1=I1 )
    for( i2=I2 ) 
      rx = par.rx(i1,i2,1,1);
      ry = par.rx(i1,i2,1,2);
      sx = par.rx(i1,i2,2,1);
      sy = par.rx(i1,i2,2,2);

      ur = Dr2(un,i1,i2); us = Ds2(un,i1,i2);
      vr = Dr2(vn,i1,i2); vs = Ds2(vn,i1,i2);

      ux = rx*ur + sx*us;  % save for use below in AD
      vy = ry*vr + sy*vs;
      div(i1,i2) = ux + vy;
    end
    end
	end

	divMax=max(max(abs(div(I1,I2)))); 

	maxErr(1) = pErrMax;
	maxErr(2) = uErrMax;
	maxErr(3) = vErrMax;
	maxErr(4) = divMax;

return
end