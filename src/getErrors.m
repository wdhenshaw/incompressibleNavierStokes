%
%  Compute the true solution and the errors
%

function [maxErr,perr,uerr,verr,div] = getErrors( tn,un,vn,pn,par )

  dx  = par.dx;
  dy  = par.dy;  

	[I1,I2] = getIndex( par.gid );

  Dzx = @(u,I1,I2) ( u(I1+1,I2) -u(I1-1,I2) )*(1./(2.*dx));   % u.x
  Dzy = @(u,I1,I2) ( u(I1,I2+1) -u(I1,I2-1) )*(1./(2.*dy));   % u.y  	

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
	div(I1,I2) = Dzx(un,I1,I2)+Dzy(vn,I1,I2);
	divMax=max(max(abs(div(I1,I2)))); 

	maxErr(1) = pErrMax;
	maxErr(2) = uErrMax;
	maxErr(3) = vErrMax;
	maxErr(4) = divMax;

return
end