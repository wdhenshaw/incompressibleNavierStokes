% ----------------------------------------------------------------------------------------------
% Return the matlab index's corresponding to "interior" points where the equation is applied.
% 
% Include end points for periodic or Neumann type conditions.
%
% Input:
%  component : par.pc, par.uc, par.vc 
% ----------------------------------------------------------------------------------------------
function [I1,I2] = getIndexInterior( gid,component,par )

  % *** CHECK ME ***

  if( nargin<2  ) extra=0; end

  gidi = gid; % choose interior and boundary points by default 

  for axis=1:par.nd
  for side=1:2
  	isv(1)=0; isv(2)=0;
  	isv(axis) = 1-2*(side-1);
  	is = isv(axis);

  	if( par.bc(side,axis)==par.dirichlet )
  	  % dirichlet: (u,v,p) = given
  		gidi(side,axis)=gidi(side,axis)+is;  % exclude the boundary points
  		
  	elseif( par.bc(side,axis)==par.periodic || ...
  		      par.bc(side,axis)==par.outflow )
  	  % ouflow: p.n + alpha*p = 0 , extrap (u,v)
      % include boundary points 
	  elseif( par.bc(side,axis)==par.noSlipWall || ...
	  		    par.bc(side,axis)==par.inflow )
	    % (u,v) = given
	    % p.n = given
  	  if( component ~= par.pc )
        gidi(side,axis)=gidi(side,axis)+is;  % exclude the boundary points for u,v
      end
    elseif( par.bc(side,axis)==par.slipWall )
      if(  (axis==1 && component == par.uc) || (axis==2 && component == par.vc)  )
        gidi(side,axis)=gidi(side,axis)+is;  % exclude the boundary points nv.uv 
      end        

    elseif( par.bc(side,axis)==par.pressureInflow )
    	% Give p and t.uv, extrap n.uv 
 	    if( component == par.pc || (axis==1 && component == par.vc) || (axis==2 && component == par.uc)  )
        gidi(side,axis)=gidi(side,axis)+is;  % exclude the boundary points for p and t.u
      end    	

	  else
	  	fprintf('getIndexInterior:ERROR -- finish me');
	  	pause
	  end
  end

  % fprintf('gidi=[%d,%d]x[%d,%d]\n',gidi(1,1),gidi(2,1),gidi(1,2),gidi(2,2));


  [I1,I2]=getIndex( gidi );


end


