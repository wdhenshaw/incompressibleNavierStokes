% ----------------------------------------------------------------------------------------------
% Return the matlab index's corresponding to points on a 
% given ghost line on a boundary defined by (side,axis)
%
%    ADJUST TANGENTIAL DIRECTIONS BASED ON THE BOUNDARY CONDITIONS
% ----------------------------------------------------------------------------------------------
function [I1g,I2g] = getAdjustedGhostIndex(side,axis,ghost,par)


  axis2 = mod(axis,par.nd)+1; % adjacent axis  axis=1 -> axis2=2, axis=2 -> axis2=1
  bid=zeros(2,1); % boundary index in tangential direction
  for( side2=1:2 )
    bid(side2) = par.gid(side2,axis2);
  end

  if( par.boundaryCondition(side,axis)==par.dirichlet )
  	% return extended boundary for Dirichlet BC
	for( side2=1:2 )
	  bid(side2) = par.dim(side2,axis2);
	end  	

  elseif( par.nd==2 && ...
  	      (par.boundaryCondition(side,axis)==par.characteristic || ...
  	       par.boundaryCondition(side,axis)==par.neumann ) )
    % exclude ends when adjacent side is a Dirichlet BC
    for( side2=1:2 )
    	is2 = 1-2*(side2-1); 
    	if( par.boundaryCondition(side2,axis2)==par.dirichlet )
    		bid(side2)=par.gid(side2,axis2)+is2; 
    	end
    end
  end

  if( axis==1 )
    I1g=par.gid(side,axis) - (1-2*(side-1))*ghost;
    I2g=bid(1):bid(2);
  elseif( axis==2 )
    I1g=bid(1):bid(2);
    I2g=par.gid(side,axis) - (1-2*(side-1))*ghost;
  else
    fprintf('getAdjustedGhostIndex:ERROR: invalid value for axis=%d\n',axis); pause; 
  end;  

  return;

end