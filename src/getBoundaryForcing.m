%
% Evaluate forcing functions for the boundary conditions
% 
% bd(I1,I2,1:nvc) : boundary data
function [par] = getBoundaryForcing( side,axis, I1,I2, par )

  FINISH ME 

  if( par.bc(side,axis)==par.inflow )
    par.bd(I1,I2,1) = par.uInflow;  
    par.bd(I1,I2,2) = 0;
  elseif( strcmp(par.ic,'shear') )
    if( axis==2 )
      if( side==1 )
        par.bd(I1,I2,1) = -1;
      else
        par.bd(I1,I2,1) = +1;
      end
    end
  end

  return
end