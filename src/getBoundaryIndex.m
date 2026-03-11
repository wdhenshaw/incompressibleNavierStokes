% ----------------------------------------------------------------------------------------------
% Return the matlab index's corresponding to points on a given boundary defined by (side,axis)
%
% ----------------------------------------------------------------------------------------------
function [I1b,I2b] = getBoundaryIndex(side,axis,par)

  % globalDeclarations;

  if( axis==1 )
    I1b=par.gid(side,axis);
    I2b=par.gid(1,2):par.gid(2,2);
  elseif( axis==2 )
    I1b=par.gid(1,1):par.gid(2,1);  
    I2b=par.gid(side,axis);
  else
    fprintf('getBoundaryIndex:ERROR: invalid value for axis=%d\n',axis); pause; 
  end 

end
