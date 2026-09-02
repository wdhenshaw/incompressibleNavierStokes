% ----------------------------------------------------------------------------------------------
% Return the matlab index's corresponding to points on a given boundary defined by (side,axis)
%    ADJUST TANGENTIAL DIRECTIONS BASED ON THE BOUNDARY CONDITIONS
% Note:
%   (1) This routine is for the VELOCITY components
%   (2) For the implicit matrix Dirichlet-type BC's should only be considered once at corners
% ----------------------------------------------------------------------------------------------
function [I1b,I2b] = getAdjustedBoundaryIndex(side,axis,par)

  % if( nargin<4  ) extra=0; end

  axis2 = mod(axis,par.nd)+1; % adjacent axis  axis=1 -> axis2=2, axis=2 -> axis2=1
  bid=zeros(2,1); % boundary index in tangential direction
  for( side2=1:2 )
    bid(side2) = par.gid(side2,axis2);
  end

  % Dirichlet BC's should only be considered once at corners 
  if( par.bc(side,axis)==par.dirichlet || par.bc(side,axis)==par.noSlipWall || par.bc(side,axis)==par.inflow )
  	dirichletType = 1; % full Dirichlet
  elseif( par.bc(side,axis)==par.slipWall || par.bc(side,axis)==par.pressureInflow )
  	dirichletType = 2; % partial dirichlet 
  elseif( par.bc(side,axis)==par.outflow || par.bc(side,axis)==par.traction )
  	dirichletType = 0 ;
  else
  	fprintf('getAdjustedBoundaryIndex:ERROR: unknown bc = %d for (side,axis)=(%d,%d)\n', par.bc(side,axis),side,axis);
  	error('error');
  end

  %    (1) At a corner between two D-type BC's the corner is assigned by the left or right side 
  %    (2) D-type takes precedence over non-D-type
  %    (3) Corners are excluded for two adjacent non-D-types (these must be treated separately)
  % 
  %             X                                    X                             X 
  %             X                                    X                             X 
  %   D-type    X                         non-D-type X                  non-D-type X 
  %             X                                    X                             X
  %             X Y Y Y Y Y Y Y                      Y Y Y Y Y Y Y                 ? Y Y Y Y Y Y  Y
  %                D-type, or non-D-type                  D-type                        non-D-type 
  % 
  if( (dirichletType && axis==2) || dirichletType~=1  )

     % exclude ends when adjacent side is a Dirichlet type BC
    for( side2=1:2 )
    	is2 = 1-2*(side2-1); 
    	if( ( par.bc(side2,axis2)==par.dirichlet || par.bc(side2,axis2)==par.noSlipWall || par.bc(side2,axis2)==par.inflow ) ||  ... 
    		  ( par.bc(side,axis)==par.slipWall && par.bc(side2,axis2)==par.slipWall ) )
    		bid(side2)=par.gid(side2,axis2)+is2; 
    	end
    end
  end

  if( axis==1 )
    I1b=par.gid(side,axis);
    I2b=bid(1):bid(2);
  elseif( axis==2 )
    I1b=bid(1):bid(2);
    I2b=par.gid(side,axis);
  else
    fprintf('getAdjustedBoundaryIndex:ERROR: invalid value for axis=%d\n',axis); pause; 
  end;  

  if( 1==0 ) fprintf('getAdjustedBoundaryIndex: [I1b,I2b]=[%d,%d][%d,%d]\n', I1b(1),I1b(end), I2b(1),I2b(end)); end

  return
end