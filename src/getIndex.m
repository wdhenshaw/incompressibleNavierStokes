% ----------------------------------------------------------------------------------------------
% Return the matlab index's corresponding to the gid(0:1,0:nd-1) array
%
% ----------------------------------------------------------------------------------------------
function [I1,I2] = getIndex( gid,extra )

  if( nargin<2  ) extra=0; end

  I1 = (gid(1,1)-extra):(gid(2,1)+extra);
  I2 = (gid(1,2)-extra):(gid(2,2)+extra);

end