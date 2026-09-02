%
% Evaluate the outward unit normal on the face (side,axis)
% 
function [n1,n2,rxNorm] = getBoundaryNormal( side,axis,I1b,I2b, gf,cur, par )

  is = 1-2*(side-1);

  % fprintf('getBoundaryNormal cur=%d side=%d axis=%d I1b=%d I2b=%d\n',cur,side,axis,I1b,I2b);
  n1 = -is*gf{cur}.rx(I1b,I2b,axis,1); % outward normal is (n1,n2)
  n2 = -is*gf{cur}.rx(I1b,I2b,axis,2); 
  rxNorm = sqrt( n1.^2 + n2.^2 ); 
  n1 = n1./rxNorm;
  n2 = n2./rxNorm;

  return 
end