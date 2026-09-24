% Function to compute modulus with base 1 
function [ mm ] = modBase1( m,numberOfTimeLevels )
  mm = mod( m-1,numberOfTimeLevels)+1; 
end

