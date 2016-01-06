function [ psi ] = phase_field_function( x,epsilon )
%PHASE_FIELD_FUNCTION 

psi = 1-2./(1+exp(-distfunc(x)/epsilon));

end

