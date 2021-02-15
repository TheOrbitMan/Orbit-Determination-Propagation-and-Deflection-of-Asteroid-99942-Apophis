% -----------------------------------------------------------------------------
%
%                           procedure finitediff
%
% this procedure perturbs the components of the state vector for processing
% with the finite differencing for the a matrix.
%
%  author        : david vallado                  719-573-2600   15 jan 2008
%
%  references    :
%    vallado       2007, 753-765
%
% This code is taken from:
% https://celestrak.com/software/vallado-sw.php-------->Computer software in Matlab
% 
% Author: David Vallado
% --------------------------------------------------------------------------- */

function [deltaamt, xnomp] = finitediff(pertelem, percentchg, xnom)

          xnomp=xnom;
          deltaamt = xnom(pertelem) * percentchg;
          xnomp(pertelem, 1) = xnom(pertelem, 1) + deltaamt;
          
end  

        
