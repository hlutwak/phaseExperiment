function err = fitPsychometricFunction(p,results,functionName)

%err = fitPsychometricFunction(p,results,functionName)
%
%Calculates maximum likelihood fit of a psychometric function to
%psychophysical results.
%
%Inputs
%   p            structure containing parameters for the function
%   results      structure containing fields 'intensity' and 'response',
%                which are vectors for intensity values shown for each trial
%                and the corresponding binary response (1 = correct)
%   functionName 'Weibull' (requires p.b and p.t)
%                 'Normal'  (requires p.u and p.s)
%
%Outputs
%   err          negative of maximum likelihood probability (so that small is good)

%11/13/2007     gmb wrote it.
%3/7/2008       gmb fixed 'hack' to keep p.t>0 only for Weibull function

if ~exist('functionName')
    functionName = 'Weibull';
end

if isstr(functionName)
    evalStr = sprintf('w = %s(p,results.intensity);',functionName);
    eval(evalStr);
else    
    w= functionName(p,results.intensity);
end


%hack!  pull the values of w away from 0 and 1

w = w*.99+.005;

err = -sum(results.response.*log(w) + (1-results.response).*log(1-w));

if strcmp(functionName,'Weibull')
    err = err + (-min(p.t,0))^2;
end
