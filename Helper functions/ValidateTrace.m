function flag = ValidateTrace(Settings,Trace)
% Flag = ValidateTrace(Settings,Trace) returns 1 if the Trace is a valid 
% Trace or 0 if not, based on the requirements that the trace is not:
% -	smaller than the minimum length
%%
flag = 1;
%% If the trace is not long enough
if length(Trace) < Settings.minimum_traclength
    flag = 0;
end
 