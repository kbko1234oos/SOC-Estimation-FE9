% SOC and Vc calculations
% NOTE: I(t), Cbat, Cc, Rc, Vc, values all change w/ respect to time, so all
% the variable values in the below code is with respect to a single time point.
% Regardless, the equations will stay the same.

%------------------------------------------------------------------------------


% I(t) = (function definition);
% Cbat = (value);
% Cc = (value);
% Rc = (value);
% Vc = (value);

SOC = (-1)*I(t) / Cbat;
Vc = (1/Cc) * I(t) - ((1/(Cc * Rc)) * Vc);
