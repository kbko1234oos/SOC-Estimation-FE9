
%--------------------------------------------------------------------------
% At time k-1
%--------------------------------------------------------------------------

% k = (time value);
% T = (value);
% Cc = (value);
% Rc = (value);

% xk1_k1 = (value);
% uk1 = (value);

% Pk1_k1 = (value);
% Aprime_k1 = (matrix);
% Eprime_k1 = (matrix);
% Qk1 = (value or matrix);

%--------------------------------------------------------------------------
% At time k
%--------------------------------------------------------------------------

k = k+1;

% fk(xk1_k1, uk1, 0) = (function definition);

% Cprime_k = (matrix);
% Fprime_k = (matrix);
% Rk = (value);

% yk = (value or matrix);
% uk = (value);
% hk(xk_k1, uk, 0) = (function definition);


xk_k1 = fk(xk1_k1, uk1, 0);
Pk_k1 = Aprime_k1 * Pk1_k1 * transpose(Aprime_k1) + Eprime_k1 * Qk1 * transpose(Eprime_k1);
Lk = Pk_k1 * transpose(Cprime_k) * inv(Cprime_k * Pk_k1 * transpose(Cprime_k) + Fprime_k * Rk * transpose(Fprime_k));
xk_k = xk_k1 + Lk * (yk - hk(xk_k1, uk, 0));
Pk_k = Pk_k1 - Lk * Cprime_k * Pk_k1;
