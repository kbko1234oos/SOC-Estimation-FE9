%% Discretization example for FSAE folks 
clc; clear; close all;

%% Analytical solution

dt = 0.01;
td = 0:dt:10;
ta = 0:0.01:10;

x = zeros(length(td), 1); % x(1) = 0
v = zeros(length(td), 1); % v(1) = 0
v(1) = 1;

for i = 1:length(td)
    if i ~= length(td)
        v(i+1) = velocity_discrete(v(i), dt);
        x(i+1) = position_discrete(x(i), v(i+1), dt);
    end
    error_x(i) = x(i) - position_analytical(td(i));
    error_v(i) = v(i) - velocity_analytical(td(i));
end

figure
subplot(2, 1, 1)
plot(td, x)
hold on
plot(ta, position_analytical(ta), '--');
xlabel("Time (s)")
ylabel("Position (m)")
legend("Discrete", "Analytical")
subplot(2, 1, 2)
plot(td, v)
hold on
plot(ta, velocity_analytical(ta), '--');
xlabel("Time (s)")
ylabel("Velocity (m/s)")
legend("Discrete", "Analytical")

figure
plot(td, error_v)
hold on
plot(td, error_x)

%%
clear x

for i = 1:length(td)-1
   x(:, i+1) = forward_euler(@falling_body_equations, x(i), dt);
end

figure
subplot(2, 1, 1)
plot(td, x(1, :))
hold on
plot(ta, position_analytical(ta), '--');
xlabel("Time (s)")
ylabel("Position (m)")
legend("Discrete", "Analytical")
subplot(2, 1, 2)
plot(td, v(2, :))
hold on
plot(ta, velocity_analytical(ta), '--');
xlabel("Time (s)")
ylabel("Velocity (m/s)")
legend("Discrete", "Analytical")

%% Discrete solution to d/dt(x) = v; d/dt(v) = 2 - v;
% Using finite difference (x_{k+1} - x_{k})/dt approx d/dt(x)

function xkp1 = position_discrete(xk, vkp1, dt)
    xkp1 = xk + vkp1 * dt;
end

function vkp1 = velocity_discrete(vk, dt)
    vkp1 = vk + dt * (2 - vk);
end

%% Solution to d/dt(x) = v; d/dt(v) = 2 - v;

function x = position_analytical(t)
    x = exp(-t) + 2*t - 1;
end

function v = velocity_analytical(t)
    v = -exp(-t) + 2;
end

%% Falling body equations
function xdot = falling_body_equations(x)
    % x is a vector, x(1) = x, x(2) = v
    v = x(2);
    xdot(1) = x(2);
    xdot(2) = 2 - x(2);
end

%% Forward euler solver
function xkp1 = forward_euler(xdot, xk, dt)
    xkp1 = xdot(xk) * dt + xk;
end