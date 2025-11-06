% The 2nd order ODE y''+(omega^2)*y = 0
% can be written as two coupled first order equations
% y(n+1) = y(n) + h*yd(n) and
% yd(n+1) = yd(n) - h*(omega^2)*y(n)
% As Part C does not mention noise, assume noise free signal

function y = euler(h, omega, t, ic_y, ic_yd)
    % Euler algorithm to evaluate y''+(omega^2)*y = 0
    % Given initial conditons y(0) = sin(phi), yd(0) = omega*cos(phi)
    % Note domain should be a column vector
    
    % Initialise matrices for y and y'
    y = zeros(1, length(t));
    yd = zeros(1, length(t));
    
    % State ICs
    y(1) = ic_y;
    yd(1) = ic_yd;

    % Implement Euler's algorithm
    for n = 1:(length(t)-1)
        y(n+1) = y(n) + h*yd(n);
        yd(n+1) = yd(n) - h*(omega^2)*y(n);
    end
end

function U = rk4(F, t, U, h)
    % RK4 algorithm to evaluate y''+(omega^2)*y = 0
    % Given initial conditons y(0) = sin(phi), yd(0) = omega*cos(phi)

    % Implement RK4 algorithm
    for n = 1:(length(t)-1)
        k1 = h*F(t, U(:, n));
        k2 = h*F(t + 0.5*h, U(:, n) + 0.5*k1);
        k3 = h*F(t + 0.5*h, U(:, n) + 0.5*k2);
        k4 = h*F(t + h, U(:, n) + k3);

        U(:, n+1) = U(:, n) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
    end
end

% Parameters for synthetic signal to try
T = 2*pi; % Total time interval
H = [0.2, 0.1, 0.05, 0.025, 0.0125]; % Sequence of step sizes
omega = 2.0; % Frequency
phi = pi/6; % Phase
ic_y = sin(phi);
ic_yd = omega*cos(phi);

% Initialise array to store errors
euler_error_percentage = zeros(1, length(H));
rk4_error_percentage = zeros(1, length(H));

figure
tiledlayout(5, 1)
i = 1;
for h = H
    t = [0:h:2*pi]';
    exact_soln = sin(omega*t + phi);

    % Necessary for RK4
    % For RK4, need to pass a vector field F as a function of the state u
    % Let u = [y y']'. Then based on the coupled equations, let F = du/dt
    % so F = [y' y'']' = [y' -(omega^2)*y]'. This is defined like so:
    F = @(t, u) [u(2), -(omega^2)*u(1)]';

    % Create state vector U
    U = zeros(2, length(t));
    U(:, 1) = [ic_y ic_yd]';
    % Necessary for RK4

    euler_result = euler(h, omega, t, ic_y, ic_yd);
    U = rk4(F, t, U, h);
    rk4_result = U(1, :);

    euler_error_percentage(i) = abs(euler_result(end) - exact_soln(end))/exact_soln(end)*100;
    rk4_error_percentage(i) = abs(rk4_result(end) - exact_soln(end))/exact_soln(end)*100;

    nexttile  
    plot(t, exact_soln, '-rx', t, euler_result, '-go', t, rk4_result, '-bo')
    xlabel('time')
    ylabel('y')
    legend('Exact', 'Euler', 'RK4')

    i = i + 1;
end

figure
loglog(H, euler_error_percentage, '-rx', H, rk4_error_percentage, '-bx')
xlabel('step size h')
ylabel('percentage error')
legend('Euler', 'RK4')
slope_euler = sprintf('Euler slope = %0.1f', log10(euler_error_percentage(3)/euler_error_percentage(end))/log10(H(3)/H(end)));
slope_rk4 = sprintf('RK4 slope = %0.1f', log10(rk4_error_percentage(1)/rk4_error_percentage(end))/log10(H(1)/H(end)));
text(0.05, 1e-5, slope_euler)
text(0.05, 1e-4, slope_rk4)