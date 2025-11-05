% Parameters for synthetic signal to try
T = 2*pi; % Total time interval
dt = 0.001; % Base sampling step
omega = 2.0; % Frequency
phi = pi/6; % Phase
sigma = 0.01; % Noise standard deviation
K = 2.^(0:7); % Step size multipliers

% Generate synthetic data y = sin(omega*t + phi)
t = 0:dt:T;
y = sin(omega*t + phi);
y_prime_true = omega*cos(omega*t + phi);

% Add noise
noise = normrnd(0, sigma, [1, length(t)]);
y_noise = y + noise;

% Arrays to store L2 errors for different schemes, for each step size
forward_diff = zeros(1, length(K));
forward_diff_noise = zeros(1, length(K));
backward_diff = zeros(1, length(K));
backward_diff_noise = zeros(1, length(K));
central_diff = zeros(1, length(K));
central_diff_noise = zeros(1, length(K));

% Forward difference scheme y'_n = [f((n+1)h) - f(nh)]/h
function [forward_diff_error, forward_diff_error_noise] = forward(h, t_sampled, y_sampled, y_noise_sampled, y_prime_true_sampled)
    t_sampled_for_deriv = t_sampled(1:end-1);
    y_prime = zeros(1, length(t_sampled_for_deriv));
    y_prime_noise = zeros(1, length(t_sampled_for_deriv));

    % Compute derivative at all points
    % Simultaneously, compute squared error at each point
    squared_error = zeros(1, length(t_sampled_for_deriv));
    squared_error_noise = zeros(1, length(t_sampled_for_deriv));
    for n = 1:length(y_prime_noise)
        y_prime(n) = (y_sampled(n+1) - y_sampled(n))./h;
        squared_error(n) = (y_prime(n) - y_prime_true_sampled(n)).^2;
        y_prime_noise(n) = (y_noise_sampled(n+1) - y_noise_sampled(n))./h;
        squared_error_noise(n) = (y_prime_noise(n) - y_prime_true_sampled(n)).^2;
    end

    % Compute L2/RMS error
    sum_of_squares = sum(squared_error);
    mean_squared = sum_of_squares/length(t_sampled_for_deriv);
    forward_diff_error = sqrt(mean_squared);
    sum_of_squares_noise = sum(squared_error_noise);
    mean_squared_noise = sum_of_squares_noise/length(t_sampled_for_deriv);
    forward_diff_error_noise = sqrt(mean_squared_noise);
end

% Backward difference scheme y'_n = [f(nh) - f((n-1)h)]/h
function [backward_diff_error, backward_diff_error_noise] = backward(h, t_sampled, y_sampled, y_noise_sampled, y_prime_true_sampled)
    t_sampled_for_deriv = t_sampled(2:end);
    y_prime = zeros(1, length(t_sampled_for_deriv));
    y_prime_noise = zeros(1, length(t_sampled_for_deriv));

    % Compute derivative at all points
    % Simultaneously, compute squared error at each point
    squared_error = zeros(1, length(t_sampled_for_deriv));
    squared_error_noise = zeros(1, length(t_sampled_for_deriv));
    for n = 2:length(y_prime_noise)+1
        y_prime(n-1) = (y_sampled(n) - y_sampled(n-1))./h;
        squared_error(n-1) = (y_prime(n-1) - y_prime_true_sampled(n)).^2;
        y_prime_noise(n-1) = (y_noise_sampled(n) - y_noise_sampled(n-1))./h;
        squared_error_noise(n-1) = (y_prime_noise(n-1) - y_prime_true_sampled(n)).^2;
    end

    % Compute L2/RMS error
    sum_of_squares = sum(squared_error);
    mean_squared = sum_of_squares/length(t_sampled_for_deriv);
    backward_diff_error = sqrt(mean_squared);
    sum_of_squares_noise = sum(squared_error_noise);
    mean_squared_noise = sum_of_squares_noise/length(t_sampled_for_deriv);
    backward_diff_error_noise = sqrt(mean_squared_noise);
end

% Central difference scheme y'_n = [f((n+1)h) - f((n-1)h)]/2h
function [central_diff_error, central_diff_error_noise] = central(h, t_sampled, y_sampled, y_noise_sampled, y_prime_true_sampled)
    t_sampled_for_deriv = t_sampled(2:end-1);
    y_prime = zeros(1, length(t_sampled_for_deriv));
    y_prime_noise = zeros(1, length(t_sampled_for_deriv));

    % Compute derivative at all points
    % Simultaneously, compute squared error at each point
    squared_error = zeros(1, length(t_sampled_for_deriv));
    squared_error_noise = zeros(1, length(t_sampled_for_deriv));
    for n = 2:length(y_prime_noise)+1
        y_prime(n-1) = (y_sampled(n+1) - y_sampled(n-1))./(2*h);
        squared_error(n-1) = (y_prime(n-1) - y_prime_true_sampled(n)).^2;
        y_prime_noise(n-1) = (y_noise_sampled(n+1) - y_noise_sampled(n-1))./(2*h);
        squared_error_noise(n-1) = (y_prime_noise(n-1) - y_prime_true_sampled(n)).^2;
    end
    
    % Compute L2/RMS error
    sum_of_squares = sum(squared_error);
    mean_squared = sum_of_squares/length(t_sampled_for_deriv);
    central_diff_error = sqrt(mean_squared);
    sum_of_squares_noise = sum(squared_error_noise);
    mean_squared_noise = sum_of_squares_noise/length(t_sampled_for_deriv);
    central_diff_error_noise = sqrt(mean_squared_noise);
end

% Loop for all step size multipliers
for i = 1:length(K)
    k = K(i);
    h = k*dt;
    t_sampled = 0:h:T;
    y_sampled = y(1:k:length(y));
    y_noise_sampled = y_noise(1:k:length(y));
    y_prime_true_sampled = y_prime_true(1:k:length(y));

    [forward_diff(i), forward_diff_noise(i)] = forward(h, t_sampled, y_sampled, y_noise_sampled, y_prime_true_sampled);
    [backward_diff(i), backward_diff_noise(i)] = backward(h, t_sampled, y_sampled, y_noise_sampled, y_prime_true_sampled);
    [central_diff(i), central_diff_noise(i)] = central(h, t_sampled, y_sampled, y_noise_sampled, y_prime_true_sampled);
end

figure
% tiledlayout(2, 1)
% nexttile
% loglog(K*dt, forward_diff, '-o', K*dt, backward_diff, '-x', K*dt, central_diff, '-o')
% legend('forward', 'backward', 'central')
% xlabel('h')
% ylabel('difference without noise')
nexttile
loglog(K*dt, forward_diff_noise, '-o', K*dt, backward_diff_noise, '-x', K*dt, central_diff_noise, '-o')
legend('forward', 'backward', 'central')
xlabel('h')
ylabel('difference noise')