% QN 1
% Parameters for synthetic signal to try
T = 2*pi; % Total time interval
dt = 0.001; % Base sampling step
omega = 2.0; % Frequency
phi = pi/6; % Phase
sigma = 0.01; % Noise standard deviation

% Generate synthetic data y = sin(omega*t + phi)
t = 0:dt:T;
y = sin(omega*t + phi);
y_prime_true = omega*cos(omega*t + phi);

% Add noise
noise = normrnd(0, sigma, [1, length(t)]);
y_noise = y + noise;
% QN 1

% QN 2
K = 2.^(0:7);
forward_diff = zeros(1, 8);
forward_diff_noise = zeros(1, 8);
backward_diff = zeros(1, 8);
backward_diff_noise = zeros(1, 8);
central_diff = zeros(1, 8);
central_diff_noise = zeros(1, 8);
i = 1;
% Forward difference scheme y'_n = [f((n+1)h) - f(nh)]/h
for k = K
    h = k*dt;
    t_sampled = 0:h:T;
    y_noise_sampled = y_noise(1:k:length(y));
    y_sampled = y(1:k:length(y));
    y_prime_true_sampled = y_prime_true(1:k:length(y));

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
    fwd_rms = sqrt(mean_squared);
    sum_of_squares_noise = sum(squared_error_noise);
    mean_squared_noise = sum_of_squares_noise/length(t_sampled_for_deriv);
    fwd_rms_noise = sqrt(mean_squared_noise);

    i = i + 1;
end

% Backward difference scheme y'_n = [f(nh) - f((n-1)h)]/h
i = 1;
for k = K
    h = k*dt;
    t_sampled = 0:h:T;
    y_noise_sampled = y_noise(1:k:length(y));
    y_sampled = y(1:k:length(y));
    y_prime_true_sampled = y_prime_true(1:k:length(y));

    t_sampled_for_deriv = t_sampled(2:end);
    y_prime = zeros(1, length(t_sampled_for_deriv));
    y_prime_noise = zeros(1, length(t_sampled_for_deriv));

    % Compute derivative at all points
    % Simultaneously, compute squared error at each point
    squared_error = zeros(1, length(t_sampled_for_deriv));
    squared_error_noise = zeros(1, length(t_sampled_for_deriv));
    for n = 2:length(y_prime_noise)+1
        y_prime(n-1) = (y_sampled(n) - y_sampled(n-1))./h;
        squared_error(n-1) = (y_prime(n-1) - y_prime_true_sampled(n-1)).^2;
        y_prime_noise(n-1) = (y_noise_sampled(n) - y_noise_sampled(n-1))./h;
        squared_error_noise(n-1) = (y_prime_noise(n-1) - y_prime_true_sampled(n-1)).^2;
    end

    % Compute L2/RMS error
    sum_of_squares = sum(squared_error);
    mean_squared = sum_of_squares/length(t_sampled_for_deriv);
    bwd_rms = sqrt(mean_squared);
    sum_of_squares_noise = sum(squared_error_noise);
    mean_squared_noise = sum_of_squares_noise/length(t_sampled_for_deriv);
    bwd_rms_noise = sqrt(mean_squared_noise);

    i = i + 1;
end

% Central difference scheme y'_n = [f((n+1)h) - f((n-1)h)]/2h
i = 1;
for k = K
    h = k*dt;
    t_sampled = 0:h:T;
    y_noise_sampled = y_noise(1:k:length(y));
    y_sampled = y(1:k:length(y));
    y_prime_true_sampled = y_prime_true(1:k:length(y));

    t_sampled_for_deriv = t_sampled(2:end-1);
    y_prime = zeros(1, length(t_sampled_for_deriv));
    y_prime_noise = zeros(1, length(t_sampled_for_deriv));

    % Compute derivative at all points
    % Simultaneously, compute squared error at each point
    squared_error = zeros(1, length(t_sampled_for_deriv));
    squared_error_noise = zeros(1, length(t_sampled_for_deriv));
    for n = 2:length(y_prime_noise)+1
        y_prime(n-1) = (y_sampled(n+1) - y_sampled(n-1))./(2*h);
        squared_error(n-1) = (y_prime(n-1) - y_prime_true_sampled(n-1)).^2;
        y_prime_noise(n-1) = (y_noise_sampled(n+1) - y_noise_sampled(n-1))./(2*h);
        squared_error_noise(n-1) = (y_prime_noise(n-1) - y_prime_true_sampled(n-1)).^2;
    end
    
    % Compute L2/RMS error
    sum_of_squares = sum(squared_error);
    mean_squared = sum_of_squares/length(t_sampled_for_deriv);
    cen_rms = sqrt(mean_squared);
    sum_of_squares_noise = sum(squared_error_noise);
    mean_squared_noise = sum_of_squares_noise/length(t_sampled_for_deriv);
    cen_rms_noise = sqrt(mean_squared_noise);

    i = i + 1;
end

figure
tiledlayout(2, 1)
nexttile
loglog(K*dt, fwd_rms, '-o', K*dt, bwd_rms, '-x', K*dt, cen_rms, '-o')
legend('forward', 'backward', 'central')
xlabel('h')
ylabel('difference')
nexttile
loglog(K*dt, fwd_rms_noise, '-o', K*dt, bwd_rms_noise, '-x', K*dt, cen_rms_noise, '-o')
legend('forward', 'backward', 'central')
xlabel('h')
ylabel('difference noise')
% QN 2