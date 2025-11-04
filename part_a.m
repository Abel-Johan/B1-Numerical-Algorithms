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

    t_sampled_for_deriv = t_sampled(1:end-1);
    y_prime = zeros(1, length(t_sampled_for_deriv));
    y_prime_noise = zeros(1, length(t_sampled_for_deriv));

    % tiledlayout(4,1)
    for n = 1:length(y_prime_noise)
        y_prime(n) = (y_sampled(n+1) - y_sampled(n))./h;
        y_prime_noise(n) = (y_noise_sampled(n+1) - y_noise_sampled(n))./h;
    end

    % Compute error on 129th data point because its accessible for all k's
    exact_deriv = y_prime_true(129);
    sim_deriv = y_prime(((129-1)/k)+1);
    sim_noise_deriv = y_prime_noise(((129-1)/k)+1);
    forward_diff(i) = abs(exact_deriv - sim_deriv);
    forward_diff_noise(i) = abs(exact_deriv - sim_noise_deriv);

    % nexttile
    % plot(t, y)
    % ylabel("sin(\omegat + \phi)")
    % 
    % nexttile
    % plot(t_sampled, y_noise_sampled)
    % ylabel("noisy sin(\omegat + \phi)")
    % 
    % nexttile
    % plot(t_sampled_for_deriv, y_prime)
    % ylabel("\omegacos(\omegat + \phi)")
    % 
    % nexttile
    % plot(t_sampled_for_deriv, y_prime_noise)
    % ylabel("noisy \omegacos(\omegat + \phi)")

    i = i + 1;
end

% Backward difference scheme y'_n = [f(nh) - f((n-1)h)]/h
i = 1;
for k = K
    h = k*dt;
    t_sampled = 0:h:T;
    y_noise_sampled = y_noise(1:k:length(y));
    y_sampled = y(1:k:length(y));

    t_sampled_for_deriv = t_sampled(2:end);
    y_prime = zeros(1, length(t_sampled_for_deriv));
    y_prime_noise = zeros(1, length(t_sampled_for_deriv));

    % tiledlayout(4,1)
    for n = 2:length(y_prime_noise)+1
        y_prime(n-1) = (y_sampled(n) - y_sampled(n-1))./h;
        y_prime_noise(n-1) = (y_noise_sampled(n) - y_noise_sampled(n-1))./h;
    end

    % Compute error on 129th data point because its accessible for all k's
    exact_deriv = y_prime_true(129);
    sim_deriv = y_prime(((129-1)/k));
    sim_noise_deriv = y_prime_noise(((129-1)/k));
    backward_diff(i) = abs(exact_deriv - sim_deriv);
    backward_diff_noise(i) = abs(exact_deriv - sim_noise_deriv);

    % nexttile
    % plot(t, y)
    % ylabel("sin(\omegat + \phi)")
    % 
    % nexttile
    % plot(t_sampled, y_noise_sampled)
    % ylabel("noisy sin(\omegat + \phi)")
    % 
    % nexttile
    % plot(t_sampled_for_deriv, y_prime)
    % ylabel("\omegacos(\omegat + \phi)")
    % 
    % nexttile
    % plot(t_sampled_for_deriv, y_prime_noise)
    % ylabel("noisy \omegacos(\omegat + \phi)")

    i = i + 1;
end

% Central difference scheme y'_n = [f((n+1)h) - f((n-1)h)]/2h
i = 1;
for k = K
    h = k*dt;
    t_sampled = 0:h:T;
    y_noise_sampled = y_noise(1:k:length(y));
    y_sampled = y(1:k:length(y));

    t_sampled_for_deriv = t_sampled(2:end-1);
    y_prime = zeros(1, length(t_sampled_for_deriv));
    y_prime_noise = zeros(1, length(t_sampled_for_deriv));

    % tiledlayout(4,1)
    for n = 2:length(y_prime_noise)+1
        y_prime(n-1) = (y_sampled(n+1) - y_sampled(n-1))./(2*h);
        y_prime_noise(n-1) = (y_noise_sampled(n+1) - y_noise_sampled(n-1))./(2*h);
    end
    
    % Compute error on 129th data point because its accessible for all k's
    exact_deriv = y_prime_true(129);
    sim_deriv = y_prime(((129-1)/k));
    sim_noise_deriv = y_prime_noise(((129-1)/k));
    central_diff(i) = abs(exact_deriv - sim_deriv);
    central_diff_noise(i) = abs(exact_deriv - sim_noise_deriv);

    % nexttile
    % plot(t, y)
    % ylabel("sin(\omegat + \phi)")
    % 
    % nexttile
    % plot(t_sampled, y_noise_sampled)
    % ylabel("noisy sin(\omegat + \phi)")
    % 
    % nexttile
    % plot(t_sampled_for_deriv, y_prime)
    % ylabel("\omegacos(\omegat + \phi)")
    % 
    % nexttile
    % plot(t_sampled_for_deriv, y_prime_noise)
    % ylabel("noisy \omegacos(\omegat + \phi)")

    i = i + 1;
end

figure
tiledlayout(2, 1)
nexttile
loglog(K*dt, forward_diff, '-o', K*dt, backward_diff, '-x', K*dt, central_diff, '-o')
legend('forward', 'backward', 'central')
xlabel('h')
ylabel('difference')
nexttile
loglog(K*dt, forward_diff_noise, '-o', K*dt, backward_diff_noise, '-x', K*dt, central_diff_noise, '-o')
legend('forward', 'backward', 'central')
xlabel('h')
ylabel('difference noise')
% QN 2