clear;
clc;

function [fit, condition, alpha_fed] = regress(x, y, M, condition)
    % Fits and plots Mth order polynomial to data points (x, y)

    % There is a MATLAB function vander(v) which automatically returns a
    % Vandermonde matrix based on x-points 'v' but that is only for square
    % Vandermondes (thus order of polynomial is fixed to length(x) - 1).
    % Here we allow the user to choose the order of polynomial themselves.

    V = zeros(length(x), M + 1);
    % Populate Vandermonde matrix by taking integer powers of x along the
    % columns
    for i = 0:M
        V(:, i+1) = x.^i;
    end

    V_inv = inv(V'*V)*V';   % Pseudoinverse of Vandermonde

    condition(M) = cond(V'*V);  % Record the condition number for this order

    % Compute polynomial coefficients
    % y is by default a row vector. Turn it into column vector
    alpha = V_inv*y;

    % Find the y-values of the fitted polynomial
    fit = V*alpha;

    % Return the polynomial coefficients in a consistent format
    % Leftmost coefficient refers to highest order term
    alpha_fed = zeros(1, 11);
    alpha_fed(1, end-M:end) = alpha;
end

function [fit, condition] = plt(x, alpha)
    % Plots a polynomial using Vandermonde given the coefficients of the
    % polynomial. Leftmost coefficient is highest order term.

    M = length(alpha) - 1;
    V = zeros(length(x), M + 1);
    % Populate Vandermonde matrix by taking integer powers of x along the
    % columns
    for i = 0:M
        V(:, i+1) = x.^i;
    end

    condition = cond(V'*V);  % Record the condition number for this order

    % Find the y-values of the fitted polynomial
    fit = V*alpha;
end

% Central difference scheme y'_n = [f((n+1)h) - f((n-1)h)]/2h
function central_diff_error_noise = central(h, t, y_noise, y_prime)
    t_sampled_for_deriv = t(2:end-1);
    y_prime_noise = zeros(1, length(t_sampled_for_deriv));

    % Compute derivative at all points
    % Simultaneously, compute squared error at each point
    squared_error_noise = zeros(1, length(t_sampled_for_deriv));
    for n = 2:length(y_prime_noise)+1
        y_prime_noise(n-1) = (y_noise(n+1) - y_noise(n-1))./(2*h);
        squared_error_noise(n-1) = (y_prime_noise(n-1) - y_prime(n)).^2;
    end

    % Compare the central difference scheme with the analytical solution
    figure
    plot(t_sampled_for_deriv, y_prime_noise, 'r', t, y_prime, 'b')
    legend('central diff', 'analytical')
    xlabel('time')
    ylabel("y'")
    
    % Compute L2/RMS error
    central_diff_error_noise = rms_error(squared_error_noise, t_sampled_for_deriv);
end

function error = rms_error(squared_error, samples)
    % Computes L2/RMS error given the squared error
    sum_of_squares = sum(squared_error);
    mean_squared = sum_of_squares/length(samples);
    error= sqrt(mean_squared);
end

% Generate synthetic data y = sin(omega*t + phi)
t = [-pi:0.1:pi]';
omega = 2.0;
phi = 0.6;
sigma = 0.01;
y = sin(omega*t + phi);
y_prime = omega*cos(omega*t + phi);
noise = normrnd(0, sigma, [1, length(t)]);
y_noise = y + noise';

% Initialise array to store polynomial fit data
M = 10;     % max polynomial order
fits = zeros(M, length(t));
condition = zeros(1, M);
alphas = zeros(M, M+1);

figure
tiledlayout(5, 2)
% Get fitted polynomial, condition number, and polynomial coefficients for
% each proposed polynomial degree
for m = 1:M
    [fits(m, :), condition, alphas(m, :)] = regress(t, y_noise, m, condition);

    % Compare the fitted polynomial with the original noisy y curve
    nexttile
    plot(t, fits(m,:), 'r', t , y_noise, 'b')
    legend('fitted poly', 'original noisy y')
    xlabel("t")
    ylabel("y")
end

% Plot condition number vs polynomial degree, log vertical axis
figure
semilogy(1:M, condition, '-x')
xlabel("polynomial degree m")
ylabel("condition number (log-axis)")

% Initialise arrays for finding 
diff_alphas = zeros(M, M);
diff_fits = zeros(M, length(t));
diff_condition = zeros(1, M);
squared_error_fitted = zeros(M, length(t));
diff_fitted_poly_error = zeros(1, M);

figure
tiledlayout(5, 2)

for m = 1:M
    % Differentiate the fitted polynomials
    diff_alphas(:, m) = alphas(:, m) .* (10-m+1);
    
    % Turn the alpha matrix into a form usable for V*alpha above
    use_diff_alphas = nonzeros(diff_alphas(:, m));
    
    % Find the fitted polynomial derivatives and their condition numbers
    [diff_fits(m, :), diff_condition(m)] = plt(t, use_diff_alphas);

    % Find the RMS/L2 Error
    for i = 1:length(t)
        squared_error_fitted(m, i) = (diff_fits(m, i) - y_prime(i)).^2;
    end
    diff_fitted_poly_error(1, m) = rms_error(squared_error_fitted(m, 15:end-15), t(15:end-15));

    % Compare derivative retrieved by polynomial fitting and the analytical
    % derivative
    nexttile
    plot(t, diff_fits(m, :), 'r', t, y_prime, 'b')
    legend('via poly fitting', 'analytical soln')
    xlabel("t")
    ylabel("y'")
end

central_diff_error = central(0.1, t, y_noise, y_prime);

figure
plot(1:M, diff_fitted_poly_error, '-x')
xlabel("polynomial degree m")
ylabel("RMS Error")
