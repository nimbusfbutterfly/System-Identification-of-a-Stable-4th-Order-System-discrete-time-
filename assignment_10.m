clc; clear all; close all;

a0_true = 3.207; a1_true = -3.853; a2_true = 2.055; a3_true = -0.4107;
b1_true = 0.04269; b2_true = -0.0332; b3_true = -0.03096; b4_true = 0.02414;

max = 1000;
y = zeros(max, 1);
Noise = 0.05 * wgn(max, 1, 10);
duration = 1000.0;
sampling_rate = 1;
t1 = 0:1/sampling_rate:(duration - 1/sampling_rate); 
% u = sind(2 * pi * t1) + sind(5 * pi * t1) + sind(7 * pi * t1) + sind(3 * pi * t1); % Sinusoidal input 
% u = zeros(1,max); u(5) =1; % Impulse input 
% u=randn(1,max);
% u=ones(1,max);%u(1)=0; %step input
u(1:1000, 1) = ones(1000, 1);
lambda = 1;
P = 100000 * eye(8);
theta = zeros(8, 1);
a0_values = zeros(max, 1);
a1_values = zeros(max, 1);
a2_values = zeros(max, 1);
a3_values = zeros(max, 1);
b1_values = zeros(max, 1);
b2_values = zeros(max, 1);
b3_values = zeros(max, 1);
b4_values = zeros(max, 1);
reset_interval = 200;

for t = 8:max
    if mod(t, reset_interval) == 0
        P = 100000 * eye(8);
    end

    y(t) = a0_true * y(t-1) + a1_true * y(t-2) + a2_true * y(t-3) + a3_true * y(t-4) ...
        + b1_true * u(t-1) + b2_true * u(t-2) + b3_true * u(t-3) + b4_true * u(t-4)+ Noise(t);
    phiint = [y(t-1); y(t-2); y(t-3); y(t-4); u(t-1); u(t-2); u(t-3); u(t-4)];
    K = P * phiint / (lambda + phiint' * P * phiint);
    theta = theta + K * (y(t) - phiint' * theta);
    P = (eye(8) - K * phiint') * P / lambda;
    a0_values(t) = theta(1);
    a1_values(t) = theta(2);
    a2_values(t) = theta(3);
    a3_values(t) = theta(4);
    b1_values(t) = theta(5);
    b2_values(t) = theta(6);
    b3_values(t) = theta(7);
    b4_values(t) = theta(8);
end

disp(['a0: ', num2str(a0_values(max))]);
disp(['a1: ', num2str(a1_values(max))]);
disp(['a2: ', num2str(a2_values(max))]);
disp(['a3: ', num2str(a3_values(max))]);
disp(['b1: ', num2str(b1_values(max))]);
disp(['b2: ', num2str(b2_values(max))]);
disp(['b3: ', num2str(b3_values(max))]);
disp(['b4: ', num2str(b4_values(max))]);


figure;
subplot(4, 1, 1)
plot(a0_values, 'LineWidth', 2);
hold on;
plot(1:max, a0_true * ones(size(1:max)), 'LineWidth', 2);  
xlabel('Time');
ylabel('A0');
legend('Estimated', 'True');
grid on;

subplot(4, 1, 2)
plot(a1_values, 'LineWidth', 2);
hold on;
plot(1:max, a1_true * ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel('A1');
grid on;

subplot(4, 1, 3)  
plot(a2_values, 'LineWidth', 2);
hold on;
plot(1:max, a2_true * ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel('A2');
grid on;

subplot(4, 1, 4)
plot(a3_values, 'LineWidth', 2);
hold on;
plot(1:max, a3_true * ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel('A3');
grid on;

figure;
subplot(4, 1, 1)
plot(b1_values, 'LineWidth', 2);
hold on;
plot(1:max, b1_true * ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel('B1');
legend('Estimated', 'True');
grid on;

subplot(4, 1, 2)
plot(b2_values, 'LineWidth', 2);
hold on;
plot(1:max, b2_true * ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel('B2');
grid on;

subplot(4, 1, 3)
plot(b3_values, 'LineWidth', 2);
hold on;
plot(1:max, b3_true * ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel('B3');
grid on;

subplot(4, 1, 4)
plot(b4_values, 'LineWidth', 2);
hold on;
plot(1:max, b4_true * ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel('B4');
grid on;
