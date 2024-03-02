clc; clear all; close all;

% a0 = 3.901483; a1 = -5.707852; a2 = 3.711206; a3 = -0.904837;
% b1 = 0.009853; b2 = -0.028870; b3 = 0.028204; b4 = -0.009187;
a0=3.207; a1=-3.853; a2=2.055; a3=-0.4107;
b1=0.04269; b2=-0.0332; b3=-0.03096; b4=0.02414;
max = 10000;
Noise=0.7*wgn(max,1,1);
y = zeros(max, 1);
duration = 1000.0;
sampling_rate =1;
t1 = 0:1/sampling_rate:duration;
u = sind(2 * pi * t1)+sind(5 * pi * t1)+sind(7 * pi * t1)+sind(3 * pi *t1); %sin wave input 
% u = zeros(1,max); u(5) =1; % Impulse input 
% u=randn(1,max);
% u=ones(1,max);%u(1)=0; %step input
% u(1:50, 1) = ones(50, 1);
T = 100;
for i = 101:max
    u(i) = u(i - T);
end

for t = 8:max
        y(t) = a0 * y(t-1) + a1 * y(t - 2) + a2 * y(t - 3) + a3 * y(t - 4) ...
         + b1 * u(t - 1) + b2 * u(t - 2) + b3 * u(t - 3) + b4 * u(t - 4)+ Noise(t);%-0.5 * Noise(t - 1);
end
mean(Noise);
P = 100 * eye(8);
theta = zeros(8, 1);
a0_values(1) =theta(1);
a1_values(1) =theta(2);
a2_values(1) = theta(3);
a3_values(1) =theta(4);
b1_values(1) =theta(5);
b2_values(1) =theta(6);
b3_values(1) =theta(7);
b4_values(1) =theta(8);


forgetting_factors = [1, 0.999, 0.99,1.1, 0.95];
landa=forgetting_factors(1);
 
    % RLS Algorithm with forgetting factor landa
   for i = 8:max
        phiint = [y(i - 1);y(i - 2);y(i - 3);y(i - 4); u(i - 1);u(i - 2);u(i - 3);u(i - 4)];%for ELS +Noise(t-1)
        K = P * phiint * inv(landa * eye(1) + phiint' * P * phiint);
        theta = theta + K * (y(i) - phiint' * theta);
        P = (eye(8) - K * phiint') * P / landa;
        traceP(t)=trace(P);
        a0_values(i) = theta(1);
        a1_values(i) = theta(2);
        a2_values(i) = theta(3);
        a3_values(i) = theta(4);
        b1_values(i) = theta(5);
        b2_values(i) = theta(6);
        b3_values(i) = theta(7);
        b4_values(i) = theta(8);
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
subplot(4,1,1)
plot(a0_values, 'LineWidth', 2);
hold on;
plot(1:max, a0*ones(size(1:max)), 'LineWidth', 2);  
xlabel('Time');
ylabel(['A0']);
legend("estimated","real")
grid on;

subplot(4,1,2)
plot(a1_values, 'LineWidth', 2);
hold on;
plot(1:max, a1*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel(['A1']);
grid on;

subplot(4,1,3)  
plot(a2_values, 'LineWidth', 2);
hold on;
plot(1:max, a2*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel(['A2']);
grid on;

subplot(4,1,4)
plot(a3_values, 'LineWidth', 2);
hold on;
plot(1:max, a3*ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel(['A3']);
grid on;

figure;
subplot(4,1,1)
plot(b1_values, 'LineWidth', 2);
hold on;
plot(1:max, b1*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel('B1');
legend('Estimated', 'Real');
grid on;


subplot(4,1,2)
plot(b2_values, 'LineWidth', 2);
hold on;
plot(1:max, b2*ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel(['B2']);
grid on;

subplot(4,1,3)
plot(b3_values, 'LineWidth', 2);
hold on;
plot(1:max, b3*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel(['B3']);
grid on;

subplot(4,1,4)
plot(b4_values, 'LineWidth', 2);
hold on;
plot(1:max, b4*ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel(['B4']);
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear P;
gama=1;
alpha=0.01;
for i=8:max
    phiint=[y(i-1);y(i-2);y(i-3);y(i-4);u(i-1);u(i-2);u(i-3);u(i-4)];
    P=gama/(alpha+transpose(phiint)*phiint);
%   P=inv(phiint'*phiint);
    theta=theta+P*phiint*(y(i)-transpose(phiint)*theta);
    a0e(i)=theta(1);
    a1e(i)=theta(2);
    a2e(i)=theta(3);
    a3e(i)=theta(4);
    b1e(i)=theta(5);
    b2e(i)=theta(6);
    b3e(i)=theta(7);
    b4e(i)=theta(8);
end
    disp("PA:")
    disp(['a0: ', num2str(a0e(max))]);
    disp(['a1: ', num2str(a1e(max))]);
    disp(['a2: ', num2str(a2e(max))]);
    disp(['a3: ', num2str(a3e(max))]);
    disp(['b1: ', num2str(b1e(max))]);
    disp(['b2: ', num2str(b2e(max))]);
    disp(['b3: ', num2str(b3e(max))]);
    disp(['b4: ', num2str(b4e(max))]);

figure;
subplot(4,1,1)
plot(a0e, 'LineWidth', 2);
hold on;
plot(1:max, a0*ones(size(1:max)), 'LineWidth', 2);  
xlabel('Time');
ylabel(['A0']);
legend("estimated","real")
grid on;

subplot(4,1,2)
plot(a1e, 'LineWidth', 2);
hold on;
plot(1:max, a1*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel(['A1']);
grid on;

subplot(4,1,3)  
plot(a2e, 'LineWidth', 2);
hold on;
plot(1:max, a2*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel(['A2']);
grid on;

subplot(4,1,4)
plot(a3e, 'LineWidth', 2);
hold on;
plot(1:max, a3*ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel(['A3']);
grid on;

figure;
subplot(4,1,1)
plot(b1e, 'LineWidth', 2);
hold on;
plot(1:max, b1*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel('B1');
legend('Estimated', 'Real');
grid on;


subplot(4,1,2)
plot(b2e, 'LineWidth', 2);
hold on;
plot(1:max, b2*ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel(['B2']);
grid on;

subplot(4,1,3)
plot(b3e, 'LineWidth', 2);
hold on;
plot(1:max, b3*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel(['B3']);
grid on;

subplot(4,1,4)
plot(b4e, 'LineWidth', 2);
hold on;
plot(1:max, b4*ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel(['B4']);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%LMS

clear P;
gama=1e-5;
P=gama;
for i=8:max
    phiint=[y(i-1);y(i-2);y(i-3);y(i-4);u(i-1);u(i-2);u(i-3);u(i-4)];
    P=gama;
    theta=theta+P*phiint*(y(i)-phiint'*theta);
    a0e(i)=theta(1);
    a1e(i)=theta(2);
    a2e(i)=theta(3);
    a3e(i)=theta(4);
    b1e(i)=theta(5);
    b2e(i)=theta(6);
    b3e(i)=theta(7);
    b4e(i)=theta(8);
end
    disp("LMS:")
    disp(['a0: ', num2str(a0e(max))]);
    disp(['a1: ', num2str(a1e(max))]);
    disp(['a2: ', num2str(a2e(max))]);
    disp(['a3: ', num2str(a3e(max))]);
    disp(['b1: ', num2str(b1e(max))]);
    disp(['b2: ', num2str(b2e(max))]);
    disp(['b3: ', num2str(b3e(max))]);
    disp(['b4: ', num2str(b4e(max))]);

figure;
subplot(4,1,1)
plot(a0e, 'LineWidth', 2);
hold on;
plot(1:max, a0*ones(size(1:max)), 'LineWidth', 2);  
xlabel('Time');
ylabel(['A0']);
legend("estimated","real")
grid on;

subplot(4,1,2)
plot(a1e, 'LineWidth', 2);
hold on;
plot(1:max, a1*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel(['A1']);
grid on;

subplot(4,1,3)  
plot(a2e, 'LineWidth', 2);
hold on;
plot(1:max, a2*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel(['A2']);
grid on;

subplot(4,1,4)
plot(a3e, 'LineWidth', 2);
hold on;
plot(1:max, a3*ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel(['A3']);
grid on;

figure;
subplot(4,1,1)
plot(b1e, 'LineWidth', 2);
hold on;
plot(1:max, b1*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel('B1');
legend('Estimated', 'Real');
grid on;


subplot(4,1,2)
plot(b2e, 'LineWidth', 2);
hold on;
plot(1:max, b2*ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel(['B2']);
grid on;

subplot(4,1,3)
plot(b3e, 'LineWidth', 2);
hold on;
plot(1:max, b3*ones(size(1:max)), 'LineWidth', 2); 
xlabel('Time');
ylabel(['B3']);
grid on;

subplot(4,1,4)
plot(b4e, 'LineWidth', 2);
hold on;
plot(1:max, b4*ones(size(1:max)), 'LineWidth', 2);
xlabel('Time');
ylabel(['B4']);
grid on;

