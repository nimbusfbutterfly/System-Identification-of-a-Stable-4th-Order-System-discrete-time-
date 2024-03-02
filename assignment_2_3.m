clc;clear all;close all;
a0 = 3.901483; a1 = -5.707852; a2 = 3.711206; a3 = -0.904837;
b1 = 0.009853; b2 = -0.028870; b3 = 0.028204; b4 = -0.009187;
A = [3.901483, -5.707852, 3.711206, -0.904837];
B = [0.009853, -0.028870, 0.028204, -0.009187];
% sigm_e=0.1;
N=1000;

y=zeros(1,N);
t = 1:N;
duration = 1000.0;
sampling_rate =1;
t1 = 0:1/sampling_rate:duration;
u = sind(2 * pi * t1)+sind(5 * pi * t1)+sind(7 * pi * t1)+sind(3 * pi * t1); %sin wave input
% u = zeros(1,N); u(5) =1; % Impulse input
% u= randn(1,N);
% u=ones(1,N);u(1:4)=0; %step input
% ep=wgn(1,N,0);
% sigm_ep=1/(N)*sum(ep.^2);
% e=sqrt(sigm_e/(sigm_ep))*ep;
order = 4;
for t = order + 1:N
    y(t) = a0 * y(t-1) + a1 * y(t - 2) + a2 * y(t - 3) + a3 * y(t - 4) ...
         + b1 * u(t - 1) + b2 * u(t - 2) + b3 * u(t - 3) + b4 * u(t - 4);
end
y1=y;
figure;
subplot(2,1,1); stairs(1:1:N,"b","LineWidth",2);grid on;ylabel("y");
subplot(2,1,2);stairs(1:N,u(1:N),"b","LineWidth",2);grid on;ylabel("u");xlabel("Time");

%Least Squares
phi=zeros(N,2*order);
Y=zeros(N,1);

a_0 = zeros(1, N);
a_1 = zeros(1, N);
a_2 = zeros(1, N);
a_3 = zeros(1, N);
b_1 = zeros(1, N);
b_2 = zeros(1, N);
b_3 = zeros(1, N);
b_4 = zeros(1, N);

for t = order + 1:N
    phi(t,:) = [y(t-1), y(t-2), y(t-3), y(t-4), u(t-1), u(t-2), u(t-3), u(t-4)];
    Y(t) = y(t);
    theta = (phi(1:t, :)' * phi(1:t, :)) \ (phi(1:t, :)' * Y(1:t));
    a_0(t) = theta(1);
    a_1(t) = theta(2);
    a_2(t) = theta(3);
    a_3(t) = theta(4);
    b_1(t) = theta(5);
    b_2(t) = theta(6);
    b_3(t) = theta(7);
    b_4(t) = theta(8);
end

y_hat = phi* theta;
j = 1/2*sum((Y - y_hat).^2);
disp(['Cost function value: ' num2str(j)]);
disp("a0:");disp(theta(1));
disp("a1:");disp(theta(2));
disp("a2:");disp(theta(3));
disp("a3:");disp(theta(4));
disp("b1:");disp(theta(5));
disp("b2:");disp(theta(6));
disp("b3:");disp(theta(7));
disp("b4:");disp(theta(8));

figure;
subplot(4,1,1)
plot(a_0, 'LineWidth', 2);
hold on;
plot(1:N, A(1)*ones(size(1:N)), 'LineWidth', 2);  
xlabel('Time');
ylabel(['A0']);
legend("estimated","real")
axis([1 N -5 5])
grid on;

subplot(4,1,2)
plot(a_1, 'LineWidth', 2);
hold on;
plot(1:N, A(2)*ones(size(1:N)), 'LineWidth', 2); 
xlabel('Time');
ylabel(['A1']);
axis([1 N -7 7])
grid on;

subplot(4,1,3)  
plot(a_2, 'LineWidth', 2);
hold on;
plot(1:N, A(3)*ones(size(1:N)), 'LineWidth', 2);  
xlabel('Time');
ylabel(['A2']);
axis([1 N -5 5])
grid on;

subplot(4,1,4)
plot(a_3, 'LineWidth', 2);
hold on;
plot(1:N, A(4)*ones(size(1:N)), 'LineWidth', 2);
xlabel('Time');
ylabel(['A3']);
axis([1 N -1 1])
ylim([-1 1])
grid on;

figure;
subplot(4,1,1)
plot(b_1, 'LineWidth', 2);
hold on;
plot(1:N, B(1)*ones(size(1:N)), 'LineWidth', 2); 
xlabel('Time');
ylabel('B1');
legend('Estimated', 'Real');
axis([1 N 0 0.02]);
grid on;


subplot(4,1,2)
plot(b_2, 'LineWidth', 2);
hold on;
plot(1:N, B(2)*ones(size(1:N)), 'LineWidth', 2);
xlabel('Time');
ylabel(['B2']);
axis([1 N -0.05 0.05])
grid on;

subplot(4,1,3)
plot(b_3, 'LineWidth', 2);
hold on;
plot(1:N, B(3)*ones(size(1:N)), 'LineWidth', 2); 
xlabel('Time');
ylabel(['B3']);
axis([1 N -0.05 0.05])
grid on;

subplot(4,1,4)
plot(b_4, 'LineWidth', 2);
hold on;
plot(1:N, B(4)*ones(size(1:N)), 'LineWidth', 2);
xlabel('Time');
ylabel(['B4']);
axis([1 N -0.05 0.05])
grid on;