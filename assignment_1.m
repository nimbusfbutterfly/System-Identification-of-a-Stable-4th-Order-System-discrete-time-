clc;clear all;close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = 999;
num=[0 1 7 24 24];
den=conv(conv(conv([1 2],[1 3]),[1 4]),[1 5]);
sys = tf(num, den);
sys_d=c2d(sys,0.01,'zoh')
figure;
rlocus(sys);grid on
figure;
pzmap(sys_d);grid on
disp("zeroes:");disp(zero(sys_d));disp("poles:");disp(pole(sys_d));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
num_d=sys_d.Numerator;
den_d=sys_d.Denominator;
sys = tf(num_d, den_d, 1, 'Variable', 'z');
[num, den] = tfdata(sys, 'v');
a = den(2:end);
b = num;
n = length(a);
m = length(b) - 1;
fprintf('y(t) = ');
for i = 1:n
    fprintf('%f*y(t-%d) ', -a(i), i);
    if i < n
        fprintf('+ ');
    end
end
for i = 1:m+1
    fprintf('%f*u(t-%d) ', b(i), i-1);
    if i < m+1
        fprintf('+ ');
    end
end
fprintf('\n');