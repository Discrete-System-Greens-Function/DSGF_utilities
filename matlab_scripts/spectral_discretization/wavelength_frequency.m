clear;
c = 3e8;
h = 6.6261e-34;
electron = 1.60217e-19;
%lambda = 0.4e-6;
%f = c/lambda

omega = 9.75e13;
f = omega/(2*pi);
%f = 1.8e16 
lambda = c/f *10^6 %[um]
E = h*c/(lambda*10^(-6)) %[J]
eV = E/electron

T = 2898/lambda 