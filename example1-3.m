clear, close all, clc

% First we generate our signal
low_freq = 20;
high_freq = 100;
n = 50*high_freq; % Total length of the signal
t = linspace(-1,1,n);
x = sin(2*pi*low_freq*t)+2*sin(2*pi*high_freq*t);
x = x';
%plot(t,x);
%xlim([-0.1,0.1]);

m = 100; % Number of samples
r = randperm(length(t), m);
Phi = zeros(m,n);
for ii = 1:m
    ii;
    ek = zeros(1,n);
    ek(r(ii)) = 1;
    ek;
    Phi(ii,:) = ek;
end

% We take the samples
y = Phi*x;

% Generating the sensing matrix
Theta = zeros(m,n);
for ii = 1:n
    ii;
    ek = zeros(1,n);
    ek(ii) = 1;
    psi = idct(ek)';
    Theta(:,ii) = Phi*psi;
end

s2 = pinv(Theta)*y; % As approximate solution we calculate the minimum 
% energy solution using the pseudoinverse of Theta

% 
s1 = l1eq_pd(s2,Theta,Theta',y,1e-3, 25); % L1-magic toolbox

x1 = zeros(n,1);
for ii = 1:n
    ii;
    ek = zeros(1,n);
    ek(ii) = 1;
    psi = idct(ek)';
    x1 = x1+psi*s1(ii);
end

x0=10;
y0=10;
width=650;
height=400;
f1 = figure('name', 'Original signal and samples');
figure(f1)
set(gcf,'position',[x0,y0,width,height]);
plot(x, 'Linewidth', 1), hold on
scatter(r,y','x', 'Linewidth',1.5),
xlim([2250,2750]),
xticks([2250, 2375, 2500, 2625, 2750])
xticklabels({'-0.1','-0.05','0', '0.05', '0.1'});
xlabel('time (s)')

f2 = figure('name', 'Reconstructed signal');
figure(f2)
set(gcf,'position',[x0,y0,width,height]);
plot(x1','Linewidth', 1)
xlim([2250,2750]),
ylim([-3,3])
xticks([2250, 2375, 2500, 2625, 2750])
xticklabels({'-0.1','-0.05','0', '0.05', '0.1'});
xlabel('time (s)')

f3 = figure('name', 'DCT signal');
figure(f3)
set(gcf,'position',[x0,y0,width,height]);
plot(dct(x),'Linewidth', 1)
xlim([0,500])
