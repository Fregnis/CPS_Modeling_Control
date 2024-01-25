%% TASK 4
clear all
close all

load ("task4data.mat")
load ("task4_sensors_positions.mat")

j = 4; %target
p = 100;
q = 25; % sensors
h = 2; % sensors under attack

k = 2;

I = eye(size(D,1));

G = [D I];

G = normalize(G);

S_x_hat = randperm(p,j);

x_reale=zeros(p,1);
for i=1:length(S_x_hat)
    x_reale(S_x_hat(i))=1;
end

S_a_hat = randperm(q,h);

a_reale=zeros(q,1);

for i=1:length(S_a_hat)

    a_reale(S_a_hat(i)) = 100;
end

noise = 0.2 *randn * ones(q,1);

T = 100;

e=1e-8;

tau = norm(G)^-2 - e;

Lambda = zeros(p+q,1);

Lambda(1:p) = 10;

Lambda(p+1:p+q) = 20;
  

  

z_prev=zeros(p+q,1);
x_hat=zeros(p,1);
a_hat=zeros(q,1);
figure(1)

axis([0 1000 0 1000])

for t=1:T

    
    y= D * x_reale + noise ;%+ a_reale;
    for i=1:length(S_a_hat)

    y(S_a_hat(i)) = 1/2*y(S_a_hat(i));
    end
    z_int = shrinkage_thresholding(z_prev + tau * G' * (y - G * z_prev),tau * Lambda,p+q);
    x_hat = A * z_int(1:p);
    a_hat = z_int(p+1:p+q);
    z_prev=[x_hat' a_hat']';
    x_reale = A * x_reale;
       

figure(1)
max = maxk(x_hat,4);

index=zeros(1,4);


plot(sensors_pos(1,:),sensors_pos(2,:),'.g')
hold on
for i=1:length(S_a_hat)
    plot(sensors_pos(1,S_a_hat(i)),sensors_pos(2,S_a_hat(i)),'xr')
end
hold on

max_a = maxk(a_hat,2);
M = min(max_a);
idx_a = find(a_hat >= M);

plot(sensors_pos(1,idx_a(1)),sensors_pos(2,idx_a(1)),'ob')

plot(sensors_pos(1,idx_a(2)),sensors_pos(2,idx_a(2)),'ob')
hold off

hold on
for i=1:length(S_x_hat)

  scatter(mod(S_x_hat(i)-1,10)*100+50,fix((S_x_hat(i)-1)/10)*100+50,'sy','filled')

end

hold on
max_x = maxk(x_hat,4);
M = min(max_x);
idx_x = find(x_hat >= M);

for i=1:length(max_x)

    scatter(mod(idx_x(i)-1,10)*100+50,fix((idx_x(i)-1)/10)*100+50,'sb')
end
hold off
grid on

end