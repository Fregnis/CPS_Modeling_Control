close all
clear all
clc

%% TASK 5
clear all
close all
load ("stochasticmatrices.mat")
correct_est_a = zeros(1,20);

tmp=0;
norms = zeros(1,20);
for w=1:20
 n = 10;
 q = 20;
 h = 2; 
 C = randn(q,n);
 x_tilde = randn(n,1);

 S_a = randperm(q,h);

 %unaware attack

 val = unifrnd([-2 1], [-1 2]);

 a = zeros(q,1);

 for i = 1:length(val)
     a(S_a(i)) = val(i);
 end


noise =1e-2*randn * ones(q,1);

tau = 0.03;

lambda = 2e-4 / tau;

tau_Lambda = [zeros(n,1); tau*lambda*ones(q,1)];

x_prev = zeros(q,n);

T = 1e5;

delta = 1e-7;

y = C*x_tilde + noise + a;

I = eye(q);
G = [C I];
z_hat = zeros(q,n+q);
z_prev = zeros(q,n+q);

for k = 1:T

   for i= 1:q
    sum2 = 0;
      for j=1:q
          sum2 = sum2+ Q4(i,j) .* z_prev(j,:)';
      end

        z_hat(i,:) = shrinkage_thresholding(sum2 + tau * G(i,:)' * (y(i) - G(i,:)*z_prev(i,:)'),tau_Lambda,n+q); 
   end
   
   somma = 0;
  for l = 1:q
    
      somma = somma + vecnorm(z_hat(l) - z_prev(l));
     
  end

  if somma < delta
      break;
  end

   z_prev = z_hat;
end

x_hat = z_hat(:,1:n);
means = mean(x_hat);

norm_x = norm(x_tilde - means')^2;


a_hat = z_hat(:,n+1:n+q);

for i = 1:q
    for j=1:20

        if abs(a_hat(i,j)) < 0.2
            a_hat(i,j) = 0;
        end
    end
end



means_a = mean(a_hat);

idx_a = find(means_a);
S_a = sort(S_a);

if length(idx_a) == h
correct_est_a(w) = idx_a(1) == S_a(1) && idx_a(2) == S_a(2);
end
tmp = tmp+k;
norms(w) = norm_x;
end 

sum(correct_est_a,"all")/20
tmp/20