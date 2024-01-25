clear all
close all
clc

n = 10;
q = 20;
h = 2;
C = rand(q,n);

Tmax = 150;

counter = 0;
norms = zeros(20,1);
for t=1:20
x_reale = rand(n,1);

% unaware attack
 S_a =sort(randperm(q,h));

 
 a_reale = zeros(q,1);
 
  for i=1:length(S_a)
    
        a_reale(S_a(i)) = unifrnd(-1,1);
    
        if a_reale(S_a(i)) > 0
            a_reale(S_a(i)) = a_reale(S_a(i)) + 1;
        else
            a_reale(S_a(i)) = a_reale(S_a(i)) -1;
        end
    
  end

  tau_lambda = 2e-3;
  tau_Lambda = tau_lambda*[zeros(n,1);ones(q,1)];

  noise1 = zeros(q,1);

  sigma = 1e-2;
  noise2 = sigma*randn(q,1);

  y = C * x_reale + noise1;%+ a_reale;

  z_prev = zeros(n+q,1);
  z_hat = zeros(n+q,1);

  

  G = [C eye(q)];
 
  tau = norm(G)^-2 -1e-8;

  for i=1:length(S_a)

      y(S_a(i)) = y(S_a(i))/2;
  end

for i = 1:Tmax
    
    
        z_prev = z_hat;
        z_hat = shrinkage_thresholding(z_prev + tau * G' * (y - G * z_prev),tau_Lambda,n+q); 
    
        if norm(z_hat(1:n) - z_prev(1:n)) < 1e-12
            break;
        end
   
end

   a_hat = z_hat(n+1:n+q);
   x_hat = z_hat(1:n);

   idx_a = find(abs(a_hat)>0.2);

   if isequal(idx_a', S_a)
       counter = counter +1;
   end

    norms(t) = norm(x_reale-x_hat)^2;
end
    
% unaware -> due sensori sotto attacco è molto influente infatti le nomrme
% hanno valori elevati quindi l'algoritmo non converge mentre è resiliente
% al rumore in quanto le stime si avvicinano a quelle reali.

%con gli aware attack la stima non è buona ma è migliore rispetto al caso
%unaware in quanto influiscono con un valore più piccolo rispetto al caso
%precedente
