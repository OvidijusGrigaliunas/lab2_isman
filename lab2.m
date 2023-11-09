close all
clear
clc
% 1. Duomenų paruošimas
X = 0:0.05:1;
% for H = 1:19
% for K = 1:19
Y = -sin(4*pi*X) + 6 * cos(pi*X);
step = 0.1;
n = 8; % neuronu skaicius
w = randn(n,2);
b = randn(n,2);
v = zeros(n,2);
y = zeros(n,1);
delta = zeros(6,2);
y_apskaiciuota = zeros(1,length(X));
for k = 1:1000
    for j = 1:length(X)
        % 1 sluosknis
        for i = 1:n
            v(i,1) = X(j) * w(i,1) + b(i,1);
            %y(i) = 1/(1+exp(-v(i,1)));
            y(i) = (exp(v(i,1))-exp(-v(i,1)))/(exp(v(i,1))+exp(-v(i,1)));
        end
        % 2 sluosknis
        v(1,2) = sum(w(:,2)'.*y') + b(1,2);
        y_apskaiciuota(j) = v(1,2);
        e = Y(j) - y_apskaiciuota(j);
        % Koeficientu atnaujinimas
        delta_out = e;
        b(1,2) = b(1,2) + step*delta_out;
        for i = 1:n
            %delta = y(i)*(1 - y(i))*(delta_out*w(i,2));
            delta = (1-y(i)^2)*(delta_out*w(i,2));
            w(i,2) = w(i,2) + step*delta_out*y(i);
            w(i,1) = w(i,1) + step*delta*X(j);
            b(i,1) = b(i,1) + step*delta;
        end
    end
end

X2 = 0:0.01:1;

for j = 1:length(X2)
    v(1,2) = 0;
    for i = 1:n
        v(i,1) = X2(j) * w(i,1) + b(i,1);
        y(i) = (exp(v(i,1))-exp(-v(i,1)))/(exp(v(i,1))+exp(-v(i,1)));
    end
    v(1,2) = sum(w(:,2)'.*y') + b(1,2);
    y_apskaiciuota(j) = v(1,2);
end
figure
plot(X,Y, '*',X2, y_apskaiciuota, 'rx')
legend('Tikroji',"Apskaičiuota")
