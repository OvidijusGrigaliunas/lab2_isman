clear all
clc
close all

% 1. Duomenų paruošimas
X1 = 0:0.05:1;
X2 = 0:0.05:1;
for i=1:length(X1)
    for j=1:length(X2)
        Y(i,j) = -sin(4*pi*X1(i)) + 2 * sin(pi*X2(j));
    end
end
surf(X1,X2,Y)

step = 0.1;
n = 32; % neuronu skaicius
w_1 = randn(n,2); w_2 = randn(n,1); b = randn(n,2);
%load('w_1.mat'); load('w_2.mat'); load('b.mat');

v_1 = zeros(n,2);
v_2 = 0;
y = zeros(n,1);
y_apskaiciuota = zeros(1,1);
for k = 1:30000
    for j = 1:length(X1)
         for l = 1:length(X2)
            % 1 sluosknis
            for i = 1:n
                v_1 = X1(j) * w_1(i,1) + X2(l) * w_1(i,2) + b(i,1);
                %y(i) = 1/(1+exp(-v(i,1)));
                y(i) = tanh(v_1);
            end
            % 2 sluosknis
            v_2 = sum(w_2(:).*y) + b(1,2);
            e = Y(j, l) - v_2;
            % Koeficientu atnaujinimas
            delta_out = e;
            b(1,2) = b(1,2) + step*delta_out;
            for i = 1:n
                %delta = y(i)*(1 - y(i))*(delta_out*w(i,2));
                delta = (1-y(i)^2)*(delta_out*w_2(i));
                step_delta = step*delta;
                w_2(i) = w_2(i) + step_delta*y(i);
                w_1(i,1) = w_1(i,1) + step_delta*X1(j);
                w_1(i,2) = w_1(i,2) + step_delta*X2(l);
                b(i,1) = b(i,1) + step_delta;
            end
        end
    end
    
end
X3 = 0:0.01:1;
X4 = 0:0.01:1;
for j = 1:length(X3)
     for l = 1:length(X4)
        for i = 1:n
            v_1 = X3(j) * w_1(i,1) + X4(l) * w_1(i,2) + b(i,1);
            y(i) = tanh(v_1);
        end
        v_2 = sum(w_2(:).*y) + b(1,2);
        y_apskaiciuota(j, l) = v_2;
     end
end
figure;
surf(X3,X4,y_apskaiciuota)

% Kadangi reikia daug treniruoti:
save("w_1.mat", "w_1")
save("w_2.mat", "w_2")
save("b.mat", "b")
