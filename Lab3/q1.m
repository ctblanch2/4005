clear all;
close all;
clc;

lambda = 5
U = xlsread('100K_uniform.xlsx'); %import unfiform random numbers

X = ( -1 / lambda ) * log(1-U);


figure();
histogram(U,77,'Normalization','probability');
title('Uniform Random Variable Input');

figure();
histogram(X,77,'Normalization','probability');
title('Exponential Random Variable');

V = -0.1 + (0.2).*U;
V2 = -0.01 + (0.02).*U;


i = 2;
CorrX = zeros(1, 100000);
Corr2 = zeros(1, 100000);

CorrX(1) = U(1);
Corr2(1) = U(1);
while i <= 99999
    CorrX(i) = mod((CorrX(i-1) + V(i)), 1);
    Corr2(i) = mod((Corr2(i-1) + V2(i)), 1);
    i = i + 1;
end

Stitched = zeros(1,100000);
Stitched2 = zeros(1,100000);

i = 1;
while i <= 100000
    if CorrX(i) <= 0.7
        Stitched(i) = CorrX(i)/0.7;
    else
        Stitched(i) = (1 - CorrX(i)) / (1-0.7);
    end
    
    if Corr2(i) <= 0.7
        Stitched2(i) = Corr2(i)/0.7;
    else
        Stitched2(i) = (1 - Corr2(i)) / (1-0.7);
    end
    
    i = i + 1;
end

corrExponential = ( -1 / lambda ) * log(1-Stitched);
corrExponential2 = ( -1 / lambda ) * log(1-Stitched2);

figure();
histogram(corrExponential,77,'Normalization','probability');
title('TES Exponential');


figure();
histogram(corrExponential2,77,'Normalization','probability');
title('TES Exponential 2');
    