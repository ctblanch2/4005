clear all;
close all;
clc;

lambda = 5
U = xlsread('100K_uniform.xlsx'); %import unfiform random numbers

X = ( -1 / lambda ) * log(1-U);


%figure();
%histogram(U,77,'Normalization','probability');
%title('Uniform Random Variable Input');

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

figure();
autocorr(X, 20);
title('Uncorrelated');

figure();
autocorr(corrExponential, 20);
title('a=-b=0.1');

figure();
autocorr(corrExponential2, 20);
title('a=-b=0.01');

edges = 0:0.005:max(X);
edges2 = 0:0.005:max(corrExponential);
edges3 = 0:0.005:max(corrExponential2);

lambda = 1 / mean(X);
lambda2 = 1 / mean(corrExponential);
lambda3 = 1 / mean(corrExponential2);
expectedFrequencies = zeros(1, length(edges) - 1);
expectedFrequencies2 = zeros(1, length(edges2) - 1);
expectedFrequencies3 = zeros(1, length(edges3) - 1);
i = 2;
while i < length(edges)
    expectedFrequencies(i-1) = length(X) * (expcdf(edges(i), 1 / lambda) - expcdf(edges(i-1), 1/lambda));
    i = i + 1;
end

i = 2;
while i < length(edges2)
    expectedFrequencies2(i-1) = length(corrExponential) * (expcdf(edges2(i), 1 / lambda2) - expcdf(edges2(i-1), 1/lambda2));
    i = i + 1;
end

i = 2;
while i < length(edges3)
    expectedFrequencies3(i-1) = length(corrExponential2) * (expcdf(edges3(i), 1 / lambda3) - expcdf(edges3(i-1), 1/lambda3));
    i = i + 1;
end


expectedFrequencies(1:10)
expectedFrequencies2(1:10)
expectedFrequencies3(1:10)


[newEdges, newValues] = combineBuckets(edges, expectedFrequencies);
[newEdges2, newValues2] = combineBuckets(edges2, expectedFrequencies2);
[newEdges3, newValues3] = combineBuckets(edges3, expectedFrequencies3);

"New Edges"
newEdges(1:10)
[observedFrequencies, E] = histcounts(X, newEdges);
[o2, E2] = histcounts(corrExponential, newEdges2);
[o3, E3] = histcounts(corrExponential2, newEdges3);

observedFrequencies(1:10)
o2(1:10)
o3(1:10)

chisq = zeros(1, length(newValues));
chisq2 = zeros(1, length(newValues2));
chisq3 = zeros(1, length(newValues3));

i = 1;
while i <= length(observedFrequencies)
    chisq(i) = ((observedFrequencies(i) - newValues(i))^2)/newValues(i);
    i = i + 1;
end

i = 1;
while i <= length(o2)
    chisq2(i) = ((o2(i) - newValues(i))^2)/newValues(i);
    i = i + 1;
end

i = 1;
while i <= length(o3)
    chisq3(i) = ((o3(i) - newValues3(i))^2)/newValues3(i);

    i = i + 1;
end
"Chi square 1: " + sum(chisq)
"Chi square 2: " + sum(chisq2)
"Chi square 3: " + sum(chisq3)
chisq(1:10);
chisq2(1:10);
chisq3(1:10);
chi2inv(0.95, length(observedFrequencies))
chi2inv(0.95, length(o2))
chi2inv(0.95, length(o3))



function [edges, values] = combineBuckets(currentEdges, currentValues)
    edges = zeros(1, length(currentEdges));
    values = zeros(1, length(currentValues));
    i = 1;
    j = 1;
    edges(1) = -0.005;
    while i < length(currentEdges)
        if values(j) >= 5
           j = j + 1;
           edges(j) = edges(j-1);
        end
        edges(j) = edges(j) + 0.005;
        values(j) = values(j) + currentValues(i);
        i = i + 1;
    end
    if values(j) < 5
        j = j-1;
        values(j) = values(j) + values(j+1);
        edges(j) = edges(j) + 0.005;
    end
    edges = edges(1:j);
    values = values(1:j);
end

