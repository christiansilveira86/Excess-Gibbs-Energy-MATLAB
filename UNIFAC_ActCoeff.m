%% ABOUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program uses the UNIFAC model [1] to compute the activity
% coefficient of mixtures.
%
% This program is written by:
% Christian Silveira, D. Sc. (christiansilveira86@gmail.com)
% Chemical Engineering Department - Universidade Federal de Santa Maria
%
% For suggestions, error reports or comments, please send us an e-mail,
% we will be glad to hear from you.
%
% [1] Fredenslund, A. , Jones, R. L. and Prausnitz, J. M. (1975)
% Group?contribution estimation of activity coefficients in nonideal liquid
% mixtures. AIChE J., 21: 1086-1099. doi:10.1002/aic.690210607
%%
%% INSTRUCTIONS
% The UNIFAC routine is nested in this file, the user may only provide the
% following information:
% NC - Number of components of the mixture
% NG - Number of functional groups on the mixture
% v  - The functional groups frequency in each component
% Rk - Functional group volumetric parameter
% Qk - Functional group area parameter
% x  - Molar fractions array
% T  - Temperature
% a  - Matrix of functional-group energy interaction parameters
% The values for Qk, Rk, and a can be found at:
% http://www.ddbst.com/published-parameters-unifac.html
%
% NC, NG, and T are scalars (one value only).
% v is a matrix of size (NC,NG).
% Rk and Qk are arrays of NG elements.
% Rk and Qk
% x is an array with NC elements.
% a is a (NG,NG) matrix 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function main
% This example computes the ativity coefficient of a mixture of
% water (1) + ethanol (2) at 333 K and molar fractions of
% 0.9 and 0.1 of water and ethanol, respectively.
%% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NC = 2;
% Functional groups:
% Water: 1 H2O
% Ethanol: 1 CH3, 1 CH2, 1 OH
% Sequence: CH3, CH2, OH, H2O
NG = 4;
v = [0 0 0 1;
     1 1 1 0];
Rk = [0.9011 0.6744 1.0000 0.9200];
Qk = [0.8480 0.5400 1.2000 1.4000];
% 1 1 5 7
a  = [0.00000 0.00000 986.50 1318.00
      0.00000 0.00000 986.50 1318.00
      156.400 156.400 0.0000 353.500
      300.000 300.000 -229.1 0.0000];
x = [0.9 0.1];
T = 333;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

gamma = UNIFAC(x, T, NC, NG, v, Rk, Qk, a);

disp('Activity Coefficients:')
for i = 1:NC
    disp(['Component ', num2str(i), ': ', num2str(gamma(i))]);
end

function gamma = UNIFAC(x, T, NC, NG, v, Rk, Qk, a)

% X is the mole fraction of group m in the mixture
for m = 1:NG
    X(m) = sum(v(:,m).*x')/sum(x*v);
end
% Xi is the mole fraction of group m in component i
for i = 1:NC
    for m = 1:NG
        Xi(i,m) = v(i,m)*x(i)/sum(v(i,:).*x(i));
    end
end
% THETA is the area fraction of group m in the mixture
for m = 1:NG
    THETA(m) = Qk(m)*X(m)/sum(Qk.*X);
end
% THETAi is the area fraction of group m in component i
for i = 1:NC
    for m = 1:NG
        THETAi(i,m) = Qk(m)*Xi(i,m)/sum(Qk.*Xi(i,:));
    end
end
% Group interaction parameters psi
psi = zeros(NG,NG);
for m = 1:NG
    for n = 1:NG
        psi(m,n) = exp(-a(m,n)/T);
    end
end
% Functional-group activity coefficient lnGAMMA in a mixture
for k = 1:NG
    sumA = 0;
    sumB = 0;
    for m = 1:NG
        sumden = 0;
        for n = 1:NG
            sumden = THETA(n)*psi(n,m) + sumden;
        end
        sumA = THETA(m)*psi(k,m)/sumden + sumA;
        sumB = THETA(m)*psi(m,k) + sumB;
    end
    lnGAMMA(k) = Qk(k)*(1 - log(sumB) - sumA);
end

% Functional-group activity coefficient lnGAMMAi in component i
clear sumA sumB sumden
for i = 1:NC
    for k = 1:NG
        sumA = 0;
        sumB = 0;
        for m = 1:NG
            sumden = 0;
            for n = 1:NG
                sumden = THETAi(i,n)*psi(n,m) + sumden;
            end
            sumA = THETAi(i,m)*psi(k,m)/sumden + sumA;
            sumB = THETAi(i,m)*psi(m,k) + sumB;
        end
        lnGAMMAi(i,k) = Qk(k)*(1 - log(sumB) - sumA);
    end
end

% Residual contribution
sumres = zeros(NC);
for i = 1:NC
    for k = 1:NG
        sumres(i) = v(i,k)*(lnGAMMA(k) - lnGAMMAi(i,k)) + sumres(i);
    end
    lngammaR(i) = sumres(i);
end

% r and q
for i = 1:NC
    r(i) = sum(v(i,:).*Rk);
    q(i) = sum(v(i,:).*Qk);
end

% Li
z = 10;
for i = 1:NC
    L(i) = (z/2)*(r(i)-q(i)) - (r(i)-1);
end
% Area and segment fractions, theta and phi
for i = 1:NC
    theta(i) = q(i)*x(i)/sum(q.*x);
    phi(i)   = r(i)*x(i)/sum(r.*x);
end

% Combinatorial contribution
for i = 1:NC
    lngammaC(i) = log(phi(i)/x(i)) + (z/2)*q(i)*log(theta(i)/phi(i)) + L(i) - phi(i)/x(i)*sum(L.*x);
end

for i = 1:NC
    lngamma(i) = lngammaC(i) + lngammaR(i);
end
gamma = exp(lngamma);