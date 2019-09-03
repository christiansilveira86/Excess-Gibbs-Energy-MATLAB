%% ABOUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program uses the NRTL model [1] to compute the activity
% coefficient of mixtures.
%
% This program is written by:
% Christian Silveira, D. Sc. (christiansilveira86@gmail.com)
% Chemical Engineering Department - Universidade Federal de Santa Maria
%
% For suggestions, error reports or comments, please send us an e-mail,
% we will be glad to hear from you.
%
% [1] Renon H., Prausnitz J. M., "Local Compositions in Thermodynamic 
% Excess Functions for Liquid Mixtures", AIChE J., 14(1), S.135–144, 1968
%%
%% INSTRUCTIONS
% The NRTL routine is nested in this file, the user may only provide the
% following information:
% NC - Number of components of the mixture
% x  - Molar fractions array
% T  - Temperature
% g  - Binary interaction parameters matrix
% alpha - Non-randomness parameter
%
% NC and T are scalars (one value only).
% x is an array with NC elements.
% g and alpha are (NC,NC) matrices, where the main diagonal elements are
% zeros and the element (i,j) is the interaction parameter between
% components i and j.
%   g(i,j)   ~=   g(j,i)
% alpha(i,j) == alpha(j,i)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
function main
% In this example, the activity coefficients of a ternary mixture
% containing Acetone (1), Methanol (2), and Water (3) are calculated at
% 324.44 K with molar fractions of [0.1 0.1 0.8].
%% INPUTS
NC = 3;             % Number of components
T = 342.44;         % Temperature [K]
R    = 1.9872;      % Universal gas constant [cal/(mol.K)]
x = [0.1 0.1 0.8];  % Molar fractions
% Binary interaction parameters matrix
g = [0.00000 1332.55 189.11;
     -519.89 0.00000 -94.23;
     1581.27 348.010 0.0000];
% Non-randomness parameters matrix
alpha = [0.000 0.284 0.301;
         0.284 0.000 0.309;
         0.301 0.309 0.000];

[gamma] = NRTL(x, T, NC, g, alpha, R)
disp('Activity Coefficients:')
for i = 1:NC
    disp(['Component ', num2str(i), ': ', num2str(gamma(i))]);
end

function gamma = NRTL(x, T, NC, g, alpha, R)
tau = g./(R*T);

for i = 1:NC
    for j = 1:NC
        G(i,j) = exp(-alpha(i,j)*tau(i,j));
    end
end

for i = 1:NC
    num1(i) = 0;
    den1(i) = 0;
    num2(i) = 0;
    
    for j = 1:NC
        num1(i) = tau(j,i)*G(j,i)*x(j) + num1(i);
        den1(i) = G(j,i)*x(j) + den1(i);
        den2(j) = 0;
        num3(j) = 0;
        den3(j) = 0;
        for k = 1:NC
            den2(j) = x(k)*G(k,j) + den2(j);
            num3(j) = x(k)*tau(k,j)*G(k,j) + num3(j);
            den3(j) = x(k)*G(k,j) + den3(j);
        end
        num2(i) = (x(j)*G(i,j))/(den2(j))*(tau(i,j) - num3(j)/den3(j)) + num2(i);
    end
    lngamma(i) = num1(i)/den1(i) + num2(i);
end

gamma = exp(lngamma);