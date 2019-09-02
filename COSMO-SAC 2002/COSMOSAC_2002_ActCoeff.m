%% ABOUT
% This program uses the sigma profiles of two pure components to calculate 
%   the liquid-phase activity coefficients in a solution.  This is the first 
%   step in predicting VLE for mixtures.  
%
%	This program uses the COSMO-SAC model as published (Lin, S.T., 
%	S.I. Sandler, Ind. Eng. Chem. Res. 41, (2002), 899-913).
%
%	THIS PROGRAM WAS ORIGINALLY WRITTEN BY:
%	  RICHARD OLDLAND (roldland@vt.edu)    MIKE ZWOLAK (zwolak@caltech.edu)
%	  DEPARTMENT OF CHEMICAL ENGINEERING   PHYSICS DEPARTMENT
%	  VIRGINIA TECH                        CALIFORNIA INSTITUTE OF TECHNOLOGY
%	  BLACKSBURG, VA 24060                 PASADENA, CA 91125
%
%%
%%   MATLAB IMPLEMENTATION AND MODIFICATIONS BY
%   Christian Luiz da Silveira (christiansilveira86@gmail.com)
%   Department of Chemical Engineering
%   Universidade Federal de Santa Maria
%   Santa Maria - RS, Brazil
%
%   Please, give us your feedback, we want to hear your opinion and
%   suggestions.
%%
%% INSTRUCTIONS
% The COSMO-SAC routine is nested in this file, the user may only provide
% the following information:
% x  - Molar fractions array of size (1,NC)
% T  - Temperature (scalar)
% NC - Number of components (scalar)
% comp - String array of size (1,NC)
%%
% Additional information:
% The sigma profiles used in this version are available at https://www.design.che.vt.edu/VT-Databases.html
% Some of the most common sigma-profiles will be added alongside this
% routine, for more components the user should visit the Database made
% available by the Liu research group (Virginia Tech).
%%
function main
% In this example, the activity coefficients of a mixture containing
% ethanol and water are computed at 303.15K and molar fractions of
% 0.25 for ethanol and 0.75 for water.
comp=["ethanol","water"];
NC = 2;
x=[0.25 0.75];
T=303.15;
gamma = COSMOSAC(x, T, NC, comp);

disp('Activity Coefficients:')
for i = 1:NC
    disp(['Component ', num2str(i), ': ', num2str(gamma(i))]);
end

function [gamma] = COSMOSAC(x, T, NC, comp)

fclose all;                 % Closes any previously opened file
Vcosmo = VCOSMODB(comp);    % Reads the COSMO volume of the components
%% Parameters
e0 = 2.395e-4;      % Permittivity of free space [(e^2.mol)/(kcal.Ang)]
aeff = 7.5;         % Effective area [Angs^2]
Rgas = 0.00198587;  % Ideal gas constant [kcal/(K.mol)]
Vnorm = 66.69;      % Standard volume [Angs^3]
Anorm = 79.53;      % Standard surface area [Angs^2]

compseg = 51;       % Number of semgments
coord = 10;         % Coordination number
SIGMAHB = 0.0084;   % Hydrogen-bonding cutoff
cHB = 85580.0;      % Hydrogen-bonding constant
fpol = 0.6916;      % Polarization factor
alpha = (0.3*(aeff^1.5))/e0;
alphaprime = fpol*alpha;    % Constant of the misfit energy

%% Reads the sigma-profile of the components
for i = 1:NC
    ff(i)   = join([comp(i) '.txt'],'');
    file(i) = fopen(ff(i), 'r');
    FF(2*i-1:2*i,:)   = fscanf(file(i), '%f %f', [2 Inf]);
end
FF = FF';

for i = 1:NC
    SIGMA(:,i) = FF(:,2*i);
    Acosmo(i) = sum(SIGMA(:,i));
end
SIGMA_mn = FF(:,1); % e/Angs^2

%% Calculate the mixture sigma-profile
for j = 1:compseg
    numer(j) = sum(x.*SIGMA(j,:));
    denom(j) = sum(x.*Acosmo);
end
profile = numer./denom;

%% Computes the exchange energy matrix deltaW
for j = 1:compseg
    for k = 1:compseg
        if (SIGMA_mn(j) >= SIGMA_mn(k))
            SIGMAACC = SIGMA_mn(j);
            SIGMADON = SIGMA_mn(k);
        end
        if (SIGMA_mn(j) < SIGMA_mn(k))
            SIGMADON = SIGMA_mn(j);
            SIGMAACC = SIGMA_mn(k);
        end
        % Misfit energy contribution
        EMF = (alphaprime/2)*(SIGMA_mn(j)+SIGMA_mn(k))^2;
        % Hydrogen bonding contribution
        EHB = cHB*max(0, (SIGMAACC - SIGMAHB))*min(0, (SIGMADON + SIGMAHB));
        deltaW(j,k) = EMF + EHB;
    end
end

%% Iteration for segment activity coeff (mixture)
SEGGAMMA = ones(compseg,1);
CONVERG = ones(compseg,1);
while max(CONVERG) > 1e-6
    SEGGAMMAOLD = SEGGAMMA;
    for j = 1:compseg
        SUM1 = 0;
        for k = 1:compseg
            SUM1 = SUM1 + profile(k)*SEGGAMMAOLD(k)*exp(-deltaW(j,k)/(Rgas*T));
        end
        SEGGAMMA(j) = exp(-log(SUM1));
        SEGGAMMA(j) = (SEGGAMMA(j) + SEGGAMMAOLD(j))/2;
    end
    for j = 1:compseg
        CONVERG(j) = abs((SEGGAMMA(j) - SEGGAMMAOLD(j))/SEGGAMMAOLD(j));
    end
end

%% Iteration for segment activity coefficient (pure species)
SEGGAMMAPR = ones(compseg,NC);
SEGGAMMAOLDPR = ones(compseg,NC);
for L = 1:NC
    SEGGAMMAPR(:,L) = 1;
    CONPR = ones(compseg,1);
    while max(CONPR) > 1e-6
        SEGGAMMAOLDPR(:,L) = SEGGAMMAPR(:,L);
        for j = 1:compseg
            SUM2 = 0;
            for k = 1:compseg
                SUM2 = SUM2 + (SIGMA(k,L)/Acosmo(L))*SEGGAMMAOLDPR(k,L)*exp(-deltaW(j,k)/(Rgas*T));
            end
            SEGGAMMAPR(j,L) = exp(-log(SUM2));
            SEGGAMMAPR(j,L) = (SEGGAMMAPR(j,L) + SEGGAMMAOLDPR(j,L))/2; 
        end
        for j = 1:compseg
            CONPR(j,L) = abs((SEGGAMMAPR(j,L) - SEGGAMMAOLDPR(j,L))/SEGGAMMAOLDPR(j,L));
        end
    end
end
%% Staverman-Guggenheim Equation
for j = 1:NC
    Rnorm(j) = Vcosmo(j)/Vnorm;
    Qnorm(j) = Acosmo(j)/Anorm;
end
BOTTHETA = sum(x.*Qnorm);
BOTPHI   = sum(x.*Rnorm);
for j = 1:NC
    THETA(j) = x(j)*Qnorm(j)/BOTTHETA;
    PHI(j)   = x(j)*Rnorm(j)/BOTPHI;
    L(j)     = (coord/2)*(Rnorm(j)-Qnorm(j)) - (Rnorm(j) - 1);
end
for j = 1:NC
    LNGAMMASG(j) = log(PHI(j)/x(j))+(coord/2)*Qnorm(j)*log(THETA(j)/PHI(j))+L(j)-(PHI(j)/x(j))*(sum(x.*L));
end

%% Calculation of the activity coefficients
SUMGAMMA = zeros(NC,1);
for j = 1:compseg
    for i = 1:NC
        SUMGAMMA(i) = SUMGAMMA(i) + ((SIGMA(j,i)/aeff)*(log(SEGGAMMA(j)/SEGGAMMAPR(j,i))));
    end
end

for j = 1:NC
    gamma(j) = exp(SUMGAMMA(j) + LNGAMMASG(j));
end
fclose all;

%% REFERENCES
%  LITERATURE CITED:
%  Klamt, A. Conductor-like Screening Model for Real Solvents: A New Approach to the
%       Quantitative Calculation of Solvation Phenomena. J. Phys. Chem 1995, 99, 2224.
%  Klamt, A.; Jonas, V.; Burger, T.; Lohrenz, J. Refinement and Parameterization of 
%	COSMO-RS. J. Phys. Chem A 1998, 102, 5074.
%  Klamt, A.; Eckert, F.; COSMO-RS: A Novel and Efficient Method for the a Priori 
%	Prediction of Thermophysical Data of Liquids.  Fluid Phase Equilibria 2000, 
%	172, 43.
%  Lin, S.T.; Sandler, S. A Priori Phase Equilibrium Prediction from a Segment 
%       Contribution Solvation Model. Ind. Eng. Chem. Res, 2002, 41, 899 
%  Lin, S.T.;  Quantum Mechanical Approaches to the Prediction of Phase Equilibria: 
%	Solvation Thermodynamics and Group Contribution Methods, PhD. Dissertation, 
%	University of Delaware, Newark, DE, 2000