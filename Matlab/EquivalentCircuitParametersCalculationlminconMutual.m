%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written by Juan A. Vasquez-Peralvo
%Luxembourg-Luxembourg
%1/1/2024
%Modified 1/1/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Description%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code allows to calculate the parameters of an equivalent circuit.
%This code is intended to be used in a bandpass FSS in the case of a band
%stop the code will be the same but the forumla inverted and the C and L
%are going to be interchanged.

%%%%%%%%%%%%%%%%%%%%%%%%Initialize Matlab%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Add Required Paths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------------------------------------------------
%                            Main Script Logic
% ------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( 1, '------------------------------------------------------------\n' );
fprintf( 1, '------ Equivalent circuit element calculation  --------\n' );
fprintf( 1, '------------------------------------------------------------\n' );
fprintf( 1, '------ Author(s):                                   --------\n' );
fprintf( 1, '------          Juan Andrés Vásquez Peralvo         --------\n' );
fprintf( 1, '------------------------------------------------------------\n' );
fprintf( 1, '------ Date:      January 2024                         --------\n' );
fprintf( 1, '------------------------------------------------------------\n' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start Timer
tini = clock;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% Initiate constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z0 = 377;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Define folders %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('./Functions/')
%%%%%%%%%%%%%%%%%%%%%%%%%Import Data from CST%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CSTProject = '../CST/UnitCellMetal.cst';
CST = TCSTInterface();
CST.OpenProject(CSTProject);
TreeItem = '1D Results\S-Parameters\SZmax(1),Zmax(1)'; 
[frequency,S11,Zref,RunIDs,Info] = CST.Get1DResultFromTreeItem(TreeItem);
TreeItem = '1D Results\S-Parameters\SZmin(1),Zmax(1)'; 
[frequency,S21,Zref,RunIDs,Info] = CST.Get1DResultFromTreeItem(TreeItem);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% Define imported varaibles %%%%%%%%%%%%%%%%%%%%%%%%%%%
w=frequency*2*pi*1e9; % define w
S21=abs(S21(1:end,1)); % define S21
S11=abs(S11(1:end,1)); % define S11
zFssS = S21*z0./(2*(1-S21)); % Calculate the impedance
[A I]=max(zFssS); %extract minimum of Z
[A II]=min(zFssS(I:end)); %extract maximum of Z avoiding the first one
II = II + I - 1; %locate the exact position of the maximum
S21 = S21(1:II-10); % cut the values of Z to avoid wrong calculations
w = w(1:II-10); % limit values of w to avoid wrong calculations
w1=2*pi*frequency(I(1))*1e9; %% Z Pole
w2=2*pi*frequency(II(1))*1e9; %% Z cero
zFssS = zFssS(1:II-10); % limit the values of Z to avoid errors
yFssS = zFssS.^(-1); % Calculate the admitance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% plot the Z of the imported result %%%%%%%%%%%%%%%%%%%%%%%%
plot(w, zFssS)
title('Simulated Impedance')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Curve Fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Guesses for C
initialGuess = [1e-16]; % a initial value for that 

% Curve Fitting
% Optimization
% The following values will allow the optimizer to find a suitable value
% eventhough that the values are extremelly small.
options = optimset('PlotFcns',@optimplotfval, 'TolX', 1e-80, 'TolFun', 1e-80);
C = fminsearch(@(x) costFunctionimpedanceC(x, w1, w2, w, ...
    yFssS), initialGuess, options); % lunch the optimizer

% Display the results
    % Display the optimized value of C
    fprintf('Optimized C: %e Farads\n', C);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%% plot the results of the first optimization %%%%%%%%%%
figure
plot(w, yFssS); 
yFssC = abs((((((1./(C*((w1)^2)))+((1./(C*((w1)^2)))/(((w2)^2/(w1)^2-1))))./ ...
    (C)-w.^2*(1./(C*((w1)^2)))*((1./(C*((w1)^2)))/(((w2)^2/(w1)^2-1))))./ ...
    (1j*w*(1./(C*((w1)^2)))+1./(1j*w*C))).^(-1)));
hold on
plot(w, yFssC)
legend('Simulated', 'Calculated')
title('First Optimization Results')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Second Optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Curve Fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Guesses for R, L, and C
initialGuess = [1e-9 1e-9]; % Example values: R = 1 uOhm, L = 1 nH, C = 1 pF
lb = [ 0 0 ];
ub = [ 1e-6 1e-6];
% Curve Fittin
    options = optimset('Display', 'iter','TolX', 1e-40, 'TolFun', 1e-40 );
    x = fminsearch(@(x) costFunctionimpedanceLL1(x, C, w, yFssS), initialGuess,  options);

% options = optimoptions('lsqcurvefit', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt', ...
%     'FiniteDifferenceType','central',...
%     'Display','iter');
%objective = @(x)impedanceC(x, w, w1, w2);
%x = lsqcurvefit(@(x,w)impedanceLL1(x, C/1e-15, w./1e12), initialGuess, w(1:end-30), yFssS(1:end-30), lb, [], options);
L = x(1)
L1 = x(2)

figure
plot(w, yFssS)
yFssC = abs((((L+L1)./(C)-w.^2*L*L1)./(1j*w*L+1./(1j*w*C))).^(-1));
hold on
plot(w, yFssC)
legend('Simulated', 'Calculated')

figure
s21C = 20*log10(abs((2*((L+L1)./(C)-w.^2*L*L1)./(j*w*L+1./(j*w*C)))./(z0+2*((L+L1)./(C)-w.^2*L*L1)./(j*w*L+1./(j*w*C)))));
plot(w, 20*log10(S21));
hold on
plot(w, s21C)
legend('Simulated', 'Calculated')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Stop Timer
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
duration=etime(clock, tini);
dminutes=floor(duration/60);
dseconds=duration-dminutes*60;
dhours=floor(dminutes/60);
dminutes=dminutes-dhours*60;
fprintf( 1, '------------------------------------------------------------\n' );
disp(['    Total Time = '  num2str(dhours) 'h ' num2str(dminutes) 'min ' num2str(dseconds) 'sec']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf( 1, '------------------------------------------------------------\n' );
fprintf( 1, '--------------------    Thank you        -------------------\n' );
fprintf( 1, '--------------------       END           -------------------\n' );
fprintf( 1, '------------------------------------------------------------\n' );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%