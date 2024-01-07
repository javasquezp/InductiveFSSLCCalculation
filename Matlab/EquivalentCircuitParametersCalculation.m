%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written by Juan A. Vasquez-Peralvo
%Luxembourg-Luxembourg
%1/1/2024
%Modified 1/1/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Description%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code allows to calculate the parameters of an equivalent circuit

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
tini=clock;
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
yFssS = zFssS.^(-1); % Calculate the impedance
[A I]=max(zFssS); %extract minimum of S11
[A II]=min(zFssS); %extract maximum of S11
%S21=S21(1:II); % cut the values of S21
%w=w(1:II); % limit values of w
w1=2*pi*frequency(I(1))*1e9; %% S11 Pole
w2=2*pi*frequency(II(1))*1e9; %% S11 cero
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% plot the Z of the imported result %%%%%%%%%%%%%%%%%%%%%%%%
plot(frequency, zFssS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Curve Fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Guesses for R, L, and C
initialGuess = [1e-15]; % Example values: R = 1 uOhm, L = 1 nH, C = 1 pF

% Curve Fitting
options = optimoptions('lsqcurvefit', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt', 'TolFun', 1e-40, 'TolX', 1e-40, 'DiffMaxChange', 1e-30, 'DiffMinChange', 1e-40 );
%objective = @(x)impedanceC(x, w, w1, w2);
x = lsqcurvefit(@(C,w)impedanceC(C, w1, w2, w), initialGuess, w(1:end-30), yFssS(1:end-30), [], [], options);

% Display the results

C= x(1);

fprintf('Fitted C: %e Farads\n', C);
figure
plot(frequency(1:end-30), yFssS(1:end-30))
yFssC = abs((((((1./(C*((w1)^2)))+((1./(C*((w1)^2)))/(((w2)^2/(w1)^2-1))))./(C)-w.^2*(1./(C*((w1)^2)))*((1./(C*((w1)^2)))/(((w2)^2/(w1)^2-1))))./(1j*w*(1./(C*((w1)^2)))+1./(1j*w*C))).^(-1)));
hold on
plot(frequency(1:end-30), yFssC(1:end-30))
legend('Simulated', 'Calculated')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Second Optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Curve Fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Guesses for R, L, and C
initialGuess = [1e-10/1e-9 1e-10/1e-9]; % Example values: R = 1 uOhm, L = 1 nH, C = 1 pF
lb = [ 0 0 ];
% Curve Fittin
options = optimoptions('lsqcurvefit', 'Display', 'iter', 'Algorithm', 'levenberg-marquardt', ...
    'FiniteDifferenceType','central',...
    'Display','iter');
%objective = @(x)impedanceC(x, w, w1, w2);
x = lsqcurvefit(@(x,w)impedanceLL1(x, C/1e-15, w./1e12), initialGuess, w(1:end-30), yFssS(1:end-30), lb, [], options);
L = x(1)
L1 = x(2)

figure
plot(frequency(1:end-30), yFssS(1:end-30))
yFssC = abs((((L+L1)./(C)-w.^2*L*L1)./(1j*w*L+1./(1j*w*C))).^(-1));
hold on
plot(frequency(1:end-30), yFssC(1:end-30))
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