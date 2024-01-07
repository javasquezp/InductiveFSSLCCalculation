%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Written by Juan A. Vasquez-Peralvo
%Luxembourg-Luxembourg
%1/1/2024
%Modified 1/1/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Description%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code allows to calculate the parameters of an equivalent circuit for
%a LC FSS """"""INDUCTIVE"""""". This codes needs as the input the cst file 
% where the FSS is located. As an output it gives the calculated LC that values.

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
[A I]=max(zFssS); %extract minimum of Z
[A II]=min(zFssS(I:end)); %extract maximum of Z avoiding the first one
II = II + I - 1; %locate the exact position of the maximum
S21 = S21(1:II-10); % cut the values of Z to avoid wrong calculations
w = w(1:II-10); % limit values of w to avoid wrong calculations
w1=2*pi*frequency(I(1))*1e9; %% Z Pole
w2=2*pi*frequency(II(1))*1e9; %% Z cero
zFssS = zFssS(1:II-10); % limit the values of Z to avoid errors
yFssS = zFssS.^(-1); % Calculate the admitance
frequency = frequency(1:II-10);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% plot the Z of the imported result %%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(w, zFssS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%% Second Optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Curve Fit %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Guesses for R, L, and C
initialGuess = [20 5]; % Example values: R = 1 uOhm, L = 1 nH, C = 1 pF
lb = [ 0 0 ];
ub = [ 1e-7 1e-15];
% Curve Fittin
    options = optimset('Display', 'iter');
    x = fmincon(@(x) costFunctionimpedanceLC(x, w, 20*log10(S21)), initialGuess, [], [], [], [], lb, [], [], options);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Show Optimized Results %%%%%%%%%%%%%%%%%%%%%%%%%%
disp(['Calculated L value:' num2str(x(1)*1e-12) 'H']);
disp(['Calculated C value:' num2str(x(2)*1e-12) 'F'])
L = x(1)*1e-12;
C = x(2)*1e-15;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot(w, zFssS)
yFssC = abs(  (  (1j.*w*L./(1j*w*C)) ./  (1j.*w*L +1./(1j*w*C)) ) );
hold on
plot(w, yFssC, 'LineWidth', 2)
legend('Simulated', 'Calculated')
title('Impedance')
 xlabel('Frequency (GHz)')
 ylabel('|Z| (\ohm)')
%Calculate the S21
s21C = 20*log10(abs((2*(  (1j.*w*L./(1j*w*C)) ./  (1j.*w*L +1./(1j*w*C)) ))./ ...
    (z0+2*(  (1j.*w*L./(1j*w*C)) ./  (1j.*w*L +1./(1j*w*C)) ))));
figure
plot(frequency, s21C)
hold on
plot(frequency, 20*log10(S21), 'LineWidth',2)
legend('Simulated', 'Calculated')
xlabel('Frequency (GHz)')
ylabel('|S21| (dB)')
title('|S21|')

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