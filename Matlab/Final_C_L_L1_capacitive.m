clear all
close all
clc
addpath ./Functions
CSTProject = '../CST/Unit_cell_2.cst';
CST = TCSTInterface();
CST.OpenProject(CSTProject);
TreeItem = '1D Results\S-Parameters\SZmax(1),Zmax(1)'; 
[frequency,S11,Zref,RunIDs,Info] = CST.Get1DResultFromTreeItem(TreeItem);
TreeItem = '1D Results\S-Parameters\SZmin(1),Zmax(1)'; 
[frequency,S21,Zref,RunIDs,Info] = CST.Get1DResultFromTreeItem(TreeItem);
w=frequency*2*pi*1e9;
S21=abs(S21(1:end,1));
S11=abs(S11(1:end,1));
w=w(1:end);
 [A I]=min(S21);
 [A II]=max(S21);
 S21=abs(S21(1:II(1)));
S11=abs(S11(1:II(1)));
w=w(1:II(1));
 w1=2*pi*frequency(I(1))*1e9; %% S11 Pole
 w2=2*pi*frequency(II(1))*1e9; %% S11 cero
[fit_results gfo]=C_Calculation_cap(w,20*(log10(S21)))
C=fit_results.C;
%C=5.514475873546713e-13;
 L=(1./(C*((w1)^2)));
 C1=C/(((w2)^2/(w1)^2-1));
[fit_results gfo]=C1_L_Calculation(w,20*log10((S21)), C)


bwdown=find(20*log10(abs(S11(:,1)))<-10,1);
bwup=find(20*log10(abs(S11(bwdown:end,1)))>-10,1)+bwdown;
fr=frequency(I(1));
FBW=(frequency(bwup)-frequency(bwdown))/(frequency(I(1)));

%C=(6*FBW)./(377*2*pi*fr*1e9);
%L=1./((2*pi*fr*1e9).^2.*C);

% 
% 
% j*(((w.^2*L*C-1))./w.*((-w.^2*L*C*C1+C1+C)))
% %(j*(w.^2*(1/(w2^2*C))*C-1))./(w.*(-w.^2*(1/(w2^2*C))*C*((C)/((w1^2/w2^2)-1))+((C)/((w1^2/w2^2)-1))+C))
% (j*(((w.^2*(1/(w2^2*C))*C-1))./(w.*((-w.^2*(1/(w2^2*C))*C*((C)/((w1^2/w2^2)-1))+((C)/((w1^2/w2^2)-1))+C)))))
% ((j*w*L+1./(j*w*C)).*1./(j*w*C1))./(j*w*L+1./(j*w*C)+1./(j*w*C1))
% (((j*w*(1/(w2^2*C))+1./(j*w*C)).*1./(j*w*((C)/((w1^2/w2^2)-1))))./(j*w*(1/(w2^2*C))+1./(j*w*C)+1./(j*w*((C)/((w1^2/w2^2)-1)))))
% 
% %L=(1/(w2^2*C))
% %C1=((C)/((w1^2/w2^2)-1))
% ((j.*w.*(1/(w2^2*C))+1./(j.*w.*C))./(((1/(w2^2*C))+((C)/((w1^2/w2^2)-1)))/C-w.^2*(1/(w2^2*C))*((C)/((w1^2/w2^2)-1))))