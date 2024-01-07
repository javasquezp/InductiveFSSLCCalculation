clear all
close all
clc

CSTProject = '../CST/UnitcellMetal.cst';
CST = TCSTInterface();
CST.OpenProject(CSTProject);
TreeItem = '1D Results\S-Parameters\SZmax(1),Zmax(1)'; 
[frequency,S11,Zref,RunIDs,Info] = CST.Get1DResultFromTreeItem(TreeItem);
TreeItem = '1D Results\S-Parameters\SZmin(1),Zmax(1)'; 
[frequency,S21,Zref,RunIDs,Info] = CST.Get1DResultFromTreeItem(TreeItem);
w=frequency*2*pi*1e9;

 S21=abs(S21(1:end,1));
 S11=abs(S11(1:end,1));
 [A I]=min(S11);
 [A II]=max(S11(I:end));
  II = II + I - 1; %locate the exact position of the maximum
  S21=S21(1:II)
w=w(1:II);
 w1=2*pi*frequency(I(1))*1e9; %% S11 Pole
 w2=2*pi*frequency(II(1))*1e9; %% S11 cero
[fit_results gfo]=C_calculation(w,20*log10(S21))
C=fit_results.C;
%C=5.514475873546713e-13;
 L=(1./(C*((w1)^2)));
 L1=((1./(C*((w1)^2)))/(((w2)^2/(w1)^2-1)));
[fit_results gfo]=L_L1_calculations(w,20*log10((S21)), C)


bwdown=find(20*log10(abs(S11(:,1)))<-10,1);
bwup=find(20*log10(abs(S11(bwdown:end,1)))>-10,1)+bwdown;
fr=frequency(I(1));
FBW=(frequency(bwup)-frequency(bwdown))/(frequency(I(1)));

%C=(6*FBW)./(377*2*pi*fr*1e9);
%L=1./((2*pi*fr*1e9).^2.*C);
