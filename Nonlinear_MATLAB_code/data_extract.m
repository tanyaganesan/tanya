clc
clear all
clc;
X = [0 0 0 0;0 2 3 0;0 4 5 0;0 0 0 0] %example matrix
% x2 = conv2(x,ones(2),'valid') %windowed mean of 3x3 grid
% x(2:end-1,2:end-1) = x2 %insert

  FWIN = 'hamming';
      F = [13 8];
      N = 4;
      Pnoise = 0.30;
      PNaNs  = 0.20;
%       X = peaks(N);                                     % original
      Y = X + ((rand(size(X))-0.5)*2)*max(X(:))*Pnoise; % add noise
      Y(round(1 + (N^2-1).*rand(N^2*PNaNs,1))) = NaN;   % add NaNs
      [Z0,W] = ndnanfilter(Y,FWIN,F)                  % filters
      Z1 = Z0; Z2 = Y; inan = isnan(Y);
      Z1(inan) = NaN;
      Z2(inan) = Z0(inan);
      subplot(231), imagesc(X), clim = caxis; axis equal tight
                    title('Original data')
      subplot(232), imagesc(Y),  caxis(clim), axis equal tight
                    title('Data + NOISE + NaNs')
      subplot(234), imagesc(Z0), caxis(clim), axis equal tight
                    title('FILTERS + NaNs interpolation')
      subplot(235), imagesc(Z1), caxis(clim), axis equal tight
                    title('FILTERS ignoring NaNs')
      subplot(236), imagesc(Z2), caxis(clim), axis equal tight
                    title('GAP-filling with interpolated NaNs')
      subplot(233), imagesc(-F(1):F(1),-F(2):F(2),W), axis equal tight,
                     title([upper(FWIN) ' 2D window']), view(2)


break

for J = 1:1000
count = num2str(J);
ext = '.txt';
filename = ['disp' count ext];
filename1 = ['mom' count ext];
filename2 = ['sf' count ext];
% 
A = load(filename);
B = load(filename1);
C = load(filename2);
% 
DISP(:,J) = A;
MOM(:,J) = B;
SHE(:,J) = C;
end
for I = 1:51
    DISP_(I,1)   = mean(DISP(I,:));
    DISP_std(I,1) = std(DISP(I,:));
end
for I = 1:51
    MOM_(I,1)   = mean(MOM(I,:));
    Mom_std(I,1) = std(MOM(I,:));
end
for I = 1:51
    SHE_(I,1)   = mean(SHE(I,:));
    SF_std(I,1) = std(SHE(I,:));
end
% DISP_
% [I1,J1] = find(DISP_ == min(DISP_));
% 
% [I2,J2] = find(MOM_ == max(MOM_));
% 
% [I3,J3] = find(SHE_ == max(SHE_));
% 
% [DISP_std(I1);Mom_std(I2);SF_std(I3)]

X0 = linspace(0,15,51);
figure1 =  figure;
axes1 = axes('FontName','Times New Roman','FontSize',12,'Parent',figure1);
grid(axes1,'on');
box(axes1,'on');
hold(axes1,'all');
xlabel('Beam length (m)','FontName','Times New Roman','FontSize',12)
ylabel('Deflection (m)','FontName','Times New Roman','FontSize',12)
set(plot(X0,DISP),'Color',[0.5 0.5 0.5],'LineWidth',1);
set(plot(X0,DISP_),'Color',[1 0 0],'LineWidth',3);

figure2 =  figure;
axes1 = axes('FontName','Times New Roman','FontSize',12,'Parent',figure2);
grid(axes1,'on');
box(axes1,'on');
hold(axes1,'all');
xlabel('Beam length (m)','FontName','Times New Roman','FontSize',12)
ylabel('Bending Moment (kNm)','FontName','Times New Roman','FontSize',12)
set(plot(X0,MOM),'Color',[0.5 0.5 0.5],'LineWidth',1);
set(plot(X0,MOM_),'Color',[1 0 0],'LineWidth',3);

figure3 =  figure;
axes1 = axes('FontName','Times New Roman','FontSize',12,'Parent',figure3);
grid(axes1,'on');
box(axes1,'on');
hold(axes1,'all');
xlabel('Beam length (m)','FontName','Times New Roman','FontSize',12)
ylabel('Shear force (kN)','FontName','Times New Roman','FontSize',12)
hold(axes1,'all');
set(plot(X0,SHE),'Color',[0.5 0.5 0.5],'LineWidth',1);
set(plot(X0,SHE_),'Color',[1 0 0],'LineWidth',3);
