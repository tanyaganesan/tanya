function post_processor_plot(X,Y,z)
% This function 
% LX = linspace(0,Lx,m+1);
% 
% LY = linspace(0,Ly,n+1);
% for I = 1:n
%     zi(I,:) = interp1(X(I,:),z(I,:),LX,'linear','extrap');
% end
% % % 
% z1 = zi;
% for J = 1:m+1
%     zi(n+1,J) = interp1(Y(:,1)',z1(:,J)',LY(end),'linear','extrap');
% end
% % 
% [x,y]  = meshgrid(LX,LY);
% 
surface(X,Y,z)
% shading  interp
shading FLAT
colorbar
% colormap gray
axis ij
return