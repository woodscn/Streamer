clear;
exact=load('dat1D_Rie_Test_10_X_exact.dat');
numx=load('dat1D_Rie_Test_10_X.dat');
x=numx(18,:);
dx=numx(2,:);
ux=numx(3,:);
px=numx(1,:);
% numy=load('dat1D_Rie_Test_12_Y.dat');
% dy=numy(2,:);
% numz=load('dat1D_Rie_Test_12_Z.dat');
% dz=numz(2,:);
pe=exact(1,:);
de=exact(2,:);
ue=exact(3,:);
plot(x,de,'b',x,dx,'b.',x,ue,'g',x,ux,'g.',x,pe,'r',x,px,'r.')
   