clear;close all;
load transonic_duct.mat;
nX=size(inputs,2);
nY=size(inputs,3);
X=reshape(inputs(18,:,:),nX*nY,1);
Y=reshape(inputs(19,:,:),nX*nY,1);
Z=reshape(inputs(20,:,:),nX*nY,1);
cd STLRead;
out=stlread('../out2.stl');
cd ..
trisurf(out.faces,out.vertices(:,1),out.vertices(:,2),out.vertices(:,3),out.faces*0+.5)
hold on
scatter3(X,Y,Z,'b.')
