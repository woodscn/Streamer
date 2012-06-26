close all;
nodes = [NodeA;NodeB;NodeC;NodeA];
plot3(nodes(:,1),nodes(:,2),nodes(:,3))
hold on
plot3(Point(1),Point(2),Point(3),'*')
Point2 = V2+NodeA;
plot3(Point2(1),Point2(2),Point2(3),'x')
axis equal
view(gca,Normal)
fprintf('Normalcy of Normal vector = %f \n\n',...
        (Normal*(NodeA-NodeB)')^2+(Normal*(NodeA-NodeC)')^2)
Tangential = (.3*(NodeC-NodeA)+.3*(NodeB-NodeA));
fprintf('Normal * Tangential = %f\n\n',Normal*Tangential')
%pause
view(gca,Tangential)
axis equal
load round_duct.mat;
nX=size(inputs,2);
nY=size(inputs,3);
X=reshape(inputs(18,:,:),nX*nY,1);
Y=reshape(inputs(19,:,:),nX*nY,1);
Z=reshape(inputs(20,:,:),nX*nY,1);
cd STLRead;
out=stlread('../out.stl');
cd ..
close all
fprintf('Don''t try to use the index value\n returned from Fortran here; \n the indices are different.\n\n')
figure
colors=out.faces*0+.5;
% colors(Inda,:)=[1,0,0];
trisurf(out.faces,out.vertices(:,1),out.vertices(:,2),...
    out.vertices(:,3),colors)
hold on
scatter3(X,Y,Z,'b.')
plot3(Point(1),Point(2),Point(3),'r*')
plot3(nodes(:,1),nodes(:,2),nodes(:,3),'k*','LineWidth',.5)
NodeAvg=(NodeA+NodeB+NodeC)/3.;
% quiver3(NodeAvg(1),NodeAvg(2),NodeAvg(3),Normal(1),Normal(2),Normal(3))