function prim_cons_check()
test = rand(21,1);
test(21)=Jacobian(test);
test
primtocons(test)
end

function out = primtocons(in)
J = Jacobian(in);
e = .5*sum(in(3:5).^2)+in(1)/(in(2)*.4);
out = in;
out(1) = J*in(2);
out(2) = J*in(2)*in(3);
out(3) = J*in(2)*in(4);
out(4) = J*in(2)*in(5);
out(5) = J*in(2)*e;
end

function out = constoprim(in)
J = Jacobian(in);
out = in;
out(2) = in(1)/J;
out(3) = in(2)/in(1);
out(4) = in(3)/in(1);
out(5) = in(4)/in(1);
out(1) = (in(5)/in(1)-.5*sum(out(3:5).^2))*out(2)*.4;
end

function out = Jacobian(in)
a=in( 6); b=in( 7); c=in( 8);
l=in( 9); m=in(10); n=in(11);
p=in(12); q=in(13); r=in(14);
out=a*m*r+b*n*p+c*l*q...
    -a*n*q-b*l*r-c*m*p;
end