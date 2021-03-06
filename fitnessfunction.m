function f= fitnessfunction(x)
r1 = x(1);
r2=  x(2);
x1 = x(3);
x2=  x(4);
xm=x(5);
zr=r2+x2*i;
zs=r1+x1*i;
z= zs + ((zr*xm)/(zr+xm));
vph = 400;
psi=cos(angle(vph/zs));
s=0.04;
ws=25;
Tn=506.7196;
Tk=421.9163;
Td=1203.0386;
f1 = (((vph*r2/s))/( ws*(r1+r2/s)^2 + (x2+x1)^2)) - Tn;
f2 = ((vph^2*r2)/ (ws*((r1+r2)^2 + (x1+x2)^2)))  - Tk;
f3 = ((vph^2)/ (2*ws*(r1+sqrt(r1^2+(x1+x2)^2))))  - Td;
f4= psi - 0.4272;
f= f1^2 + f2^2 +f3^2+f4^2;
return
