clc;close all;clear all;
M=20;
N=8000;
noise=0.1;
A=randn(M,N);
I1=1;
I2=2;
lambda=0.1;

b=A(:,I1)+A(:,I2)+noise*randn(M,1);

AA= [A,-1*ones(size(A,1),1),zeros(size(A,1),1);...
    -A,-1*ones(size(A,1),1),zeros(size(A,1),1)];
Aeq=[ones(1,size(A,2)),0,-1];
beq=0;
bb=[b;-b];
lb=[zeros(size(A,2),1);0;0];
hb=[1*ones(size(A,2),1);inf;inf];
ff=[zeros(size(A,2),1);1;lambda];
    
x = linprog(ff,AA,bb,Aeq,beq,lb,hb);
stem(x(1:N));

x0=zeros(N+2,1);
x0(1)=1;
x0(2)=1;
x0(end)=2;
x0(end-1)=2;
%%
xc=nan(N,1);
for I=1:N,
    cf=corrcoef(b,A(:,I));
    xc(I)=abs(cf(1,2));
end