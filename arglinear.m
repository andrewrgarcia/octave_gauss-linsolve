%arglinear.m Octave code
%Developer: Andrew Garcia, 2014

clear all
close all
clc

%define values of A b (these are the values we determined previously):

A=[1 0 0 1 0 0;-1 1 1 0 0 0;0 1 0 0 0 -1;420.3728 0 -2199.61 2689.53 0 0;-420.3728 327.1985 0 -2234.84 2689.081 0;0 327.1985 0 0 2236.64 -2599.346];
b=[22680 0 4536 2394973.98 0 1004737.1544]';
n=length(A);

%Triangular Echelon Form
for k=1:n-1
 for i=k+1:n
  m(i,k)=A(i,k)/A(k,k);
  b(i)=b(i)-m(i,k)*b(k);
  for j=1:n
   A(i,j)=A(i,j)-m(i,k)*A(k,j);
  end
 end
end
  
%Gauss-Jordan solution
x(n)=b(n)/A(n,n);
for i=n-1:-1:1
 sum=0;
 for j=i+1:n
  sum=sum+A(i,j)*x(j);
 end
 x(i)=(b(i)-sum)/A(i,i);
end

A=eye(length(A))

disp("solution vector b")
disp(x')
