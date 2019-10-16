% Reference: http://read.pudn.com/downloads446/sourcecode/math/1879345/GMRES.m__.htm

function [X,r,t]=GMRESfromWeb(A,m,B,s) 
%����Global full orthogonalization method�ⷽ��AX=B 
% ������ 
% A�� N  x N ��ʵ������ 
% m�� Krylov�ռ�ά�� 
% B�� N  x s ��ʵ������ 
% s:  ������������B������������ 
% ����� 
% r��  ���� 
% t��  �������� 
tic                %ͳ�Ƹó�������ʱ�� 
[n,~]=size(A); 
%%������ʼ��             
X0=zeros(n,s);    % X0: ��ʼ��  
epsi=1e-8;     % epsi:���ȿ��Ʋ��� 
N=50;             % N: ���������� 
t=0; 
e=1; 
I=eye(s); 
R=[]; 
R0=B-A*X0; 
 
 while t<N && e>epsi 
      
    V=zeros(n,(m+1)*s); 
    H=zeros(m+1,m); 
    V(:,1:s)=R0/norm(R0,'fro'); 
    for j=1:m 
        W=A*V(:,(j-1)*s+1:j*s); 
        for i=1:j 
           H(i,j)=trace(W'*V(:,(i-1)*s+1:i*s)); 
           W=W-H(i,j)*V(:,(i-1)*s+1:i*s); 
        end 
        H(j+1,j)=norm(W,'fro'); 
        if H(j+1,j)==0 
            break 
        end 
        V(:,j*s+1:(j+1)*s)=W/H(j+1,j); 
    end 
    y=H\(norm(R0,'fro')*speye(m+1,1)); 
    V1=V(:,1:m*s); 
    X=X0+V1*kron(y,I); 
    X0=X; 
    R0=B-A*X0; 
    r=norm(R0,'fro'); 
    e=r; 
%     R=[R,log10(r)];      %�˴����Ա���ͼ 
    t=t+1; 
     
end 
% t1=1:t; 
% plot(t1, R, 'r*--'); 

toc 
end 
