function [f1, f2]=fitDiff(x,y,o,fast)
% o - offset
    if nargin<3
       o = 0; 
    end
    s=size(x);
    N=max(s);
    if min(s)>1
        disp('Incorrect data format');
    end
    if (N==5 || N==11) && o==0 && nargin>3
       [f1, f2] = fastDiff(y,N);
       return
    end
    Dx=0.5*(x(N)-x(1));
    if o<0
       x=1-(x(end+o)-x)/Dx;
    else
       x=(x-x(1+o))/Dx-1;
    end
    
    M=3; % number of terms in fit
    c=4;
    w=exp(-c*x.^2);
    A=zeros(M);
    %Q=zeros(M,N);
    b=zeros(M,1);
    for n=1:M
        %Q(n,:) = w.*x.^(n-1);
        b(n)=sum(w.*x.^(n-1).*y);
        A(n,n)=sum(w.*x.^(n+n-2));
        for k=n+1:M
            A(n,k)=sum(w.*x.^(k+n-2));
            A(k,n)=A(n,k);
        end
    end
    %A
    %Q
    a=A\b;
    f1=a(2)/Dx;
    f2=2*a(3)/Dx^2;
end