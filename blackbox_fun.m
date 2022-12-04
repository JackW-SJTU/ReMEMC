function [Y,STD]=blackbox_fun(X,funtype,noiselevelX,noiselevelY,N)
if ~exist('noiselevelX','var')
    noiselevelX=0;
end

if ~exist('noiselevelY','var')
    noiselevelY=0;
end

if ~exist('N','var')
    N=1;
end

X=X+max(abs(X),1).*randn(size(X)).*noiselevelX;
Y = zeros(size(X,1),1);
for i = 1:size(X,1)
    p = X(i,:);
    switch funtype
        case 'ackley'
            Y(i) = 20 + exp(1) ...
                - 20 * exp( -0.2 * sqrt( (1/length(p)) * sum(p .^ 2)))...
                - exp((1/length(p)) * sum(cos(2*pi .* p)));
        case 'rosenbrock'
            temp=0;
            for ii=1:length(p)-1
                temp = temp + 100*( p(ii+1) - p(ii).^2 ).^2 + (p(ii)-1).^2;
            end
            Y(i) = temp;
        case 'hill'
            p_eff=zeros(size(p));
            p_adv=p_eff;
            a=linspace(-0.9,0.9,length(p));
            for ii=1:length(p)
                temp=p;
                temp(ii)=[];
                p_eff(ii)=p(ii).*prod(1+a(ii).*temp./(0.5+temp));
                p_adv(ii)=p(ii).*prod(1+a(ii).*temp./(1+temp));
            end
            Y(i)=prod(1./(1+p_eff.^2))-prod(1./(1+p_adv.^3));  
    end
end
Y=Y+max(abs(Y),1).*randn(size(Y,1),N).*noiselevelY;
STD=std(Y,[],2);
Y=mean(Y,2);
end
