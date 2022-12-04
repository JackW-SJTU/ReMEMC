% This function use linear regression only, other models should be
% transformed into linear model in advance
% STD_S=standard sample deviation
% N=number of replicates
function [Coefs,Rs]=MC_Regression(X,Y,Model,STD_S,noiselevelX,N,IterN,Method)

%default IterN=1000
if ~exist('IterN','var')
    IterN=1000;
end

% default Method='normal_Stepwise'
if ~exist('Method','var')
    Method='normal_Stepwise';
end

%Pre-allocate memory space
Coefs=zeros(IterN,size(x2fx(X,Model),2));
Rs=zeros(IterN,1);
% start regression (use parallel computering by default)
parfor i=1:IterN     
%     disp(i)
    y=Y+randn(length(Y),1).*STD_S./sqrt(N); % generate a set of possible Y values 
    x=X+randn(size(X)).*noiselevelX; % generate a set of possible X values 
    x=x2fx(x,Model);
    x(:,1)=[];
    switch Method % choose regression method
        case 'CV_Stepwise'
            [terms, fitinfo]=CVStepwise(x,y,'linear','loocv');
            inmodel=terms>0;
            temp=zeros(size(inmodel));
            temp(inmodel)=fitinfo.Coefficients';
            Coefs(i,:)=temp;
            Rs(i)=fitinfo.RValue;
        case 'normal_Stepwise'
            mdl=stepwiselm(x,y,'linear','upper','linear','criterion','adjrsquared','premove',-0.0001,'verbose',0);
            inmodel=mdl.Formula.InModel;
            temp=zeros(size(inmodel));
            temp(logical([1 inmodel(1:end-1)]))=mdl.Coefficients.Estimate;
            Coefs(i,:)=temp;
            Rs(i)=mdl.Rsquared.Ordinary;
        case 'lasso_Regression'
            [b,f]=lasso(x,y,'Alpha',1);
            A=find(f.DF<=length(Y).*0.8,1);
            Coefs(i,:)=[f.Intercept(A); b(:,A)]';
            R=corrcoef(y,f.Intercept(A)+x*b(:,A));
            Rs(i)=R(1,2);
    end
end
end
            
