function output = opti_process(initialDOE,n_iter,method,funtype,noiselevelX,noiselevelY,N)
% n_comb combinations per iteration
% X range [-xmin,xmax]
[n_comb,n_dim]=size(initialDOE);
switch funtype
    case 'hill'
        xmax=2;
        xmin=0;
    case 'deltah'
        xmax=5;
        xmin=1;
    case 'contact angle'
        xmax=10;
        xmin=1;
    otherwise
        xmax=2;
        xmin=-2;
end

warning off all
output=ones(n_iter,1)*999;
In={};
Out={};
In{1}=initialDOE; % random input
output(1)=min(blackbox_fun(cell2mat(In'),funtype)); % first true output
switch method
    case 'random'
        for i=2:n_iter
            In{i}=(xmax-xmin)*rand(n_comb,n_dim)+xmin; % random input of round i
            output(i)=min(blackbox_fun(cell2mat(In'),funtype)); % true output of round i
        end
        
    case 'regression'
        Out{1}=blackbox_fun(cell2mat(In'),funtype,noiselevelX,noiselevelY,N); % experiment
        for i=2:n_iter
            if n_comb*(i-1)<2*n_dim
                mdl=stepwiselm(cell2mat(In'),cell2mat(Out'),'criterion',...
                    'AdjRsquared','upper','linear','verbose',0); % generate model
            elseif n_comb*(i-1)<4*n_dim
                mdl=stepwiselm(cell2mat(In'),cell2mat(Out'),'criterion',...
                    'AdjRsquared','upper','interactions','verbose',0);
            else
                mdl=stepwiselm(cell2mat(In'),cell2mat(Out'),'criterion',...
                    'AdjRsquared','upper','quadratic','verbose',0);
            end
            
            temp=[];
            for k=1:n_comb
                x=fminsearch(@(x)regmodel(x,mdl),ones(1,n_dim).*((k-1).*xmax+(n_comb-k).*xmin)./(n_comb-1),optimset('Display','off')); % prediction
                temp=[temp; x];
            end
            temp(temp<xmin)=xmin;temp(temp>xmax)=xmax;
            temp=unique(temp,'rows');
            if size(temp,1)<n_comb
                temp=[temp;(xmax-xmin)*rand(n_comb-size(temp,1),n_dim)+xmin];
            end   
            
            In{i}=temp;
            Out{i}=blackbox_fun(In{i},funtype,noiselevelX,noiselevelY,N); % new experiment output
            output(i)=min(blackbox_fun(cell2mat(In'),funtype)); % record true output 
        end
        
    case 'lasso'
        Out{1}=blackbox_fun(cell2mat(In'),funtype,noiselevelX,noiselevelY,N); % experiment
        for i=2:n_iter
            if n_comb*(i-1)<2*n_dim
                model='linear';
            elseif n_comb*(i-1)<4*n_dim
                model='interactions';
            else
                model='quadratic';
            end
            X=x2fx(cell2mat(In'),model);
            X=X(:,2:end);
            
            [b,f]=lasso(X,cell2mat(Out'),'Alpha',1);
            A=find(f.DF<=length(cell2mat(Out'))-2*i,1);
            Intercept=f.Intercept(A); 
            Coefs=b(:,A);
            
            temp=[];
            for k=1:n_comb
                x=fminsearch(@(x)lassoLinearModel(Coefs,Intercept,model,x),ones(1,n_dim).*((k-1).*xmax+(n_comb-k).*xmin)./(n_comb-1),optimset('Display','off')); % prediction
                temp=[temp; x];
            end
            temp(temp<xmin)=xmin;temp(temp>xmax)=xmax;
            temp=unique(temp,'rows');
            if size(temp,1)<n_comb
                temp=[temp;(xmax-xmin)*rand(n_comb-size(temp,1),n_dim)+xmin];
            end   
            
            In{i}=temp;
            Out{i}=blackbox_fun(In{i},funtype,noiselevelX,noiselevelY,N); % new experiment
            output(i)=min(blackbox_fun(cell2mat(In'),funtype)); % record fmin         
        end
        
    case 'DE'
        Out{1}=blackbox_fun(cell2mat(In'),funtype,noiselevelX,noiselevelY,N); % iter1
        for i=2:n_iter
            temp=In{i-1}-0.75*(In{i-1}-[In{i-1}(2:end,:); In{i-1}(1,:)]); % mutation
            temp(temp<xmin)=xmin;temp(temp>xmax)=xmax;
            pp=rand(n_comb,n_dim); % crossover
            temp(pp>0.75)=In{i-1}(pp>0.75);
            In{i}=temp; % iter2
            Out{i}=blackbox_fun(In{i},funtype,noiselevelX,noiselevelY,N); % experiment
            In{i}(Out{i-1}<Out{i},:)=In{i-1}(Out{i-1}<Out{i},:); % competition
            Out{i}(Out{i-1}<Out{i},:)=Out{i-1}(Out{i-1}<Out{i},:);
            
            output(i)=min(blackbox_fun(cell2mat(In'),funtype)); % record fmin
        end
        
    case 'STRICT'
        Out{1}=blackbox_fun(cell2mat(In'),funtype,noiselevelX,noiselevelY,N); % experiment
        for i=2:n_iter
            % regression
            if n_comb*(i-1)<2*n_dim
                model='linear';
            elseif n_comb*(i-1)<4*n_dim
                model='interactions';
            else
                model='quadratic';
            end
            X=x2fx(cell2mat(In'),model);
            X=X(:,2:end);
            
            [b,f]=lasso(X,cell2mat(Out'),'Alpha',0.5);
            A=find(f.DF<=length(cell2mat(Out'))-2,1);
            Intercept=f.Intercept(A); 
            Coefs=b(:,A);

            x=fminsearch(@(x)lassoLinearModel(Coefs,Intercept,model,x),zeros(1,n_dim),optimset('Display','off')); % prediction          
            temp=[x;repmat(x,[2 1])+0.5*(rand(2,n_dim)-0.5)];
            temp(temp<xmin)=xmin;temp(temp>xmax)=xmax;

            % DE
            temp2=[cell2mat(Out') cell2mat(In')];
            temp2=sortrows(temp2,1);
            temp2=temp2(1:n_comb-2,2:end);
            new=temp2(1:n_comb-3,:)-0.75*(temp2(1:n_comb-3,:)-temp2(2:n_comb-2,:)); % mutation
            new(new<xmin)=xmin;new(new>xmax)=xmax;
            temp2=temp2(1:n_comb-3,:);
            pp=rand(n_comb-3,n_dim); % crossover
            new(pp>0.75)=temp2(pp>0.75);
            
            In{i}=[temp; new]; % iter2            
                       
            Out{i}=blackbox_fun(In{i},funtype,noiselevelX,noiselevelY,N); % new experiment         
            output(i)=min(blackbox_fun(cell2mat(In'),funtype)); % record fmin 
        end
        
        
    case 'Monte Carlo'
        [Out{1}, STD{1}]=blackbox_fun(cell2mat(In'),funtype,noiselevelX,noiselevelY,N); % experiment
        for i=2:n_iter % start with 2nd round
            if n_comb*(i-1)<2*n_dim
                model='linear';
            elseif n_comb*(i-1)<4*n_dim
                model='interactions';
            else
                model='quadratic';
            end
            X=x2fx(cell2mat(In'),model);
            X=X(:,2:end); 
            
            Coefs=MC_Regression(X,cell2mat(Out'),'linear',cell2mat(STD').*N./(N-1),noiselevelX,N,100,'lasso_Regression');
            
            temp=[];
            
            for k=1:n_comb
                x=fminsearch(@(x)multiLinearModel(Coefs,model,x),ones(1,n_dim).*((k-1).*xmax+(n_comb-k).*xmin)./(n_comb-1),optimset('Display','off')); % prediction
                temp=[temp; x];
            end
            temp(temp<xmin)=xmin;temp(temp>xmax)=xmax;
            temp=unique(temp,'rows');
            if size(temp,1)<n_comb
                temp=[temp;(xmax-xmin)*rand(n_comb-size(temp,1),n_dim)+xmin];
            end            
            In{i}=temp;
            [Out{i},STD{i}]=blackbox_fun(In{i},funtype,noiselevelX,noiselevelY,N); % new experiment
            output(i)=min(blackbox_fun(cell2mat(In'),funtype)); % record fmin         
        end
        
    case 'GP'
        Out{1}=blackbox_fun(cell2mat(In'),funtype,noiselevelX,noiselevelY,N); % experiment
        for i=2:n_iter % start with 2nd round
            GPMDL=fitrgp(cell2mat(In'),cell2mat(Out'));
            
            temp=[];
            for k=1:n_comb
                x=fminsearch(@(x)GPPredict(GPMDL,x),ones(1,n_dim).*((k-1).*xmax+(n_comb-k).*xmin)./(n_comb-1),optimset('Display','off')); % prediction
                temp=[temp; x];
            end
            temp(temp<xmin)=xmin;temp(temp>xmax)=xmax;
            temp=unique(temp,'rows');
            if size(temp,1)<n_comb
                temp=[temp;(xmax-xmin)*rand(n_comb-size(temp,1),n_dim)+xmin];
            end            
            In{i}=temp;
            Out{i}=blackbox_fun(In{i},funtype,noiselevelX,noiselevelY,N); % new experiment
            output(i)=min(blackbox_fun(cell2mat(In'),funtype)); % record fmin
        end
               
end


end

function Y=regmodel(x,mdl)
if sum(x>2)+sum(x<-2)>0
    Y=9999;
else
    Y=predict(mdl,x);
end
end

function y=linearModel(beta,X)
y=X*beta(2:end)+beta(1);
end

function y=lassoLinearModel(Coefs,Intercept,model,X)
X=x2fx(X,model);
y=X*[Intercept; Coefs];
end

function UCB=multiLinearModel(Coefs,model,X)
X4P=x2fx(X,model);
YP=X4P*Coefs';
UCB=abs(prctile(YP,75,2)-prctile(YP,25,2))+prctile(YP,50,2);
end

function y=GPPredict(MDL,X)
    [t1,t2]=predict(MDL,X);
    y=t1+t2;
end