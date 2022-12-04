% run this code to initiate in silico benchmark
clear
method={'Monte Carlo','lasso','regression'};
funtype={'hill','ackley','rosenbrock'};
n_dim=[5 10 15];
n_iter=5;
n_exp=3;
n_comb=48;
noiselevelX=0.05;
noiselevelY=0.1;
result=cell(length(method),length(funtype));
for k=n_dim
    for ii=1:length(funtype)
        switch funtype{ii}
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
        initialDOE=cellfun(@(x) (xmax-xmin)*rand(n_comb,k)+xmin,cell(50,1),'UniformOutput',false);
        for i=1:length(method)
            output=ones(n_iter,50)*999;
            disp([i ii])
            parfor iii=1:50
                output(:,iii)=opti_process(initialDOE{iii},n_iter,method{i},funtype{ii},noiselevelX,noiselevelY,n_exp); %(n_dim,n_iter,method,funtype,noiselevelX,noiselevelY,N,xmin,xmax)
            end
            result{i,ii}=output;
        end
    end
    save(['dim_' num2str(k) '_0.1.mat'])
end


%%
clear
method={'Monte Carlo','lasso','regression'};
funtype={'hill','ackley','rosenbrock'};
n_dim=[5 10 15];
n_iter=5;
n_exp=3;
n_comb=48;
noiselevelX=0.05;
noiselevelY=0.3;
result=cell(length(method),length(funtype));
for k=n_dim
    for ii=1:length(funtype)
        switch funtype{ii}
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
        initialDOE=cellfun(@(x) (xmax-xmin)*rand(n_comb,k)+xmin,cell(50,1),'UniformOutput',false);
        for i=1:length(method)
            output=ones(n_iter,50)*999;
            disp([i ii])
            parfor iii=1:50
                output(:,iii)=opti_process(initialDOE{iii},n_iter,method{i},funtype{ii},noiselevelX,noiselevelY,n_exp); %(n_dim,n_iter,method,funtype,noiselevelX,noiselevelY,N,xmin,xmax)
            end
            result{i,ii}=output;
        end
    end  
    save(['dim_' num2str(k) '_0.3.mat'])
end

%mailto('wangboqian_723@163.com')

%%
clear
method={'Monte Carlo','lasso','regression'};
funtype={'hill','ackley','rosenbrock'};
n_dim=[5 10 15];
n_iter=5;
n_exp=3;
n_comb=48;
noiselevelX=0.05;
noiselevelY=0.5;
result=cell(length(method),length(funtype));

for k=n_dim
    for ii=1:length(funtype)
        switch funtype{ii}
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
        initialDOE=cellfun(@(x) (xmax-xmin)*rand(n_comb,k)+xmin,cell(50,1),'UniformOutput',false);
        for i=1:length(method)
            output=ones(n_iter,50)*999;
            disp([i ii])
            parfor iii=1:50
                output(:,iii)=opti_process(initialDOE{iii},n_iter,method{i},funtype{ii},noiselevelX,noiselevelY,n_exp); %(n_dim,n_iter,method,funtype,noiselevelX,noiselevelY,N,xmin,xmax)
            end
            result{i,ii}=output;
        end
    end
    save(['dim_' num2str(k) '_0.5.mat'])
end
        
        
