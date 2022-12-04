clear
for dim=[5 10 15]
    for nl=[0.1 0.3 0.5]
        load(['dim_' num2str(dim) '_' num2str(nl) '.mat'])
        %c=prism;
        funtype={'Hill','Ackley','Schwefel 1.2','Schwefel 2.2','Rosenbrock'};
        method={'Monte Carlo','GP','Lasso','Stepwise'};
        lines=lines;
        funtype=funtype([1 2 5]);
        result=result(:,[1 2 5]);
        method=method([1 3 4]);
        result=result([1 3 4],:);
        close all
        figure()
        for ii=1:length(funtype)
            subplot(1,length(funtype),ii)
            box on
            for i=1:length(method)
                plot([1:0.1:5],interp1(1:0.5:5,interp1(1:5,prctile(result{i,ii}',50),1:0.5:5,'linear'),1:0.1:5,'makima'),'color',lines(i,:),'linewidth',1);     

                xlim([1 5])
                set(gca,'ticklength',[0 0])
                %'color',c(i,:));
                hold on
                % if ii==1
                    ylabel('Function Value','fontsize',12)
                % elseif ii==5
                    xlabel('Iterations','fontsize',12)
                % end
                %plot(prctile(result{i,ii}',25),':','linewidth',1,'color',c(i,:));
                %plot(prctile(result{i,ii}',75),':','linewidth',1,'color',c(i,:));
            end
            for i=1:length(method)
                p25sm=interp1(1:0.5:5,interp1(1:5,prctile(result{i,ii}',25),1:0.5:5,'linear'),1:0.1:5,'makima');
                p75sm=interp1(1:0.5:5,interp1(1:5,prctile(result{i,ii}',75),1:0.5:5,'linear'),1:0.1:5,'makima');
        %        patch([1:5 5:-1:1],[smoothdata(prctile(result{i,ii}',25),'movmean',1) flip(smoothdata(prctile(result{i,ii}',75),'movmean',1))],lines(i,:),...
        %            'edgecolor',lines(i,:),'linestyle',':')
               patch([1:0.1:5 5:-0.1:1],[p25sm flip(p75sm)],lines(i,:),...
                   'edgecolor',lines(i,:),'linestyle',':')
                alpha(0.3)

            end


            hold off
            title(funtype{ii},'fontsize',12)
            legend(method)
        end
        set(gcf,'unit','inches')
        set(gcf,'position',[1 3 10 2.8])
        print(['dim_' num2str(dim) '_' num2str(nl) '.png'],'-dpng','-r600')
        close all
    end
end
