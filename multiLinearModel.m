function y=multiLinearModel(Coefs,model,X)
X4P=x2fx(X,model);
YP=X4P*Coefs';
y=sqrt((prctile(YP,75,2)-prctile(YP,25,2)).^2+prctile(YP,50,2).^2)./sqrt(2);
end