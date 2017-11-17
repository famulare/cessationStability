function h=plotWithCI(x,p,ci,col,linestyleStr,markersize,solidMarkerFlag,dottedCIFlag)
if nargin<4
    col='b';
end
if nargin<5
    linestyleStr='-';
end
if nargin<6
    markersize=10;
end
if nargin<7
    solidMarkerFlag=false;
end
if nargin<8
    dottedCIFlag=true;
end
    
hold on;
if solidMarkerFlag
    plot(x,p,linestyleStr,'linewidth',1.5,'color',col,'markerfacecolor',col)
else
    plot(x,p,linestyleStr,'linewidth',1.5,'color',col)
end

if dottedCIFlag
    for n=1:2
       plot(x,ci(:,n),'--','color',col,'markersize',markersize);
    end
else
    for n=1:length(p)
       plot(x(n)*[1,1],ci(n,:),regexprep(linestyleStr,'[ox+*sdv^<>ph.]',''),'color',col,'markersize',markersize);
    end
end
end