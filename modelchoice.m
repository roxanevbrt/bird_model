function [L_f,L_tbt,L_tmt,rK1,rA,rP,L_h,name, masse, rK2, L_fibula, dav,dar, Lg, Tg] = modelchoice(i)
    T=readtable('bdd_article.csv','ReadVariableNames',false);
    L_f=table2array(T(i,3));
    L_tbt=table2array(T(i,4));
    L_tmt=table2array(T(i,5));
    rK1=table2array(T(i,6));
    rA=table2array(T(i,7));
    rP=table2array(T(i,8));
    L_h=table2array(T(i,9));
    name=char(table2cell(T(i,2)));
    masse=table2array(T(i,12));
    L_fibula=table2array(T(i,13));
    rK2=table2array(T(i,14));
    dav=table2array(T(i,10));
    dar=table2array(T(i,11));
    Lg=table2array(T(i,15));
    Tg=table2array(T(i,16));
end

