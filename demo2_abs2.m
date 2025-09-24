clearvars; close all; clc
%%% Startup
set(groot,'DefaultFigurePosition', [200 100 1000 700]);
set(groot,'defaultlinelinewidth',2)
set(groot,'defaultlinemarkersize',4)
set(groot,'defaultaxesfontsize',18)
list_factory = fieldnames(get(groot,'factory'));index_interpreter = find(contains(list_factory,'Interpreter'));for i = 1:length(index_interpreter);set(groot, strrep(list_factory{index_interpreter(i)},'factory','default'),'latex');end
%%% 
CAS = '1b';
r   = 6;
%%%
[pts,val,data]  = zol.example(CAS);
[la,mu,W,V]     = zol.example2data(pts,val,data);
opt             = [];
opt.target      = r;
[h4,info]       = zol.loewner(la,mu,W,V,opt);
%%%
[x,idx] = sort(pts);
val     = val(idx);
for i = 1:length(x)
    h4_x(i)     = h4(x(i));
    h4opt_x(i)  = data.z4x{r}(x(i));
end
figure, hold on
plot(x,val,'.'), grid on
plot(x,h4opt_x,'--'),
plot(x,h4_x),
xlabel('$x$'), ylabel('$h(x)$'), 
legend({'Original','Optimal','Loewner'})

[hnum,hden,M]   = zol.get_numden(h4);
HNUM            = [0;hnum];
HDEN            = hden;
[data.z4{r}{1} HNUM]
[data.z4{r}{2} HDEN]
