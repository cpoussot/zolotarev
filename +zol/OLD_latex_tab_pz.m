function str = latex_tab_pz(ZER,POL,name)

format shortE
n       = size(ZER,1);
name    = regexprep(name, 'e\+?(-?\d+)', ' \\cdot 10^{$1}');
str_top_add1 = [];
for i = 1:numel(name)
    str_top_add1 = [str_top_add1 [' & \multicolumn{2}{c}{\text{' name{i} '}} ']];
end
str_top_add2 = repmat(' & \text{zeros} & \text{poles} ',1,size(ZER,2)-1);
str_top1 = [ str_top_add1 '  \\ \hline '];
str_top2 = ['\text{\#} & \text{zeros} & \text{poles} ' str_top_add2 ' \\ \hline '];
APPROX = [];
for i = 1:size(ZER,2)
    APPROX = [APPROX [ZER(:,i) POL(:,i)]];
    %APPROX = [APPROX [real(approxN(:,i)); real(approxD(:,i))] [imag(approxN(:,i)); imag(approxD(:,i))]];
end
str = latex(vpa([(1:n).' APPROX],2));
%
str = regexprep(str, 'e\+?(-?\d+)', ' \\cdot 10^{$1}');
%str = strrep(str,['s^' num2str(n_num)],[' \hline s^' num2str(n_num)]);
%str = strrep(str,['s^' num2str(n_den)],[' \hline s^' num2str(n_den)]);
str = strrep(str,'\left(','');
str = strrep(str,'\right)','');
str = strrep(str,'\begin{array}{ccc','\begin{array}{c||c>{\columncolor{colhl}}c|');
str = strrep(str,'\end{array}',' \\ \hline \end{array}');
str = strrep(str,'cc','c>{\columncolor{colhl}}c|');
str = strrep(str,'c>{\columncolor{colhl}}c|}',['c>{\columncolor{colhl}}c|} ' str_top1  str_top2]);
str = strrep(str,'c|c>','c||c>');
str = strrep(str,'\','\\');
str = strrep(str,'\mathrm{i}','\imath');

%%% OLD
%str = strrep(str,'\begin{array}{ccc','\begin{array}{c||c:>{\columncolor{colhl}}c|');
%str = strrep(str,'\end{array}',' \\ \hline \end{array}');
%str = strrep(str,'cc','c:>{\columncolor{colhl}}c|');
%str = strrep(str,'c:>{\columncolor{colhl}}c|}',['c:>{\columncolor{colhl}}c|} ' str_top1  str_top2]);
%str = strrep(str,'c|c:>','c||c:>');


