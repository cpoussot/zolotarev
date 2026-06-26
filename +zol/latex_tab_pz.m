function str = latex_tab_pz2(ZER,name)
    format shortE
    n       = size(ZER,1);
    name    = regexprep(name, 'e\+?(-?\d+)', ' \\cdot 10^{$1}');
    str_top_add1 = [];
    for i = 1:numel(name)
        str_top_add1 = [str_top_add1 [' & \multicolumn{1}{c}{\text{' name{i} '}} ']];
    end
    %str_top_add2 = repmat(' & \text{zeros} ',1,size(ZER,2)-1);
    str_top1 = [ str_top_add1 '  \\ \hline '];
    %str_top2 = ['\text{\#} & \text{zeros} ' str_top_add2 ' \\ \hline '];
    APPROX = [];
    %ZER = round(ZER,3);
    for i = 1:size(ZER,2)
        APPROX = [APPROX ZER(:,i)];
    end
    %APPROX = vpa(APPROX,2)
    str = latex(vpa([(1:n).' APPROX],3));
    %str = latex([(1:n).' APPROX])
    str = regexprep(str, 'e\+?(-?\d+)', ' \\cdot 10^{$1}');
    str = strrep(str,'\left(','');
    str = strrep(str,'\right)','');
    str = strrep(str,'\begin{array}{ccc','\begin{array}{c||cc');
    str = strrep(str,'cc}',['c>{\columncolor{colhl}}c} ' str_top1  ]);
    str = strrep(str,'cc','c|c|');
    str = strrep(str,'\end{array}',' \\ \hline \end{array}');
    str = strrep(str,'\','\\');
    str = strrep(str,'\mathrm{i}','\imath');
end