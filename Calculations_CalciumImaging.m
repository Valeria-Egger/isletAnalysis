function Calculations_CalciumImaging
[file, path] = uigetfile("*.csv", "Select a csv file");
dataset = readtable(fullfile(path, file));
[~, baseName, ~] = fileparts(file);

T = dataset;
columns = startsWith(T.Properties.VariableNames, "Mean");
selected = T{:, columns};

min_Values = min(selected, [], 1);
max_Values = max(selected, [], 1);

scaled = (selected-min_Values)./(max_Values-min_Values);

T{:, columns} = scaled;

longT = stack(T, T.Properties.VariableNames(columns), ... 
    'NewDataVariableName', 'value', ... 
    'IndexVariableName', 'trace');
unique(longT.trace)

figure
hold on
traces = categories(longT.trace);

for k = 1:numel(traces)
    idx = longT.trace == traces{k};
    y = longT.value(idx);
    x = 1:numel(y);   % force each trace to start at x=1
    plot(x, y)
end

title('Scaled traces')
hold off

if ~exist("scaled", "dir")
    mkdir("scaled");
end

output = T(:, ['Var1', T.Properties.VariableNames(columns)]);
outputName = baseName + "_scaled.csv";
writetable(output, fullfile("scaled", outputName))
