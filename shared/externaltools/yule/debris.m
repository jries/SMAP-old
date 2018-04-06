function debrisi=debris(data)

means = mean(data);
stds  = std(data);

std1 = ~inrange(data, [means-stds*2; means+stds*2]);
std2 = ~inrange(data, [means-stds*3; means+stds*3]);

% if any feature is more than three std away from mean or two features that
% are two stds away
debrisi = (sum(std2, 2) | sum(std1, 2) >= 2);

% figure;
% hold on;
% scatter(data(:, 1) , data(:, 2), '.b');
% scatter(data(debrisi, 1) , data(debrisi, 2), 150, '.r');
end