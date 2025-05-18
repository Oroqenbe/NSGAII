clc; clear; close all;

%% ========== 加载数据 ===========
demand_data = readtable('模拟参数数据_含单位.xlsx', 'Sheet', '需求点数据', 'VariableNamingRule', 'preserve');
candidate_data = readtable('模拟参数数据_含单位.xlsx', 'Sheet', '候选站点数据', 'VariableNamingRule', 'preserve');

num_demand = height(demand_data);
num_sites = height(candidate_data);

h_j = demand_data{:, "运输需求量 h_j (件)"}; 
coord_demand = demand_data{:, ["经度 (°)", "纬度 (°)"]};

F_i = candidate_data{:, "固定建设成本 F_i (万元)"}; 
s_i_base = candidate_data{:, "最大配送能力 s_i (件)"}; 
delta_i_base = candidate_data{:, "禁飞区 δ_i (0/1)"}; 
gamma_i = candidate_data{:, "气象风险系数 γ_i"}; 
w_air = candidate_data{:, "空域风险权重 w^airspace"}; 
w_weather = candidate_data{:, "气象风险权重 w^weather"}; 
coord_site = candidate_data{:, ["经度 (°)", "纬度 (°)"]}; 
site_names = string(candidate_data{:,1});

%% ========== 参数设置 ==========
num_gen = 200;
pop_size = 100;
min_sites = 4;
max_sites = 7;
max_radius_km = 10;

p_cross = 0.9;
p_mut = 0.1;

% 并行常量变量
max_radius_km_constant = parallel.pool.Constant(max_radius_km);

%% ========== 初始化 & 基础优化运行 ========== 
pop = generateInitialPopulation(pop_size, num_sites, min_sites, max_sites, delta_i_base);
D_full = computeDistanceMatrix(coord_site, coord_demand);

objectives = zeros(pop_size, 3);
for gen = 1:num_gen
    valid_count = zeros(pop_size,1);
    parfor i = 1:pop_size
        y = pop(i, :);
        [obj, valid] = evaluateSolution(y, D_full, F_i, s_i_base, h_j, delta_i_base, gamma_i, w_air, w_weather, min_sites, max_sites, max_radius_km_constant.Value);
        objectives(i,:) = obj;
        if valid
            valid_count(i) = 1;
        end
    end
    fprintf("基础优化第 %d 代：有效解数量 = %d\n", gen, sum(valid_count));
    rank = paretoRank(objectives);
    pop = geneticOperators(pop, rank, p_cross, p_mut, num_sites);
end

visualizeParetoFront(objectives);
visualizeTopThreeSolutions_balanced(pop, objectives, D_full, coord_demand, coord_site, min_sites, max_sites);
exportBalancedSolutionsToExcel(pop, objectives, num_sites);

%% ========== 敏感性分析主循环 ==========
delta_set = [0, 1];
mu_set = [300, 500, 800];

base_output = 'output/sensitivity';
if ~exist(base_output, 'dir'); mkdir(base_output); end
summary_all = [];

for delta = delta_set
    for mu = mu_set
        fprintf('运行 δ=%d, μ=%d 的敏感性实验...\n', delta, mu);

        delta_i = delta_i_base;
        s_i = min(s_i_base, mu);

        pop = generateInitialPopulation(pop_size, num_sites, min_sites, max_sites, delta_i);
        D = D_full;  % 距离矩阵不变，重用即可

        objectives = zeros(pop_size, 3);
        for gen = 1:num_gen
            valid_count = zeros(pop_size,1);
            parfor i = 1:pop_size
                y = pop(i, :);
                [obj, valid] = evaluateSolution(y, D, F_i, s_i, h_j, delta_i, gamma_i, w_air, w_weather, min_sites, max_sites, max_radius_km_constant.Value);
                objectives(i,:) = obj;
                if valid
                    valid_count(i) = 1;
                end
            end
            rank = paretoRank(objectives);
            pop = geneticOperators(pop, rank, p_cross, p_mut, num_sites);
        end

        folder = fullfile(base_output, sprintf('delta%d_mu%d', delta, mu));
        if ~exist(folder, 'dir'); mkdir(folder); end

        rank = paretoRank(objectives);
        pareto_mask = (rank == 1);

        if sum(pareto_mask) == 0
            fprintf(' δ=%d, μ=%d 无可行解，跳过输出。\n', delta, mu);
            continue;
        end

        pareto_objectives = objectives(pareto_mask, :);

        % === 图像保存 ===
        figure;
        scatter3(objectives(:,1), objectives(:,2), objectives(:,3), 30, 'b', 'filled'); hold on;
        scatter3(pareto_objectives(:,1), pareto_objectives(:,2), pareto_objectives(:,3), 60, 'r', 'filled');
        xlabel('总成本（万元）'); ylabel('运输量距离（件·公里）'); zlabel('综合风险');
        title(sprintf('帕累托前沿图 δ=%d, μ=%d', delta, mu)); grid on;
        legend('所有解', '帕累托最优解');
        saveas(gcf, fullfile(folder, sprintf('pareto_delta%d_mu%d.png', delta, mu))); close;

        % === CSV 导出 ===
        T = table;
        T.delta = repmat(delta, size(objectives, 1), 1);
        T.mu = repmat(mu, size(objectives, 1), 1);
        T.cost = objectives(:, 1);
        T.distance = objectives(:, 2);
        T.risk = objectives(:, 3);
        for i = 1:size(pop, 2)
            T.(sprintf('site_%d', i)) = pop(:, i);
        end
        csv_filename = fullfile(folder, sprintf('results_delta%d_mu%d.csv', delta, mu));
        writetable(T, csv_filename);
        fprintf('已导出 CSV：%s\n', csv_filename);

        % === 地图输出 ===
        selected_idx = 1;
        selected = find(pop(selected_idx, :) == 1);
        [~, assign_idx] = min(D(selected, :), [], 1);
        figure;
        geoscatter(coord_demand(:,2), coord_demand(:,1), 40, 'b', 'filled'); hold on;
        geoscatter(coord_site(:,2), coord_site(:,1), 40, 'k', '^');
        geoscatter(coord_site(selected,2), coord_site(selected,1), 80, 'r', 'filled');
        for j = 1:length(assign_idx)
            site_id = selected(assign_idx(j));
            geoplot([coord_demand(j,2), coord_site(site_id,2)], [coord_demand(j,1), coord_site(site_id,1)], '--', 'LineWidth', 1, 'Color', [0.6 0.6 0.6]);
        end
        for k = 1:length(selected)
            lon = coord_site(selected(k), 1) + 0.001;
            lat = coord_site(selected(k), 2);
            text(lon, lat, sprintf('%d', selected(k)), 'FontSize', 8, 'Color', 'm');
        end
        title(sprintf('敏感性分析方案地图 δ=%d, μ=%d', delta, mu));
        geobasemap('streets');
        legend({'需求点','候选站点','选定站点'}, 'Location','bestoutside');
        saveas(gcf, fullfile(folder, sprintf('map_detailed_delta%d_mu%d.png', delta, mu))); close;

        % === Excel 导出 ===
        pareto_solutions = pop(pareto_mask, :);
        obj_min = min(pareto_objectives, [], 1);
        obj_max = max(pareto_objectives, [], 1);
        norm_obj = (pareto_objectives - obj_min) ./ (obj_max - obj_min + eps);
        dist_to_ideal = sqrt(sum(norm_obj.^2, 2));
        [~, top_idx] = mink(dist_to_ideal, 3);

        solutionTags = ["均衡方案一"; "均衡方案二"; "均衡方案三"];
        top_solutions = pareto_solutions(top_idx, :);
        top_objectives = pareto_objectives(top_idx, :);

        site_flags = array2table(top_solutions, 'VariableNames', strcat('站点', string(1:num_sites)));
        site_flags.TotalCost_Yuan = top_objectives(:,1);
        site_flags.TotalDistance_ItemKm = top_objectives(:,2);
        site_flags.RiskValue = top_objectives(:,3);
        site_flags.SolutionID = strcat("方案", string(1:height(site_flags)))';
        site_flags.Tag = solutionTags;

        writetable(site_flags, fullfile(folder, 'top3_solutions.xlsx'));
        fprintf('已保存 Excel：%s\n', fullfile(folder, 'top3_solutions.xlsx'));

        % === 汇总添加 ===
        summary_all = [summary_all; {delta, mu, top_objectives(1,1), top_objectives(1,2), top_objectives(1,3)}];
    end
end

summary_table = cell2table(summary_all, 'VariableNames', {'delta', 'mu', 'BestCost', 'BestDist', 'BestRisk'});
writetable(summary_table, fullfile(base_output, 'summary_all.xlsx'));

fprintf('📊 所有敏感性实验汇总已保存 summary_all.xlsx\n');

disp('敏感性分析全部完成');
%% ========== 函数部分 ==========
function pop = generateInitialPopulation(pop_size, num_sites, min_sites, max_sites, delta_i)
    pop = zeros(pop_size, num_sites);
    valid = find(delta_i == 0);  % 排除禁飞区
    for i = 1:pop_size
        k = randi([min_sites, max_sites]);
        idx = randperm(length(valid), min(k, length(valid)));
        pop(i, valid(idx)) = 1;
    end
end



function new_pop = geneticOperators(pop, rank, p_cross, p_mut, num_sites)
    pop_size = size(pop,1);
    new_pop = zeros(size(pop));

    % 选择：锦标赛选择
    for i = 1:pop_size
        c = randperm(pop_size, 2);
        [~, idx] = min([rank(c(1)), rank(c(2))]);
        new_pop(i,:) = pop(c(idx),:);
    end

    % 交叉：单点交叉
    for i = 1:2:pop_size-1
        if rand < p_cross
            cp = randi([1, num_sites-1]);
            temp = new_pop(i,cp+1:end);
            new_pop(i,cp+1:end) = new_pop(i+1,cp+1:end);
            new_pop(i+1,cp+1:end) = temp;
        end
    end

    % 变异：按位翻转
    for i = 1:pop_size
        for j = 1:num_sites
            if rand < p_mut
                new_pop(i,j) = 1 - new_pop(i,j);
            end
        end
    end
end

function rank = paretoRank(objs)
    n = size(objs,1);
    rank = zeros(n,1);
    dominated = cell(n,1);
    dominationCount = zeros(n,1);

    for i = 1:n
        for j = 1:n
            if i ~= j
                if all(objs(i,:) <= objs(j,:)) && any(objs(i,:) < objs(j,:))
                    dominated{i} = [dominated{i}, j];
                    dominationCount(j) = dominationCount(j) + 1;
                end
            end
        end
    end

    currentFront = find(dominationCount == 0);
    currentRank = 1;
    while ~isempty(currentFront)
        rank(currentFront) = currentRank;
        nextFront = [];
        for i = 1:length(currentFront)
            idx = currentFront(i);
            for j = dominated{idx}
                dominationCount(j) = dominationCount(j) - 1;
                if dominationCount(j) == 0
                    nextFront = [nextFront, j];
                end
            end
        end
        currentFront = nextFront;
        currentRank = currentRank + 1;
    end
end

function D = computeDistanceMatrix(coord_site, coord_demand)
    R = 6371;
    lat1 = deg2rad(coord_site(:,2));
    lon1 = deg2rad(coord_site(:,1));
    lat2 = deg2rad(coord_demand(:,2));
    lon2 = deg2rad(coord_demand(:,1));
    D = zeros(size(coord_site,1), size(coord_demand,1));
    for i = 1:size(coord_site,1)
        dlat = lat2 - lat1(i);
        dlon = lon2 - lon1(i);
        a = sin(dlat/2).^2 + cos(lat1(i)) .* cos(lat2) .* sin(dlon/2).^2;
        c = 2 * atan2(sqrt(a), sqrt(1 - a));
        D(i,:) = R * c;
    end
end


function visualizeTopThreeSolutions_balanced(pop, objectives, D, coord_demand, coord_site, min_sites, max_sites)
    rank = paretoRank(objectives);
    pareto_mask = (rank == 1);
    pareto_solutions_all = pop(pareto_mask, :);
    pareto_objectives_all = objectives(pareto_mask, :);

    %  筛选满足站点数 >= min_sites 的帕累托解
    valid_idx = [];
    for i = 1:size(pareto_solutions_all,1)
        num_selected = sum(pareto_solutions_all(i,:) == 1);
        if num_selected >= min_sites
            valid_idx(end+1) = i;
        end
    end

    % 没有合法解报错退出
    if isempty(valid_idx)
        error('没有满足站点数限制的帕累托解！');
    end

    pareto_solutions = pareto_solutions_all(valid_idx, :);
    pareto_objectives = pareto_objectives_all(valid_idx, :);


    % 归一化目标值
    obj_min = min(pareto_objectives, [], 1);
    obj_max = max(pareto_objectives, [], 1);
    norm_obj = (pareto_objectives - obj_min) ./ (obj_max - obj_min + eps);

    % 不用 vecnorm，改为手动计算欧氏距离
    dist_to_ideal = sqrt(sum(norm_obj.^2, 2));
    [~, top_idx] = mink(dist_to_ideal, 3);

    rep_labels = ["均衡方案一", "均衡方案二", "均衡方案三"];
    for idx = 1:3
        sol_idx = top_idx(idx);
        selected = find(pareto_solutions(sol_idx,:));
        [~, assign_idx] = min(D(selected, :), [], 1);

        figure;
        geoscatter(coord_demand(:,2), coord_demand(:,1), 40, 'b', 'filled'); hold on;
        geoscatter(coord_site(:,2), coord_site(:,1), 40, 'k', '^');
        geoscatter(coord_site(selected,2), coord_site(selected,1), 80, 'r', 'filled');

        for j = 1:length(assign_idx)
            site_id = selected(assign_idx(j));
            geoplot([coord_demand(j,2), coord_site(site_id,2)], ...
                    [coord_demand(j,1), coord_site(site_id,1)], '--', 'LineWidth', 1, 'Color', [0.6 0.6 0.6]);
        end

        for k = 1:length(selected)
            lon = coord_site(selected(k), 1) + 0.001;
            lat = coord_site(selected(k), 2);
            text(lon, lat, sprintf('%d', selected(k)), 'FontSize', 8, 'Color', 'm');
        end

        total_cost = pareto_objectives(sol_idx, 1);
        total_dist = pareto_objectives(sol_idx, 2);
        risk_val = pareto_objectives(sol_idx, 3);
        title_str = sprintf('%s：成本=%.2f万元 距离=%.2f件·公里 风险=%.4f', rep_labels(idx), total_cost, total_dist, risk_val);
        title(title_str);
        geobasemap('streets');
        legend({'需求点','候选站点','选定站点'}, 'Location','bestoutside');

        saveas(gcf, sprintf('Balanced_Solution_%d.png', idx));
    end
end

function exportBalancedSolutionsToExcel(pop, objectives, num_sites)
    rank = paretoRank(objectives);
    pareto_mask = (rank == 1);
    pareto_solutions = pop(pareto_mask, :);
    pareto_objectives = objectives(pareto_mask, :);

    obj_min = min(pareto_objectives, [], 1);
    obj_max = max(pareto_objectives, [], 1);
    norm_obj = (pareto_objectives - obj_min) ./ (obj_max - obj_min + eps);
    dist_to_ideal = sqrt(sum(norm_obj.^2, 2));  % 🚫 替代 vecnorm
    [~, top_idx] = mink(dist_to_ideal, 3);

    solutionTags = ["均衡方案一"; "均衡方案二"; "均衡方案三"];
    top_solutions = pareto_solutions(top_idx, :);
    top_objectives = pareto_objectives(top_idx, :);

    site_flags = array2table(top_solutions, 'VariableNames', strcat('站点', string(1:num_sites)));
    site_flags.TotalCost_Yuan = top_objectives(:,1);
    site_flags.TotalDistance_ItemKm = top_objectives(:,2);
    site_flags.RiskValue = top_objectives(:,3);
    site_flags.SolutionID = strcat("方案", string(1:height(site_flags)))';
    site_flags.Tag = solutionTags;

    writetable(site_flags, '均衡最优选址方案.xlsx');
    disp("三个均衡帕累托解已导出到 Excel（均衡最优选址方案.xlsx）");
end

function visualizeParetoFront(objectives)
    % 计算每个解的帕累托前沿排名
    rank = paretoRank(objectives);
    pareto_mask = (rank == 1);  % 提取帕累托前沿解
    pareto_objectives = objectives(pareto_mask, :);

    % 创建一个 3D 散点图，显示所有解和帕累托最优解
    figure;
    scatter3(objectives(:,1), objectives(:,2), objectives(:,3), 30, 'b', 'filled'); hold on;
    scatter3(pareto_objectives(:,1), pareto_objectives(:,2), pareto_objectives(:,3), 60, 'r', 'filled');
    xlabel('总成本（万元）');
    ylabel('运输量距离（件·公里）');
    zlabel('综合风险');
    title('帕累托前沿图');
    grid on;
    legend('所有解', '帕累托最优解');

    % Save the plot as an image
    saveas(gcf, 'Pareto_Front_Plot.png');
%     close;  % Close the figure after saving
end


function [obj, valid] = evaluateSolution(y, D, F_i, s_i, h_j, delta_i, gamma_i, w_air, w_weather, min_sites, max_sites, max_radius_km)
    obj = [Inf, Inf, Inf];
    valid = false;
    selected = find(y);
    
    % 检查站点数是否超过最大限制
    if numel(selected) < min_sites || numel(selected) > max_sites
        return;
    end

    fix_cost = sum(F_i(selected));
    [min_d, assign_idx] = min(D(selected, :), [], 1);

    % 检查最大服务半径约束
    if any(min_d .* h_j' > max_radius_km * max(h_j))  % 模拟约束 x_ij * d_ij <= D_max
        return;
    end
    base_dist = sum(min_d .* h_j');  % 修改为运输距离总量（件数×距离）

    assign_count = zeros(1, numel(y));
    for j = 1:length(h_j)
        site_id = selected(assign_idx(j));
        assign_count(site_id) = assign_count(site_id) + h_j(j);
    end

    overload_val = sum(max(assign_count(selected) - s_i(selected)', 0));
    risk_penalty = sum(w_air(selected) .* delta_i(selected));
    uncovered = (numel(unique(assign_idx)) ~= length(h_j));

    % ===== 单位一致化惩罚 =====
    penalty_cost = overload_val * 0.01;             % 万元
    penalty_dist = 1000 * double(uncovered);        % 公里

    total_cost = fix_cost + penalty_cost;           % 万元
    total_dist = base_dist + penalty_dist;          % 公里
    risk_val = mean(w_weather(selected) .* gamma_i(selected)) + mean(w_air(selected) .* delta_i(selected));

    obj = [total_cost, total_dist, risk_val];
    valid = true;
end
