
saveDir = 'testEcnLgg';
if ~exist(saveDir, 'dir')
    mkdir(saveDir)
end

ubMax = 20;
nSamples = 100;

% load data
ec = load('ECN_20200729.mat');
ec = ec.model;
lg = load('LrGG_20200729.mat');
lg = lg.model;
data = load('metabolomicData.mat');
[mappedMets, metsFoldChangeEc, metsFoldChangeLg, metsFoldChangeHeader] = ...
    deal(data.mappedMets, data.metsFoldChangeEc, data.metsFoldChangeLg, data.metsFoldChangeHeader);
clear data

% other media components not in the exometabolome data that the microbes need for growth
otherMediaComponentsEc = {'EX_h_e0';'EX_fe3_e0';'EX_mn2_e0';'EX_co2_e0';'EX_zn2_e0';'EX_mg2_e0';'EX_h2o_e0';'EX_ca2_e0';'EX_k_e0';'EX_ni2_e0';'EX_cu2_e0';'EX_na1_e0';'EX_fe2_e0';'EX_cobalt2_e0';'EX_sel_e0';'EX_mobd_e0';'EX_cl_e0'};
otherMediaComponentsLg = {'EX_h2o(e)';'EX_co2(e)';'EX_mn2(e)';'EX_zn2(e)';'EX_cu2(e)';'EX_ca2(e)';'EX_h(e)';'EX_cl(e)';'EX_cobalt2(e)';'EX_k(e)';'EX_mg2(e)';'EX_fol(e)';'EX_na1(e)';'EX_fe3(e)';'EX_fe2(e)'};

exEc = sum(ec.S ~= 0, 1) == 1;
ec.lb(exEc) = 0;
ec = changeRxnBounds(ec, 'EX_o2_e0',  -0.1, 'l');
ec = changeRxnBounds(ec, metsFoldChangeEc(:, 2), -1000, 'l');
ec = changeRxnBounds(ec, 'ATPM',  1, 'l');
ec = changeRxnBounds(ec, otherMediaComponentsEc,  -1000, 'l');

exLg = sum(lg.S ~= 0, 1) == 1;
lg.lb(exLg) = 0;
lg = changeRxnBounds(lg, 'EX_o2(e)',  -0.1, 'l');
lg = changeRxnBounds(lg, metsFoldChangeLg(:, 2), -1000, 'l');
lg = changeRxnBounds(lg, 'ATPM',  1, 'l');
lg = changeRxnBounds(lg, otherMediaComponentsLg,  -1000, 'l');

%% construct the optimization model
[mE, nE] = size(ec.S);
[mL, nL] = size(lg.S);
mExMapped = numel(mappedMets);

LP = struct();
% variables: [ECN fluxes, LGG fluxes, concentrations of media components]
info = struct();
info.colNames = [ec.rxns; lg.rxns; strcat(mappedMets, '_media')];
info.colID.ecRxns = 1:nE;
info.colID.lgRxns = (nE + 1):(nE + nL);
info.colID.media = (nE + nL + 1):(nE + nL + mExMapped);
nCol = nE + nL + mExMapped;

info.rowNames = strcat([ec.mets; lg.mets], '_massBalance');
info.rowID.ecMets = 1:mE;
info.rowID.lgMets = (mE + 1):(mE + mL);
nRow = mE + mL;

LP.A = [ec.S, sparse(mE, nL + mExMapped); ...  % Sv = 0 for ECN
    sparse(mL, nE), lg.S, sparse(mL, mExMapped)];  % Sv = 0 for LGG
LP.csense = repmat('E', nRow, 1);

eps0 = 0.1;
eps = eps0;
medMin = 1e-3;

% add relative uptake/production constraints based on exo-metabolomes

for j = 1:size(metsFoldChangeEc, 1)
    idEx = findRxnIDs(ec, metsFoldChangeEc{j, 2});
    idMappedMets = find(strcmp(mappedMets, metsFoldChangeEc{j, 1}));
    fc = metsFoldChangeEc{j, 3};  % fold change
    % a few inconsistencies that need a large error range for feasible solutions
    if any(strcmp({'spermidine'; 'pantothenate'; 'thiamin (Vitamin B1)'; 'pyridoxine (Vitamin B6)'}, mappedMets{idMappedMets}))
        eps = 1;
    end
    if fc > 1 && strcmp(metsFoldChangeEc{j, 4}, 'significant')
        % UB: v_ex(i) - (FC - 1)(1 + eps)c_i <= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.ecRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, -(fc - 1) * (1 + eps)], 1, nCol);
        info.rowNames(end + 1, :) = {['ECN_produce_' metsFoldChangeEc{j, 1}, '_UB_significant']};
        % LB: v_ex(i) - (FC - 1)(1 - eps)c_i >= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.ecRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, -(fc - 1) * (1 - eps)], 1, nCol);
        info.rowNames(end + 1, :) = {['ECN_produce_' metsFoldChangeEc{j, 1}, '_LB_significant']};
        LP.csense = [LP.csense; 'L'; 'G'];
    elseif fc < 1 && strcmp(metsFoldChangeEc{j, 4}, 'significant')
        % UB: v_ex(i) - (FC - 1)(1 - eps)c_i <= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.ecRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, -(fc - 1) * (1 - eps)], 1, nCol);
        info.rowNames(end + 1, :) = {['ECN_consume_' metsFoldChangeEc{j, 1}, '_UB_significant']};
        % LB: v_ex(i) - (FC - 1)(1 + eps)c_i >= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.ecRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, -(fc - 1) * (1 + eps)], 1, nCol);
        info.rowNames(end + 1, :) = {['ECN_consume_' metsFoldChangeEc{j, 1}, '_LB_significant']};
        LP.csense = [LP.csense; 'L'; 'G'];
    elseif strcmp(metsFoldChangeEc{j, 4}, 'insignificant')
        % UB: v_ex(i) - |FC - 1|c_i <= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.ecRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, -abs(fc - 1)], 1, nCol);
        info.rowNames(end + 1, :) = {['ECN_produce_' metsFoldChangeEc{j, 1}, '_UB_insignificant']};
        % LB: v_ex(i) + |FC - 1|c_i >= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.ecRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, abs(fc - 1)], 1, nCol);
        info.rowNames(end + 1, :) = {['ECN_consume_' metsFoldChangeEc{j, 1}, '_LB_insignificant']};
        LP.csense = [LP.csense; 'L'; 'G'];
    else
        error('Bad category!')
    end
    eps = eps0;
end

info.rowID.ecExchange = (nRow + 1):(nRow + size(metsFoldChangeEc, 1) * 2);
nRow = nRow + size(metsFoldChangeEc, 1) * 2;

for j = 1:size(metsFoldChangeLg, 1)
    idEx = findRxnIDs(lg, metsFoldChangeLg{j, 2});
    idMappedMets = find(strcmp(mappedMets, metsFoldChangeLg{j, 1}));
    fc = metsFoldChangeLg{j, 3};  % fold change
    if any(strcmp({'spermidine'; 'pantothenate'; 'thiamin (Vitamin B1)'; 'pyridoxine (Vitamin B6)'}, mappedMets{idMappedMets}))
        eps = 1;
    end
    if fc > 1 && strcmp(metsFoldChangeLg{j, 4}, 'significant')
        % UB: v_ex(i) - (FC - 1)(1 + eps)c_i <= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.lgRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, -(fc - 1) * (1 + eps)], 1, nCol);
        info.rowNames(end + 1, :) = {['LGG_produce_' metsFoldChangeLg{j, 1}, '_UB_significant']};
        % LB: v_ex(i) - (FC - 1)(1 - eps)c_i >= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.lgRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, -(fc - 1) * (1 - eps)], 1, nCol);
        info.rowNames(end + 1, :) = {['LGG_produce_' metsFoldChangeLg{j, 1}, '_LB_significant']};
        LP.csense = [LP.csense; 'L'; 'G'];
    elseif fc < 1 && strcmp(metsFoldChangeLg{j, 4}, 'significant')
        % UB: v_ex(i) - (FC - 1)(1 - eps)c_i <= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.lgRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, -(fc - 1) * (1 - eps)], 1, nCol);
        info.rowNames(end + 1, :) = {['LGG_consume_' metsFoldChangeLg{j, 1}, '_UB_significant']};
        % LB: v_ex(i) - (FC - 1)(1 + eps)c_i >= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.lgRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, -(fc - 1) * (1 + eps)], 1, nCol);
        info.rowNames(end + 1, :) = {['LGG_consume_' metsFoldChangeLg{j, 1}, '_LB_significant']};
        LP.csense = [LP.csense; 'L'; 'G'];
    elseif strcmp(metsFoldChangeLg{j, 4}, 'insignificant')
        % UB: v_ex(i) - |FC - 1|c_i <= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.lgRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, -abs(fc - 1)], 1, nCol);
        info.rowNames(end + 1, :) = {['LGG_produce_' metsFoldChangeLg{j, 1}, '_UB_insignificant']};
        % LB: v_ex(i) + |FC - 1|c_i >= 0
        LP.A(end + 1, :) = sparse([1; 1], [info.colID.lgRxns(idEx); info.colID.media(idMappedMets)], ...
            [1, abs(fc - 1)], 1, nCol);
        info.rowNames(end + 1, :) = {['LGG_consume_' metsFoldChangeLg{j, 1}, '_LB_insignificant']};
        LP.csense = [LP.csense; 'L'; 'G'];
    else
        error('Bad category!')
    end
    eps = eps0;
end

info.rowID.lgExchange = (nRow + 1):(nRow + size(metsFoldChangeLg, 1) * 2);
nRow = nRow + size(metsFoldChangeLg, 1) * 2;

% constrain acetate production less than or equal to lactate production in
% LGG. Constraint-based models have been reported to be unable to predict
% lower-yield higher-rate metabolism such as lactate production in lactic acid bacteria 
% (compared to acetate production) since constraint-based models maximize yield. 
% Add an empirical constraint to simulate a more likely metabolism for LGG
LP.A(end + 1, info.colID.lgRxns(findRxnIDs(lg, {'EX_lac__L(e)'; 'EX_ac(e)'}))) = [-1, 1];
nRow = nRow + 1;
LP.csense(end + 1) = 'L';
info.rowNames(end + 1) = {'LactateToAcetateMinRatio'};

% equal growth constraints. From the experimental data, the final OD of the
% two probiotics are approximately the same.
LP.A(end + 1, [info.colID.ecRxns(findRxnIDs(ec, 'bio1')); info.colID.lgRxns(findRxnIDs(lg, 'bio1'))]) = [1, -1];
nRow = nRow + 1;
LP.csense(end + 1) = 'E';

info.rowNames(end + 1) = {'EqualGrowth'};
info.rowID.others = (nRow + 1):(nRow + 2);

% Add absolute variables |v|
LP.A = [LP.A, sparse(nRow, nE + nL); ...
    spdiag(1:(nE + nL)), sparse(nE + nL, mExMapped), -spdiag(1:(nE + nL)); ... v - |v| <= 0
    -spdiag(1:(nE + nL)), sparse(nE + nL, mExMapped), -spdiag(1:(nE + nL))]; ... -v - |v| <= 0
LP.csense((end + 1):(end + (nE + nL)*2)) = 'L';
info.rowNames((end + 1):(end + nE + nL)) = strcat([ec.rxns; lg.rxns], '_abs1');
info.rowNames((end + 1):(end + nE + nL)) = strcat([ec.rxns; lg.rxns], '_abs2');
info.rowID.rxnAbs1 = (nRow + 1):(nRow + nE + nL);
nRow = nRow + nE + nL;
info.rowID.rxnAbs2 = (nRow + 1):(nRow + nE + nL);
nRow = nRow + nE + nL;

info.colNames((end + 1):(end + nE + nL)) = strcat([ec.rxns; lg.rxns], '_abs');
info.colID.rxnAbs = (nCol + 1):(nCol + nE + nL);
nCol = nCol + nE + nL;

rxnBiomass = [info.colID.ecRxns(findRxnIDs(ec, 'bio1')); info.colID.lgRxns(findRxnIDs(lg, 'bio1'))];
LP.b = zeros(nRow, 1);
LP.c = zeros(nCol, 1);
LP.lb = [ec.lb; lg.lb; 0 * ones(mExMapped, 1); zeros(nE + nL, 1)];
LP.ub = [ec.ub; lg.ub; 1000 * ones(mExMapped, 1); 1000 * ones(nE + nL, 1)];

%% sample uptake bounds

LP.osense = 1;

% determine the minimum values for the media variables that sustain minimum growth
LP.lb(info.colID.media) = medMin;
LP.lb(rxnBiomass) = 0.1;
LP.c(info.colID.media) = 1;
sol = solveCobraLP(LP);
medUbMin = sol.full(info.colID.media);
medRandomUB = ((ubMax - medUbMin) * ones(1, nSamples)) .* rand(mExMapped, nSamples) + medUbMin * ones(1, nSamples);
    
%% pFBA simulation with randomly sampled uptake bounds

for j = 1:nSamples
    
    LP.lb(info.colID.media) = medMin;
    LP.ub(info.colID.media) = round(medRandomUB(:, j), 6);
    LP.lb(rxnBiomass) = 0;
    LP.ub(rxnBiomass) = 1000;
    
    % maximize biomass
    LP.c(:) = 0;
    LP.c(rxnBiomass) = 1;
    LP.osense = -1;
    
    sol = solveCobraLP(LP);
    
    % take the shadow price with maximum biomass as objective function
    shadowPriceEc = sol.dual(info.rowID.ecExchange);
    shadowPriceLg = sol.dual(info.rowID.lgExchange);
    
    % fix biomass, minimize total sum of absolute fluxes
    LP.lb(rxnBiomass) = round(sol.full(rxnBiomass), 6) - 1e-6;
    LP.ub(rxnBiomass) = round(sol.full(rxnBiomass), 6) + 1e-6;
    
    LP.c(:) = 0;
    LP.c(info.colID.rxnAbs) = 1;
    LP.osense = 1;
    
    sol2 = solveCobraLP(LP);
    
    % take the pFBA flux distributions 
    fluxEc = sol2.full(info.colID.ecRxns);
    fluxLg = sol2.full(info.colID.lgRxns);
    
    save(sprintf('%s%ssample%05d.mat', saveDir, filesep, j), 'fluxEc', 'fluxLg', 'shadowPriceEc', 'shadowPriceLg', 'sol', 'sol2'); 
    fprintf('Finished %d. %04d-%02d-%02d %02d:%02d:%02.0f\n', j, clock)
end

%% analyze shadow prices

metExShadowPrice = NaN(mExMapped, nSamples, 2);
growthPromotingProducts = false(mExMapped, nSamples, 2);
growthCompetingProducts = false(mExMapped, nSamples, 2);
growthPromotingSubstrates = false(mExMapped, nSamples, 2);
growthCompetingSubstrates = false(mExMapped, nSamples, 2);
[yn, idEcExToMetEx] = ismember(metsFoldChangeEc(:, 1), mappedMets);
assert(all(yn))
[yn, idLgExToMetEx] = ismember(metsFoldChangeLg(:, 1), mappedMets);
assert(all(yn))
for j = 1:nSamples
    data = load(sprintf('%s%ssample%05d.mat', saveDir, filesep, j), 'fluxEc', 'fluxLg', 'shadowPriceEc', 'shadowPriceLg', 'sol', 'sol2'); 
    shadowPriceEc = data.shadowPriceEc(1:2:end) + data.shadowPriceEc(2:2:end);
    for k = 1:size(metsFoldChangeEc, 1)
        if metsFoldChangeEc{k, 3} > 1 && strcmp(metsFoldChangeEc{k, 4}, 'significant')
            % significantly produced
            % positive shadow price (in the lower bound constraint)
            % means growth-competing product
            % negative shadow price (in the upper bound constraint) means 
            % growth-promoting product
            metExShadowPrice(idEcExToMetEx(k), j, 1) = shadowPriceEc(k);
            growthCompetingProducts(idEcExToMetEx(k), j, 1) = shadowPriceEc(k) > 0;
            growthPromotingProducts(idEcExToMetEx(k), j, 1) = shadowPriceEc(k) < 0;
        elseif metsFoldChangeEc{k, 3} < 1 && strcmp(metsFoldChangeEc{k, 4}, 'significant')
            % significantly consumed
            % positive shadow price (in the lower bound constraint)
            % means growth-promoting substrate
            % negative shadow price (in the upper bound constraint) means 
            % growth-competing substrate
            metExShadowPrice(idEcExToMetEx(k), j, 1) = shadowPriceEc(k);
            growthCompetingSubstrates(idEcExToMetEx(k), j, 1) = shadowPriceEc(k) < 0;
            growthPromotingSubstrates(idEcExToMetEx(k), j, 1) = shadowPriceEc(k) > 0;
        end
    end
    
    shadowPriceLg = data.shadowPriceLg(1:2:end) + data.shadowPriceLg(2:2:end);
    for k = 1:size(metsFoldChangeLg, 1)
        if metsFoldChangeLg{k, 3} > 1 && strcmp(metsFoldChangeLg{k, 4}, 'significant')
            % significantly produced
            % positive shadow price (in the lower bound constraint)
            % means growth-competing product
            % negative shadow price (in the upper bound constraint) means 
            % growth-promoting product
            metExShadowPrice(idLgExToMetEx(k), j, 2) = shadowPriceLg(k);
            growthCompetingProducts(idLgExToMetEx(k), j, 2) = shadowPriceLg(k) > 0;
            growthPromotingProducts(idLgExToMetEx(k), j, 2) = shadowPriceLg(k) < 0;
        elseif metsFoldChangeLg{k, 3} < 1 && strcmp(metsFoldChangeLg{k, 4}, 'significant')
            % significantly consumed
            % positive shadow price (in the lower bound constraint)
            % means growth-promoting substrate
            % negative shadow price (in the upper bound constraint) means 
            % growth-competing substrate
            metExShadowPrice(idLgExToMetEx(k), j, 2) = shadowPriceLg(k);
            growthCompetingSubstrates(idLgExToMetEx(k), j, 2) = shadowPriceLg(k) < 0;
            growthPromotingSubstrates(idLgExToMetEx(k), j, 2) = shadowPriceLg(k) > 0;
        end
    end
end

% probability of metabolites being growth-promoting/-competing substrates/products
growthPromotingProductsProb = permute(sum(growthPromotingProducts, 2), [1 3 2]) / nSamples;
growthCompetingProductsProb = permute(sum(growthCompetingProducts, 2), [1 3 2]) / nSamples;
growthPromotingSubstratesProb = permute(sum(growthPromotingSubstrates, 2), [1 3 2]) / nSamples;
growthCompetingSubstratesProb = permute(sum(growthCompetingSubstrates, 2), [1 3 2]) / nSamples;

% mean and s.d. of shadow price
metExShadowPriceMean = permute(mean(metExShadowPrice, 2), [1 3 2]);
metExShadowPriceSD = permute(std(metExShadowPrice, 0, 2), [1 3 2]);