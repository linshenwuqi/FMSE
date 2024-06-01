clear; clc

load('rngState.mat', 'rngState'); % Load the saved random seed.
rng(rngState); % Restore the random seed state.

datasetList = {'Seeds', 'Ecoli'};
methodList = {'FMSE_a','FMSE_v'};

methodFunctions = {
    'Method1', @FMSE_a;
    'Method2', @FMSE_v;
};

M = 20;
cntTimes = 10; % How many times will be run.
resultACC = NaN(length(datasetList), length(methodList));
resultNMI = NaN(length(datasetList), length(methodList));
resultARI = NaN(length(datasetList), length(methodList));
resultStdACC = NaN(length(datasetList), length(methodList));
resultStdNMI = NaN(length(datasetList), length(methodList));
resultStdARI = NaN(length(datasetList), length(methodList));

for methodIdx = 1:length(methodList)
    methodName = methodList{methodIdx};
    methodFunction = methodFunctions{methodIdx, 2};
    
    for datasetIdx = 1:length(datasetList)
        dataName = datasetList{datasetIdx};
        load(['data_', dataName, '.mat']);
        disp('*************************************************************');
        disp(['DataName:',dataName]);

        [N, ~] = size(fea);
        clsNums = length(unique(gt));
        [poolSize, ~] = size(FCE_Cell);
        bcIdx = zeros(cntTimes, M);
        result_all = zeros(cntTimes, 4);
        for i = 1:cntTimes
            tmp = randperm(poolSize);
            bcIdx(i,:) = tmp(1:M);
        end

        for runIdx = 1:cntTimes % Average results over cntTimes runs
            disp('*************************************************************');
            disp(['methodName:',methodName, ' -->  ','dataSet:' dataName,' -->  ',...
                ' run:', num2str(runIdx)]);

            FCE_CSC = [];
            clsArr = zeros(1,M);
            for i = 1:M
                FCE = FCE_Cell{bcIdx(runIdx,i)};
                clsArr(i) = size(FCE,2);
                FCE_CSC = [FCE_CSC, FCE];
            end

            % Execute current method function
            tic;
            results = methodFunction(FCE_CSC, clsArr, clsNums);
            elapsedTime = toc;
            disp([methodName, ': ',num2str(elapsedTime),' sec.']);

            % Compute clustering metrics
            res = ClusteringMeasure(gt, results);
            ACC = res(1);
            NMI = res(2);
            ARI = compute_ARI(gt, results);

            % Store results and runtime
            result_all(runIdx, 1) = ACC;
            result_all(runIdx, 2) = NMI;
            result_all(runIdx, 3) = ARI;
            result_all(runIdx, 4) = elapsedTime;

        end
        
        % Compute average results
        resultACC(datasetIdx, methodIdx) = mean(result_all(:, 1));
        resultNMI(datasetIdx, methodIdx) = mean(result_all(:, 2));
        resultARI(datasetIdx, methodIdx) = mean(result_all(:, 3));
        resultStdACC(datasetIdx, methodIdx) = std(result_all(:, 1));
        resultStdNMI(datasetIdx, methodIdx) = std(result_all(:, 2));
        resultStdARI(datasetIdx, methodIdx) = std(result_all(:, 3));
    end
end

save('results\resultACC.mat', 'resultACC');
save('results\resultNMI.mat', 'resultNMI');
save('results\resultARI.mat', 'resultARI');
