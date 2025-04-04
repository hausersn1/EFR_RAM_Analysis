function [f, P1_env, P1_tfs, PLV_env, PLV_tfs, T_tfs, T_env] = getSpectAverage(pos,neg, fs, subset, k_iters)
    
    trials_p = size(pos,2);
    trials_n = size(neg,2);

    
    P1_env = zeros(floor(length(pos)/2)+1,k_iters);
    P1_tfs = P1_env;

    PLV_env = P1_env;
    PLV_tfs = P1_env;
    
    for i=1:k_iters

        % Get random subset
        pos_sub = pos(:,randperm(trials_p,subset));
        neg_sub = neg(:,randperm(trials_n,subset));

        % Frequency domain calculations
        [f,P1_env(:,i),P1_tfs(:,i),PLV_env(:,i), PLV_tfs(:,i)] = helper.getSpectMag(pos_sub,neg_sub,fs);

        % Time domain calculations
       tfs_t(:,i) = mean((pos_sub - neg_sub)/2,2); 
       env_t(:,i) = mean((pos_sub + neg_sub)/2,2);

    end

PLV_env = mean(PLV_env,2);
PLV_tfs = mean(PLV_tfs,2);

P1_env = mean(P1_env,2);
P1_tfs = mean(P1_tfs,2);

T_tfs = mean(tfs_t,2);
T_env = mean(env_t,2);

end

