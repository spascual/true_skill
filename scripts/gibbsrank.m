clear
load tennis_data

randn('seed',27); % set the pseudo-random number generator seed

M = size(W,1);            % number of players
N = size(G,1);            % number of games in 2011 season 

pv = 0.5*ones(M,1);           % prior skill variance 
prior_mu = zeros(M,1);     
w = prior_mu;                 % set skills to prior mean

iterations = 1000;

skills = []; %Collects skills of 3 random players
% skillBs = [];
% skillCs = [];
ind = randi(M,3,1); 

%for seed = [12,27,48]
% 
%     skillA = []; 
%     skillB = []; 
%     skillC = []; 
for i = 1:iterations 

  % First, sample performance differences given the skills and outcomes
  
  t = nan(N,1); % contains a t_g variable for each game
  for g = 1:N   % loop over games
    s = w(G(g,1))-w(G(g,2));  % difference in skills
    t(g) = randn()+s;         % performace difference sample
    while t(g) < 0  % rejection sampling: only positive perf diffs accepted
      t(g) = randn()+s; % if rejected, sample again
    end
  end 
 
  
  % Second, jointly sample skills given the performance differences
  
  m = nan(M,1);  % container for the mean of the conditional 
                 % skill distribution given the t_g samples
  for p = 1:M
      m(p) = prior_mu(p)/pv(p) + sum(((G(:,1) == p) - (G(:,2) == p)).* t);
  end
  
  iS = zeros(M,M); % container for the sum of precision matrices contributed
                   % by all the games (likelihood terms)

  for g = 1:N 
      iS(G(g,1),G(g,1)) = iS(G(g,1),G(g,1)) + 1; 
      iS(G(g,2),G(g,2)) = iS(G(g,2),G(g,2)) + 1;
      iS(G(g,1),G(g,2)) = iS(G(g,1),G(g,2)) - 1; 
      iS(G(g,2),G(g,1)) = iS(G(g,1),G(g,2));  
  end

  iSS = diag(1./pv) + iS; % posterior precision matrix
  % prepare to sample from a multivariate Gaussian
  % Note: inv(M)*z = R\(R'\z) where R = chol(M);
  iR = chol(iSS);  % Cholesky decomposition of the posterior precision matrix
  mu = iR\(iR'\m); % equivalent to inv(iSS)*m but more efficient
    
  % sample from N(mu, inv(iSS))
  w = mu + iR\randn(M,1);

%   skillA = [skillA, w(ind(1))]; %Collects skills of 3 random players
%   skillB = [skillB, w(ind(2))];
%   skillC = [skillC, w(ind(3))];

  skills = [skills, w];
  i
end
%% Exercise c) 
%   skillAs = [skillAs; skillA]; %Collects skills of 3 random players
%   skillBs = [skillBs; skillB];
%   skillCs = [skillCs; skillC];
%   
%   skillAs = skillAs';
%   skills = [skillAs,skillBs',skillCs'];
%   
%   figure(1)
%   plot(skills(200:8:1000,:))
%   xlabel('Iterations')
%   ylabel('Skill sample produced')
%   legend('Player A', 'Player B', 'Player C')
%   xlim([0 100])
% end
% figure(1), clf
% plot(skillAs')
% xlabel('Iterations')
% ylabel('Skill sample produced')
% legend('Seed 12', 'Seed 27', 'Seed 48')
% 
% figure(2), clf
% plot(skillBs')
% xlabel('Iterations')
% ylabel('Skill sample produced')
% legend('Seed 12', 'Seed 27', 'Seed 48')
% 
% figure(3), clf
% plot(skillCs')
% xlabel('Iterations')
% ylabel('Skill sample produced')
% legend('Seed 12', 'Seed 27', 'Seed 48')
% 


% figure(2), clf
% histogram(skills(20:100,1), 10)
% ylabel('Frequencies')
% xlabel('Skill sample produced')
% figure(3), clf
% %plot(xcov(skills(:,1))) 
% stem(autocor(skillAs(200:8:1000,1),8))
% xlabel('Lag')
% ylabel('Autocorretalation \theta_{t}')
% ylim([-0.2,0.2])



%% Exercise e)

% for player = 1:M
%     win_prob = zeros(hyp_games,1);
%     RandomRivals = randi(M,hyp_games,1);
%         for g = 1:hyp_games
%             win_prob(g) = normcdf(w(player) - w(RandomRivals(g)));
%         end
%     Pnew(player) = mean(win_prob); 
% end

%% Exercise f)

% P = zeros(4, 4); skill_selection = skills(:,200:8:1000);
% players = [16, 1 ,5, 11];
% i = 1; 
% for player = players
%     j = 1
%     for player2 = players
%         win_prob = [];
%         for sample = 1:101
%             ws = skill_selection(:,sample);
%             if player2 ~= player
%             win_prob = [win_prob; normcdf(ws(player) - ws(player2))];
%             end
%         end
%         P(i, j) = mean(win_prob);
%         j = j+1
%     end
%    i = i+1
% end % Then use cw2.m

%% Exercise j) 
P = zeros(4, 4); skill_selection = skills(:,200:8:1000);
players = [16, 1 ,5, 11];
i = 1; 
for player = players
    j = 1
    for player2 = players
        win_prob = [];
        for sample = 1:101
            ws = skill_selection(:,sample);
            if player2 ~= player
            win_prob = [win_prob; ws(player) > ws(player2)];
            end
        end
        P(i, j) = mean(win_prob);
        j = j+1
    end
   i = i+1
end % Then use cw2.m



