load tennis_data

M = size(W,1);            % number of players
N = size(G,1);            % number of games in 2011 season 

psi = inline('normpdf(x)./normcdf(x)');
lambda = inline('(normpdf(x)./normcdf(x)).*( (normpdf(x)./normcdf(x)) + x)');

pv = 0.5;            % prior skill variance (prior mean is always 0)

% initialize matrices of skill marginals - means and precisions
Ms = nan(M,1); 
Ps = nan(M,1);

% initialize matrices of game to skill messages - means and precisions
Mgs = zeros(N,2); 
Pgs = zeros(N,2);

% allocate matrices of skill to game messages - means and precisions
Msg = nan(N,2); 
Psg = nan(N,2);

ind = randi(M,3,1);


%% Parameters for messages and n
ms = [];
ps = [];
iter = 0 ;


for iter= 1:3
    
    
    
    
  % (1) compute marginal skills 
  for p=1:M
    % precision first because it is needed for the mean update
    Ps(p) = 1/pv + sum(Pgs(G==p)); 
    Ms(p) = sum(Pgs(G==p).*Mgs(G==p))./Ps(p);
  end
    ps = [ps; Ps'];
    ms = [ms; Ms'];
  % (2) compute skill to game messages
  % precision first because it is needed for the mean update
  Psg = Ps(G) - Pgs;
  Msg = (Ps(G).*Ms(G) - Pgs.*Mgs)./Psg;
    
  % (3) compute game to performance messages
  vgt = 1 + sum(1./Psg, 2);
  mgt = Msg(:,1) - Msg(:,2); % player 1 always wins the way we store data
   
  % (4) approximate the marginal on performance differences
  Mt = mgt + sqrt(vgt).*psi(mgt./sqrt(vgt));
  Pt = 1./( vgt.*( 1-lambda(mgt./sqrt(vgt)) ) );
    
  % (5) compute performance to game messages
  ptg = Pt - 1./vgt;
  mtg = (Mt.*Pt - mgt./vgt)./ptg;   
    
  % (6) compute game to skills messages
  Pgs = 1./(1 + repmat(1./ptg,1,2) + 1./Psg(:,[2 1]));
  Mgs = [mtg, -mtg] + Msg(:,[2 1]);
end

while (norm(ms(iter) - ms(iter-2))> 0.001)*(norm(ps(iter) - ps(iter-2))> 0.001)
   
  % (1) compute marginal skills 
  for p=1:M
    % precision first because it is needed for the mean update
    Ps(p) = 1/pv + sum(Pgs(G==p)); 
    Ms(p) = sum(Pgs(G==p).*Mgs(G==p))./Ps(p);
  end
    ps = [ps; Ps'];
    ms = [ms; Ms'];
  % (2) compute skill to game messages
  % precision first because it is needed for the mean update
  Psg = Ps(G) - Pgs;
  Msg = (Ps(G).*Ms(G) - Pgs.*Mgs)./Psg;
    
  % (3) compute game to performance messages
  vgt = 1 + sum(1./Psg, 2);
  mgt = Msg(:,1) - Msg(:,2); % player 1 always wins the way we store data
   
  % (4) approximate the marginal on performance differences
  Mt = mgt + sqrt(vgt).*psi(mgt./sqrt(vgt));
  Pt = 1./( vgt.*( 1-lambda(mgt./sqrt(vgt)) ) );
    
  % (5) compute performance to game messages
  ptg = Pt - 1./vgt;
  mtg = (Mt.*Pt - mgt./vgt)./ptg;   
    
  % (6) compute game to skills messages
  Pgs = 1./(1 + repmat(1./ptg,1,2) + 1./Psg(:,[2 1]));
  Mgs = [mtg, -mtg] + Msg(:,[2 1]);
  
  iter = iter + 1;
  
end

figure(1),clf
plot(ms(:, ind))
xlabel('Number of iterations')
ylabel('Marginal Skill Means')
legend('Player A', 'Player B', 'Player C')
xlim([2 iter])
figure(2), clf
plot(ps(:, ind))
xlabel('Number of iterations')
ylabel('Marginal Skill Precisions')
legend('Player A', 'Player B', 'Player C')
xlim([2 iter])
Vars = 1./Ps;

%% Exercise h) 

% P = nan(M,1);

% for player = 1:M
%     win_prob = [];
%     for player2 = 1:M
%         
%             if player ~= player2
%             win_prob = [win_prob; normcdf((Ms(player) - Ms(player2))/(sqrt(1 + Vars(player) + Vars(player2))))];
%             end
%         P(player) = mean(win_prob);
%     end
% end % Then use cw2.m
% figure(3)
% [kk, ii] = sort(P, 'descend');
% 
% np = 107;
% barh(kk(np:-1:1))
% set(gca,'YTickLabel',W(ii(np:-1:1)),'YTick',1:np,'FontSize',8)
% axis([0 1 0.5 np+0.5])

%% Exercise i)
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
            win_prob = [win_prob; normcdf((Ms(player) - Ms(player2))/(sqrt(Vars(player) + Vars(player2))))];
            end
        end
        P(i, j) = mean(win_prob);
        j = j+1
    end
   i = i+1
end % Then use cw2.m




