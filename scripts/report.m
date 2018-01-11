clear
load tennis_data 
% P = zeros(np, 1); np = 107;
% for i = 1:np
%    Tot_games = sum([G(:,1)== i; G(:,2)== i]) ; 
%    Tot_wins = sum(G(:,1)== i);
%    P(i) = Tot_wins/Tot_games;  
% end
% %Then use cw2 for P
% [kk, ii] = sort(P, 'descend'); 
% barh(kk(np:-1:1))
% 
% set(gca,'YTickLabel',W(ii(np:-1:1)),'YTick',1:np,'FontSize',7.5)
% axis([0 1 0.5 np+0.5])



randn('seed',27); % set the pseudo-random number generator seed

M = 20;
M = size(W,1);            % number of players
N = size(G,1);            % number of games in 2011 season 

pv = 0.5*ones(M,1);           % prior skill variance 

w = zeros(M,1);               % set skills to prior mean

m = nan(M,1);  % container for the mean of the conditional 
              % skill distribution given the t_g samples
iS1 = zeros(M,M);
for l = 1:M
      for k = 1:M
          if l == k 
              iS1(l,k) = sum((G(:,1) == l) + (G(:,2) == l));
          else
              iS1(l,k) = - sum((G(:,1) == l).*(G(:,2) == k) + (G(:,1) == k).*(G(:,2) == l));
          end
      end
end
iS = zeros(M,M);
for g = 1:N 
    iS(G(g,1),G(g,1)) = iS(G(g,1),G(g,1)) + 1; 
    iS(G(g,2),G(g,2)) = iS(G(g,2),G(g,2)) + 1;
    iS(G(g,1),G(g,2)) = iS(G(g,1),G(g,2)) - 1; 
    iS(G(g,2),G(g,1)) = iS(G(g,1),G(g,2));  
end





% 
% for i = 1:100
% 
%   % First, sample performance differences given the skills and outcomes
%   
%   t = nan(N,1); % contains a t_g variable for each game
%   for g = 1:N   % loop over games
%     s = w(G(g,1))-w(G(g,2));  % difference in skills
%     t(g) = randn(1)+s;         % performace difference sample
%     while t(g) < 0  % rejection sampling: only positive perf diffs accepted
%       t(g) = randn(1)+s; % if rejected, sample again
%     end
%   end 
%  
%   
%   % Second, jointly sample skills given the performance differences
%   
%   m = nan(M,1);  % container for the mean of the conditional 
%   f = nan(M,1);               % skill distribution given the t_g samples
%   for p = 1:M
%         
%       m(p) = w(p)/pv(p) + sum(((G(:,1) == p) - (G(:,2) == p)));
%      
%   end
%   
% 
% 
% 
% for p = 1:M
%     win_index = find(G(:,1)==p,1802);
%     loss_index = find(G(:,2)==p,1802);
%     f(p) = sum(t(win_index)) - sum(t(loss_index));
%   end
% end
% %   
% %   
% %   iS = zeros(M,M); % container for the sum of precision matrices contributed
% %                    % by all the games (likelihood terms)
% %   for l = 1:M
% %       for k = 1:M
% %           if l == k 
% %               iS(l,k) = sum((G(:,1) == l) + (G(:,2) == l));
% %           else
% %               iS(l,k) =  - sum((G(:,1) == l).*(G(:,2) == k) + (G(:,1) == k).*(G(:,2) == l));
% %           end
% %       end
% %   end
% %   
% % end