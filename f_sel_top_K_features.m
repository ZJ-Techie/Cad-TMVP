function idx_selected = f_sel_top_K_features(B, top_k_selected)

% B: p x c

if nargin < 2 
	top_k_selected = 10;
end

b_weight = sum(abs(B), 2);
[b_weigth_sorted, b_weigth_sorted_idx] = sort(b_weight, 'descend');

idx_selected = b_weigth_sorted_idx(1 : max(top_k_selected));

%%%% end of f_select_feauture %%%%