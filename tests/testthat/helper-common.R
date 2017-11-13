# Shorthand: map 1->1, 2->2, 3->3, 4->1, ..., 149->2, 150->3, ... (if base==3)
I = function(i, base)
	(i-1) %% base + 1
