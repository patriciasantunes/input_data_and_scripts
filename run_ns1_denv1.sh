cd /home/patricia/Apta-MCTS/src
python apta_mcts.py \
-i '/home/patricia/aptamers_design/apta-mcts_AIchemy_RNA2_zdock/ns1_denv1_protein/apta-mcts/consensus_ns1_denv1.fasta' \
-k 10 \
-bp 70 \
-n 1000 \
-s 'score_functions/rf-ictf-li2014/mcc0.484-ppv1.000-acc0.822-sn0.290-sp1.000-npv0.809-yd0.290-77trees' \
-e '/home/patricia/aptamers_design/apta-mcts_AIchemy_RNA2_zdock/ns1_denv1_protein/apta-mcts/except.fasta' \
-o '/home/patricia/aptamers_design/apta-mcts_AIchemy_RNA2_zdock/ns1_denv1_protein/apta-mcts' 
