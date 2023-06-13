cuffnorm -p 8 -o cuffnorm --library-type fr-firststrand --labels Col0,cyp71,se $gtf $Col0_rep1,$Col0_rep2,$Col0_rep3 $cyp71_rep1,$cyp71_rep2,$cyp71_rep3 $se_rep1,$se_rep2,$se_rep3
cuffdiff -p 8 -o cuffnorm --library-type fr-firststrand --labels Col0,cyp71,se $gtf $Col0_rep1,$Col0_rep2,$Col0_rep3 $cyp71_rep1,$cyp71_rep2,$cyp71_rep3 $se_rep1,$se_rep2,$se_rep3
