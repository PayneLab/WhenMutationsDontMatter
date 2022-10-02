genes=(PTEN PIK3CA KRAS EGFR)
for gene in ${genes[@]}; do 
	echo "Replacing $gene"
	rm -r $gene/*Effect_output;
	rm $gene/full_analysis.ipynb;
	cp TP53/full_analysis.ipynb $gene/full_analysis.ipynb;
done;
