for gene in ${genes[@]}; do 
	rm -r $gene/*Effect_output;
	rm $gene/full_analysis.ipynb;
	cp TP53/full_analysis.ipynb $gene/full_analysis.ipynb;
done;
