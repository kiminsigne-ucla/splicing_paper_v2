# TODO: add commands to process sequencing data. This will run everything for
# splicemod, from processing sequencing, formatting reference, and making figures

all: format_ref rescore_ref clean_data qc_figs

format_ref: ../../ref/splicemod/splicemod_ref_formatted_converted.txt

rescore_ref: ../../ref/splicemod/splicemod_ref_rescored.txt

clean_data: ../../processed_data/splicemod/splicemod_data_clean.txt

qc_figs: ../../figs/splicemod/smn1/splicemod_smn1_replicates.png

category_figs: ../../figs/splicemod/nat_index_legend.png

figs: qc_figs category_figs


../../ref/splicemod/splicemod_ref_formatted_converted.txt: ../../ref/splicemod/splicemod_ref.txt
	@ echo "Formatting reference..."
	@ Rscript splicemod_format_ref.R


../../ref/splicemod/splicemod_ref_rescored.txt: ../../ref/splicemod/splicemod_ref_formatted_converted.txt
	@ echo "Re-scoring reference file..."
	@ cd ../splicemod_src && \
	python ../score_lib.py $< ./data/motifs/Ke2011/ \
	../../ref/Rosenberg_2015/exonic_mer6_scores.series 30 $@


../../processed_data/splicemod/splicemod_data_clean.txt: ../../ref/splicemod/splicemod_ref_formatted_converted.txt
	@ echo "Cleaning up data..."
	@ Rscript splicemod_data_clean.R


# these scripts make multiple figures, list the first figure as a target so the script
# will run anytime it is older than script that makes it. This will make all the
# figures which is fine since it's quick
../../figs/splicemod/smn1/splicemod_smn1_replicates.png: splicemod_graph_qc.R ../../processed_data/splicemod/splicemod_data_clean.txt
	@ echo "Making quality control figures..."
	@ Rscript splicemod_graph_qc.R


../../figs/splicemod/nat_index_legend.png: splicemod_graph_categories.R ../../processed_data/splicemod/splicemod_data_clean.txt
	@ echo "Making figures for splicemod categories..."
	@ Rscript splicemod_graph_categories.R




