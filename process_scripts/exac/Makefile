all: format_ref rescore_ref clean_data create_fig3 create_fig4

format_ref: ../../ref/exac/exac_ref_formatted_converted_flipped.txt

rescore_ref: ../../ref/exac/exac_ref_rescored.txt

clean_data: ../../processed_data/exac/exac_data_clean.txt

run_HAL: ../../processed_data/exac/exac_HAL_scores.txt

run_SPANR: ../../processed_data/exac/exac_SPANR_scores.txt

create_fig3: ../../figs/exac/exac_snv_density.png

create_fig4: ../../figs/exac/exac_fig3D_pLI_enrichment.png


../../ref/exac/exac_ref_formatted_converted_flipped.txt: ../../ref/exac/exac_ref_all.txt
	@ echo "Formatting reference..."
	@ Rscript exac_format_ref.R


../../ref/exac/exac_ref_rescored.txt: ../../ref/exac/exac_ref_formatted_converted_flipped.txt
	@ echo "Re-scoring reference file..."
	@ cd ../splicemod_src && \
	python ../score_lib.py $< --rev ./data/motifs/Ke2011/ \
	../../ref/Rosenberg_2015/exonic_mer6_scores.series 30 $@


../../processed_data/exac/exac_data_clean.txt: ../../ref/exac/exac_ref_formatted_converted_flipped.txt exac_data_clean.R
	@ echo "Cleaning up data..."
	@ Rscript exac_data_clean.R


../../processed_data/exac/exac_HAL_scores.txt: ../../processed_data/exac/exac_data_clean.txt exac_run_HAL.py
	@ echo "Calculating HAL scores..."
	@ python exac_run_HAL.py $< $@


../../processed_data/exac/exac_SPANR_scores.txt: ../../processed_data/exac/exac_data_clean.txt exac_run_SPANR.py
	@ echo "Calculating SPANR scores..."
	@ python exac_run_SPANR.py $< $@


../../figs/exac/exac_snv_density.png: ../../processed_data/exac/exac_data_clean.txt exac_graph_fig3_lib_stats.R
	@ echo "Creating graphs for Figure 3..."
	@ Rscript exac_graph_fig3_lib_stats.R


../../figs/exac/exac_fig3D_pLI_enrichment.png: ../../processed_data/exac/exac_func_annot.txt exac_graph_fig4_func_annot.R
	@ echo "Creating graphs for Figure 4..."
	@ Rscript exac_graph_fig4_func_annot.R
