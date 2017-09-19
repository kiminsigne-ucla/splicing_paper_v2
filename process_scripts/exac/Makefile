all: format_ref rescore_ref clean_data add_annot create_fig3 create_fig4 create_curves graph_curves

format_ref: ../../ref/exac/exac_ref_formatted_converted_flipped.txt

rescore_ref: ../../ref/exac/exac_ref_rescored.txt

clean_data: ../../processed_data/exac/exac_data_clean.txt

add_annot: ../../processed_data/exac/exac_func_annot.txt

run_HAL: ../../processed_data/exac/exac_HAL_scores.txt

run_SPANR: ../../processed_data/exac/exac_SPANR_scores.txt

create_fig3: ../../figs/exac/exac_snv_density.png

create_fig4: ../../figs/exac/exac_fig3D_pLI_enrichment.png

create_curves: ../../processed_data/exac/exac_models_roc_curves_intron.txt

graph_curves: ../../figs/exac/exac_fig4E_exac_roc_curves_type_legend.tiff


../../ref/exac/exac_ref_formatted_converted_original_seq.txt: ../../ref/exac/exac_ref_all.txt
	@ echo "Formatting reference..."
	@ Rscript exac_format_ref.R


../../ref/exac/exac_ref_rescored.txt: ../../ref/exac/exac_ref_formatted_converted_original_seq.txt
	@ echo "Re-scoring reference file..."
	@ cd ../splicemod_src && \
	python ../score_lib.py $< original_seq ./data/motifs/Ke2011/ \
	../../ref/Rosenberg_2015/exonic_mer6_scores.series 30 $@


../../processed_data/exac/exac_data_clean.txt: ../../ref/exac/exac_ref_formatted_converted_original_seq.txt exac_data_clean.R
	@ echo "Cleaning up data..."
	@ Rscript exac_data_clean.R


../../processed_data/exac/exac_func_annot.txt: ./../processed_data/exac/exac_data_clean.txt exac_add_annot.R
	@ echo "Adding functional annotation..."
	@ Rscript exac_add_annot.R


../../processed_data/exac/exac_HAL_scores.txt: ../../processed_data/exac/exac_data_clean.txt exac_run_HAL.py
	@ echo "Calculating HAL scores..."
	@ python exac_run_HAL.py $< $@


../../processed_data/exac/exac_SPANR_scores.txt: ../../processed_data/exac/exac_data_clean.txt exac_run_SPANR.py
	@ echo "Calculating SPANR scores..."
	@ python exac_run_SPANR.py $< $@


../../figs/exac/exac_snv_density.svg: ../../processed_data/exac/exac_data_clean.txt exac_graph_fig3_lib_stats.R
	@ echo "Creating graphs for Figure 3..."
	@ Rscript exac_graph_fig3_lib_stats.R


../../figs/exac/exac_fig4D_allele_frequency_binned.tiff: ../../processed_data/exac/exac_func_annot.txt exac_graph_fig4_func_annot.R
	@ echo "Creating graphs for Figure 4..."
	@ Rscript exac_graph_fig4_func_annot.R


../../processed_data/exac/exac_models_roc_curves_intron.txt: ../../processed_data/exac/exac_func_annot.txt \
	../../processed_data/exac/exac_SPANR_scores_capped.txt \
	../../processed_data/exac/exac_HAL_scores.txt
	@ echo "Creating PR and ROC curves..."
	@ python exac_models.py


../../figs/exac/exac_fig4E_exac_roc_curves_type_legend.tiff: ../../processed_data/exac/exac_models_roc_curves_intron.txt
	@ echo "Graphing curves..."
	@ Rscript exac_graph_models.R
