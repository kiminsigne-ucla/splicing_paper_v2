all: format_ref rescore_ref clean_data

format_ref: ../../ref/exac/exac_ref_formatted_converted.txt

rescore_ref: ../../ref/exac/exac_ref_rescored.txt

clean_data: ../../processed_data/exac/exac_data_clean.txt


../../ref/exac/exac_ref_formatted_converted.txt: ../../ref/exac/exac_ref_all.txt
	@ echo "Formatting reference..."
	@ Rscript exac_format_ref.R


../../ref/exac/exac_ref_rescored.txt: ../../ref/exac/exac_ref_formatted_converted.txt
	@ echo "Re-scoring reference file..."
	@ cd ../splicemod_src && \
	python ../score_lib.py $< --rev ./data/motifs/Ke2011/ \
	../../ref/Rosenberg_2015/exonic_mer6_scores.series 30 $@


../../processed_data/exac/exac_data_clean.txt: ../../ref/exac/exac_ref_formatted_converted.txt
	@ echo "Cleaning up data..."
	@ Rscript exac_data_clean.R