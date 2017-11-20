all: results/subcellular_location.swissprot.tsv

.PHONY: clean
clean:
	rm -rfv results/*

results/subcellular_location.swissprot.tsv: data/uniprot_sprot.xml getSubcellularLocation.pl
	./getSubcellularLocation.pl $< > $@
