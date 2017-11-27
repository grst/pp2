from collections import defaultdict

# read id mappings, required for deeploc_data.fasta
file = open("./data/biomart_hgnc_uniprot.tsv", "r")
id_mapping = {}
for line in file:
    cols = line[:-1].split("\t")
    if (len(cols) != 2) | (not cols[1]) | (not cols[0]):
        continue
    else:
        id_mapping[cols[1]] = cols[0]

# training and test data used by DeepLoc, available at http://www.cbs.dtu.dk/services/DeepLoc/data.php
# save list of IDs for each DeepLoc location
file = open("./data/deeploc_data.fasta", "r")
deeploc_location_id = defaultdict(list)
for line in file:
    if line.startswith('>'):
        cols = line.split(" ")
        uniprot_id = cols[0][1:]
        if uniprot_id in id_mapping:
            hgnc = id_mapping[uniprot_id]
            location = cols[1].split("-")[0].replace('.', ' ')
            deeploc_location_id[location].append(hgnc)

# swissprot location for each ID
file = open("./results/swissprot_filtered.tsv", "r")
swissprot_id_location = {}
for line in file:
    cols = line.split("\t")
    swissprot_id_location[cols[2]] = cols[3]

# hpa location for each ID
file = open("./results/hpa_filtered.tsv", "r")
hpa_id_location = {}
for line in file:
    cols = line.split("\t")
    hpa_id_location[cols[0]] = cols[2][:-1]

# for every DeepLoc location we look up the location of every ID assigned to it
# thus a location (Swissprot / HPA) can be mapped to multiple DeepLoc locations
deeploc_swissprot_mapping = defaultdict(set)
deeploc_hpa_mapping = defaultdict(set)
for location in deeploc_location_id:
    for hgnc in deeploc_location_id[location]:
        if hgnc in swissprot_id_location:
            deeploc_swissprot_mapping[location].add(swissprot_id_location[hgnc])
        elif hgnc in hpa_id_location:
            deeploc_hpa_mapping[location].add(hpa_id_location[hgnc])

# write output file
file = open("./results/location_mapping.tsv", "w")
file.write("deeploc_location\ttarget_location\ttarget_source\n")
for deeploc_location in deeploc_swissprot_mapping:
    for swissprot_location in deeploc_swissprot_mapping[deeploc_location]:
        file.write(deeploc_location + "\t" + swissprot_location + "\tswissprot\n")
for deeploc_location in deeploc_hpa_mapping:
    for hpa_location in deeploc_hpa_mapping[deeploc_location]:
        file.write(deeploc_location + "\t" + hpa_location + "\thpa\n")

