<?php

    $fh = fopen('../data/merged_locations.csv', 'r');

    $next_line = function ($fh) { return fgetcsv($fh, 0, "\t"); };

    $header = $next_line($fh);

    $proteins = [];

    while ($line = $next_line($fh)) {
        $data = array_combine($header, $line);
        $ID   = $data['UniProt_ID'];
        unset($data['UniProt_ID']);
        $proteins[$ID] = $data;
    }
    fclose($fh);

    echo count($proteins) . "\n";
    print_r(current($proteins));

    $fh        = fopen('data/uniprot-human.fasta', 'r');
    $fh        = fopen('data/human_proteome_cdhit.fasta', 'r');
    $skip      = true;
    $ID        = null;
    $sequences = [];
    while ($line = fgets($fh)) {
        if (strpos($line, ">") === 0) {
            # header line
            $ID = explode("|", $line)[1];

            $skip = array_key_exists($ID, $proteins) === false;
        }
        elseif (! $skip && $ID) {
            if (! array_key_exists($ID, $sequences))
                $sequences[$ID] = "";

            $sequences[$ID] .= trim($line);
        }
    }
    fclose($fh);

    echo count($sequences) . "\n";
    print_r(current($sequences));


    $fh        = fopen('../results/multiclass_data.csv', 'w+');
    $writeline = function ($fh, $array) { fputcsv($fh, $array, "\t"); };

    $writeline($fh, array_merge($header, ["sequence"]));
    $hist = array_fill(0, 29, 0);
    foreach (array_intersect_key($proteins, $sequences) as $ID => $data) {
        $writeline($fh, array_merge([$ID], $data, [$sequences[$ID]]));
        $hist[count(array_filter($data))]++;
    }
    fclose($fh);

    $fh = fopen('../results/multiclass_data_num_comp.csv', 'w+');
    $writeline($fh, ["num_compartments", "prot_count"]);
    foreach ($hist as $num => $count)
        $writeline($fh, [$num, $count]);


    #print_r(array_slice(array_keys(array_diff_key($proteins, $sequences)), 0, 20));