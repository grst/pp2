<?php

    $fh = fopen('../results/multiclass_data.csv', 'r');

    $next_line = function ($fh) { return fgetcsv($fh, 0, "\t"); };

    $header = $next_line($fh);

    $proteins = [];

    while ($line = $next_line($fh)) {
        $data          = array_combine($header, $line);
        $ID            = $data['UniProt_ID'];
        $proteins[$ID] = $data;
    }
    fclose($fh);


    $fh = fopen('../data/protVec_100d_3grams.csv', 'r');

    $next_line = function ($fh) { return fgetcsv($fh, 0, "\t"); };

    $protVec = [];

    while ($line = $next_line($fh)) {
        $protVec[$line[0]] = array_slice($line, 1);
    }
    fclose($fh);

    echo "#prot " . count($proteins) . "\n";
    echo "#vec " . count($protVec) . "\n";

    $fh = fopen('../results/multiclass_data_prot_vec.csv', 'w+');

    $writeline = function ($fh, $array) { fputcsv($fh, $array, "\t"); };
    $unknown   = 0;
    foreach ($proteins as $ID => $data) {
        $vec      = array_fill(0, 100, 0);
        $sequence = $data['sequence'];
        for ($i = 0; $i < (strlen($sequence) - 2); $i++) {
            $three_mer = substr($sequence, $i, 3);
            if (! array_key_exists($three_mer, $protVec)) {
                $three_mer = "<unk>";
                $unknown++;
            }

            foreach ($protVec[$three_mer] as $k => $value)
                $vec[$k] += $value;

        }
        $writeline($fh, array_merge([$ID], $vec));
    }
    fclose($fh);
    echo "#unkown " . $unknown . "\n";