<?php

    #region parse
    $dom = new DomDocument;
    $dom->loadXml(file_get_contents('data/swissprot-locations-all.rdf'));
    $xph = new DOMXPath($dom);
    $xph->registerNamespace('', "http://www.w3.org/2005/Atom");
    $xph->registerNamespace('rdf', "http://www.w3.org/1999/02/22-rdf-syntax-ns#");

    $locations  = [];
    $hierarchy  = [];
    $root_nodes = [];
    foreach ($xph->query('//rdf:Description') as $location) {
        /** @var DOMNode $location */
        $about = $location->attributes->getNamedItem('about')->nodeValue;

        #no location => skip
        if (preg_match('!/locations/(\d+)!', $about, $matches) === 0)
            continue;

        #insert into location hashmap
        $id             = $matches[1];
        $locations[$id] = $location;

        #fill hierarchy by partOf and subClassOf nodes
        $is_root = true;
        foreach ($location->childNodes as $child_node) {
            /** @var DOMNode $child_node */
            if ($child_node->nodeName == "partOf" || $child_node->nodeName == "rdfs:subClassOf") {
                foreach ($child_node->attributes as $attribute)
                    /**@var DOMNode $attribute */
                    if (preg_match('!/locations/(\d+)!', $attribute->nodeValue, $matches) === 1) {
                        $hierarchy[$matches[1]][] = $id;
                        $is_root                  = false;
                    }
            }
        }

        #is root node?
        if ($is_root)
            $root_nodes[] = $id;
    }

    #remove root nodes having no children
    $filterted_root_nodes = [];
    foreach ($root_nodes as $root_node_id) {
        if (! array_key_exists($root_node_id, $hierarchy))
            continue;
        $filterted_root_nodes[] = $root_node_id;
    }
    $root_nodes = $filterted_root_nodes;
    #endregion


    //    foreach ($root_nodes as $root_node_id) {
    //        echo "$root_node_id\t" . getNodeName($locations[$root_node_id]) . "\n";
    //    }
    //    echo "------------------------\n";
    //    foreach (getAllChildren(204, $hierarchy) as $child_id) {
    //        echo "$child_id\t" . getNodeName($locations[$child_id]) . "\n";
    //    }

    printf("%s\t%s\t%s\n", 'deeploc_location', 'target_location', 'target_source');

    #region print-swissprot-mapping
    $deeploc_to_swissprot = [
        'Cell membrane'         => [39],
        'Nucleus'               => [191],
        'Cytoplasm'             => [86],
        'Extracellular'         => [243],
        'Mitochondrion'         => [173],
        'Endoplasmic reticulum' => [95],
        'Golgi apparatus'       => [132],
        'Lysosome/Vacuole'      => [158, 272],
        'Peroxisome'            => [202, 203],
    ];
    foreach ($deeploc_to_swissprot as $subcell_loc => $swissprot_ids) {
        foreach ($swissprot_ids as $swissprot_id) {
            foreach (getAllChildren($swissprot_id, $hierarchy) as $child_id) {
                printf("%s\t%s\t%s\n", $subcell_loc, getNodeName($locations[$child_id]), 'swissprot');
            }
        }
    }
    #endregion

    #region print-hpa-mapping

    $deeploc_to_hpa = [
        'Cell membrane'         => ['Plasma membrane'],
        'Nucleus'               => ['Nucleoplasm', 'Nuclear bodies', 'Nucleus', 'Nucleoli', 'Nuclear speckles', 'Nuclear membrane', 'Nucleoli fibrillar center'],
        'Cytoplasm'             => ['Cytosol', 'Cytoplasmic bodies', 'Actin filaments',
                                    'Centrosome', 'Microtubule organizing center',
                                    'Microtubules', 'Cytokinetic bridge', 'Intermediate filaments',
                                    'Mitotic spindle', 'Midbody', 'Microtubule ends', 'Aggresome',],
        'Extracellular'         => [],
        'Mitochondrion'         => ['Mitochondria'],
        'Endoplasmic reticulum' => ['Endoplasmic reticulum'],
        'Golgi apparatus'       => ['Golgi apparatus'],
        'Lysosome/Vacuole'      => ['Lysosomes'],
        'Peroxisome'            => ['Peroxisomes'],
    ];

    /**
     * UNMAPPED HPA LOCATIONS
     *
     * **********************
     *
     * Vesicles
     * Endosomes
     * Focal adhesion sites
     * Cell Junctions
     * Lipid droplets
     * Rods & Rings
     * Midbody ring
     */
    foreach ($deeploc_to_hpa as $subcell_loc => $hpa_locs) {
        foreach ($hpa_locs as $hpa_loc) {
            printf("%s\t%s\t%s\n", $subcell_loc, $hpa_loc, 'hpa');
        }
    }
    #endregion


    #region aux-functions
    /**
     * @param int $id
     * @param [] $hierarchy
     *
     * @return array
     */
    function getAllChildren($id, $hierarchy) {
        $children = [$id];
        if (! array_key_exists($id, $hierarchy))
            return $children;
        foreach ($hierarchy[$id] as $child) {
            $children = array_merge($children, getAllChildren($child, $hierarchy));
        }
        return $children;
    }

    function getNodeName($location_node) {
        foreach ($location_node->childNodes as $child_node) {
            if ($child_node->nodeName == "skos:prefLabel") {
                return $child_node->nodeValue;
            }
        }
        return "";
    }
    #endregion
