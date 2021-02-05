#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import time
import math
from granatum_sdk import Granatum
import networkx as nx

def parse_gmt(gmt_str):
    lines = gmt_str.split("\n")
    gsets = []
    for line in lines:
        if line == "":
            continue

        name, url, *gene_ids = line.split("\t")
        gsets.append({"name": name, "url": url, "gene_ids": gene_ids})
    return gsets


gset_group_id_to_filename = {
    "kegg": "c2.cp.kegg.v7.2.symbols.gmt",
    "go": "c5.all.v7.2.symbols.gmt",
}


def load_gsets(gset_group_id):
    data_fn = gset_group_id_to_filename[gset_group_id]
    gsets = parse_gmt(open(data_fn, 'r').read())

    return gsets

def main():
    tic = time.perf_counter()

    gn = Granatum()

    clustersvsgenes = gn.pandas_from_assay(gn.get_import('clustersvsgenes'))
    gset_group_id = gn.get_arg('gset_group_id')
    min_zscore = gn.get_arg('min_zscore')

    clustercomparisonstotest = list(clustersvsgenes.index)

    # Load all gene sets
    gsets = load_gsets(gset_group_id)

    G = nx.Graph()
    node
    clusternames = list(clustersvsgenes.T.columns)
    individualclusters = [n[:n.index(" vs rest")] for n in clusternames if n.endswith("vs rest")]
    print(individualclusters, flush=True)

    # {pathway : {"cluster1":score1, "cluster2":score2}, pathway2 : {}}
    resultsmap = {}
    for gset in gsets:
        for cluster in clustercomparisonstotest:
            try:
                resultdf = clustersvsgenes.loc[cluster, gset["gene_ids"]]
                score = np.nanmean(resultdf)
                if score >= min_zscore:
                    print("Score = {}".format(score))
                    olddict = resultsmap.get(gset["name"], {})
                    olddict[cluster] = score
                    resultsmap[gset["name"]] = olddict
                    print("Done adding score", flush=True)
            except:
                print("Key error with {}".format(gset["name"]), flush=True)

    print(resultsmap, flush=True)
    time.sleep(5)
    
    # gn.export(return_df.T.to_csv(), 'differential_gene_sets.csv', kind='raw', meta=None, raw=True)

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished differential expression sets step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
