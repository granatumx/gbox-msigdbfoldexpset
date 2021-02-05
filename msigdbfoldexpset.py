#!/usr/bin/env python

from base64 import b64encode
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import time
import math
from granatum_sdk import Granatum
import networkx as nx
import re
from networkx.drawing.nx_agraph import write_dot
import os

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

    G = nx.MultiDiGraph()
    clusternames = list(clustersvsgenes.T.columns)
    individualclusters = [n[:n.index(" vs rest")] for n in clusternames if n.endswith("vs rest")]
    print(individualclusters, flush=True)
    for cl in individualclusters:
        G.add_node(cl)

    # {pathway : {"cluster1":score1, "cluster2":score2}, pathway2 : {}}
    resultsmap = {}
    relabels = {}
    keys = {}
    currentkeyindex = 0
    for gset in gsets:
        for cluster in clustercomparisonstotest:
            try:
                resultdf = clustersvsgenes.loc[cluster, gset["gene_ids"]]
                score = np.nanmean(resultdf)
                if score >= min_zscore:
                    keys[gset["name"]] = keys.get(gset["name"], currentkeyindex)
                    print("Score = {}".format(score))
                    olddict = resultsmap.get(gset["name"], {})
                    olddict[cluster] = score
                    resultsmap[gset["name"]] = olddict
                    from_to = re.split(' vs ', cluster)
                    if from_to[1] != 'rest':
                        G.add_weighted_edges_from([(from_to[0], from_to[1], score*10.0)], label=str(currentkeyindex), weight=score*10.0)
                    else:
                        relabel_dict = relabels.get(from_to[0], "")
                        if relabel_dict == "":
                            relabel_dict = from_to[0] + ": " + currentkeyindex
                        else:
                            relabel_dict = relabel_dict + ", " + currentkeyindex
                        relabels[from_to[0]] = relabel_dict
                    currentkeyindex = currentkeyindex + 1
            except Exception as inst:
                print("Key error with {}".format(gset["name"]), flush=True)
                print("Exception: {}".format(inst), flush=True)

    print("Relabels {}".format(relabels), flush=True)
    G = nx.relabel_nodes(G, relabels)
    pos = nx.spring_layout(G)
    edge_labels = nx.get_edge_attributes(G, 'label')
    write_dot(G, 'plot.dot')
    os.system("dot plot.dot -Tpng -Kcirco -Gdpi=600 > plot.png")
    with open('plot.png', "rb") as f:
        image_b64 = b64encode(f.read()).decode("utf-8")

    gn.results.append(
        {
                "type": "png",
                "width": 650,
                "height": 480,
                "description": 'Network of clusters based on expression',
                "data": image_b64,
            })

    for k, v in sorted(keys.items(), key=lambda item: item[1])
        gn.add_result("{}: {}".format(v, k), "markdown")

    timing = "* Finished differential expression sets step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")
    # gn.export(return_df.T.to_csv(), 'differential_gene_sets.csv', kind='raw', meta=None, raw=True)

    toc = time.perf_counter()
    time_passed = round(toc - tic, 2)

    timing = "* Finished differential expression sets step in {} seconds*".format(time_passed)
    gn.add_result(timing, "markdown")

    gn.commit()


if __name__ == '__main__':
    main()
