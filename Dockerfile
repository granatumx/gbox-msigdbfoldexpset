FROM granatumx/gbox-py-sdk:1.0.0

RUN pip install networkx

RUN apt-get install -y graphviz

RUN wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/c2.cp.kegg.v7.2.symbols.gmt
RUN wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.2/c5.all.v7.2.symbols.gmt

COPY . .

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
