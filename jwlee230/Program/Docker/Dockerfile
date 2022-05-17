FROM quay.io/qiime2/core:2021.4
LABEL maintainer="jwlee230@unist.ac.kr"

RUN ln -sf /usr/share/zoneinfo/Asia/Seoul /etc/localtime

RUN apt-get update --allow-releaseinfo-change && apt-get upgrade -y && apt-get install -y make gcc g++ zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev unzip libcurl4-openssl-dev libfontconfig1 r-base r-base-core r-recommended r-base-dev cabextract && wget "https://www.freedesktop.org/software/fontconfig/webfonts/webfonts.tar.gz" && tar -zxvf webfonts.tar.gz && cd msfonts && cabextract *.exe && mkdir -p /opt/conda/envs/qiime2-2021.4/lib/python3.8/site-packages/matplotlib/mpl-data/fonts/ttf && cp *.ttf *.TTF /opt/conda/envs/qiime2-2021.4/lib/python3.8/site-packages/matplotlib/mpl-data/fonts/ttf

# PIP
ADD requirements.txt /
RUN pip install --requirement /requirements.txt && pip install scikit-bio

# R Packages
RUN mkdir -p $HOME/.R && echo "MAKEFLAGS = -j" > $HOME/.R/Makevars && Rscript -e "install.packages(c('optparse', 'BiocManager', 'data.table', 'e1071', 'colorRamps', 'devtools', 'remotes', 'R.utils'), repos='http://cran.us.r-project.org')" && Rscript -e "BiocManager::install(c('DESeq2', 'limma'))"

CMD ["/bin/bash"]
