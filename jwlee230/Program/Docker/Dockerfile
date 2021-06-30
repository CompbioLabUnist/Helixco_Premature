FROM qiime2/core:2020.6
LABEL maintainer="jwlee230@unist.ac.kr"

RUN ln -sf /usr/share/zoneinfo/Asia/Seoul /etc/localtime

RUN apt-get update && apt-get upgrade -y && apt-get install -y make gcc g++ zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev unzip libcurl4-openssl-dev libfontconfig1

RUN mkdir /Tools
ADD Tools/* /Tools/

# PIP
ADD requirements.txt /
RUN pip install --requirement /requirements.txt && pip install scikit-bio

# LEfSe
RUN cd /Tools && bash Miniconda.sh -b
ENV PATH=/home/qiime2/anaconda2/bin:$PATH
RUN conda install -c biobakery lefse

CMD ["/bin/bash"]
