FROM qiime2/core:2020.6
LABEL maintainer="jwlee230@unist.ac.kr"

RUN ln -sf /usr/share/zoneinfo/Asia/Seoul /etc/localtime

RUN apt-get update && apt-get dist-upgrade -y && apt-get install -y make gcc g++ zlib1g-dev libbz2-dev liblzma-dev libncurses5-dev unzip libcurl4-openssl-dev libfontconfig1

RUN conda install -c cduvallet q2_perc_norm

# PIP
ADD requirements.txt /
RUN pip install --requirement /requirements.txt && pip install scikit-bio

CMD ["/bin/bash"]
