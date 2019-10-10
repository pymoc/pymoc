FROM ubuntu:18.04

# ADD . /pymoc
WORKDIR /pymoc

RUN apt-get update -y && \
    apt-get install -y python3.6-dev python3-pip git curl
RUN pip3 install numpy==1.13 && \
    pip3 install scipy==1.3 && \
    pip3 install matplotlib && \
    pip3 install pytest==5.1.1 && \
    pip3 install codecov && \
    pip3 install pytest-cov && \
    pip3 install sphinx && \
    pip3 install yapf && \
    pip3 install futures && \
    pip3 install recommonmark && \
    pip3 install sphinx-rtd-theme && \
    pip3 install jupyter-sphinx-theme && \
    pip3 install funcsigs
RUN ln -s /usr/bin/python3.6 /usr/bin/python

CMD bash
