FROM ubuntu:18.04

# ADD . /pymoc

ARG NB_USER=nbuser
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
# WORKDIR /pymoc

RUN apt-get update -y && \
    apt-get install -y python3.6-dev python3-pip git curl
RUN pip3 install --no-cache --upgrade pip && \
    pip3 install numpy==1.13 && \
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
    pip3 install funcsigs && \
    pip3 install jupyter && \
    pip3 install --no-cache notebook
RUN ln -s /usr/bin/python3.6 /usr/bin/python

USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

CMD bash
