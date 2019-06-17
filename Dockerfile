FROM ubuntu:18.04

# ADD . /pymoc
WORKDIR /pymoc

RUN apt-get update -y && \
    apt-get install -y python3.6-dev python3-pip git
RUN pip3 install numpy scipy matplotlib pytest codecov
RUN ln -s /usr/bin/python3.6 /usr/bin/python

CMD bash