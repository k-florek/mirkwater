FROM python:3

# metadata
LABEL base.image="python:3"
LABEL dockerfile.version="1"
LABEL maintainer="Kelsey Florek"
LABEL maintainer.email="Kelsey.florek@slh.wisc.edu"

WORKDIR /usr/src/app

COPY requirements.txt ./
RUN pip install --no-cache-dir -r requirements.txt

COPY partial_matching.py ./
COPY ./constellations/constellations/definitions ./definitions
