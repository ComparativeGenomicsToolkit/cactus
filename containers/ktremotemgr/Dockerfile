FROM ubuntu:16.04

RUN apt-get update && apt-get install -y python kyototycoon


# Set WORKDIR to /data -- predefined mount location.
RUN mkdir /data
WORKDIR /data

ENV PATH /usr/bin:/bin:/usr/sbin:/sbin:/usr/local/bin:/opt/bin

# And clean up
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["ktremotemgr"]
