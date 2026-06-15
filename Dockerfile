FROM ensemblorg/ensembl-vep:release_109.3 AS spython-base
ADD /mnt/hsf/AC_tests/morfal2/netstart-1.0c /opt/netstart-1.0c
ADD /mnt/hsf/AC_tests/morfal2/MORFAL_publi.py /opt/morfal/MORFAL.py
ENV PATH=/opt/python3.8/bin:/opt/netstart-1.0c:$PATH
ENV LC_ALL=C
ENV LANG=C
RUN export DEBIAN_FRONTEND=noninteractive
RUN sed -i 's:archive.ubuntu.com:fr.archive.ubuntu.com:gi' /etc/apt/sources.list || true
RUN apt-get update && apt-get install -y \
wget \
build-essential \
zlib1g-dev \
libncurses5-dev \
libgdbm-dev \
libnss3-dev \
libssl-dev \
libreadline-dev \
libffi-dev \
libsqlite3-dev \
libbz2-dev \
liblzma-dev \
tk-dev \
uuid-dev \
samtools \
tcsh \
ca-certificates \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*
RUN mkdir -p /opt/morfal /data /out /refs /tmp/python-build
RUN chmod -R 755 /opt/netstart-1.0c
RUN cd /tmp/python-build
RUN wget https://www.python.org/ftp/python/3.8.20/Python-3.8.20.tgz --no-check-certificate
RUN tar -xzf Python-3.8.20.tgz
RUN cd Python-3.8.20
RUN ./configure --prefix=/opt/python3.8 --enable-optimizations
RUN make -j"$(nproc)"
RUN make install
RUN /opt/python3.8/bin/python3.8 --version
RUN chmod 755 /opt/morfal/MORFAL.py
RUN cat > /usr/local/bin/morfal << 'EOF'
RUN exec /opt/python3.8/bin/python3.8 /opt/morfal/MORFAL.py "$@"
RUN EOF
RUN chmod 755 /usr/local/bin/morfal
RUN ln -sf /opt/python3.8/bin/python3.8 /usr/local/bin/python3
RUN ln -sf /opt/python3.8/bin/python3.8 /usr/local/bin/python
RUN rm -rf /tmp/python-build
CMD exec /usr/local/bin/morfal "$@"
