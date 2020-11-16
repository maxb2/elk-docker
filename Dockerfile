FROM ubuntu:20.04 as base
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qq update && apt-get -qq install build-essential liblapack-dev libblas-dev gfortran libomp-dev libopenmpi-dev python


FROM base as dep-builder
WORKDIR /app
COPY libxc-5.0.0 /app/libxc-5.0.0
COPY wannier90-3.1.0 /app/wannier90-3.1.0

RUN cd wannier90-3.1.0 && make && make lib && make test-serial && PREFIX='/app/usr' make install

RUN cd libxc-5.0.0 && ./configure --prefix='/app/usr' && make && make check && make install


FROM base as final

COPY elk-6.8.4 /app/elk-6.8.4

WORKDIR /app/elk-6.8.4
COPY --from=dep-builder /app/libxc-5.0.0/src/libxcf90.f90  /app/elk-6.8.4/src/libxcf90.f90
COPY --from=dep-builder /app/usr/lib/libxcf90.a            /app/elk-6.8.4/src/libxcf90.a
COPY --from=dep-builder /app/usr/lib/libxc.a               /app/elk-6.8.4/src/libxc.a
COPY --from=dep-builder /app/usr/lib/libwannier.a          /app/elk-6.8.4/src/libwannier.a
RUN make

CMD ['/bin/bash']

FROM final as test

RUN make test
