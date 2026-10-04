# Container for the SWEVO reproducibility package (1.0.0-rc).
# Build and run the full verification from the unpacked package folder:
#   docker build -t swevo-repro .
#   docker run --rm -v "$PWD":/pkg -w /pkg swevo-repro sh regenerate_all.sh /tmp/regen
# Optional PyWake cross-checks:  docker build --build-arg WITH_PYWAKE=1 -t swevo-repro-pywake .
# The verification environment was CPython 3.11.15 on Linux x86_64 (Ubuntu 24.04) with the NumPy/SciPy wheels
# below (bundled OpenBLAS 0.3.31). If the python:3.11.15-slim tag is no longer available, use python:3.11-slim and
# check `python --version`; all outputs except MS-SLSQP runs (SciPy SLSQP is platform dependent) are expected to
# reproduce bit for bit on x86_64 with the pinned wheels.
FROM python:3.11.15-slim
ARG WITH_PYWAKE=0
ENV PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
COPY requirements.txt requirements-pywake.txt /tmp/
RUN python -m pip install --no-cache-dir -r /tmp/requirements.txt && \
    if [ "$WITH_PYWAKE" = "1" ]; then python -m pip install --no-cache-dir -r /tmp/requirements-pywake.txt; fi
WORKDIR /pkg
CMD ["sh", "regenerate_all.sh", "/tmp/regen"]
