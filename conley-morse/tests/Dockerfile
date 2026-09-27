FROM python:3.11-slim
RUN pip install --no-cache-dir pytest==8.4.1 pytest-json-ctrf==0.3.5 numpy==2.1.3 scipy==1.14.1 networkx==3.4.2
RUN mkdir -p /app /logs/verifier /logs/artifacts
COPY . /tests/
RUN chown -R root:root /tests \
 && chmod 0755 /tests \
 && chmod 0700 /tests/heldout /tests/oracles \
 && find /tests/heldout /tests/oracles -type f -exec chmod 0600 {} +
USER root
